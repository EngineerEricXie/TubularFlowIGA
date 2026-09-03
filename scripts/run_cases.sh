#!/usr/bin/env bash

set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
config_file=$repo_dir/execution.conf
requested_cases=()
built_cpu_3d=false
built_cuda_3d=false
built_cpu_1d=false

Usage()
{
	cat <<USAGE
usage: $0 [--config FILE] [CASE ...]

Validate and optionally prepare and solve one or more 1D or 3D cases. CASE
values are directory names directly below CASE_ROOT. With no CASE arguments,
every immediate subdirectory of CASE_ROOT is run.

Examples:
  $0 MyCase
  $0 --config execution.conf CaseA CaseB
USAGE
}

Die()
{
	printf 'error: %s\n' "$*" >&2
	exit 2
}

IsTruthy()
{
	case ${1,,} in
		1|true|yes|on) return 0 ;;
		0|false|no|off) return 1 ;;
		*) Die "expected a boolean value, got: $1" ;;
	esac
}

Trim()
{
	local value=$1
	value="${value#"${value%%[![:space:]]*}"}"
	value="${value%"${value##*[![:space:]]}"}"
	printf '%s' "$value"
}

ResolvePath()
{
	if [[ $1 == /* ]]; then realpath -m "$1"; else realpath -m "$repo_dir/$1"; fi
}

ParseArgs()
{
	while (( $# > 0 )); do
		case $1 in
			-h|--help) Usage; exit 0 ;;
			--config)
				(( $# >= 2 )) || Die "--config requires a file"
				config_file=$2
				shift 2
				;;
			--config=*) config_file=${1#*=}; shift ;;
			--*) Die "unknown option: $1" ;;
			*) requested_cases+=("$1"); shift ;;
		esac
	done
}

LoadConfig()
{
	config_file=$(realpath -m "$config_file")
	[[ -f $config_file ]] || Die "configuration file not found: $config_file"

	CASE_ROOT=Input
	OUTPUT_ROOT=
	RANKS=2
	BACKEND=cpu
	BUILD_SOLVERS=1
	CLEAN=0
	RUN_MESH_CHECK=1
	RUN_SOLVER=1
	RUN_VALIDATION=1
	SOLVER=auto
	SYSTEM=
	MPIEXEC=mpiexec
	OMP_NUM_THREADS=2
	SOLVER_ARGS=
	DRY_RUN=0
	PETSC_DIR=${PETSC_DIR:-}
	PETSC_ARCH=${PETSC_ARCH:-}
	HDF5_CFLAGS=${HDF5_CFLAGS:-}
	HDF5_LIBS=${HDF5_LIBS:-}

	local line key value line_number=0
	while IFS= read -r line || [[ -n $line ]]; do
		((line_number += 1))
		line=$(Trim "$line")
		[[ -z $line || ${line:0:1} == '#' ]] && continue
		[[ $line == *=* ]] || Die "$config_file:$line_number: expected KEY=VALUE"
		key=$(Trim "${line%%=*}")
		value=$(Trim "${line#*=}")
		[[ $key =~ ^[A-Z][A-Z0-9_]*$ ]] \
			|| Die "$config_file:$line_number: invalid key: $key"
		case $key in
			CASE_ROOT|OUTPUT_ROOT|RANKS|BACKEND|BUILD_SOLVERS|CLEAN|RUN_MESH_CHECK|RUN_SOLVER|RUN_VALIDATION|SOLVER|SYSTEM|MPIEXEC|OMP_NUM_THREADS|SOLVER_ARGS|DRY_RUN|PETSC_DIR|PETSC_ARCH|HDF5_CFLAGS|HDF5_LIBS) ;;
			*) Die "$config_file:$line_number: unknown setting: $key" ;;
		esac
		if (( ${#value} >= 2 )); then
			if [[ ${value:0:1} == '"' && ${value: -1} == '"' ]] \
				|| [[ ${value:0:1} == "'" && ${value: -1} == "'" ]]; then
				value=${value:1:${#value}-2}
			fi
		fi
		printf -v "$key" '%s' "$value"
	done < "$config_file"
}

ValidateConfig()
{
	[[ $RANKS =~ ^[0-9]+$ ]] && (( RANKS >= 1 )) \
		|| Die "RANKS must be a positive integer"
	[[ $OMP_NUM_THREADS =~ ^[0-9]+$ ]] && (( OMP_NUM_THREADS >= 1 )) \
		|| Die "OMP_NUM_THREADS must be a positive integer"
	case $BACKEND in cpu|cuda) ;; *) Die "BACKEND must be cpu or cuda; got: $BACKEND" ;; esac
	case $SOLVER in
		auto|navier_stokes|transport) ;;
		*) Die "SOLVER must be auto, navier_stokes, or transport; got: $SOLVER" ;;
	esac
	IsTruthy "$BUILD_SOLVERS" || true
	IsTruthy "$CLEAN" || true
	IsTruthy "$RUN_MESH_CHECK" || true
	IsTruthy "$RUN_SOLVER" || true
	IsTruthy "$RUN_VALIDATION" || true
	IsTruthy "$DRY_RUN" || true
}

PrintCommand()
{
	printf '+'
	printf ' %q' "$@"
	printf '\n'
}

Run()
{
	PrintCommand "$@"
	if ! IsTruthy "$DRY_RUN"; then "$@"; fi
}

RunMpi()
{
	local launcher=()
	read -r -a launcher <<< "$MPIEXEC"
	(( ${#launcher[@]} > 0 )) || Die "MPIEXEC cannot be empty"
	Run "${launcher[@]}" -np "$RANKS" "$@"
}

SelectCases()
{
	local case_root=$1 case_name
	case_dirs=()
	if (( ${#requested_cases[@]} > 0 )); then
		for case_name in "${requested_cases[@]}"; do
			[[ $case_name =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ ]] \
				|| Die "invalid case directory name: $case_name"
			[[ -d $case_root/$case_name ]] \
				|| Die "case directory not found: $case_root/$case_name"
			case_dirs+=("$(realpath "$case_root/$case_name")")
		done
	else
		while IFS= read -r -d '' case_name; do
			case_dirs+=("$(realpath "$case_name")")
		done < <(find "$case_root" -mindepth 1 -maxdepth 1 -type d ! -name '.*' -print0 | sort -z)
	fi
	(( ${#case_dirs[@]} > 0 )) || Die "no case directories found under $case_root"
}

EnsureConfigChecker()
{
	local checker=$repo_dir/solvers/cpu/iga_config_check
	[[ -x $checker ]] && return 0
	if IsTruthy "$DRY_RUN"; then
		Die "iga_config_check is unavailable; run 'make cpu' before DRY_RUN=1"
	fi
	Run make -C "$repo_dir/solvers/cpu" iga_config_check
}

InspectExecutionPlan()
{
	local simulation_config=$1 checker=$repo_dir/solvers/cpu/iga_config_check
	local plan line name kind
	if ! plan=$("$checker" "$simulation_config" --execution-plan); then
		Die "invalid simulation configuration: $simulation_config"
	fi
	plan_dimension=
	plan_requires_mesh=
	plan_coupling_mode=
	plan_backends=()
	plan_navier_systems=()
	plan_transport_systems=()
	plan_one_d_flow_systems=()
	while IFS= read -r line; do
		case $line in
			dimension=*) plan_dimension=${line#*=} ;;
			requires_mesh=*) plan_requires_mesh=${line#*=} ;;
			coupling_mode=*) plan_coupling_mode=${line#*=} ;;
			backend=*) plan_backends+=("${line#*=}") ;;
			system=*)
				if [[ $line =~ ^system=([^[:space:]]+)[[:space:]]+kind=([^[:space:]]+) ]]; then
					name=${BASH_REMATCH[1]}
					kind=${BASH_REMATCH[2]}
					case $kind in
						navier_stokes) plan_navier_systems+=("$name") ;;
						linear_transport) plan_transport_systems+=("$name") ;;
						network_flow_1d) plan_one_d_flow_systems+=("$name") ;;
					esac
				fi
				;;
		esac
	done <<< "$plan"
	[[ $plan_dimension == 1d || $plan_dimension == 3d ]] \
		|| Die "execution plan has unsupported dimension: $plan_dimension"
	local supported=false backend
	for backend in "${plan_backends[@]}"; do
		[[ $backend == "$BACKEND" ]] && supported=true
	done
	[[ $supported == true ]] || Die \
		"BACKEND=$BACKEND is unsupported for dimension=$plan_dimension coupling=$plan_coupling_mode"
}

SelectSystem()
{
	selected_solver=
	selected_system=$SYSTEM
	if [[ $plan_dimension == 1d ]]; then
		[[ $SOLVER == auto ]] \
			|| Die "1D cases use the coupled iga_1d flow driver and require SOLVER=auto"
		if [[ -n $selected_system ]]; then
			local found=false name
			for name in "${plan_one_d_flow_systems[@]}"; do
				[[ $name == "$selected_system" ]] && found=true
			done
			[[ $found == true ]] || Die "1D flow SYSTEM '$selected_system' was not found"
		else
			(( ${#plan_one_d_flow_systems[@]} == 1 )) \
				|| Die "1D cases with multiple flow systems require SYSTEM"
			selected_system=${plan_one_d_flow_systems[0]}
		fi
		selected_solver=one_d
		return 0
	fi

	selected_solver=$SOLVER
	local name kind=
	if [[ -n $selected_system ]]; then
		for name in "${plan_navier_systems[@]}"; do
			[[ $name == "$selected_system" ]] && kind=navier_stokes
		done
		for name in "${plan_transport_systems[@]}"; do
			[[ $name == "$selected_system" ]] && kind=transport
		done
		[[ -n $kind ]] || Die "SYSTEM '$selected_system' was not found"
		if [[ $selected_solver == auto ]]; then selected_solver=$kind; fi
		[[ $selected_solver == "$kind" ]] \
			|| Die "SYSTEM '$selected_system' does not match SOLVER=$selected_solver"
	fi
	if [[ $selected_solver == auto ]]; then
		if (( ${#plan_navier_systems[@]} == 1 )); then
			selected_solver=navier_stokes
			selected_system=${plan_navier_systems[0]}
		elif (( ${#plan_navier_systems[@]} > 1 )); then
			Die "multiple Navier-Stokes systems are not supported by the current solver"
		elif (( ${#plan_transport_systems[@]} == 1 )); then
			selected_solver=transport
			selected_system=${plan_transport_systems[0]}
		elif (( ${#plan_transport_systems[@]} > 1 )); then
			Die "multiple transport systems require SYSTEM"
		else
			Die "configuration has no supported 3D equation system"
		fi
	elif [[ $selected_solver == navier_stokes ]]; then
		if [[ -z $selected_system ]]; then
			(( ${#plan_navier_systems[@]} == 1 )) \
				|| Die "SOLVER=navier_stokes requires SYSTEM when the choice is not unique"
			selected_system=${plan_navier_systems[0]}
		fi
	else
		if [[ -z $selected_system ]]; then
			(( ${#plan_transport_systems[@]} == 1 )) \
				|| Die "SOLVER=transport requires SYSTEM when the choice is not unique"
			selected_system=${plan_transport_systems[0]}
		fi
	fi
}

BuildRuntime()
{
	if ! IsTruthy "$RUN_SOLVER" && ! IsTruthy "$RUN_MESH_CHECK"; then return 0; fi
	IsTruthy "$BUILD_SOLVERS" || return 0
	local make_args=()
	[[ -z $HDF5_CFLAGS ]] || make_args+=("HDF5_CFLAGS=$HDF5_CFLAGS")
	[[ -z $HDF5_LIBS ]] || make_args+=("HDF5_LIBS=$HDF5_LIBS")
	if [[ $plan_dimension == 1d ]]; then
		[[ $built_cpu_1d == false ]] || return 0
		[[ -z $PETSC_DIR ]] || make_args+=("PETSC_DIR=$PETSC_DIR" "PETSC_ARCH=$PETSC_ARCH")
		Run make -C "$repo_dir/solvers/one_d" petsc "${make_args[@]}"
		built_cpu_1d=true
	elif [[ $BACKEND == cpu ]]; then
		[[ $built_cpu_3d == false ]] || return 0
		if [[ -z $PETSC_DIR ]] && ! IsTruthy "$DRY_RUN"; then
			Die "PETSC_DIR must be set in $config_file or exported when building the 3D CPU solver"
		fi
		[[ -z $PETSC_DIR ]] || make_args+=("PETSC_DIR=$PETSC_DIR" "PETSC_ARCH=$PETSC_ARCH")
		Run make -C "$repo_dir/solvers/cpu" petsc "${make_args[@]}"
		built_cpu_3d=true
	else
		[[ $built_cuda_3d == false ]] || return 0
		Run make -C "$repo_dir/solvers/cuda" "${make_args[@]}"
		built_cuda_3d=true
	fi
}

GeneratedDirectory()
{
	local source_dir=$1 output_root=$2
	local candidate
	if [[ -n $output_root ]]; then
		candidate=$output_root/${source_dir##*/}
	else
		candidate=$source_dir/generated
	fi
	[[ ! -L $candidate ]] || Die "refusing to replace a symbolic-link output: $candidate"
	realpath -m "$candidate"
}

ValidateGeneratedTarget()
{
	local source_dir=$1 generated_dir=$2 case_name=$3 manifest_case= manifest_path=
	[[ $generated_dir != / && $generated_dir != "$repo_dir" && $generated_dir != "$source_dir" ]] \
		|| Die "unsafe generated-output directory: $generated_dir"
	if IsTruthy "$CLEAN" && [[ -d $generated_dir ]] \
		&& [[ -n $(find "$generated_dir" -mindepth 1 -maxdepth 1 -print -quit) ]]; then
		if [[ -f $generated_dir/manifest.json ]]; then
			manifest_path=$generated_dir/manifest.json
		elif [[ -f $generated_dir/run_manifest.json ]]; then
			manifest_path=$generated_dir/run_manifest.json
		else
			Die "refusing to clean an unrecognized directory without a manifest: $generated_dir"
		fi
		manifest_case=$(sed -n 's/^[[:space:]]*"case_name":[[:space:]]*"\([^"]*\)".*$/\1/p' \
			"$manifest_path")
		[[ $manifest_case == "$case_name" ]] \
			|| Die "refusing to clean output whose manifest belongs to '$manifest_case': $generated_dir"
	fi
}

JsonEscape()
{
	local value=$1
	value=${value//\\/\\\\}
	value=${value//\"/\\\"}
	value=${value//$'\n'/\\n}
	value=${value//$'\r'/\\r}
	value=${value//$'\t'/\\t}
	printf '%s' "$value"
}

WriteRunManifest()
{
	local path=$1 case_name=$2 result=$3 validation=$4
	{
		printf '{\n'
		printf '  "schema_version": 1,\n'
		printf '  "case_name": "%s",\n' "$(JsonEscape "$case_name")"
		printf '  "dimension": "%s",\n' "$plan_dimension"
		printf '  "backend": "%s",\n' "$BACKEND"
		printf '  "ranks": %s,\n' "$RANKS"
		printf '  "primary_system": "%s",\n' "$(JsonEscape "$selected_system")"
		printf '  "solver": "%s",\n' "$selected_solver"
		printf '  "coupling_mode": "%s",\n' "$plan_coupling_mode"
		printf '  "result": "%s",\n' "$(JsonEscape "$result")"
		printf '  "validation": "%s",\n' "$validation"
		printf '  "completed": true\n'
		printf '}\n'
	} > "$path"
}

RunOneDCase()
{
	local source_dir=$1 output_root=$2
	local case_name=${source_dir##*/}
	local generated_dir
	generated_dir=$(GeneratedDirectory "$source_dir" "$output_root")
	if IsTruthy "$RUN_SOLVER"; then ValidateGeneratedTarget "$source_dir" "$generated_dir" "$case_name"; fi
	if IsTruthy "$RUN_MESH_CHECK"; then
		RunMpi "$repo_dir/solvers/one_d/iga_1d" "$source_dir" \
			--system "$selected_system" --check
	fi
	IsTruthy "$RUN_SOLVER" || return 0

	local extra_args=()
	if [[ -n $SOLVER_ARGS ]]; then read -r -a extra_args <<< "$SOLVER_ARGS"; fi
	local final_results=$generated_dir/results/$selected_system
	if IsTruthy "$DRY_RUN"; then
		RunMpi "$repo_dir/solvers/one_d/iga_1d" "$source_dir" \
			--system "$selected_system" --output-dir "$final_results" "${extra_args[@]}"
		return 0
	fi
	if [[ -d $generated_dir && -n $(find "$generated_dir" -mindepth 1 -maxdepth 1 -print -quit) ]] \
		&& ! IsTruthy "$CLEAN"; then
		Die "generated-output directory is not empty; set CLEAN=1 to replace it: $generated_dir"
	fi

	local output_parent=${generated_dir%/*} output_base=${generated_dir##*/}
	local staging_dir=$output_parent/.${output_base}.staging.$BASHPID
	local backup_dir=$output_parent/.${output_base}.previous.$BASHPID
	[[ ! -e $staging_dir && ! -e $backup_dir ]] \
		|| Die "temporary output path already exists for $generated_dir"
	(
		set -e
		trap 'rm -rf -- "$staging_dir"' EXIT
		mkdir -p "$staging_dir/results/$selected_system"
		RunMpi "$repo_dir/solvers/one_d/iga_1d" "$source_dir" \
			--system "$selected_system" \
			--output-dir "$staging_dir/results/$selected_system" "${extra_args[@]}"
		validation=not_requested
		if IsTruthy "$RUN_VALIDATION"; then
			[[ -s $staging_dir/results/$selected_system/summary.json ]] \
				|| Die "1D solver did not create summary.json"
			validation=summary_present
		fi
		WriteRunManifest "$staging_dir/run_manifest.json" "$case_name" \
			"results/$selected_system" "$validation"
		cp "$source_dir/simulation_config.json" "$staging_dir/"
		if [[ -e $generated_dir ]]; then mv "$generated_dir" "$backup_dir"; fi
		if mv "$staging_dir" "$generated_dir"; then
			[[ ! -e $backup_dir ]] || rm -rf -- "$backup_dir"
		else
			[[ ! -e $backup_dir ]] || mv "$backup_dir" "$generated_dir"
			exit 1
		fi
		trap - EXIT
	)
	printf '1D case completed\ncase: %s\ngenerated root: %s\nrun manifest: %s\n' \
		"$source_dir" "$generated_dir" "$generated_dir/run_manifest.json"
}

RunThreeDCase()
{
	local source_dir=$1 output_root=$2
	local case_name=${source_dir##*/}
	(( RANKS >= 2 )) || Die "3D generation requires RANKS of at least 2 for mpmetis"
	local generated_dir
	generated_dir=$(GeneratedDirectory "$source_dir" "$output_root")
	ValidateGeneratedTarget "$source_dir" "$generated_dir" "$case_name"
	local database=$generated_dir/database/$case_name-$RANKS.ntiga
	local runtime_case=$generated_dir/preprocessing
	local results_dir=$generated_dir/results/$selected_system
	local generate_args=("$source_dir" --output "$generated_dir" --ranks "$RANKS")
	if IsTruthy "$CLEAN"; then generate_args+=(--clean); fi
	Run env "OMP_NUM_THREADS=$OMP_NUM_THREADS" \
		"$repo_dir/scripts/generate_case.sh" "${generate_args[@]}"
	if ! IsTruthy "$DRY_RUN"; then
		[[ -s $database ]] || Die "generated database not found: $database"
		mkdir -p "$results_dir"
	fi
	if IsTruthy "$RUN_MESH_CHECK"; then
		if [[ $BACKEND == cpu ]]; then
			RunMpi "$repo_dir/solvers/cpu/iga_mesh_check" "$database"
		else
			Run "$repo_dir/solvers/cuda/iga_cuda" mesh-check "$database"
		fi
	fi
	IsTruthy "$RUN_SOLVER" || return 0

	local extra_args=()
	if [[ -n $SOLVER_ARGS ]]; then read -r -a extra_args <<< "$SOLVER_ARGS"; fi
	local result=$results_dir/${selected_solver}-${BACKEND}.txt
	case $BACKEND:$selected_solver in
		cpu:navier_stokes)
			RunMpi "$repo_dir/solvers/cpu/iga_navier_stokes" \
				"$database" "$runtime_case" --system "$selected_system" \
				--output "$result" "${extra_args[@]}"
			;;
		cpu:transport)
			RunMpi "$repo_dir/solvers/cpu/iga_solve" \
				"$database" "$runtime_case" --system "$selected_system" \
				--output "$result" "${extra_args[@]}"
			;;
		cuda:navier_stokes)
			Run "$repo_dir/solvers/cuda/iga_cuda" navier-stokes \
				"$database" "$runtime_case" --system "$selected_system" \
				--output "$result" "${extra_args[@]}"
			;;
		cuda:transport)
			Run "$repo_dir/solvers/cuda/iga_cuda" solve \
				"$database" "$runtime_case" --system "$selected_system" \
				--output "$result" "${extra_args[@]}"
			;;
	esac
	local validation=not_available
	if IsTruthy "$RUN_VALIDATION" && [[ $selected_solver == navier_stokes ]]; then
		Run "$repo_dir/solvers/cpu/iga_flow_validate" "$database" "$result"
		validation=flow_balance
	elif ! IsTruthy "$RUN_VALIDATION"; then
		validation=not_requested
	fi
	if ! IsTruthy "$DRY_RUN"; then
		WriteRunManifest "$results_dir/run_manifest.json" "$case_name" \
			"${result#$generated_dir/}" "$validation"
	fi
}

RunCase()
{
	local source_dir=$1 output_root=$2
	local case_name=${source_dir##*/}
	[[ -f $source_dir/simulation_config.json ]] \
		|| Die "$source_dir is missing simulation_config.json"
	InspectExecutionPlan "$source_dir/simulation_config.json"
	SelectSystem
	BuildRuntime
	printf '\n== Case: %s (%s, %s, coupling=%s) ==\n' \
		"$case_name" "$plan_dimension" "$BACKEND" "$plan_coupling_mode"
	if [[ $plan_dimension == 1d ]]; then
		RunOneDCase "$source_dir" "$output_root"
	else
		RunThreeDCase "$source_dir" "$output_root"
	fi
}

Main()
{
	ParseArgs "$@"
	LoadConfig
	ValidateConfig
	local case_root output_root= case_dir
	case_root=$(ResolvePath "$CASE_ROOT")
	[[ -d $case_root ]] || Die "CASE_ROOT does not exist: $case_root"
	if [[ -n $OUTPUT_ROOT ]]; then
		output_root=$(ResolvePath "$OUTPUT_ROOT")
		[[ $output_root != / && $output_root != "$repo_dir" ]] \
			|| Die "unsafe OUTPUT_ROOT: $output_root"
		if ! IsTruthy "$DRY_RUN"; then mkdir -p "$output_root"; fi
	fi
	SelectCases "$case_root"
	EnsureConfigChecker
	printf 'execution profile: %s\ncase root:         %s\noutput mode:       %s\nbackend:           %s\nranks:             %s\ncases:             %s\n' \
		"$config_file" "$case_root" "${output_root:-CASE/generated}" \
		"$BACKEND" "$RANKS" "${#case_dirs[@]}"
	for case_dir in "${case_dirs[@]}"; do RunCase "$case_dir" "$output_root"; done
	printf '\ncompleted %s case(s)\n' "${#case_dirs[@]}"
}

Main "$@"
