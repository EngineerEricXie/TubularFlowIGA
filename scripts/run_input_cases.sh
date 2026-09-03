#!/usr/bin/env bash

set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
config_file=$repo_dir/TFI_input.inp
requested_cases=()

Usage()
{
	cat <<USAGE
usage: $0 [--config FILE] [CASE ...]

Prepare and run case directories found directly under Input/.
With no CASE arguments, every immediate subdirectory is run. CASE arguments
are directory names such as AGH50Bifurcation.

Examples:
  $0
  $0 AGH50Bifurcation
  $0 --config /path/to/TFI_input.inp CaseA CaseB
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

ParseArgs()
{
	while (( $# > 0 )); do
		case $1 in
			-h|--help)
				Usage
				exit 0
				;;
			--config)
				(( $# >= 2 )) || Die "--config requires a file"
				config_file=$2
				shift 2
				;;
			--config=*)
				config_file=${1#*=}
				shift
				;;
			--*) Die "unknown option: $1" ;;
			*)
				requested_cases+=("$1")
				shift
				;;
		esac
	done
}

LoadConfig()
{
	[[ -f $config_file ]] || Die "configuration file not found: $config_file"

	PROCESSORS=2
	INPUT_DIR=Input
	OUTPUT_DIR=Output
	BUILD_SOLVERS=1
	CLEAN_OUTPUT=0
	RUN_MESH_CHECK=1
	RUN_SOLVER=1
	SOLVER=auto
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
			PROCESSORS|INPUT_DIR|OUTPUT_DIR|BUILD_SOLVERS|CLEAN_OUTPUT|RUN_MESH_CHECK|RUN_SOLVER|SOLVER|MPIEXEC|OMP_NUM_THREADS|SOLVER_ARGS|DRY_RUN|PETSC_DIR|PETSC_ARCH|HDF5_CFLAGS|HDF5_LIBS) ;;
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

ResolvePath()
{
	if [[ $1 == /* ]]; then
		printf '%s\n' "$1"
	else
		printf '%s/%s\n' "$repo_dir" "$1"
	fi
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
	if ! IsTruthy "$DRY_RUN"; then
		"$@"
	fi
}

BuildSolvers()
{
	IsTruthy "$BUILD_SOLVERS" || return 0
	IsTruthy "$RUN_SOLVER" || return 0
	[[ -n $PETSC_DIR ]] || Die "PETSC_DIR must be set in $config_file or exported"

	local make_args=("PETSC_DIR=$PETSC_DIR" "PETSC_ARCH=$PETSC_ARCH")
	[[ -z $HDF5_CFLAGS ]] || make_args+=("HDF5_CFLAGS=$HDF5_CFLAGS")
	[[ -z $HDF5_LIBS ]] || make_args+=("HDF5_LIBS=$HDF5_LIBS")
	Run make -C "$repo_dir/solvers/cpu" petsc "${make_args[@]}"
}

SelectCases()
{
	local input_root=$1 case_name
	case_dirs=()
	if (( ${#requested_cases[@]} > 0 )); then
		for case_name in "${requested_cases[@]}"; do
			[[ $case_name =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ ]] \
				|| Die "invalid case directory name: $case_name"
			[[ -d $input_root/$case_name ]] \
				|| Die "case directory not found: $input_root/$case_name"
			case_dirs+=("$input_root/$case_name")
		done
	else
		while IFS= read -r -d '' case_name; do
			case_dirs+=("$case_name")
		done < <(find "$input_root" -mindepth 1 -maxdepth 1 -type d ! -name '.*' -print0 | sort -z)
	fi
	(( ${#case_dirs[@]} > 0 )) || Die "no case directories found under $input_root"
}

RunMpi()
{
	local launcher=()
	read -r -a launcher <<< "$MPIEXEC"
	(( ${#launcher[@]} > 0 )) || Die "MPIEXEC cannot be empty"
	Run "${launcher[@]}" -np "$PROCESSORS" "$@"
}

PrepareCase()
{
	local source_dir=$1 output_root=$2
	local case_name=${source_dir##*/}
	local generated_dir=$output_root/$case_name
	local database=$generated_dir/database/$case_name-$PROCESSORS.ntiga
	local runtime_case=$generated_dir/preprocessing
	local results_dir=$generated_dir/results
	local generate_args=("$source_dir" --output "$generated_dir" --ranks "$PROCESSORS")

	[[ -f $source_dir/simulation_config.json ]] \
		|| Die "$source_dir is missing simulation_config.json"
	if IsTruthy "$CLEAN_OUTPUT"; then
		generate_args+=(--clean)
	fi

	printf '\n== Case: %s ==\n' "$case_name"
	Run env "OMP_NUM_THREADS=$OMP_NUM_THREADS" \
		"$repo_dir/scripts/generate_case.sh" "${generate_args[@]}"

	IsTruthy "$RUN_SOLVER" || return 0
	if ! IsTruthy "$DRY_RUN"; then
		[[ -s $database ]] || Die "generated database not found: $database"
		mkdir -p "$results_dir"
	fi

	if IsTruthy "$RUN_MESH_CHECK"; then
		RunMpi "$repo_dir/solvers/cpu/iga_mesh_check" "$database"
	fi

	local selected_solver=$SOLVER
	if [[ $selected_solver == auto ]]; then
		if grep -Eq '"kind"[[:space:]]*:[[:space:]]*"navier_stokes"' \
			"$source_dir/simulation_config.json"; then
			selected_solver=navier_stokes
		else
			selected_solver=transport
		fi
	fi

	local extra_args=()
	if [[ -n $SOLVER_ARGS ]]; then
		read -r -a extra_args <<< "$SOLVER_ARGS"
	fi
	case $selected_solver in
		navier_stokes)
			RunMpi "$repo_dir/solvers/cpu/iga_navier_stokes" \
				"$database" "$runtime_case" \
				--output "$results_dir/velocity-cpu.txt" "${extra_args[@]}"
			Run "$repo_dir/solvers/cpu/iga_flow_validate" \
				"$database" "$results_dir/velocity-cpu.txt"
			;;
		transport)
			RunMpi "$repo_dir/solvers/cpu/iga_solve" \
				"$database" "$runtime_case" \
				--output "$results_dir/transport-cpu.txt" "${extra_args[@]}"
			;;
		*) Die "SOLVER must be auto, navier_stokes, or transport; got: $selected_solver" ;;
	esac
}

Main()
{
	ParseArgs "$@"
	LoadConfig

	[[ $PROCESSORS =~ ^[0-9]+$ ]] && (( PROCESSORS >= 2 )) \
		|| Die "PROCESSORS must be an integer of at least 2"
	[[ $OMP_NUM_THREADS =~ ^[0-9]+$ ]] && (( OMP_NUM_THREADS >= 1 )) \
		|| Die "OMP_NUM_THREADS must be a positive integer"
	case $SOLVER in
		auto|navier_stokes|transport) ;;
		*) Die "SOLVER must be auto, navier_stokes, or transport; got: $SOLVER" ;;
	esac

	local input_root output_root case_dir
	input_root=$(ResolvePath "$INPUT_DIR")
	output_root=$(ResolvePath "$OUTPUT_DIR")
	[[ -d $input_root ]] || Die "INPUT_DIR does not exist: $input_root"
	if ! IsTruthy "$DRY_RUN"; then
		mkdir -p "$output_root"
	fi

	SelectCases "$input_root"
	printf 'configuration: %s\ninput root:   %s\noutput root:  %s\nprocessors:   %s\ncases:        %s\n' \
		"$config_file" "$input_root" "$output_root" "$PROCESSORS" "${#case_dirs[@]}"

	BuildSolvers
	for case_dir in "${case_dirs[@]}"; do
		PrepareCase "$case_dir" "$output_root"
	done

	printf '\ncompleted %s case(s) under %s\n' "${#case_dirs[@]}" "$output_root"
}

Main "$@"
