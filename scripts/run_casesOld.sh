#!/usr/bin/env bash

set -euo pipefail

SECONDS=0

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
default_petsc_dir=/Users/yamnesh/softwares/petsc/opt
ranks=${RANKS:-}
backend=${BACKEND:-cpu}
build=${BUILD:-1}
run_root=${RUN_ROOT:-}

Usage()
{
	cat >&2 <<EOF
usage: $0 [-o DIR|--output-dir DIR] [-r N|--ranks N] [CASE ...]

Run prepared TubularFlowIGA examples.

Default cases:
  vascular_flow/straight_tube
  one_d/rigid_straight

Supported CASE forms:
  vascular_flow/CASE
  neuron_transport/CASE
  one_d/CASE

Environment:
  BACKEND=cpu|cuda|prepare-only   default: cpu
  RANKS=N                         default: 1 on macOS, 2 elsewhere
  RUN_ROOT=DIR                    default: fresh temporary directory
  BUILD=0|1                       default: 1, build PETSc/CUDA solver targets
  PETSC_DIR=DIR                   default: /Users/yamnesh/softwares/petsc/opt
  PETSC_ARCH=ARCH                 optional PETSc arch
  EIGEN_DIR=DIR                   directory containing Eigen/
  MPIEXEC=COMMAND                 default: mpiexec
  CPU_FLOW_SOLVER_ARGS="ARGS"     extra args for CPU iga_navier_stokes
  CUDA_ARCHS="ARCH ..."           passed to make cuda, e.g. "70 80 89"

Examples:
  $0
  $0 -o results -r 4 vascular_flow/straight_tube
  $0 --output-dir /tmp/tfi-results --ranks 4 vascular_flow/straight_tube
  $0 vascular_flow/straight_tube neuron_transport/straight_neurite
  BACKEND=prepare-only $0 vascular_flow/y_bifurcation
  RANKS=4 RUN_ROOT=/tmp/tfi-run $0 vascular_flow/bent_tube
EOF
}

Die()
{
	printf '%s\n' "$*" >&2
	exit 2
}

IsTruthy()
{
	case ${1:-} in
		1|yes|true|on) return 0 ;;
		*) return 1 ;;
	esac
}

ConfigureHostDefaults()
{
	if [[ -z $ranks ]]; then
		if [[ $(uname -s) == Darwin ]]; then
			ranks=1
		else
			ranks=2
		fi
	fi

	if [[ -z ${PETSC_DIR:-} && -d $default_petsc_dir ]]; then
		export PETSC_DIR=$default_petsc_dir
	fi
	if [[ -n ${PETSC_DIR:-} ]]; then
		local petsc_prefix=$PETSC_DIR
		if [[ -n ${PETSC_ARCH:-} ]]; then
			petsc_prefix=$petsc_prefix/$PETSC_ARCH
		fi
		if [[ ! -f $petsc_prefix/lib/petsc/conf/petscvariables &&
			-f $default_petsc_dir/lib/petsc/conf/petscvariables ]]; then
			export PETSC_DIR=$default_petsc_dir
			unset PETSC_ARCH
		fi
	fi

	if [[ -z ${EIGEN_DIR:-} ]]; then
		if [[ -d /opt/homebrew/include/eigen3/Eigen ]]; then
			export EIGEN_DIR=/opt/homebrew/include/eigen3
		elif [[ -d /usr/local/include/eigen3/Eigen ]]; then
			export EIGEN_DIR=/usr/local/include/eigen3
		fi
	fi

	if [[ $(uname -s) == Darwin && -z ${LDLIBS+x} ]]; then
		export LDLIBS=
	fi
}

ParseArgs()
{
	cases=()
	while [[ $# -gt 0 ]]; do
		case $1 in
			--help|-h)
				Usage
				exit 0
				;;
			--output-dir|--run-root|-o)
				if [[ $# -lt 2 || -z $2 ]]; then
					Die "$1 requires a directory argument"
				fi
				run_root=$2
				shift 2
				;;
			--output-dir=*|--run-root=*)
				run_root=${1#*=}
				shift
				;;
			-o*)
				run_root=${1#-o}
				if [[ -z $run_root ]]; then
					Die "-o requires a directory argument"
				fi
				shift
				;;
			--ranks|--processors|--cpus|-r)
				if [[ $# -lt 2 || -z $2 ]]; then
					Die "$1 requires an integer argument"
				fi
				ranks=$2
				shift 2
				;;
			--ranks=*|--processors=*|--cpus=*)
				ranks=${1#*=}
				shift
				;;
			-r*)
				ranks=${1#-r}
				if [[ -z $ranks ]]; then
					Die "-r requires an integer argument"
				fi
				shift
				;;
			--*)
				Die "unknown option: $1"
				;;
			*)
				cases+=("$1")
				shift
				;;
		esac
	done
}

EnsureRunRoot()
{
	if [[ -z $run_root ]]; then
		run_root=$(mktemp -d "${TMPDIR:-/tmp}/tubularflowiga-run.XXXXXX")
	else
		mkdir -p "$run_root"
	fi
	printf 'run root: %s\n' "$run_root"
}

BuildSolversIfNeeded()
{
	local need_3d_cpu=$1
	local need_1d=$2
	local need_cuda=$3
	local make_args=()

	if [[ -n ${PETSC_DIR:-} ]]; then
		make_args+=(PETSC_DIR="$PETSC_DIR")
	fi
	if [[ -n ${PETSC_ARCH:-} ]]; then
		make_args+=(PETSC_ARCH="$PETSC_ARCH")
	fi

	if ! IsTruthy "$build"; then
		return
	fi

	if [[ $need_cuda == 1 ]]; then
		make -C "$repo_dir" cuda
	fi
	if [[ $need_3d_cpu == 1 ]]; then
		if [[ -z ${PETSC_DIR:-} ]]; then
			Die "PETSC_DIR is required to build 3D CPU solvers; expected default at $default_petsc_dir or set PETSC_DIR."
		fi
		make -C "$repo_dir/solvers/cpu" petsc "${make_args[@]}"
	fi
	if [[ $need_1d == 1 ]]; then
		make -C "$repo_dir/solvers/one_d" petsc "${make_args[@]}"
	fi
}

RequireExecutable()
{
	local path=$1
	local hint=$2
	if [[ ! -x $path ]]; then
		Die "missing executable: $path
$hint"
	fi
}

Prepare3DCase()
{
	local example=$1
	local application=${example%%/*}
	local case_name=${example##*/}
	local work_dir=$run_root/${example//\//_}
	local prepare_ranks=$ranks

	if (( ranks == 1 )); then
		prepare_ranks=2
	fi

	mkdir -p "$work_dir"
	if [[ -n $(find "$work_dir" -mindepth 1 -maxdepth 1 -print -quit) ]]; then
		Die "work directory must be empty: $work_dir"
	fi

	printf '\n== Preparing %s ==\n' "$example"
	RANKS="$prepare_ranks" "$repo_dir/scripts/prepare_example.sh" "$example" "$work_dir"

	if (( ranks == 1 )); then
		awk '{print 0}' "$work_dir/bzmeshinfo.txt.epart.$prepare_ranks" \
			> "$work_dir/bzmeshinfo.txt.epart.1"
		"$repo_dir/solvers/cpu/iga_pack" "$work_dir" 1 "$work_dir/$case_name-1.ntiga"
		"$repo_dir/solvers/cpu/iga_inspect" "$work_dir/$case_name-1.ntiga"
		"$repo_dir/solvers/cpu/iga_case_check" "$work_dir/$case_name-1.ntiga" "$work_dir"
	fi

	local database=$work_dir/$case_name-$ranks.ntiga
	if [[ ! -s $database ]]; then
		Die "expected database was not created: $database"
	fi

	if [[ $backend == prepare-only ]]; then
		printf 'prepared %s\n  case directory: %s\n  database: %s\n' "$example" "$work_dir" "$database"
		return
	fi

	case $backend in
		cpu)
			if [[ $application == vascular_flow ]]; then
				RunVascularCpu "$database" "$work_dir"
			else
				RunNeuronCpu "$database" "$work_dir"
			fi
			;;
		cuda)
			if [[ $application == vascular_flow ]]; then
				RunVascularCuda "$database" "$work_dir"
			else
				RunNeuronCuda "$database" "$work_dir"
			fi
			;;
		*) Die "BACKEND must be cpu, cuda, or prepare-only; got: $backend" ;;
	esac
}

RunVascularCpu()
{
	local database=$1
	local work_dir=$2
	local mpiexec_cmd=${MPIEXEC:-mpiexec}
	local output=$work_dir/velocity-cpu.txt
	local solver_cmd=("$repo_dir/solvers/cpu/iga_navier_stokes" "$database" "$work_dir" --output "$output")
	if [[ -z ${CPU_FLOW_SOLVER_ARGS:-} && ${work_dir##*/} == vascular_flow_AGH50Bifurcation ]]; then
		CPU_FLOW_SOLVER_ARGS="--linear-rtol 1e-6 --linear-max-it 12000 --linear-gmres-restart 200"
		printf 'note    using AGH50Bifurcation CPU flow solver args: %s\n' "$CPU_FLOW_SOLVER_ARGS"
	fi
	if [[ -n ${CPU_FLOW_SOLVER_ARGS:-} ]]; then
		local solver_args
		read -r -a solver_args <<< "$CPU_FLOW_SOLVER_ARGS"
		solver_cmd+=("${solver_args[@]}")
	fi

	RequireExecutable "$repo_dir/solvers/cpu/iga_mesh_check" \
		"build it with: make cpu-petsc PETSC_DIR=/path/to/petsc"
	RequireExecutable "$repo_dir/solvers/cpu/iga_navier_stokes" \
		"build it with: make cpu-petsc PETSC_DIR=/path/to/petsc"

	printf '\n== Running CPU vascular flow ==\n'
	if (( ranks == 1 )); then
		"$repo_dir/solvers/cpu/iga_mesh_check" "$database"
		"${solver_cmd[@]}"
	else
		"$mpiexec_cmd" -np "$ranks" "$repo_dir/solvers/cpu/iga_mesh_check" "$database"
		"$mpiexec_cmd" -np "$ranks" "${solver_cmd[@]}"
	fi
	"$repo_dir/solvers/cpu/iga_flow_validate" "$database" "$output"
	printf 'vascular CPU result: %s\n' "$output"
}

RunNeuronCpu()
{
	local database=$1
	local work_dir=$2
	local mpiexec_cmd=${MPIEXEC:-mpiexec}
	local output=$work_dir/neuron-cpu.txt

	RequireExecutable "$repo_dir/solvers/cpu/iga_solve" \
		"build it with: make cpu-petsc PETSC_DIR=/path/to/petsc"

	printf '\n== Running CPU neuron transport ==\n'
	if (( ranks == 1 )); then
		"$repo_dir/solvers/cpu/iga_solve" \
			"$database" "$work_dir" --system neuron_transport --output "$output"
	else
		"$mpiexec_cmd" -np "$ranks" "$repo_dir/solvers/cpu/iga_solve" \
			"$database" "$work_dir" --system neuron_transport --output "$output"
	fi
	printf 'neuron CPU result: %s\n' "$output"
}

RunVascularCuda()
{
	local database=$1
	local work_dir=$2
	local output=$work_dir/velocity-cuda.txt

	RequireExecutable "$repo_dir/solvers/cuda/iga_cuda" \
		"build it with: make cuda CUDA_ARCHS=YOUR_GPU_ARCH"

	printf '\n== Running CUDA vascular flow ==\n'
	"$repo_dir/solvers/cuda/iga_cuda" mesh-check "$database"
	"$repo_dir/solvers/cuda/iga_cuda" navier-stokes \
		"$database" "$work_dir" --output "$output"
	printf 'vascular CUDA result: %s\n' "$output"
}

RunNeuronCuda()
{
	local database=$1
	local work_dir=$2
	local output=$work_dir/neuron-cuda.txt

	RequireExecutable "$repo_dir/solvers/cuda/iga_cuda" \
		"build it with: make cuda CUDA_ARCHS=YOUR_GPU_ARCH"

	printf '\n== Running CUDA neuron transport ==\n'
	"$repo_dir/solvers/cuda/iga_cuda" solve \
		"$database" "$work_dir" --system neuron_transport --output "$output"
	printf 'neuron CUDA result: %s\n' "$output"
}

RunOneDCase()
{
	local example=$1
	local case_name=${example##*/}
	local output_dir=$run_root/${example//\//_}_output

	if [[ $backend != cpu ]]; then
		printf '\n== Skipping %s ==\nNative 1D cases use the PETSc CPU solver; set BACKEND=cpu to run it.\n' "$example"
		return
	fi

	RequireExecutable "$repo_dir/solvers/one_d/iga_1d" \
		"build it with: make one-d-petsc PETSC_DIR=/path/to/petsc"

	printf '\n== Running native 1D case %s ==\n' "$example"
	mkdir -p "$output_dir"
	"$repo_dir/solvers/one_d/iga_1d" "$repo_dir/examples/$example" --check
	"$repo_dir/solvers/one_d/iga_1d" "$repo_dir/examples/$example" \
		--output-dir "$output_dir"
	printf '1D result directory: %s\n' "$output_dir"
}

Main()
{
	local cases=()
	local need_3d=0
	local need_3d_cpu=0
	local need_1d=0
	local need_cuda=0

	ParseArgs "$@"

	if [[ ${#cases[@]} -eq 0 ]]; then
		cases=(vascular_flow/straight_tube one_d/rigid_straight)
	fi

	ConfigureHostDefaults
	EnsureRunRoot

	[[ $ranks =~ ^[0-9]+$ ]] || Die "RANKS must be an integer"

	for example in "${cases[@]}"; do
		if [[ ! -d $repo_dir/examples/$example ]]; then
			Die "unknown case: $example"
		fi
		case $example in
			vascular_flow/*|neuron_transport/*)
				need_3d=1
				if [[ $backend == cpu ]]; then
					need_3d_cpu=1
				elif [[ $backend == cuda ]]; then
					need_cuda=1
				fi
				;;
			one_d/*)
				if [[ $backend == cpu ]]; then
					need_1d=1
				fi
				;;
			*)
				Die "unsupported case name: $example"
				;;
		esac
	done
	if [[ $need_3d == 1 ]] && (( ranks < 1 )); then
		Die "RANKS must be at least 1"
	fi
	case $backend in
		cpu|cuda|prepare-only) ;;
		*) Die "BACKEND must be cpu, cuda, or prepare-only; got: $backend" ;;
	esac

	BuildSolversIfNeeded "$need_3d_cpu" "$need_1d" "$need_cuda"

	for example in "${cases[@]}"; do
		case $example in
			vascular_flow/*|neuron_transport/*) Prepare3DCase "$example" ;;
			one_d/*) RunOneDCase "$example" ;;
		esac
	done

	printf '\ncompleted requested cases under: %s\n' "$run_root"
}

Main "$@"

echo "Total elapsed time: $SECONDS seconds"
