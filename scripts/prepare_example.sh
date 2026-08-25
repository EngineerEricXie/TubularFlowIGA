#!/usr/bin/env bash

set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ranks=${RANKS:-2}

if [[ $# -lt 1 || $# -gt 2 ]]; then
	printf 'usage: %s EXAMPLE [EMPTY_WORK_DIRECTORY]\n' "$0" >&2
	exit 2
fi
if ! [[ $ranks =~ ^[0-9]+$ ]] || (( ranks < 2 )); then
	printf 'RANKS must be an integer of at least 2 because mpmetis rejects one partition\n' >&2
	exit 2
fi

example=$1
if [[ ! $example =~ ^(neuron_transport|vascular_flow)/[A-Za-z0-9][A-Za-z0-9_-]*$ ]]; then
	printf 'example must be APPLICATION/CASE under examples/neuron_transport or examples/vascular_flow: %s\n' "$example" >&2
	exit 2
fi
application=${example%%/*}
case_name=${example##*/}
source_dir=$repo_dir/examples/$example
if [[ ! -d $source_dir ]]; then
	printf 'unknown example: %s\navailable examples:\n' "$example" >&2
	find "$repo_dir/examples/neuron_transport" "$repo_dir/examples/vascular_flow" \
		-mindepth 1 -maxdepth 1 -type d -printf '  %h/%f\n' \
		| sed "s|$repo_dir/examples/|  |" | sort >&2
	exit 2
fi

required=(mesh_parameter.txt simulation_config.json)
for file in "${required[@]}"; do
	if [[ ! -f $source_dir/$file ]]; then
		printf 'example is missing required input: %s/%s\n' "$source_dir" "$file" >&2
		exit 2
	fi
done
skeleton_source=
for candidate in skeleton_initial.swc skeleton_initial.obj; do
	if [[ -f $source_dir/$candidate ]]; then
		if [[ -n $skeleton_source ]]; then
			printf 'example must contain only one skeleton_initial.swc or skeleton_initial.obj: %s\n' "$source_dir" >&2
			exit 2
		fi
		skeleton_source=$candidate
	fi
done
if [[ -z $skeleton_source ]]; then
	printf 'example is missing skeleton_initial.swc or skeleton_initial.obj: %s\n' "$source_dir" >&2
	exit 2
fi

if [[ $# -eq 2 ]]; then
	work_dir=$2
	mkdir -p "$work_dir"
	if [[ -n $(find "$work_dir" -mindepth 1 -maxdepth 1 -print -quit) ]]; then
		printf 'work directory must be empty: %s\n' "$work_dir" >&2
		exit 2
	fi
else
	work_dir=$(mktemp -d "${TMPDIR:-/tmp}/tubularflowiga-${case_name}.XXXXXX")
fi

"$repo_dir/scripts/check_dependencies.sh" preprocessing
for file in "${required[@]}"; do cp "$source_dir/$file" "$work_dir/$file"; done
cp "$source_dir/$skeleton_source" "$work_dir/$skeleton_source"

make -C "$repo_dir" mesh spline cpu
"$repo_dir/preprocessing/mesh/tubular_mesh" pipeline \
	"$work_dir" "$repo_dir/meshgeneration/template"
for generated in skeleton_normalized.swc skeleton.vtp skeleton_smooth.swc controlmesh.vtk; do
	if [[ ! -s $work_dir/$generated ]]; then
		printf 'preprocessing did not create required output: %s/%s\n' "$work_dir" "$generated" >&2
		exit 1
	fi
done
OMP_NUM_THREADS=${OMP_NUM_THREADS:-2} \
	"$repo_dir/preprocessing/spline/spline" "$work_dir/" --no-legacy-text
mpmetis "$work_dir/bzmeshinfo.txt" "$ranks"
database=$work_dir/$case_name-$ranks.ntiga
"$repo_dir/solvers/cpu/iga_pack" "$work_dir" "$ranks" "$database"
"$repo_dir/solvers/cpu/iga_inspect" "$database"
"$repo_dir/solvers/cpu/iga_config_check" "$work_dir/simulation_config.json"
"$repo_dir/solvers/cpu/iga_case_check" "$database" "$work_dir"

printf '%s preparation passed\ncase directory: %s\ndatabase: %s\nranks: %s\n' \
	"$example" "$work_dir" "$database" "$ranks"
if [[ $application == neuron_transport ]]; then
	printf 'CPU neuron transport:\n  mpiexec -np %s %s/solvers/cpu/iga_solve %s %s --system neuron_transport --output %s/neuron-cpu.txt\n' \
		"$ranks" "$repo_dir" "$database" "$work_dir" "$work_dir"
	printf 'CUDA neuron transport:\n  %s/solvers/cuda/iga_cuda solve %s %s --system neuron_transport --output %s/neuron-cuda.txt\n' \
		"$repo_dir" "$database" "$work_dir" "$work_dir"
else
	printf 'CPU vascular flow:\n  mpiexec -np %s %s/solvers/cpu/iga_navier_stokes %s %s --output %s/velocity-cpu.txt\n' \
		"$ranks" "$repo_dir" "$database" "$work_dir" "$work_dir"
	printf 'CUDA vascular flow:\n  %s/solvers/cuda/iga_cuda navier-stokes %s %s --output %s/velocity-cuda.txt\n' \
		"$repo_dir" "$database" "$work_dir" "$work_dir"
fi
