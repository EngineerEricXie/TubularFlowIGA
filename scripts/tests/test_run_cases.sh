#!/usr/bin/env bash

set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
work_dir=$(mktemp -d "${TMPDIR:-/tmp}/tubularflowiga-run-cases.XXXXXX")
trap 'rm -rf -- "$work_dir"' EXIT

case_root=$work_dir/case\ root
mkdir -p "$case_root/simple" "$case_root/transport" "$case_root/one_d" "$case_root/vca_3d"
cp "$repo_dir/examples/validation/simple_vascular/simulation_config.json" \
	"$case_root/simple/simulation_config.json"
cp "$repo_dir/examples/neuron_transport/straight_neurite/simulation_config.json" \
	"$case_root/transport/simulation_config.json"
cp "$repo_dir/examples/one_d/vca_pfc_closed_loop/simulation_config.json" \
	"$case_root/one_d/simulation_config.json"
cp "$repo_dir/examples/vascular_flow/vca_bifurcation/simulation_config.json" \
	"$case_root/vca_3d/simulation_config.json"

WriteConfig()
{
	local path=$1 output_root=$2 clean=$3 backend=${4:-cpu} ranks=${5:-10}
	cat > "$path" <<CONFIG
CASE_ROOT=$case_root
OUTPUT_ROOT=$output_root
RANKS=$ranks
BACKEND=$backend
BUILD_SOLVERS=0
CLEAN=$clean
RUN_MESH_CHECK=0
RUN_SOLVER=1
RUN_VALIDATION=1
SOLVER=auto
SYSTEM=
MPIEXEC=mpiexec
OMP_NUM_THREADS=2
SOLVER_ARGS=
DRY_RUN=1
CONFIG
}

config=$work_dir/execution.conf
WriteConfig "$config" "" 0
flow_output=$("$repo_dir/scripts/run_cases.sh" --config "$config" simple)
[[ $flow_output == *"iga_navier_stokes"* ]]
[[ $flow_output == *"--system blood_flow"* ]]
[[ $flow_output == *"-np 10"* ]]
[[ $flow_output == *"simple/generated"* ]]
[[ ! -e "$case_root/simple/generated" ]]

transport_output=$("$repo_dir/scripts/run_cases.sh" --config "$config" transport)
[[ $transport_output == *"iga_solve"* ]]
[[ $transport_output == *"--system neuron_transport"* ]]
[[ ! -e "$case_root/transport/generated" ]]

WriteConfig "$config" "" 0 cuda
cuda_output=$("$repo_dir/scripts/run_cases.sh" --config "$config" simple)
[[ $cuda_output == *"iga_cuda navier-stokes"* ]]
[[ $cuda_output != *"iga_mesh_check"* ]]

WriteConfig "$config" "" 0 cpu 1
one_d_output=$("$repo_dir/scripts/run_cases.sh" --config "$config" one_d)
[[ $one_d_output == *"iga_1d"* ]]
[[ $one_d_output == *"--system vascular_flow"* ]]
[[ $one_d_output == *"results/vascular_flow"* ]]
[[ $one_d_output != *"generate_case.sh"* ]]

WriteConfig "$config" "" 0 cuda 1
if "$repo_dir/scripts/run_cases.sh" --config "$config" one_d > /dev/null 2>&1; then
	printf 'CUDA was accepted for a 1D case\n' >&2
	exit 1
fi

WriteConfig "$config" "" 0 cuda 2
if "$repo_dir/scripts/run_cases.sh" --config "$config" vca_3d > /dev/null 2>&1; then
	printf 'CUDA was accepted for a 3D VCA closed-loop case\n' >&2
	exit 1
fi

WriteConfig "$config" "" 0 cpu 2
vca_output=$("$repo_dir/scripts/run_cases.sh" --config "$config" vca_3d)
[[ $vca_output == *"coupling=vca_closed_loop"* ]]
[[ $vca_output == *"iga_navier_stokes"*"--system blood_flow"* ]]
[[ $vca_output != *"iga_solve"* ]]

WriteConfig "$config" "$repo_dir" 0
if "$repo_dir/scripts/run_cases.sh" --config "$config" simple > /dev/null 2>&1; then
	printf 'repository root was accepted as OUTPUT_ROOT\n' >&2
	exit 1
fi

output_root=$work_dir/outputs
mkdir -p "$output_root/simple"
printf 'unrelated data\n' > "$output_root/simple/keep.txt"
WriteConfig "$config" "$output_root" 1
if "$repo_dir/scripts/run_cases.sh" --config "$config" simple > /dev/null 2>&1; then
	printf 'CLEAN accepted an output directory without a generated-case manifest\n' >&2
	exit 1
fi
[[ -f $output_root/simple/keep.txt ]]

mkdir -p "$output_root/one_d"
cat > "$output_root/one_d/run_manifest.json" <<'JSON'
{
  "case_name": "one_d",
  "completed": true
}
JSON
WriteConfig "$config" "$output_root" 1 cpu 1
clean_one_d_output=$("$repo_dir/scripts/run_cases.sh" --config "$config" one_d)
[[ $clean_one_d_output == *"iga_1d"* ]]
[[ -f $output_root/one_d/run_manifest.json ]]

printf 'run_cases tests passed\n'
