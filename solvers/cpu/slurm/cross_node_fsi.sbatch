#!/bin/bash
# Supply the allocation at submission: sbatch -A YOUR_PSC_PROJECT this-file.sbatch
#SBATCH -J iga-fsi-xnode
#SBATCH -p RM
#SBATCH -N 2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH -t 03:00:00
#SBATCH -o slurm/fsi-xnode-%j.out

set -euo pipefail
repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)
output=${IGA_FSI_OUTPUT:?Set IGA_FSI_OUTPUT to a new shared output directory}
: "${PETSC_DIR:?Set PETSC_DIR to the PETSc installation used for this build}"
module load openmpi/4.0.5-gcc10.2.0

if [[ ${SLURM_JOB_NUM_NODES:-0} -lt 2 || ${SLURM_NTASKS:-0} -lt 4 ]]; then
	echo "cross_node_fsi.sbatch requires two nodes and four allocated tasks" >&2
	exit 2
fi
binary=${repo_root}/solvers/cpu/immersed_moving_distributed_fsi_runtime_test
scheduled_binary=${repo_root}/solvers/cpu/distributed_immersed_extension_test
if [[ ! -x ${binary} || ! -x ${scheduled_binary} ]]; then
	echo "build immersed_moving_distributed_fsi_runtime_test and distributed_immersed_extension_test before submission" >&2
	exit 2
fi
mkdir "${output}"
checkpoint=${output}/checkpoint-4

export OMP_NUM_THREADS=1
export OMP_THREAD_LIMIT=1
export OMP_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
python3 "${repo_root}/scripts/hpc_build_manifest.py" --component cpu \
	--dependency-status passed --output "${output}/build.json" --binary "${binary}"
python3 "${repo_root}/scripts/hpc_test_tiers.py" --tier scheduled --timeout 1800 \
	--output-dir "${output}/scheduled"

run_fsi()
{
	local ranks=$1
	local destination=$2
	local ranks_per_node=1
	if (( ranks >= 4 )); then ranks_per_node=2; fi
	shift 2
	mkdir "${destination}"
	mpiexec --map-by "ppr:${ranks_per_node}:node" --bind-to core --report-bindings -np "${ranks}" \
		python3 "${repo_root}/scripts/hpc_rank_run.py" \
		--output-dir "${destination}" --expected-ranks "${ranks}" --timeout 7200 -- \
		"${binary}" strong -immersed_moving_mat_mumps_icntl_14 100 "$@"
}

run_fsi 1 "${output}/reference-1"
run_fsi 4 "${output}/writer-4" -moving_fsi_checkpoint_root "${checkpoint}"
python3 "${repo_root}/scripts/hpc_check_strong_fsi.py" \
	--reference-run "${output}/reference-1" --run "${output}/reference-1" \
	--run "${output}/writer-4" >"${output}/strong-fsi.json"

python3 "${repo_root}/scripts/hpc_hash_tree.py" --root "${checkpoint}" \
	--output "${output}/checkpoint-input-hashes.json"
run_fsi 2 "${output}/reader-2" -moving_fsi_checkpoint_root "${checkpoint}" \
	-moving_fsi_checkpoint_read_only true -moving_fsi_checkpoint_source_ranks 4
python3 "${repo_root}/scripts/hpc_check_paired_fsi_restart.py" \
	--run "${output}/reader-2" --source-ranks 4 \
	--input-hashes "${output}/checkpoint-input-hashes.json" \
	>"${output}/paired-restart.json"

python3 "${repo_root}/scripts/hpc_scheduler_record.py" \
	--output "${output}/scheduler.json" --status passed --returncode 0 -- \
	"${binary}" strong
