#!/bin/bash
# Supply the allocation at submission: sbatch -A YOUR_PSC_PROJECT this-file.sbatch
#SBATCH -J iga-scaling
#SBATCH -p RM
#SBATCH -N 2
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH -t 04:00:00
#SBATCH -o slurm/scaling-%j.out

set -euo pipefail
repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)
source_root=${IGA_SCALING_CASE_ROOT:?Set IGA_SCALING_CASE_ROOT on the shared filesystem}
output_root=${IGA_SCALING_OUTPUT:?Set IGA_SCALING_OUTPUT to a new shared output directory}
strong_case=${IGA_STRONG_CASE:?Set IGA_STRONG_CASE relative to IGA_SCALING_CASE_ROOT}
strong_database=${IGA_STRONG_DATABASE:?Set IGA_STRONG_DATABASE, using {ranks} where needed}
weak_case=${IGA_WEAK_CASE:?Set IGA_WEAK_CASE relative to IGA_SCALING_CASE_ROOT and include {ranks}}
weak_database=${IGA_WEAK_DATABASE:?Set IGA_WEAK_DATABASE, using {ranks} where needed}
: "${PETSC_DIR:?Set PETSC_DIR to the PETSc installation used for this build}"
module load openmpi/4.0.5-gcc10.2.0

if [[ ${SLURM_JOB_NUM_NODES:-0} -lt 2 ]]; then
	echo "cross_node_scaling.sbatch requires at least two nodes" >&2
	exit 2
fi
for program in iga_navier_stokes iga_mesh_check iga_flow_validate; do
	if [[ ! -x ${repo_root}/solvers/cpu/${program} ]]; then
		echo "missing solvers/cpu/${program}; build before submission" >&2
		exit 2
	fi
done

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_THREAD_LIMIT=${SLURM_CPUS_PER_TASK}
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export IGA_STAGE_SOURCE=$(realpath "${source_root}")
srun --nodes="${SLURM_JOB_NUM_NODES}" --ntasks="${SLURM_JOB_NUM_NODES}" \
	--ntasks-per-node=1 --cpus-per-task=1 \
	bash -c 'mkdir -p "$SLURM_TMPDIR/tubularflowiga-scaling"; rsync -a --delete "$IGA_STAGE_SOURCE/" "$SLURM_TMPDIR/tubularflowiga-scaling/"'
staged_root=${SLURM_TMPDIR}/tubularflowiga-scaling

read -r -a ranks <<<"${IGA_SCALING_RANKS:-1 64 128 256}"
tasks_per_node=${SLURM_NTASKS_PER_NODE}
launcher="mpiexec --bind-to core --map-by ppr:${tasks_per_node}:node"
python3 "${repo_root}/scripts/hpc_build_manifest.py" --component cpu \
	--dependency-status passed --output "${output_root}.build.json" \
	--binary "${repo_root}/solvers/cpu/iga_navier_stokes"
python3 "${repo_root}/scripts/hpc_cross_node_scaling.py" \
	--strong-case "${staged_root}/${strong_case}" --strong-database "${strong_database}" \
	--weak-case "${staged_root}/${weak_case}" --weak-database "${weak_database}" \
	--output-dir "${output_root}" --ranks "${ranks[@]}" \
	--repetitions "${IGA_SCALING_REPETITIONS:-3}" \
	--timeout "${IGA_SCALING_TIMEOUT:-3600}" --launcher "${launcher}"
