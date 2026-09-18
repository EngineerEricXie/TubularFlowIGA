#!/bin/bash
# Submit from the repository root:
# sbatch -A mch260002p solvers/cpu/slurm/run_liver_refined2_14.sbatch
#SBATCH --job-name=liver_refined2_14
#SBATCH --partition=RM
#SBATCH --qos=rm
#SBATCH --nodes=2
#SBATCH --ntasks=24
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
# Use the RM partition default memory for each full node.
#SBATCH --time=05:00:00
#SBATCH --output=liver_refined2_14-%j.out
#SBATCH --error=liver_refined2_14-%j.err

set -euo pipefail
module load openmpi/4.0.5-gcc10.2.0
repo_root=/ocean/projects/mch260002p/agrawaly/Codes/TubularFlowIGA
case_root="$repo_root/Output/liver_refined2_14"
database="$case_root/database/liver_refined2-24.ntiga"
solver="$repo_root/solvers/cpu/iga_navier_stokes"
[[ -r "$database" ]] || { echo "Cannot read database: $database" >&2; exit 1; }
[[ -x "$solver" ]] || { echo "Missing solver executable: $solver" >&2; exit 1; }
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export PETSC_OPTIONS="-ksp_type fgmres -ksp_gmres_restart 200 -pc_type asm -pc_asm_overlap 2 -sub_pc_type ilu -sub_pc_factor_levels 1 -ksp_converged_reason -log_view"
output_dir="$case_root/results/blood_flow/retry-$SLURM_JOB_ID"
mkdir -p "$output_dir"
cd "$repo_root"
/usr/bin/time -v mpiexec -np 24 --map-by ppr:12:node --bind-to core     "$solver" "$database" "$case_root/preprocessing"     --system blood_flow --output "$output_dir/navier_stokes-cpu.txt"     --nonlinear-rtol 0.01
