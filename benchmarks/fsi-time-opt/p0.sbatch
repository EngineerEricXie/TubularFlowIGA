#!/bin/bash
#SBATCH -J fsi-opt-p0
#SBATCH -A mch260002p
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=7600M
#SBATCH -t 01:00:00
set -euo pipefail
module load anaconda3 openmpi/4.0.5-gcc10.2.0
export PETSC_DIR=/ocean/projects/mch260002p/thsieh1/TubularFlowIGA-hpc-acceptance/petsc
export PETSC_ARCH=arch-bridges2-acceptance CXX=mpicxx
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OMP_THREAD_LIMIT=$SLURM_CPUS_PER_TASK
export OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export IGA_ASSEMBLY_THREADS=$OMP_NUM_THREADS IGA_ASSEMBLY_BATCH_SIZE=${IGA_REQUESTED_BATCH_SIZE:-8} IGA_PROFILE=1
bash "$1" --root "$2" --variants "$3" --case-list "$4"
