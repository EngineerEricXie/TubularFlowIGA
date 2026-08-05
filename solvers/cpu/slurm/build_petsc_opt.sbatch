#!/bin/bash
#SBATCH -A mch260002p
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem-per-cpu=2000M
#SBATCH -t 01:00:00
#SBATCH -J petsc_opt

set -euo pipefail
module load openmpi/4.0.5-gcc10.2.0

PETSC_DIR=${PETSC_DIR:-/ocean/projects/mch260002p/thsieh1/petsc}
PETSC_ARCH=${PETSC_ARCH:-arch-linux-c-opt}
cd "$PETSC_DIR"

./configure PETSC_ARCH="$PETSC_ARCH" \
	--with-debugging=0 \
	COPTFLAGS=-O3 CXXOPTFLAGS=-O3 \
	--with-cc=mpicc --with-cxx=mpicxx --with-fc=0 \
	--with-x=0 --with-shared-libraries=1
make PETSC_ARCH="$PETSC_ARCH" -j"${SLURM_CPUS_PER_TASK}" all
