#!/bin/bash
#SBATCH --job-name=liver_refined2
#SBATCH --partition=RM
#SBATCH --nodes=2
#SBATCH --ntasks=128
#SBATCH --time=05:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

source ~/.bashrc
tfi_env
/usr/bin/time -f "Elapsed time: %e seconds" ./scripts/run_cases.sh liver_refined2 2>&1

