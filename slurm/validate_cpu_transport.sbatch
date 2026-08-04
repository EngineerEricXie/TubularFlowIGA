#!/bin/bash
#SBATCH -J iga-cuda-cpu-ref
#SBATCH -A mch260002p
#SBATCH -p RM-shared
#SBATCH -t 00:10:00
#SBATCH -N 1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH -o slurm/validate-cpu-transport-%j.out

set -euo pipefail
repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
case_root=${IGA_CASE_ROOT:-/ocean/projects/mch260002p/thsieh1/NeuronTransportIGA}
cd "${case_root}"
module load openmpi/4.0.5-gcc10.2.0

output=".nsvms_diagnostics/cuda-validation-42861850/transport-step1-cpu.txt"
/usr/bin/time -v mpiexec -np 8 "${repo_root}/../NeuronTransportIGA-CPU/iga_transport" \
  .nsvms_diagnostics/cylinder-8.ntiga example/cylinder 1 "${output}" \
  example/cylinder/velocityfield.txt

awk 'NR==FNR { a[NR]=$2; b[NR]=$3; next }
  { d0=$2-a[FNR]; d1=$3-b[FNR]; s+=d0*d0+d1*d1; r+=a[FNR]*a[FNR]+b[FNR]*b[FNR];
    if(d0<0)d0=-d0; if(d1<0)d1=-d1; if(d0>m)m=d0; if(d1>m)m=d1 }
  END { printf "transport max_abs=%.17g rel_l2=%.17g\\n",m,sqrt(s/r) }' \
  "${output}" .nsvms_diagnostics/cuda-validation-42861850/transport-step1.txt
