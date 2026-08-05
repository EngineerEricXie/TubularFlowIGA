#!/bin/bash
#SBATCH -J iga-cuda-nmo-full
#SBATCH -A mch260002p
#SBATCH -p GPU-shared
#SBATCH --gres=gpu:v100-32:1
#SBATCH -t 00:15:00
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH -o slurm/nmo-full-v100-%j.out

set -euo pipefail
repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
case_root=${IGA_CASE_ROOT:-/ocean/projects/mch260002p/thsieh1/NeuronTransportIGA}
cd "${case_root}"
module load cuda/12.4.0

database=".nsvms_diagnostics/NMO_54499_new_v2/case-16.ntiga"
case_dir=".nsvms_diagnostics/NMO_54499_new_v2"
output_dir=".nsvms_diagnostics/cuda-nmo-full-${SLURM_JOB_ID}"
mkdir -p "${output_dir}"

/usr/bin/time -v "${repo_root}/iga_cuda" navier-stokes "${database}" "${case_dir}" 8 \
  "${output_dir}/velocity.txt"
/usr/bin/time -v "${repo_root}/iga_cuda" transport "${database}" "${case_dir}" 300 \
  "${output_dir}/transport.txt" "${output_dir}/velocity.txt"

awk 'NR==FNR{for(i=1;i<=3;i++)a[NR,i]=$i;next}{for(i=1;i<=3;i++){d=$i-a[FNR,i];s+=d*d;r+=a[FNR,i]*a[FNR,i];if(d<0)d=-d;if(d>m)m=d}}END{e=sqrt(s/r);printf "velocity max_abs=%.17g rel_l2=%.17g\n",m,e;exit(e>1e-4)}' \
  "${case_dir}/velocity-v2.txt" "${output_dir}/velocity.txt"
awk 'NR==FNR{a[NR]=$1;next}{d=$1-a[FNR];s+=d*d;r+=a[FNR]*a[FNR];if(d<0)d=-d;if(d>m)m=d}END{e=sqrt(s/r);printf "pressure max_abs=%.17g rel_l2=%.17g\n",m,e;exit(e>1e-4)}' \
  "${case_dir}/velocity-v2.txt.pressure" "${output_dir}/velocity.txt.pressure"
awk 'NR==FNR{a[NR]=$2;b[NR]=$3;next}{d0=$2-a[FNR];d1=$3-b[FNR];s+=d0*d0+d1*d1;r+=a[FNR]*a[FNR]+b[FNR]*b[FNR];if(d0<0)d0=-d0;if(d1<0)d1=-d1;if(d0>m)m=d0;if(d1>m)m=d1}END{e=sqrt(s/r);printf "transport max_abs=%.17g rel_l2=%.17g\n",m,e;exit(e>1e-4)}' \
  "${case_dir}/transport-v2.txt" "${output_dir}/transport.txt"

echo "validation_output=${output_dir}"
