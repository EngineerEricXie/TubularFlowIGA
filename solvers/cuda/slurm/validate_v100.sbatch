#!/bin/bash
#SBATCH -J iga-cuda-test
#SBATCH -A mch260002p
#SBATCH -p GPU-shared
#SBATCH --gres=gpu:v100-32:1
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH -o slurm/validate-v100-%j.out

set -euo pipefail
repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
case_root=${IGA_CASE_ROOT:-/ocean/projects/mch260002p/thsieh1/NeuronTransportIGA}
cd "${case_root}"
module load cuda/12.4.0

test_root=".nsvms_diagnostics/cuda-validation-${SLURM_JOB_ID}"
mkdir -p "${test_root}"

"${repo_root}/iga_cuda" device-info
"${repo_root}/iga_cuda" mesh-check .nsvms_diagnostics/cylinder-8.ntiga
"${repo_root}/iga_cuda" transport \
  .nsvms_diagnostics/cylinder-8.ntiga example/cylinder 1 \
  "${test_root}/transport-step1.txt" example/cylinder/velocityfield.txt
"${repo_root}/iga_cuda" navier-stokes \
  .nsvms_diagnostics/cylinder-8.ntiga example/cylinder 8 \
  "${test_root}/velocity.txt"

test -s "${test_root}/transport-step1.txt.vtk"
test -s "${test_root}/velocity.txt.vtk"
grep -q '^SCALARS N0 double 1$' "${test_root}/transport-step1.txt.vtk"
grep -q '^SCALARS Nplus double 1$' "${test_root}/transport-step1.txt.vtk"
grep -q '^VECTORS Velocity double$' "${test_root}/velocity.txt.vtk"
grep -q '^SCALARS Pressure double 1$' "${test_root}/velocity.txt.vtk"


cpu_transport=".nsvms_diagnostics/cuda-validation-42861850/transport-step1-cpu.txt"
cpu_velocity=".nsvms_diagnostics/cylinder-ns-bjacobi-velocity.txt"
if [[ -f "${cpu_transport}" ]]; then
  awk 'NR==FNR{a[NR]=$2;b[NR]=$3;next}{d0=$2-a[FNR];d1=$3-b[FNR];s+=d0*d0+d1*d1;r+=a[FNR]*a[FNR]+b[FNR]*b[FNR]}END{e=sqrt(s/r);printf "transport_rel_l2=%.17g\n",e;exit(e>1e-5)}' "${cpu_transport}" "${test_root}/transport-step1.txt"
fi
if [[ -f "${cpu_velocity}" ]]; then
  awk 'NR==FNR{for(i=1;i<=3;i++)a[NR,i]=$i;next}{for(i=1;i<=3;i++){d=$i-a[FNR,i];s+=d*d;r+=a[FNR,i]*a[FNR,i]}}END{e=sqrt(s/r);printf "velocity_rel_l2=%.17g\n",e;exit(e>1e-7)}' "${cpu_velocity}" "${test_root}/velocity.txt"
  awk 'NR==FNR{a[NR]=$1;next}{d=$1-a[FNR];s+=d*d;r+=a[FNR]*a[FNR]}END{e=sqrt(s/r);printf "pressure_rel_l2=%.17g\n",e;exit(e>1e-7)}' "${cpu_velocity}.pressure" "${test_root}/velocity.txt.pressure"
fi

echo "validation_output=${test_root}"
