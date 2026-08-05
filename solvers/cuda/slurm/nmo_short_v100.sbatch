#!/bin/bash
#SBATCH -J iga-cuda-nmo-short
#SBATCH -A mch260002p
#SBATCH -p GPU-shared
#SBATCH --gres=gpu:v100-32:1
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH -o slurm/nmo-short-v100-%j.out

set -euo pipefail
repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
case_root=${IGA_CASE_ROOT:?Set IGA_CASE_ROOT to the external case-data directory}
cd "${case_root}"
module load cuda/12.4.0

database=".nsvms_diagnostics/NMO_54499_new_v2/case-16.ntiga"
case_dir=".nsvms_diagnostics/NMO_54499_new_v2"
output_dir=".nsvms_diagnostics/cuda-nmo-short-${SLURM_JOB_ID}"
mkdir -p "${output_dir}"

/usr/bin/time -v "${repo_root}/iga_cuda" mesh-check "${database}"
/usr/bin/time -v "${repo_root}/iga_cuda" transport "${database}" "${case_dir}" 1 \
  "${output_dir}/transport-step1.txt" "${case_dir}/velocity-v2.txt"

set +e
/usr/bin/time -v "${repo_root}/iga_cuda" navier-stokes "${database}" "${case_dir}" 1 \
  "${output_dir}/velocity-one-newton.txt"
ns_status=$?
set -e
if [[ ${ns_status} -ne 1 ]]; then
  echo "unexpected_navier_stokes_status=${ns_status}"
  exit 1
fi
echo "expected one-Newton nonconvergence observed"
echo "validation_output=${output_dir}"
