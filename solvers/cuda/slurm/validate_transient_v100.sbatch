#!/bin/bash
#SBATCH -J iga-transient-gpu
# Supply the allocation at submission: sbatch -A YOUR_PSC_PROJECT this-file.sbatch
#SBATCH -p GPU-shared
#SBATCH --gres=gpu:v100-32:1
#SBATCH -t 01:00:00
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH -o slurm/validate-transient-v100-%j.out

set -euo pipefail
cuda_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
case_dir=${IGA_TRANSIENT_CASE_DIR:?Set IGA_TRANSIENT_CASE_DIR}
database=${IGA_TRANSIENT_DATABASE:?Set IGA_TRANSIENT_DATABASE}
transport_system=${IGA_TRANSPORT_SYSTEM:-tracer_transport}
parity_tolerance=${IGA_PARITY_TOLERANCE:-1e-5}
mass_tolerance=${IGA_MASS_BALANCE_TOLERANCE:-1e-2}
divergence_tolerance=${IGA_DIVERGENCE_THEOREM_TOLERANCE:-1e-6}
cycle_tolerance=${IGA_CYCLE_REL_L2_TOLERANCE:-1e-3}
womersley_tolerance=${IGA_WOMERSLEY_REL_L2_TOLERANCE:-5e-2}
module load cuda/12.4.0

test_root="${case_dir}/.transient-cuda-validation-${SLURM_JOB_ID}"
mkdir -p "${test_root}"

womersley_case_dir=${IGA_WOMERSLEY_CASE_DIR:-}
womersley_manifest=${IGA_WOMERSLEY_MANIFEST:-}
womersley_config=${IGA_WOMERSLEY_CONFIG:-}
if [[ -n "${womersley_config}" && -n "${womersley_manifest}" ]]; then
    echo "Set either IGA_WOMERSLEY_CONFIG or IGA_WOMERSLEY_MANIFEST, not both" >&2
    exit 2
fi
if [[ (-n "${womersley_config}" || -n "${womersley_manifest}") \
  && -z "${IGA_CARDIAC_PERIOD:-}" ]]; then
  echo "Womersley series validation requires IGA_CARDIAC_PERIOD" >&2
  exit 2
fi

"${cuda_root}/iga_cuda" device-info
"${cuda_root}/iga_cuda" mesh-check "${database}"

"${cuda_root}/iga_cuda" navier-stokes "${database}" "${case_dir}" \
  --output "${test_root}/flow-full.txt" --output-every 1
"${cuda_root}/iga_cuda" navier-stokes "${database}" "${case_dir}" \
  --output "${test_root}/flow-interrupted.txt" \
  --checkpoint "${test_root}/flow-step1" --stop-after-step 1
"${cuda_root}/iga_cuda" navier-stokes "${database}" "${case_dir}" \
  --output "${test_root}/flow-restarted.txt" \
  --restart "${test_root}/flow-step1"

awk 'NR==FNR {for(i=1;i<=3;i++)a[NR,i]=$i;next}
  {for(i=1;i<=3;i++){d=$i-a[FNR,i];s+=d*d;r+=a[FNR,i]*a[FNR,i]}}
  END{e=(r? sqrt(s/r):(s? 1:0));printf "flow_restart_velocity_rel_l2=%.17g\n",e;exit(e>1e-12)}' \
  "${test_root}/flow-full.txt" "${test_root}/flow-restarted.txt"
awk 'NR==FNR {a[NR]=$1;next}{d=$1-a[FNR];s+=d*d;r+=a[FNR]*a[FNR]}
  END{e=(r? sqrt(s/r):(s? 1:0));printf "flow_restart_pressure_rel_l2=%.17g\n",e;exit(e>1e-12)}' \
  "${test_root}/flow-full.txt.pressure" "${test_root}/flow-restarted.txt.pressure"

flow_metrics=$("${cuda_root}/../cpu/iga_flow_validate" \
  "${database}" "${test_root}/flow-full.txt")
printf '%s\n' "${flow_metrics}"
mass_imbalance=$(awk '{for(i=1;i<=NF;i++){split($i,a,"=");if(a[1]=="relative_mass_imbalance")print a[2]}}' \
  <<<"${flow_metrics}")
divergence_error=$(awk '{for(i=1;i<=NF;i++){split($i,a,"=");if(a[1]=="relative_divergence_theorem_error")print a[2]}}' \
  <<<"${flow_metrics}")
awk -v value="${mass_imbalance}" -v gate="${mass_tolerance}" \
  'BEGIN{printf "mass_balance_gate=%.17g tolerance=%.17g\n",value,gate;exit(value>gate)}'
awk -v value="${divergence_error}" -v gate="${divergence_tolerance}" \
  'BEGIN{printf "relative_divergence_theorem_gate=%.17g tolerance=%.17g\n",value,gate;exit(value>gate)}'

if [[ -n "${IGA_CARDIAC_PERIOD:-}" ]]; then
  flow_dt=${IGA_FLOW_DT:?Set IGA_FLOW_DT with IGA_CARDIAC_PERIOD}
  flow_steps=${IGA_FLOW_STEPS:?Set IGA_FLOW_STEPS with IGA_CARDIAC_PERIOD}
  manifest="${test_root}/flow-series.csv"
  printf 'time,file\n' > "${manifest}"
  for ((step=1; step<=flow_steps; ++step)); do
    time=$(awk -v step="${step}" -v dt="${flow_dt}" 'BEGIN{printf "%.17g",step*dt}')
    printf -v snapshot 'flow-full.step%06d.txt' "${step}"
    printf '%s,%s\n' "${time}" "${snapshot}" >> "${manifest}"
  done
  cycle_metrics=$("${cuda_root}/../cpu/iga_flow_validate" "${database}" \
    --manifest "${test_root}" "$(basename "${manifest}")" "${IGA_CARDIAC_PERIOD}")
  printf '%s\n' "${cycle_metrics}"
  cycle_imbalance=$(awk '{for(i=1;i<=NF;i++){split($i,a,"=");if(a[1]=="maximum_relative_mass_imbalance")print a[2]}}' \
    <<<"${cycle_metrics}")
  cycle_error=$(awk '{for(i=1;i<=NF;i++){split($i,a,"=");if(a[1]=="maximum_cycle_relative_velocity_l2")print a[2]}}' \
    <<<"${cycle_metrics}")
  awk -v value="${cycle_imbalance}" -v gate="${mass_tolerance}" \
    'BEGIN{printf "multi_cycle_mass_gate=%.17g tolerance=%.17g\n",value,gate;exit(value>gate)}'
  awk -v value="${cycle_error}" -v gate="${cycle_tolerance}" \
    'BEGIN{printf "cycle_repeatability_gate=%.17g tolerance=%.17g\n",value,gate;exit(value>gate)}'

  if [[ -n "${womersley_config}" ]]; then
    womersley_metrics=$("${cuda_root}/../cpu/iga_flow_validate" "${database}" \
      --womersley "${womersley_config}" "${test_root}" "$(basename "${manifest}")")
    printf '%s\n' "${womersley_metrics}"
    womersley_error=$(awk '{for(i=1;i<=NF;i++){split($i,a,"=");if(a[1]=="maximum_womersley_volume_relative_l2")print a[2]}}' \
      <<<"${womersley_metrics}")
    awk -v value="${womersley_error}" -v gate="${womersley_tolerance}" \
      'BEGIN{printf "womersley_volume_gate=%.17g tolerance=%.17g\n",value,gate;exit(value>gate)}'
  elif [[ -n "${womersley_manifest}" ]]; then
    if [[ -z "${womersley_case_dir}" ]]; then
      echo "Set IGA_WOMERSLEY_CASE_DIR with IGA_WOMERSLEY_MANIFEST" >&2
      exit 2
    fi
    womersley_metrics=$("${cuda_root}/../cpu/iga_flow_validate" "${database}" \
      --compare-manifests "${womersley_case_dir}" "${womersley_manifest}" \
      "${test_root}" "$(basename "${manifest}")")
    printf '%s\n' "${womersley_metrics}"
    womersley_error=$(awk '{for(i=1;i<=NF;i++){split($i,a,"=");if(a[1]=="maximum_relative_velocity_l2")print a[2]}}' \
      <<<"${womersley_metrics}")
    awk -v value="${womersley_error}" -v gate="${womersley_tolerance}" \
      'BEGIN{printf "projected_womersley_series_gate=%.17g tolerance=%.17g\n",value,gate;exit(value>gate)}'
  fi
fi

"${cuda_root}/iga_cuda" solve "${database}" "${case_dir}" \
  --system "${transport_system}" --output "${test_root}/transport-full.txt" \
  --output-every 1
"${cuda_root}/iga_cuda" solve "${database}" "${case_dir}" \
  --system "${transport_system}" --output "${test_root}/transport-interrupted.txt" \
  --checkpoint "${test_root}/transport-step1" \
  --stop-after-step 1
"${cuda_root}/iga_cuda" solve "${database}" "${case_dir}" \
  --system "${transport_system}" --output "${test_root}/transport-restarted.txt" \
  --restart "${test_root}/transport-step1"

awk 'NR==FNR {for(i=2;i<=NF;i++)a[NR,i]=$i;next}
  {for(i=2;i<=NF;i++){d=$i-a[FNR,i];s+=d*d;r+=a[FNR,i]*a[FNR,i]}}
  END{e=(r? sqrt(s/r):(s? 1:0));printf "transport_restart_rel_l2=%.17g\n",e;exit(e>1e-12)}' \
  "${test_root}/transport-full.txt" "${test_root}/transport-restarted.txt"

if [[ -n "${IGA_CPU_FLOW_REFERENCE:-}" ]]; then
  awk -v gate="${parity_tolerance}" \
    'NR==FNR {for(i=1;i<=3;i++)a[NR,i]=$i;next}
    {for(i=1;i<=3;i++){d=$i-a[FNR,i];s+=d*d;r+=a[FNR,i]*a[FNR,i]}}
    END{e=sqrt(s/r);printf "cpu_cuda_flow_velocity_rel_l2=%.17g\n",e;exit(e>gate)}' \
    "${IGA_CPU_FLOW_REFERENCE}" "${test_root}/flow-full.txt"
fi

if [[ -n "${IGA_CPU_TRANSPORT_REFERENCE:-}" ]]; then
  awk -v gate="${parity_tolerance}" \
    'NR==FNR {for(i=2;i<=NF;i++)a[NR,i]=$i;next}
    {for(i=2;i<=NF;i++){d=$i-a[FNR,i];s+=d*d;r+=a[FNR,i]*a[FNR,i]}}
    END{e=sqrt(s/r);printf "cpu_cuda_transport_rel_l2=%.17g\n",e;exit(e>gate)}' \
    "${IGA_CPU_TRANSPORT_REFERENCE}" "${test_root}/transport-full.txt"
fi

echo "validation_output=${test_root}"
