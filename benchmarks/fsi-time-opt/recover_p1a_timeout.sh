#!/bin/bash
# Diagnostic recovery only; never use a previous allocation as an execution input.
set -u
out=/ocean/projects/mch260002p/thsieh1/TubularFlowIGA/benchmarks/fsi-time-opt/20260915T035654Z/job-46022197
printf 'recovery_job=%s node=%s LOCAL=%s\n' "$SLURM_JOB_ID" "$SLURM_JOB_NODELIST" "${LOCAL:-unset}"
found=0
for base in /local /localscratch; do
    test -d "$base" || continue
    while IFS= read -r path; do
        test -d "$path/evidence" || continue
        found=1
        mkdir -p "$out/recovered"
        cp -a "$path/evidence/." "$out/recovered/" || exit 2
        printf 'recovered_from=%s\n' "$path"
    done < <(find "$base" -maxdepth 4 -type d -name fsi-p1a-46022197 2>/dev/null)
done
printf 'evidence_found=%s\n' "$found"
