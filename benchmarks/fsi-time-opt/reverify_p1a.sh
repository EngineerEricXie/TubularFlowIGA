#!/bin/bash
# Reanalyze immutable completed output only. No solver or compiler execution.
set -euo pipefail
module load anaconda3
root=${1:?evidence root}; validation=${2:?frozen verifier}; compute=${3:?completed compute job}
: "${LOCAL:?Slurm LOCAL required}" "${SLURM_JOB_ID:?Slurm required}"
work="$LOCAL/fsi-reverify-$SLURM_JOB_ID"; out="$root/job-$SLURM_JOB_ID"
mkdir -p "$work/evidence" "$out"
finish() {
    result=$?; trap - EXIT
    printf '{"job_id":"%s","compute_job_id":"%s","verification_only":true,"status":"%s","exit_code":%s}\n' "$SLURM_JOB_ID" "$compute" "$([ "$result" -eq 0 ] && echo passed || echo inconclusive)" "$result" > "$work/evidence/status.json"
    copied=0; cp -a "$work/evidence/." "$out/" || copied=$?
    test "$copied" -eq 0 || { test "$result" -ne 0 || result=$copied; }
    exit "$result"
}
trap finish EXIT
sha256sum -c "$root/$validation-source.sha256"
mkdir -p "$work/source"
tar -xf "$root/$validation-source.tar" -C "$work/source"
cp -a "$root/job-$compute/." "$work/evidence/"
mv "$work/evidence/status.json" "$work/evidence/original-compute-status.json"
cp "$root/$validation-source.sha256" "$work/evidence/reverification-source.sha256"
sha256sum "$0" > "$work/evidence/reverification-runner.sha256"
{ date -u +%FT%TZ; python3 --version; scontrol show job "$SLURM_JOB_ID"; } > "$work/evidence/verification-environment.txt" 2>&1
# These refer to the immutable numerical configuration, not analysis core count.
export OMP_NUM_THREADS=8 IGA_ASSEMBLY_BATCH_SIZE=8
export IGA_VERIFICATION_ONLY=1 IGA_COMPUTE_JOB_ID="$compute"
export IGA_COMPUTE_NODE=$(sed -n 's/^ *NodeList=\([^ ]*\).*/\1/p' "$work/evidence/environment.txt")
python3 "$work/source/benchmarks/fsi-time-opt/verify_p1a.py" "$work/evidence" P1ABefore P1A "$root"
