#!/bin/bash
set -euo pipefail
while test $# -gt 0; do
    case "$1" in --root) root=$2; shift 2;; --variants) variant=$2; shift 2;; --case-list) cases=$2; shift 2;; *) exit 2;; esac
done
: "${LOCAL:?Slurm LOCAL required}" "${SLURM_JOB_ID:?Slurm required}"
work="$LOCAL/fsi-proof-$SLURM_JOB_ID"
out="$root/job-$SLURM_JOB_ID"
mkdir -p "$work/evidence" "$out"
stage=staging; failure_status=failed_environment
finish() {
    result=$?; trap - EXIT
    printf '{"job_id":"%s","stage":"%s","status":"%s","exit_code":%s,"ended_utc":"%s"}\n' "$SLURM_JOB_ID" "$stage" "$([ "$result" -eq 0 ] && echo passed || echo "$failure_status")" "$result" "$(date -u +%FT%TZ)" > "$work/evidence/status.json"
    copied=0; cp -a "$work/evidence/." "$out/" || copied=$?
    if test "$copied" -ne 0; then echo "copy-back failed: $work" >&2; test "$result" -ne 0 || result=$copied; fi
    exit "$result"
}
trap finish EXIT
sha256sum -c "$root/$variant-source.sha256"
mkdir -p "$work/source"
tar -xf "$root/$variant-source.tar" -C "$work/source"
src="$work/source"
cp "$cases" "$work/evidence/case.json"
cp "$root/$variant-source.sha256" "$work/evidence/source.sha256"
sha256sum "$0" "$cases" > "$work/evidence/runner-input.sha256"
{ date -u +%FT%TZ; module list; mpicxx --version; lscpu; scontrol show job "$SLURM_JOB_ID"; grep Cpus_allowed_list /proc/self/status; } > "$work/evidence/environment.txt" 2>&1
stage=build
/usr/bin/time -v -o "$work/evidence/build.time" make -f "$src/benchmarks/fsi-time-opt/fsi_case.mk" SOURCE="$src" BUILD="$work" CASE_SOURCE=test_compliant_channel_fsi.cpp CASE_NAME=probe fsi_case > "$work/evidence/build.stdout" 2> "$work/evidence/build.stderr"
cp "$work/probe" "$work/evidence/probe"
sha256sum "$work/probe" > "$work/evidence/binary.sha256"
ldd "$work/probe" > "$work/evidence/ldd.txt"
stage=FSI-probe; failure_status=failed_numerical
cd "$work/evidence"
IGA_PROFILE_DETAIL=1 /usr/bin/time -v -o timing.txt mpiexec -np 1 --map-by slot:PE="$OMP_NUM_THREADS" --bind-to core "$work/probe" --reference-output reference > stdout.txt 2> stderr.txt
stage=proof-verification
python3 "$src/benchmarks/fsi-time-opt/verify_fsi_probe.py" "$work/evidence" "$variant"
stage=complete
