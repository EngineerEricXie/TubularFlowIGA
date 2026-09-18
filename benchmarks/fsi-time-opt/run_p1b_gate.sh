#!/bin/bash
# One bounded primary job, durable completed stages, no scheduler monitoring.
set -euo pipefail
while test $# -gt 0; do
    case "$1" in --root) root=$2; shift 2;; --variants) variants=$2; shift 2;; --case-list) cases=$2; shift 2;; *) exit 2;; esac
done
: "${LOCAL:?}" "${SLURM_JOB_ID:?}"
IFS=, read -r baseline candidate validation <<< "$variants"
validation=${validation:-$candidate}
work="$LOCAL/fsi-p1b-$SLURM_JOB_ID"; out="$root/job-$SLURM_JOB_ID"
mkdir -p "$work/evidence/binaries" "$out"
stage=staging
save_stage() { printf '{"stage":"%s"}\n' "$stage" > "$work/evidence/stage.json"; cp -au "$work/evidence/." "$out/"; }
finish() {
    result=$?; trap - EXIT TERM USR1
    printf '{"job_id":"%s","stage":"%s","exit_code":%s,"status":"%s"}\n' "$SLURM_JOB_ID" "$stage" "$result" "$([ "$result" -eq 0 ] && echo passed || echo failed)" > "$work/evidence/status.json"
    cp -au "$work/evidence/." "$out/" || { test "$result" -ne 0 || result=1; }
    exit "$result"
}
trap finish EXIT; trap 'exit 124' TERM USR1
cp "$0" "$work/evidence/runner.sh"; cp "$cases" "$work/evidence/case.json"
for variant in "$baseline" "$candidate" "$validation"; do
    sha256sum -c "$root/$variant-source.sha256"
    mkdir -p "$work/$variant"; tar -xf "$root/$variant-source.tar" -C "$work/$variant"
    cp "$root/$variant-source.sha256" "$work/evidence/$variant-source.sha256"
done
tools="$work/$validation/benchmarks/fsi-time-opt"; common="$work/$validation"
{ date -u +%FT%TZ; module list; mpicxx --version; lscpu; scontrol show job "$SLURM_JOB_ID";
  grep Cpus_allowed_list /proc/self/status; env | sort | grep -E '^(OMP_|IGA_|PETSC_|SLURM_CPUS)';
} > "$work/evidence/environment.txt" 2>&1
stage=before-proof; save_stage
python3 "$tools/verify_p1b.py" --proof "$root/job-46071831" "$work/evidence"
python3 "$tools/verify_p1b.py" --prefix "$root" "$work/evidence" "$cases" "$common"
stage=build-targeted; save_stage
timeout -k 30s 600s /usr/bin/time -v -o "$work/evidence/targeted-build.time" make -f "$tools/fsi_case.mk" SOURCE="$work/$candidate" BUILD="$common" CASE_SOURCE="$common/solvers/cpu/tests/test_accepted_line_search_assembly.cpp" CASE_NAME=accepted-full fsi_case > "$work/evidence/targeted-build.stdout" 2> "$work/evidence/targeted-build.stderr"
cp "$common/accepted-full" "$work/evidence/binaries/accepted-full"
sha256sum "$common/accepted-full" > "$work/evidence/binary.sha256"
stage=targeted-regression; save_stage
mkdir -p "$work/evidence/accepted-full"; cd "$work/evidence/accepted-full"
mapfile -t test_options < <(python3 -c 'import json,sys; c=json.load(open(sys.argv[1])); print("--reuse-case\n"+c["only_case"]) if c.get("only_case") else None' "$cases")
IGA_PROFILE_DETAIL=1 timeout -k 30s 3000s /usr/bin/time -v -o timing.txt mpiexec -np 1 --map-by slot:PE="$OMP_NUM_THREADS" --bind-to core "$common/accepted-full" "${test_options[@]}" > stdout.txt 2> stderr.txt
stage=completed-targeted; save_stage
for variant in "$baseline" "$candidate"; do
    stage="build-$variant-fsi"; save_stage
    timeout -k 30s 600s /usr/bin/time -v -o "$work/evidence/$variant-build.time" make -f "$tools/fsi_case.mk" SOURCE="$work/$variant" BUILD="$work/$variant" CASE_SOURCE="$common/solvers/cpu/tests/test_compliant_channel_fsi.cpp" CASE_NAME=fsi fsi_case > "$work/evidence/$variant-build.stdout" 2> "$work/evidence/$variant-build.stderr"
    cp "$work/$variant/fsi" "$work/evidence/binaries/$variant-fsi"
    sha256sum "$work/$variant/fsi" >> "$work/evidence/binary.sha256"
    stage="small-FSI-$variant"; save_stage
    sample="$work/evidence/$variant-fsi"; mkdir -p "$sample"; cd "$sample"
    IGA_PROFILE_DETAIL=1 timeout -k 30s 1400s /usr/bin/time -v -o timing.txt mpiexec -np 1 --map-by slot:PE="$OMP_NUM_THREADS" --bind-to core "$work/$variant/fsi" --reference-output reference > stdout.txt 2> stderr.txt
    stage="completed-small-FSI-$variant"; save_stage
done
stage=verification; save_stage
python3 "$tools/verify_p1b.py" "$work/evidence" "$baseline" "$candidate"
stage=complete
