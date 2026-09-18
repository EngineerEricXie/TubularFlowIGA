#!/bin/bash
# Matched complete small-FSI screen for P3. Computation is Slurm/node-local.
set -euo pipefail
while test $# -gt 0; do
    case "$1" in --root) root=$2; shift 2;; --variants) variants=$2; shift 2;; --case-list) cases=$2; shift 2;; *) exit 2;; esac
done
: "${LOCAL:?Slurm LOCAL required}" "${SLURM_JOB_ID:?Slurm required}"
IFS=, read -r baseline candidate verifier <<< "$variants"
: "${baseline:?}" "${candidate:?}" "${verifier:?}"
work="$LOCAL/fsi-p3-gate-$SLURM_JOB_ID"; out="$root/job-$SLURM_JOB_ID"
mkdir -p "$work/evidence/binaries" "$out"
stage=staging; failure_status=failed_environment; started=$(date -u +%FT%TZ)
save_stage() {
    printf '{"job_id":"%s","stage":"%s","status":"in_progress"}\n' "$SLURM_JOB_ID" "$stage" > "$work/evidence/stage.json"
    cp -au "$work/evidence/." "$out/"
}
finish() {
    result=$?; trap - EXIT TERM USR1
    test "$result" -ne 124 || failure_status=timeout
    printf '{"job_id":"%s","stage":"%s","status":"%s","exit_code":%s,"started":"%s","ended":"%s"}\n' \
        "$SLURM_JOB_ID" "$stage" "$([ "$result" -eq 0 ] && echo passed || echo "$failure_status")" "$result" "$started" "$(date -u +%FT%TZ)" > "$work/evidence/status.json"
    copied=0; cp -au "$work/evidence/." "$out/" || copied=$?
    test "$copied" -eq 0 || { test "$result" -ne 0 || result=$copied; }
    exit "$result"
}
trap finish EXIT
trap 'failure_status=interrupted; exit 124' TERM USR1
cp "$cases" "$work/evidence/case.json"
cp "$0" "$work/evidence/runner.sh"
for variant in "$baseline" "$candidate" "$verifier"; do
    sha256sum -c "$root/$variant-source.sha256"
    mkdir -p "$work/$variant"; tar -xf "$root/$variant-source.tar" -C "$work/$variant"
    cp "$root/$variant-source.sha256" "$work/evidence/$variant-source.sha256"
done
common="$work/$verifier"; tools="$common/benchmarks/fsi-time-opt"
read -r budget_s expected_threads < <(python3 -c 'import json,sys; c=json.load(open(sys.argv[1])); print(c["per_process_timeout_seconds"],c["threads"])' "$cases")
test "$OMP_NUM_THREADS" -eq "$expected_threads"
{
    date -u +%FT%TZ; module list; mpicxx --version; lscpu; scontrol show job "$SLURM_JOB_ID"; grep Cpus_allowed_list /proc/self/status
    for name in OMP_NUM_THREADS OMP_THREAD_LIMIT OMP_DYNAMIC OMP_PROC_BIND OMP_PLACES IGA_ASSEMBLY_THREADS IGA_ASSEMBLY_BATCH_SIZE PETSC_DIR PETSC_ARCH; do printf '%s=%s\n' "$name" "${!name:-}"; done
} > "$work/evidence/environment.txt" 2>&1
cp "$PETSC_DIR/$PETSC_ARCH/lib/petsc/conf/petscvariables" "$work/evidence/petscvariables.txt"
sha256sum "$common/solvers/cpu/tests/test_compliant_channel_fsi.cpp" "$tools/verify_p3.py" "$cases" "$0" > "$work/evidence/input-runner.sha256"
for variant in "$baseline" "$candidate"; do
    stage="build-$variant-fsi"; save_stage
    timeout -k 30s 600s /usr/bin/time -v -o "$work/evidence/$variant-build.time" \
        make -f "$tools/fsi_case.mk" SOURCE="$work/$variant" BUILD="$work/$variant" \
        CASE_SOURCE="$common/solvers/cpu/tests/test_compliant_channel_fsi.cpp" CASE_NAME=fsi fsi_case \
        > "$work/evidence/$variant-build.stdout" 2> "$work/evidence/$variant-build.stderr"
    cp "$work/$variant/fsi" "$work/evidence/binaries/$variant-fsi"
    sha256sum "$work/$variant/fsi" >> "$work/evidence/binary.sha256"
    ldd "$work/$variant/fsi" > "$work/evidence/$variant-ldd.txt"
done
for variant in "$baseline" "$candidate"; do
    stage="small-FSI-$variant"; failure_status=failed_numerical; save_stage
    sample="$work/evidence/$variant-fsi"; mkdir -p "$sample"; cd "$sample"
    batch=$(python3 -c 'import json,sys; c=json.load(open(sys.argv[1])); print(c.get("batch_by_variant",{}).get(sys.argv[2],c["batch"]))' "$cases" "$variant")
    printf '%s\n' "$batch" > batch-size.txt
    IGA_ASSEMBLY_BATCH_SIZE=$batch IGA_PROFILE_DETAIL=1 timeout -k 30s "${budget_s}s" \
        /usr/bin/time -v -o timing.txt mpiexec -np 1 --map-by slot:PE="$OMP_NUM_THREADS" --bind-to core \
        "$work/$variant/fsi" --reference-output reference > stdout.txt 2> stderr.txt
    stage="completed-small-FSI-$variant"; save_stage
done
stage=verification; failure_status=inconclusive; save_stage
python3 "$tools/verify_p3.py" "$work/evidence" "$baseline" "$candidate"
stage=complete
