#!/bin/bash
# Durable, bounded gate execution. No periodic monitoring or log polling.
set -euo pipefail
while test $# -gt 0; do
    case "$1" in --root) root=$2; shift 2;; --variants) variants=$2; shift 2;; --case-list) cases=$2; shift 2;; *) exit 2;; esac
done
: "${LOCAL:?Slurm LOCAL required}" "${SLURM_JOB_ID:?Slurm required}"
IFS=, read -r baseline candidate validation <<< "$variants"
validation=${validation:-P1AValidation}
work="$LOCAL/fsi-p1a-gate-$SLURM_JOB_ID"; out="$root/job-$SLURM_JOB_ID"
mkdir -p "$work/evidence/binaries" "$out"
stage=staging; failure_status=failed_environment; started=$(date -u +%FT%TZ)
save_stage() {
    printf '{"job_id":"%s","stage":"%s","status":"in_progress","started":"%s"}\n' "$SLURM_JOB_ID" "$stage" "$started" > "$work/evidence/stage.json"
    cp -a "$work/evidence/." "$out/"
}
finish() {
    result=$?; trap - EXIT TERM USR1
    test "$result" -ne 124 || failure_status=timeout
    printf '{"job_id":"%s","stage":"%s","status":"%s","exit_code":%s,"started":"%s","ended":"%s"}\n' "$SLURM_JOB_ID" "$stage" "$([ "$result" -eq 0 ] && echo passed || echo "$failure_status")" "$result" "$started" "$(date -u +%FT%TZ)" > "$work/evidence/status.json"
    copied=0; cp -a "$work/evidence/." "$out/" || copied=$?
    if test "$copied" -ne 0; then echo "copy-back failed: $work" >&2; test "$result" -ne 0 || result=$copied; fi
    exit "$result"
}
trap finish EXIT
trap 'failure_status=interrupted; exit 124' TERM USR1
cp "$cases" "$work/evidence/case.json"
sha256sum "$0" "$cases" > "$work/evidence/runner-input.sha256"
{ date -u +%FT%TZ; module list; mpicxx --version; lscpu; scontrol show job "$SLURM_JOB_ID"; grep Cpus_allowed_list /proc/self/status;
  for name in OMP_NUM_THREADS OMP_THREAD_LIMIT OMP_DYNAMIC OMP_PROC_BIND OMP_PLACES IGA_ASSEMBLY_THREADS IGA_ASSEMBLY_BATCH_SIZE PETSC_DIR PETSC_ARCH; do printf '%s=%s\n' "$name" "${!name:-}"; done;
} > "$work/evidence/environment.txt" 2>&1
cp "$PETSC_DIR/$PETSC_ARCH/lib/petsc/conf/petscvariables" "$work/evidence/petscvariables.txt"
for variant in "$baseline" "$candidate" "$validation"; do
    sha256sum -c "$root/$variant-source.sha256"
    mkdir -p "$work/$variant"
    tar -xf "$root/$variant-source.tar" -C "$work/$variant"
    cp "$root/$variant-source.sha256" "$work/evidence/$variant-source.sha256"
done
common="$work/$validation"; tools="$common/benchmarks/fsi-time-opt"
python3 -c 'import json,pathlib,sys; c=json.load(open(sys.argv[1])); p=c.get("passed_paired_adapter_job"); assert not p or (pathlib.Path(sys.argv[2])/f"job-{p}"/"workload-test.sha256").is_file(), "prior adapter evidence missing before computation"' "$cases" "$root"
read -r matched budget_s paired w1 < <(python3 -c 'import json,sys; c=json.load(open(sys.argv[1])); print(int(c["matched_complete_fsi"]),c["per_process_timeout_seconds"],int(c.get("paired_regressions",False)),int(c.get("w1_operator",False)))' "$cases")
mapfile -t regressions < <(python3 -c 'import json,sys; sys.stdout.write("\n".join(json.load(open(sys.argv[1]))["regressions"]))' "$cases")
if test "$matched" -eq 1; then
    for variant in "$baseline" "$candidate"; do
        src="$work/$variant"; stage="build-$variant-fsi"; save_stage
        timeout -k 30s 600s /usr/bin/time -v -o "$work/evidence/$variant-build.time" make -f "$tools/fsi_case.mk" SOURCE="$src" BUILD="$work/$variant" CASE_SOURCE="$common/solvers/cpu/tests/test_compliant_channel_fsi.cpp" CASE_NAME=fsi fsi_case > "$work/evidence/$variant-build.stdout" 2> "$work/evidence/$variant-build.stderr"
        cp "$work/$variant/fsi" "$work/evidence/binaries/$variant-fsi"
        sha256sum "$work/$variant/fsi" > "$work/evidence/$variant-binary.sha256"
        ldd "$work/$variant/fsi" > "$work/evidence/$variant-ldd.txt"
    done
fi
for name in "${regressions[@]}"; do
    case "$name" in immersed) source=test_immersed_transient_flow.cpp;; adapter) source=test_moving_immersed_transient_flow_fsi_runtime.cpp;; volume) source=test_parallel_immersed_volume.cpp;; *) exit 2;; esac
    sha256sum "$common/solvers/cpu/tests/$source" >> "$work/evidence/workload-test.sha256"
    test_variants=("$candidate"); test "$paired" -eq 0 || test_variants=("$baseline" "$candidate")
    for variant in "${test_variants[@]}"; do
        label=$name; test "$variant" != "$baseline" || label=$name-before
        stage="build-$variant-$name"; save_stage
        timeout -k 30s 600s /usr/bin/time -v -o "$work/evidence/$label-build.time" make -f "$tools/fsi_case.mk" SOURCE="$work/$variant" BUILD="$work/$variant" CASE_SOURCE="$common/solvers/cpu/tests/$source" CASE_NAME="$name" fsi_case > "$work/evidence/$label-build.stdout" 2> "$work/evidence/$label-build.stderr"
        cp "$work/$variant/$name" "$work/evidence/binaries/$variant-$name"
        sha256sum "$work/$variant/$name" >> "$work/evidence/$variant-binary.sha256"
        stage="regression-$variant-$name"; failure_status=failed_numerical; save_stage
        sample="$work/evidence/$label"; mkdir -p "$sample"; cd "$sample"
        IGA_PROFILE_DETAIL=1 timeout -k 30s "${budget_s}s" /usr/bin/time -v -o timing.txt mpiexec -np 1 --map-by slot:PE="$OMP_NUM_THREADS" --bind-to core "$work/$variant/$name" > stdout.txt 2> stderr.txt
        stage="completed-$variant-$name"; save_stage
    done
done
if test "$matched" -eq 1; then
    for variant in "$baseline" "$candidate"; do
        stage="V2-$variant"; failure_status=failed_environment; save_stage
        sample="$work/evidence/$variant-fsi"; mkdir -p "$sample"; cd "$sample"
        failure_status=failed_numerical
        IGA_PROFILE_DETAIL=1 timeout -k 30s "${budget_s}s" /usr/bin/time -v -o timing.txt mpiexec -np 1 --map-by slot:PE="$OMP_NUM_THREADS" --bind-to core "$work/$variant/fsi" --reference-output reference > stdout.txt 2> stderr.txt
        stage="completed-V2-$variant"; save_stage
    done
fi
if test "$w1" -eq 1; then
    stage=W1-input-staging; failure_status=failed_environment; save_stage
    sha256sum -c "$root/surface.sha256"; sha256sum -c "$root/archived-workload.sha256"
    cp "$root/surface.txt" "$work/surface.txt"
    tar -xOf "$root/archived-workload.tar.gz" ./driver.cpp > "$tools/driver.cpp"
    sha256sum "$work/surface.txt" "$tools/driver.cpp" "$tools/benchmark.cpp" > "$work/evidence/W1-input.sha256"
    for variant in "$baseline" "$candidate"; do
        stage="build-$variant-W1"; save_stage
        timeout -k 30s 600s /usr/bin/time -v -o "$work/evidence/$variant-W1-build.time" make -f "$tools/runner.mk" SOURCE="$work/$variant" WORKLOAD_SOURCE="$common" BUILD="$work/$variant" benchmark > "$work/evidence/$variant-W1-build.stdout" 2> "$work/evidence/$variant-W1-build.stderr"
        cp "$work/$variant/benchmark" "$work/evidence/binaries/$variant-benchmark"
        sha256sum "$work/$variant/benchmark" >> "$work/evidence/$variant-binary.sha256"
    done
    for pair in 0 1 2; do
        order=("$baseline" "$candidate"); test "$pair" -ne 1 || order=("$candidate" "$baseline")
        for variant in "${order[@]}"; do
            stage="W1-$pair-$variant"; save_stage
            sample="$work/evidence/W1-pair-$pair-$variant"; mkdir -p "$sample"; cd "$sample"
            failure_status=failed_numerical
            IGA_PROFILE_DETAIL=1 timeout -k 30s 300s /usr/bin/time -v -o timing.txt mpiexec -np 1 --map-by slot:PE="$OMP_NUM_THREADS" --bind-to core "$work/$variant/benchmark" "$work/surface.txt" > stdout.txt 2> stderr.txt
            stage="completed-W1-$pair-$variant"; save_stage
        done
    done
fi
stage=verification; failure_status=inconclusive; save_stage
python3 "$tools/verify_p1a.py" "$work/evidence" "$baseline" "$candidate" "$root"
stage=complete
