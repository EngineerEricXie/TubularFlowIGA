#!/bin/bash
set -euo pipefail
while test $# -gt 0; do
    case "$1" in
        --root) root=$2; shift 2;;
        --variants) variant=$2; shift 2;;
        --case-list) cases=$2; shift 2;;
        *) exit 2;;
    esac
done
: "${LOCAL:?Slurm LOCAL required}" "${SLURM_JOB_ID:?Slurm required}" "${root:?}" "${variant:?}" "${cases:?}"
work="$LOCAL/fsi-r6-w2-$SLURM_JOB_ID"
out="$root/job-$SLURM_JOB_ID"
mkdir -p "$work/evidence" "$out" "$work/$variant"
stage=staging
failure_status=failed_environment
simulation_executed=false
started=$(date -u +%FT%TZ)
finish() {
    result=$?
    trap - EXIT TERM USR1
    test "$result" -ne 124 || failure_status=timeout
    printf '{"job_id":"%s","stage":"%s","status":"%s","exit_code":%s,"simulation_executed":%s,"started":"%s","ended":"%s"}\n' "$SLURM_JOB_ID" "$stage" "$([ "$result" -eq 0 ] && echo passed || echo "$failure_status")" "$result" "$simulation_executed" "$started" "$(date -u +%FT%TZ)" > "$work/evidence/status.json"
    cp -au "$work/evidence/." "$out/"
    exit "$result"
}
trap finish EXIT
trap 'failure_status=interrupted; exit 124' TERM USR1

read -r expected_variant selection_job < <(python3 -c 'import json,sys; c=json.load(open(sys.argv[1])); print(c["variant"],c["selection_job"])' "$cases")
test "$variant" = "$expected_variant"
threads=$(tr -d '[:space:]' < "$root/job-$selection_job/selected-w2-threads.txt")
expected=$(python3 -c 'import json,sys; print(",".join(map(str,json.load(open(sys.argv[1]))["expected_threads"])))' "$cases")
test "$threads" = "$expected"
sha256sum -c "$root/$variant-source.sha256"
tar -xf "$root/$variant-source.tar" -C "$work/$variant"
sha256sum -c "$root/surface.sha256"
cp "$root/surface.txt" "$work/surface.txt"
cp "$root/historical-heartbeat-history.csv" "$work/evidence/historical-history.csv"
cp "$root/$variant-source.sha256" "$0" "$work/evidence/"
cp "$cases" "$work/evidence/case.json"
cp "$root/job-$selection_job/selected-w2-threads.txt" "$work/evidence/selected-w2-threads.txt"
src="$work/$variant"
stage=build
timeout -k 30s 900s /usr/bin/time -v -o "$work/evidence/build.time" make -f "$src/benchmarks/fsi-time-opt/heartbeat_case.mk" SOURCE="$src" WORKLOAD_SOURCE="$src" BUILD="$src" heartbeat > "$work/evidence/build.stdout" 2> "$work/evidence/build.stderr"
cp "$src/heartbeat" "$work/evidence/heartbeat"
sha256sum "$src/heartbeat" > "$work/evidence/binary.sha256"
ldd "$src/heartbeat" > "$work/evidence/ldd.txt"
cp "$PETSC_DIR/$PETSC_ARCH/lib/petsc/conf/petscvariables" "$work/evidence/petscvariables.txt"
{
    date -u +%FT%TZ
    module list
    lscpu
    scontrol show job "$SLURM_JOB_ID"
    grep Cpus_allowed_list /proc/self/status
} > "$work/evidence/environment.txt" 2>&1
options=(-demo_steps 2 -demo_grid 16 -demo_dt .05 -demo_pressure 5
    -demo_pulse_amplitude 30 -demo_period 1 -demo_cut_depth 1
    -demo_wall_inertial_gamma 1 -demo_display_scale 100 -demo_native_reference true)
printf '%q ' "${options[@]}" > "$work/evidence/solver-options.txt"

simulation_executed=true
failure_status=failed_numerical
IFS=, read -r -a top_two <<< "$threads"
for count in "${top_two[@]}"; do
    stage="W2-threads-$count"
    sample="$work/evidence/W2-threads-$count"
    mkdir -p "$sample"
    cd "$sample"
    OMP_NUM_THREADS=$count OMP_THREAD_LIMIT=$count IGA_ASSEMBLY_THREADS=$count IGA_ASSEMBLY_BATCH_SIZE=32 IGA_PROFILE_DETAIL=1 timeout -k 30s 3200s /usr/bin/time -v -o timing.txt mpiexec -np 1 --map-by slot:PE="$count" --bind-to core "$src/heartbeat" "$work/surface.txt" "$sample/results" "${options[@]}" > stdout.txt 2> stderr.txt
done

stage=verification
python3 "$src/benchmarks/fsi-time-opt/verify_r6_w2.py" "$work/evidence"
stage=complete
