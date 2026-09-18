#!/bin/bash
# One allocation: regression gates and matched complete FSI before/after.
set -euo pipefail
while test $# -gt 0; do
    case "$1" in --root) root=$2; shift 2;; --variants) variants=$2; shift 2;; --case-list) cases=$2; shift 2;; *) exit 2;; esac
done
: "${LOCAL:?Slurm LOCAL required}" "${SLURM_JOB_ID:?Slurm required}"
IFS=, read -r baseline candidate <<< "$variants"
work="$LOCAL/fsi-p1a-$SLURM_JOB_ID"; out="$root/job-$SLURM_JOB_ID"
mkdir -p "$work/evidence" "$out"
stage=staging; failure_status=failed_environment; started=$(date -u +%FT%TZ)
finish() {
    result=$?; trap - EXIT
    printf '{"job_id":"%s","stage":"%s","status":"%s","exit_code":%s,"started":"%s","ended":"%s"}\n' "$SLURM_JOB_ID" "$stage" "$([ "$result" -eq 0 ] && echo passed || echo "$failure_status")" "$result" "$started" "$(date -u +%FT%TZ)" > "$work/evidence/status.json"
    copied=0; cp -a "$work/evidence/." "$out/" || copied=$?
    if test "$copied" -ne 0; then echo "copy-back failed: $work" >&2; test "$result" -ne 0 || result=$copied; fi
    exit "$result"
}
trap finish EXIT
cp "$cases" "$work/evidence/case.json"
sha256sum "$0" "$cases" > "$work/evidence/runner-input.sha256"
{ date -u +%FT%TZ; module list; mpicxx --version; lscpu; scontrol show job "$SLURM_JOB_ID"; grep Cpus_allowed_list /proc/self/status;
  for name in OMP_NUM_THREADS OMP_THREAD_LIMIT OMP_DYNAMIC OMP_PROC_BIND OMP_PLACES IGA_ASSEMBLY_THREADS IGA_ASSEMBLY_BATCH_SIZE PETSC_DIR PETSC_ARCH; do printf '%s=%s\n' "$name" "${!name:-}"; done;
} > "$work/evidence/environment.txt" 2>&1
cp "$PETSC_DIR/$PETSC_ARCH/lib/petsc/conf/petscvariables" "$work/evidence/petscvariables.txt"
for variant in "$baseline" "$candidate"; do
    sha256sum -c "$root/$variant-source.sha256"
    mkdir -p "$work/$variant"
    tar -xf "$root/$variant-source.tar" -C "$work/$variant"
    cp "$root/$variant-source.sha256" "$work/evidence/$variant-source.sha256"
done
common="$work/$candidate"
for variant in "$baseline" "$candidate"; do
    src="$work/$variant"; stage="build-$variant"
    /usr/bin/time -v -o "$work/evidence/$variant-build.time" make -f "$common/benchmarks/fsi-time-opt/fsi_case.mk" SOURCE="$src" BUILD="$work/$variant" CASE_SOURCE=test_compliant_channel_fsi.cpp CASE_NAME=fsi fsi_case > "$work/evidence/$variant-build.stdout" 2> "$work/evidence/$variant-build.stderr"
    cp "$work/$variant/fsi" "$work/evidence/$variant-fsi"
    sha256sum "$work/$variant/fsi" > "$work/evidence/$variant-binary.sha256"
    ldd "$work/$variant/fsi" > "$work/evidence/$variant-ldd.txt"
done
src="$common"; stage=build-regressions
for spec in immersed:test_immersed_transient_flow.cpp adapter:test_moving_immersed_transient_flow_fsi_runtime.cpp volume:test_parallel_immersed_volume.cpp; do
    name=${spec%%:*}; source=${spec#*:}
    /usr/bin/time -v -o "$work/evidence/$name-build.time" make -f "$common/benchmarks/fsi-time-opt/fsi_case.mk" SOURCE="$src" BUILD="$work/$candidate" CASE_SOURCE="$source" CASE_NAME="$name" fsi_case > "$work/evidence/$name-build.stdout" 2> "$work/evidence/$name-build.stderr"
    cp "$work/$candidate/$name" "$work/evidence/$candidate-$name"
    sha256sum "$work/$candidate/$name" >> "$work/evidence/$candidate-binary.sha256"
done
failure_status=failed_numerical
for name in immersed adapter volume; do
    stage="V0-V1-$name"; sample="$work/evidence/$name"
    mkdir -p "$sample"; cd "$sample"
    IGA_PROFILE_DETAIL=1 /usr/bin/time -v -o timing.txt mpiexec -np 1 --map-by slot:PE="$OMP_NUM_THREADS" --bind-to core "$work/$candidate/$name" > stdout.txt 2> stderr.txt
done
# Same node, inputs, flags, affinity, detail mode and unchanged native gates.
for variant in "$baseline" "$candidate"; do
    stage="V2-$variant"; sample="$work/evidence/$variant-fsi"
    mkdir -p "$sample"; cd "$sample"
    IGA_PROFILE_DETAIL=1 /usr/bin/time -v -o timing.txt mpiexec -np 1 --map-by slot:PE="$OMP_NUM_THREADS" --bind-to core "$work/$variant/fsi" --reference-output reference > stdout.txt 2> stderr.txt
done
stage=verification; failure_status=inconclusive
python3 "$common/benchmarks/fsi-time-opt/verify_p1a.py" "$work/evidence" "$baseline" "$candidate"
stage=complete
