#!/bin/bash
# Original Y two-step matched regression, not a whole-cycle contraction check.
set -euo pipefail
while test $# -gt 0; do
    case "$1" in --root) root=$2; shift 2;; --variants) variants=$2; shift 2;; --case-list) cases=$2; shift 2;; *) exit 2;; esac
done
: "${LOCAL:?Slurm LOCAL required}" "${SLURM_JOB_ID:?Slurm required}"
IFS=, read -r baseline candidate <<< "$variants"
work="$LOCAL/fsi-w2-$SLURM_JOB_ID"; out="$root/job-$SLURM_JOB_ID"
mkdir -p "$work/evidence/binaries" "$out"
stage=staging; failure_status=failed_environment; started=$(date -u +%FT%TZ)
simulation_executed=false
save_stage() {
    printf '{"job_id":"%s","stage":"%s","status":"in_progress"}\n' "$SLURM_JOB_ID" "$stage" > "$work/evidence/stage.json"
    cp -au "$work/evidence/." "$out/"
}
finish() {
    result=$?; trap - EXIT TERM USR1
    test "$result" -ne 124 || failure_status=timeout
    printf '{"job_id":"%s","stage":"%s","status":"%s","exit_code":%s,"simulation_executed":%s,"started":"%s","ended":"%s"}\n' "$SLURM_JOB_ID" "$stage" "$([ "$result" -eq 0 ] && echo passed || echo "$failure_status")" "$result" "$simulation_executed" "$started" "$(date -u +%FT%TZ)" > "$work/evidence/status.json"
    copied=0; cp -au "$work/evidence/." "$out/" || copied=$?
    test "$copied" -eq 0 || { test "$result" -ne 0 || result=$copied; }
    exit "$result"
}
trap finish EXIT
trap 'failure_status=interrupted; exit 124' TERM USR1
cp "$cases" "$work/evidence/case.json"
cp "$0" "$work/evidence/runner.sh"
read -r build_job verifier_variant < <(python3 -c 'import json,sys; c=json.load(open(sys.argv[1])); print(c["binary_build_job"],c["verifier_variant"])' "$cases")
for variant in "$baseline" "$candidate" W2Tools "$verifier_variant"; do
    sha256sum -c "$root/$variant-source.sha256"
    mkdir -p "$work/$variant"; tar -xf "$root/$variant-source.tar" -C "$work/$variant"
    cp "$root/$variant-source.sha256" "$work/evidence/$variant-source.sha256"
done
sha256sum -c "$root/surface.sha256"
cp "$root/surface.txt" "$work/surface.txt"
cp "$root/historical-heartbeat-history.csv" "$work/evidence/historical-history.csv"
sha256sum "$work/surface.txt" "$work/evidence/historical-history.csv" "$cases" "$0" > "$work/evidence/input-runner.sha256"
python3 - "$root/job-$build_job" "$work/evidence" "$baseline" "$candidate" <<'PY'
import hashlib,json,pathlib,sys
old,new=map(pathlib.Path,sys.argv[1:3]); variants=sys.argv[3:5]
assert json.loads((old/'status.json').read_text())['status']=='passed'
for variant in (*variants,'W2Tools'):
    assert (old/f'{variant}-source.sha256').read_text().split()[0]==(new/f'{variant}-source.sha256').read_text().split()[0]
for variant in variants:
    checksum=next(line.split()[0] for line in (old/'binary.sha256').read_text().splitlines() if line.split()[1].endswith(f'/{variant}/heartbeat'))
    assert hashlib.sha256((old/'binaries'/f'{variant}-heartbeat').read_bytes()).hexdigest()==checksum
assert hashlib.sha256((new/'historical-history.csv').read_bytes()).hexdigest()=='ab5f0aa39714d49274ca3442d8a1d003399ac37b0e29bd6676d5ec86a94c5025'
print('W2_binary_sources_build_gate_and_historical_input verified')
PY
for variant in "$baseline" "$candidate"; do
    cp "$root/job-$build_job/binaries/$variant-heartbeat" "$work/$variant/heartbeat"
    cp "$work/$variant/heartbeat" "$work/evidence/binaries/$variant-heartbeat"
    sha256sum "$work/$variant/heartbeat" >> "$work/evidence/binary.sha256"
    ldd "$work/$variant/heartbeat" > "$work/evidence/$variant-ldd.txt"
done
printf '%s\n' "$build_job" > "$work/evidence/binary-build-job.txt"
cp "$PETSC_DIR/$PETSC_ARCH/lib/petsc/conf/petscvariables" "$work/evidence/petscvariables.txt"
{ date -u +%FT%TZ; module list; lscpu; scontrol show job "$SLURM_JOB_ID"; grep Cpus_allowed_list /proc/self/status;
  for name in OMP_NUM_THREADS OMP_THREAD_LIMIT OMP_DYNAMIC OMP_PROC_BIND OMP_PLACES IGA_ASSEMBLY_THREADS IGA_ASSEMBLY_BATCH_SIZE PETSC_DIR PETSC_ARCH; do printf '%s=%s\n' "$name" "${!name:-}"; done;
} > "$work/evidence/environment.txt" 2>&1
options=(-demo_steps 2 -demo_grid 16 -demo_dt .05 -demo_pressure 5 -demo_pulse_amplitude 30 -demo_period 1 -demo_cut_depth 1 -demo_wall_inertial_gamma 1 -demo_display_scale 100 -demo_native_reference true)
printf '%q ' "${options[@]}" > "$work/evidence/solver-options.txt"
if python3 -c 'import json,sys; sys.exit(0 if json.load(open(sys.argv[1])).get("preflight_only") else 1)' "$cases"; then
    stage=preflight_complete; exit 0
fi
for variant in "$baseline" "$candidate"; do
    stage="W2-$variant"; failure_status=failed_environment; save_stage
    sample="$work/evidence/W2-$variant"; mkdir -p "$sample"; cd "$sample"
    failure_status=failed_numerical
    simulation_executed=true
    batch=$(python3 -c 'import json,sys; c=json.load(open(sys.argv[1])); print(c.get("batch_by_variant",{}).get(sys.argv[2],c.get("batch",8)))' "$cases" "$variant")
    printf '%s\n' "$batch" > batch-size.txt
    IGA_ASSEMBLY_BATCH_SIZE=$batch IGA_PROFILE_DETAIL=1 timeout -k 30s 4000s /usr/bin/time -v -o timing.txt mpiexec -np 1 --map-by slot:PE="$OMP_NUM_THREADS" --bind-to core "$work/$variant/heartbeat" "$work/surface.txt" "$sample/results" "${options[@]}" > stdout.txt 2> stderr.txt
    stage="completed-W2-$variant"; save_stage
done
stage=verification; failure_status=inconclusive; save_stage
python3 "$work/$verifier_variant/benchmarks/fsi-time-opt/verify_w2.py" "$work/evidence" "$baseline" "$candidate"
stage=complete
