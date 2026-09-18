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
work="$LOCAL/fsi-r7-final-$SLURM_JOB_ID"
out="$root/job-$SLURM_JOB_ID"
mkdir -p "$work/evidence/final" "$work/evidence/historical" "$out" "$work/$variant"
stage=staging
failure_status=failed_environment
simulation_executed=false
started=$(date -u +%FT%TZ)
finish() {
    result=$?
    trap - EXIT TERM USR1
    test "$result" -ne 124 || failure_status=timeout
    printf '{"job_id":"%s","stage":"%s","status":"%s","exit_code":%s,"simulation_executed":%s,"started":"%s","ended":"%s"}\n' "$SLURM_JOB_ID" "$stage" "$([ "$result" -eq 0 ] && echo passed || echo "$failure_status")" "$result" "$simulation_executed" "$started" "$(date -u +%FT%TZ)" > "$work/evidence/status.json"
    copied=0
    cp -au "$work/evidence/." "$out/" || copied=$?
    test "$copied" -eq 0 || { test "$result" -ne 0 || result=$copied; }
    exit "$result"
}
trap finish EXIT
trap 'failure_status=interrupted; exit 124' TERM USR1

read -r expected_variant historical < <(python3 -c 'import json,sys; c=json.load(open(sys.argv[1])); print(c["variant"],c["historical_reference"])' "$cases")
test "$variant" = "$expected_variant"
sha256sum -c "$root/$variant-source.sha256"
tar -xf "$root/$variant-source.tar" -C "$work/$variant"
sha256sum -c "$root/surface.sha256"
cp "$root/surface.txt" "$work/surface.txt"
cp "$root/$variant-source.sha256" "$0" "$work/evidence/"
cp "$cases" "$work/evidence/case.json"
for name in history.csv coupling.csv run.json heartbeat-verification.json heartbeat-vtk-verification.json; do cp "$historical/$name" "$work/evidence/historical/$name"; done
python3 - "$cases" "$work/evidence/historical" <<'PY'
import hashlib,json,pathlib,sys
case=json.load(open(sys.argv[1])); root=pathlib.Path(sys.argv[2])
for name,digest in case['historical_hashes'].items(): assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest
PY
src="$work/$variant"
stage=build
timeout -k 30s 900s /usr/bin/time -v -o "$work/evidence/build.time" make -f "$src/benchmarks/fsi-time-opt/heartbeat_case.mk" SOURCE="$src" WORKLOAD_SOURCE="$src" BUILD="$src" heartbeat > "$work/evidence/build.stdout" 2> "$work/evidence/build.stderr"
cp "$src/heartbeat" "$work/evidence/heartbeat"
sha256sum "$src/heartbeat" "$work/surface.txt" "$cases" > "$work/evidence/input-binary.sha256"
ldd "$src/heartbeat" > "$work/evidence/ldd.txt"
cp "$PETSC_DIR/$PETSC_ARCH/lib/petsc/conf/petscvariables" "$work/evidence/petscvariables.txt"
{
    date -u +%FT%TZ
    module list
    mpicxx --version
    lscpu
    scontrol show job "$SLURM_JOB_ID"
    grep Cpus_allowed_list /proc/self/status
    for name in OMP_NUM_THREADS OMP_THREAD_LIMIT OMP_DYNAMIC OMP_PROC_BIND OMP_PLACES IGA_ASSEMBLY_THREADS IGA_ASSEMBLY_BATCH_SIZE PETSC_DIR PETSC_ARCH; do printf '%s=%s\n' "$name" "${!name:-}"; done
} > "$work/evidence/environment.txt" 2>&1
python3 - "$work/evidence" "$root" "$variant" "$cases" "$work/surface.txt" "$src/heartbeat" "$0" "$started" <<'PY'
import hashlib,json,os,pathlib,sys
evidence=pathlib.Path(sys.argv[1]); root=pathlib.Path(sys.argv[2]); variant=sys.argv[3]
case=pathlib.Path(sys.argv[4]); surface=pathlib.Path(sys.argv[5]); binary=pathlib.Path(sys.argv[6])
runner=pathlib.Path(sys.argv[7]); started=sys.argv[8]
digest=lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
manifest={
    'schema_version':1,'phase':'R7-full-20-step-VTK-package','variant':variant,
    'job_id':os.environ['SLURM_JOB_ID'],'job_name':os.environ.get('SLURM_JOB_NAME'),
    'node':os.environ.get('SLURMD_NODENAME',os.environ.get('HOSTNAME')),
    'started_utc':started,'allocated_cpus':int(os.environ['SLURM_CPUS_PER_TASK']),
    'source_archive_sha256':(root/f'{variant}-source.sha256').read_text().split()[0],
    'binary_sha256':digest(binary),'surface_sha256':digest(surface),
    'case_sha256':digest(case),'runner_sha256':digest(runner),
    'petsc_dir':os.environ['PETSC_DIR'],'petsc_arch':os.environ['PETSC_ARCH'],
    'threads':int(os.environ['OMP_NUM_THREADS']),'batch':32,
}
(evidence/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
PY

stage=full-cycle
failure_status=failed_numerical
simulation_executed=true
sample="$work/evidence/final"
cd "$sample"
options=(-demo_visualization_only true -demo_steps 20 -demo_grid 16 -demo_dt .05
    -demo_pressure 5 -demo_pulse_amplitude 30 -demo_period 1 -demo_cut_depth 1
    -demo_wall_inertial_gamma 1 -demo_display_scale 100 -demo_native_reference false)
printf '%q ' "${options[@]}" > "$work/evidence/solver-options.txt"
IGA_ASSEMBLY_BATCH_SIZE=32 IGA_PROFILE_DETAIL=1 timeout -k 60s 15000s /usr/bin/time -v -o timing.txt mpiexec -np 1 --map-by slot:PE=8 --bind-to core "$src/heartbeat" "$work/surface.txt" "$sample/results" "${options[@]}" > stdout.txt 2> stderr.txt

stage=csv-verification
python3 "$src/examples/vascular_flow/bifurcation_fsi/verify_heartbeat.py" "$sample/results" > "$work/evidence/heartbeat-verification.stdout" 2> "$work/evidence/heartbeat-verification.stderr"
stage=vtk-verification
/jet/home/thsieh1/.conda/envs/research/bin/python "$src/examples/vascular_flow/bifurcation_fsi/verify_heartbeat_vtk.py" "$sample/results" > "$work/evidence/vtk-verification.stdout" 2> "$work/evidence/vtk-verification.stderr"
stage=package
failure_status=failed_environment
cd "$sample"
tar -czf "$work/evidence/heartbeat-paraview.tar.gz" results
gzip -t "$work/evidence/heartbeat-paraview.tar.gz"
tar -tzf "$work/evidence/heartbeat-paraview.tar.gz" > "$work/evidence/archive-files.txt"
sha256sum "$work/evidence/heartbeat-paraview.tar.gz" > "$work/evidence/archive.sha256"
stage=verification
failure_status=failed_numerical
python3 "$src/benchmarks/fsi-time-opt/verify_r7.py" "$work/evidence"
stage=complete
