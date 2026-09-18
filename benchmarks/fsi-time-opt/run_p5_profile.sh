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
: "${LOCAL:?Slurm LOCAL required}" "${SLURM_JOB_ID:?Slurm required}"
work="$LOCAL/fsi-p5-profile-$SLURM_JOB_ID"
out="$root/job-$SLURM_JOB_ID"
mkdir -p "$work/evidence" "$out"
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

cp "$cases" "$work/evidence/case.json"
cp "$0" "$work/evidence/runner.sh"
read -r build_job expected_variant < <(python3 -c 'import json,sys; c=json.load(open(sys.argv[1])); print(c["binary_build_job"],c["variant"])' "$cases")
test "$variant" = "$expected_variant"
sha256sum -c "$root/$variant-source.sha256"
sha256sum -c "$root/surface.sha256"
cp "$root/$variant-source.sha256" "$work/evidence/$variant-source.sha256"
cp "$root/surface.txt" "$work/surface.txt"
python3 - "$root/job-$build_job" "$work" "$variant" <<'PY'
import hashlib,json,pathlib,shutil,sys
old,work=map(pathlib.Path,sys.argv[1:3]); variant=sys.argv[3]
assert json.loads((old/'status.json').read_text())['status']=='passed'
checksum=next(line.split()[0] for line in (old/'binary.sha256').read_text().splitlines()
              if line.split()[1].endswith(f'/{variant}/heartbeat'))
source=old/'binaries'/f'{variant}-heartbeat'
assert hashlib.sha256(source.read_bytes()).hexdigest()==checksum
shutil.copy2(source,work/'heartbeat')
PY
sha256sum "$work/heartbeat" "$work/surface.txt" "$cases" "$0" > "$work/evidence/input.sha256"
ldd "$work/heartbeat" > "$work/evidence/ldd.txt"
cp "$PETSC_DIR/$PETSC_ARCH/lib/petsc/conf/petscvariables" "$work/evidence/petscvariables.txt"
{
    date -u +%FT%TZ
    module list
    lscpu
    scontrol show job "$SLURM_JOB_ID"
    grep Cpus_allowed_list /proc/self/status
    for name in OMP_NUM_THREADS OMP_THREAD_LIMIT OMP_DYNAMIC OMP_PROC_BIND OMP_PLACES IGA_ASSEMBLY_THREADS IGA_ASSEMBLY_BATCH_SIZE PETSC_DIR PETSC_ARCH; do
        printf '%s=%s\n' "$name" "${!name:-}"
    done
} > "$work/evidence/environment.txt" 2>&1

stage=profiling
failure_status=failed_numerical
simulation_executed=true
sample="$work/evidence/profile"
mkdir -p "$sample"
cd "$sample"
options=(-demo_steps 1 -demo_grid 16 -demo_dt .05 -demo_pressure 5
    -demo_pulse_amplitude 30 -demo_period 1 -demo_cut_depth 1
    -demo_wall_inertial_gamma 1 -demo_display_scale 100
    -demo_native_reference false -log_view :petsc-log.txt
    -immersed_transient_ksp_view)
printf '%q ' "${options[@]}" > ../solver-options.txt
IGA_ASSEMBLY_BATCH_SIZE=32 IGA_PROFILE_DETAIL=1 timeout -k 30s 2500s /usr/bin/time -v -o timing.txt mpiexec -np 1 --map-by slot:PE="$OMP_NUM_THREADS" --bind-to core "$work/heartbeat" "$work/surface.txt" results "${options[@]}" > stdout.txt 2> stderr.txt
stage=extracting
python3 - petsc-log.txt stdout.txt > ../profile.json <<'PY'
import json,re,sys
log=open(sys.argv[1],errors='replace').read()
stdout=open(sys.argv[2],errors='replace').read()
events={}
for name in ('PCSetUp','KSPSetUp','KSPSolve','MatLUFactorSym','MatLUFactorNum',
             'MatGetOrdering','MatMult','VecNorm'):
    match=re.search(r'^'+re.escape(name)+r'\s+.*$',log,re.M)
    events[name]=None if match is None else float(match.group(0).split()[3])
profile=re.findall(r'hpc_profile (\{.*\})',stdout)
print(json.dumps({
    'schema_version':1,
    'status':'passed',
    'petsc_version':'3.22.1',
    'events_time_s':events,
    'phase_profile':None if not profile else json.loads(profile[-1]),
    'ksp_view_present':'KSP Object:' in stdout,
    'pc_type_lu':'type: lu' in stdout.lower(),
    'factor_solver_type':next(iter(re.findall(r'factor solver type: ([^\n]+)',stdout,re.I)),None)
},indent=2,sort_keys=True))
PY
stage=complete
