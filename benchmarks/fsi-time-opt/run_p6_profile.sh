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
work="$LOCAL/fsi-p6-profile-$SLURM_JOB_ID"
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

read -r expected_variant < <(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["variant"])' "$cases")
test "$variant" = "$expected_variant"
sha256sum -c "$root/$variant-source.sha256"
tar -xf "$root/$variant-source.tar" -C "$work/$variant"
cp "$root/$variant-source.sha256" "$cases" "$0" "$work/evidence/"
src="$work/$variant"
stage=build
timeout -k 30s 900s /usr/bin/time -v -o "$work/evidence/build.time" make -f "$src/benchmarks/fsi-time-opt/heartbeat_case.mk" SOURCE="$src" WORKLOAD_SOURCE="$src" BUILD="$src" heartbeat > "$work/evidence/build.stdout" 2> "$work/evidence/build.stderr"
cp "$src/heartbeat" "$work/evidence/heartbeat"
sha256sum "$src/heartbeat" > "$work/evidence/binary.sha256"
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
    -demo_native_reference false)
printf '%q ' "${options[@]}" > ../solver-options.txt
IGA_ASSEMBLY_BATCH_SIZE=32 IGA_PROFILE_DETAIL=1 timeout -k 30s 2500s /usr/bin/time -v -o timing.txt mpiexec -np 1 --map-by slot:PE="$OMP_NUM_THREADS" --bind-to core "$src/heartbeat" "$root/surface.txt" results "${options[@]}" > stdout.txt 2> stderr.txt

stage=verification
python3 - "$work/evidence" <<'PY'
import json,os,pathlib,re,sys
r=pathlib.Path(sys.argv[1]); sample=r/'profile'
rows=[json.loads(line) for line in (sample/'assembly-detail.jsonl').read_text().splitlines()]
runtime=[x for x in rows if x['kind']=='runtime_build']
preallocation=[x for x in rows if x['kind']=='jacobian_preallocation']
moving=[x for x in rows if x['kind']=='moving_begin']
geometry=[x for x in rows if x['kind']=='geometry']
required={'validation_catalog_preflight','volume_basis_cache','active_layout','gauge_weights','petsc_objects_preallocation','solver_options','initial_state'}
text=(sample/'stdout.txt').read_text(); timing=(sample/'timing.txt').read_text()
profiles=re.findall(r'hpc_profile (\{.*\})',text)
checks={
    'simulation_completed':'bifurcation_fsi completed steps=1' in text and 'Exit status: 0' in timing,
    'runtime_records_complete':len(runtime)==6 and all(x['status']=='completed' for x in runtime),
    'runtime_stage_coverage':all(required==set(x['timers_s']) for x in runtime),
    'preallocation_records_complete':len(preallocation)==len(runtime) and all(x['status']=='completed' for x in preallocation),
    'moving_records_complete':len(moving)==5 and all(x['status']=='completed' for x in moving),
    'geometry_records_complete':len(geometry)==6 and all(x['status']=='completed' for x in geometry),
    'active_layout_pattern_observed':bool(runtime) and all(len(x['active_layout_pattern'])==64 for x in runtime),
    'geometry_bound_layout_identity_observed':bool(runtime) and all(len(x['layout_identity'])==64 for x in runtime),
    'phase_profile_present':len(profiles)==1
}
totals={key:sum(x['timers_s'][key] for x in runtime) for key in sorted(required)}
result={
    'schema_version':1,
    'status':'passed' if all(checks.values()) else 'failed',
    'checks':checks,
    'runtime_build_count':len(runtime),
    'moving_begin_count':len(moving),
    'geometry_build_count':len(geometry),
    'runtime_build_wall_s':sum(x['wall_s'] for x in runtime),
    'runtime_stage_totals_s':totals,
    'jacobian_preallocation_wall_s':sum(x['wall_s'] for x in preallocation),
    'active_layout_pattern_count':len({x['active_layout_pattern'] for x in runtime}),
    'layout_identity_count':len({x['layout_identity'] for x in runtime}),
    'active_nodes':sorted({int(x['active_nodes']) for x in runtime}),
    'rows':sorted({int(x['rows']) for x in runtime}),
    'moving_stage_totals_s':{key:sum(x['timers_s'].get(key,0.0) for x in moving) for key in sorted({k for x in moving for k in x['timers_s']})},
    'geometry_stage_totals_s':{key:sum(x['timers_s'].get(key,0.0) for x in geometry) for key in sorted({k for x in geometry for k in x['timers_s']})},
    'phase_profile':None if not profiles else json.loads(profiles[-1]),
    'formal_physics_acceptance':False,
    'scope':'P6 profiling only; no reuse candidate promoted'
}
(r/'profile.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
(r/'correctness.json').write_text(json.dumps({'status':result['status'],'checks':checks,'job_id':os.environ['SLURM_JOB_ID'],'formal_physics_acceptance':False},indent=2)+'\n')
(r/'summary.md').write_text(f'# P6 geometry/runtime profile: {result["status"]}\n\nRuntime records: {len(runtime)}; active-layout patterns: {result["active_layout_pattern_count"]}; preallocation: {result["jacobian_preallocation_wall_s"]:.6f} s.\nFailed: {[k for k,v in checks.items() if not v]}.\n')
sys.exit(0 if result['status']=='passed' else 1)
PY
stage=complete
