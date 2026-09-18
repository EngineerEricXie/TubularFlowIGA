#!/bin/bash
set -euo pipefail
while test $# -gt 0; do case "$1" in --root) root=$2;shift 2;;--variants) variants=$2;shift 2;;--case-list) cases=$2;shift 2;;*)exit 2;;esac;done
: "${LOCAL:?}" "${SLURM_JOB_ID:?}" "${root:?}" "${variants:?}" "${cases:?}"
IFS=, read -r baseline candidate <<< "$variants"
work="$LOCAL/fsi-p5-w1-$SLURM_JOB_ID";out="$root/job-$SLURM_JOB_ID";mkdir -p "$work/evidence/binaries" "$out";stage=staging
finish(){ result=$?;trap - EXIT TERM;printf '{"job_id":"%s","stage":"%s","status":"%s","exit_code":%s,"simulation_executed":true}\n' "$SLURM_JOB_ID" "$stage" "$([ "$result" -eq 0 ]&&echo passed||echo failed)" "$result">"$work/evidence/status.json";cp -au "$work/evidence/." "$out/"||{ test "$result" -ne 0||result=1;};exit "$result";};trap finish EXIT;trap 'exit 124' TERM
cp "$cases" "$0" "$work/evidence/"
read -r build_job expected_baseline expected_candidate < <(python3 -c 'import json,sys;c=json.load(open(sys.argv[1]));print(c["binary_build_job"],c["baseline"],c["candidate"])' "$cases")
test "$baseline" = "$expected_baseline";test "$candidate" = "$expected_candidate"
for variant in "$baseline" "$candidate"; do sha256sum -c "$root/$variant-source.sha256";cp "$root/$variant-source.sha256" "$work/evidence/";done
sha256sum -c "$root/surface.sha256";cp "$root/surface.txt" "$work/surface.txt"
python3 - "$root/job-$build_job" "$work" "$baseline" "$candidate" <<'PY'
import hashlib,json,pathlib,shutil,sys
old,work=map(pathlib.Path,sys.argv[1:3]);variants=sys.argv[3:5]
assert json.loads((old/'status.json').read_text())['status']=='passed'
for variant in variants:
    assert (old/f'{variant}-source.sha256').read_text().split()[0]==(pathlib.Path(sys.argv[1]).parent/f'{variant}-source.sha256').read_text().split()[0]
    checksum=next(line.split()[0] for line in (old/'binary.sha256').read_text().splitlines() if line.split()[1].endswith(f'/{variant}/heartbeat'))
    source=old/'binaries'/f'{variant}-heartbeat';assert hashlib.sha256(source.read_bytes()).hexdigest()==checksum
    shutil.copy2(source,work/f'{variant}-heartbeat')
PY
for variant in "$baseline" "$candidate"; do cp "$work/$variant-heartbeat" "$work/evidence/binaries/";sha256sum "$work/$variant-heartbeat" >> "$work/evidence/binary.sha256";ldd "$work/$variant-heartbeat" > "$work/evidence/$variant-ldd.txt";done
sha256sum "$work/surface.txt" "$cases" "$0" > "$work/evidence/input-runner.sha256"
{
date -u +%FT%TZ
module list
lscpu
scontrol show job "$SLURM_JOB_ID"
grep Cpus_allowed_list /proc/self/status
for name in OMP_NUM_THREADS OMP_THREAD_LIMIT OMP_DYNAMIC OMP_PROC_BIND OMP_PLACES IGA_ASSEMBLY_THREADS IGA_ASSEMBLY_BATCH_SIZE PETSC_DIR PETSC_ARCH;do printf '%s=%s\n' "$name" "${!name:-}";done
} > "$work/evidence/environment.txt" 2>&1
options=(-demo_steps 1 -demo_grid 16 -demo_dt .05 -demo_pressure 5 -demo_pulse_amplitude 30 -demo_period 1 -demo_cut_depth 1 -demo_wall_inertial_gamma 1 -demo_display_scale 100 -demo_native_reference true)
printf '%q ' "${options[@]}" > "$work/evidence/solver-options.txt"
orders=("$baseline $candidate" "$candidate $baseline" "$baseline $candidate")
for pair in 1 2 3;do
    read -r first second <<< "${orders[$((pair-1))]}"
    for variant in "$first" "$second";do
        stage="pair-$pair-$variant";sample="$work/evidence/$stage";mkdir -p "$sample";cd "$sample"
        IGA_ASSEMBLY_BATCH_SIZE=32 IGA_PROFILE_DETAIL=1 timeout -k 30s 1800s /usr/bin/time -v -o timing.txt mpiexec -np 1 --map-by slot:PE="$OMP_NUM_THREADS" --bind-to core "$work/$variant-heartbeat" "$work/surface.txt" results "${options[@]}" > stdout.txt 2> stderr.txt
    done
done
stage=verification
python3 - "$work/evidence" "$baseline" "$candidate" <<'PY'
import csv,json,math,os,pathlib,re,statistics,sys
root=pathlib.Path(sys.argv[1]);baseline,candidate=sys.argv[2:4]
physical={'fluid_velocity','fluid_pressure','membrane_displacement','membrane_velocity','surface_traction','surface_force','material_position','material_displacement','material_velocity','port_flow','port_mean_pressure','port_mean_normal_traction','port_area'}
scales={'fluid_velocity':1.0,'fluid_pressure':100.0,'membrane_displacement':1e-4,'membrane_velocity':1e-3,'surface_traction':100.0,'surface_force':1e-2,'material_position':1e-1,'material_displacement':1e-4,'material_velocity':1e-3,'port_flow':1e-6,'port_mean_pressure':100.0,'port_mean_normal_traction':100.0,'port_area':1e-4}
rtol=1e-8;checks={};samples={};elapsed={};profiles={};reuse={}
def wall(path):
    text=path.read_text();match=re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)',text);assert match
    values=[float(x) for x in match.group(1).split(':')]
    return sum(value*60**power for power,value in enumerate(reversed(values)))
def table(path):
    rows=list(csv.reader(path.open()));return rows[0],[(int(row[0]),[float(x) for x in row[1:]]) for row in rows[1:]]
for pair in (1,2,3):
  for variant in (baseline,candidate):
    sample=root/f'pair-{pair}-{variant}';samples[(pair,variant)]=sample;elapsed[(pair,variant)]=wall(sample/'timing.txt')
    text=(sample/'stdout.txt').read_text();profiles[(pair,variant)]=json.loads(text.split('hpc_profile ')[-1].splitlines()[0])
    checks[f'pair_{pair}_{variant}_complete']='bifurcation_fsi completed steps=1' in text and 'Exit status: 0' in (sample/'timing.txt').read_text()
    history=list(csv.DictReader((sample/'results/history.csv').open()));checks[f'pair_{pair}_{variant}_one_formal_step']=len(history)==1 and int(history[0]['conservation_passed'])==1
    manifest=json.loads((sample/'results/step-1/native/manifest.json').read_text());fields={f['name']:f for f in manifest['fields']}
    checks[f'pair_{pair}_{variant}_native_gate']=manifest['native_gates_passed'] is True and physical.issubset(fields)
    raw=[json.loads(line) for line in (sample/'assembly-detail.jsonl').read_text().splitlines()]
    linear=[row for row in raw if row['kind']=='linear-solve-attempt']
    rebuilds=sum(row['kind']=='linear-solve-attempt' and row.get('reuse_preconditioner')==1 and i+1<len(raw) and raw[i+1]['kind']=='linear-solve-attempt' and raw[i+1].get('reuse_preconditioner')==0 for i,row in enumerate(raw))
    reuse[(pair,variant)]={'attempts':sum(row.get('reuse_preconditioner')==1 for row in linear),'fresh':sum(row.get('reuse_preconditioner')==0 for row in linear),'rebuilds':rebuilds,'accepts':sum(row.get('reuse_preconditioner')==1 for row in linear)-rebuilds}
    if variant==candidate:checks[f'pair_{pair}_candidate_true_residual']=all(row['ksp_reason']>0 and math.isfinite(row['true_linear_relative_residual']) and row['true_linear_relative_residual']<=1e-10 for row in linear)
for field in physical:
    reference=(samples[(1,baseline)]/'results/step-1/native'/json.loads((samples[(1,baseline)]/'results/step-1/native/manifest.json').read_text())['fields'][next(i for i,x in enumerate(json.loads((samples[(1,baseline)]/'results/step-1/native/manifest.json').read_text())['fields']) if x['name']==field)]['file']).read_bytes()
    for pair in (2,3):
        manifest=json.loads((samples[(pair,baseline)]/'results/step-1/native/manifest.json').read_text());entry=next(x for x in manifest['fields'] if x['name']==field)
        checks[f'baseline_repeat_{pair}_{field}_exact']=(samples[(pair,baseline)]/'results/step-1/native'/entry['file']).read_bytes()==reference
for pair in (1,2,3):
    bm=json.loads((samples[(pair,baseline)]/'results/step-1/native/manifest.json').read_text());cm=json.loads((samples[(pair,candidate)]/'results/step-1/native/manifest.json').read_text())
    for field in physical:
        be=next(x for x in bm['fields'] if x['name']==field);ce=next(x for x in cm['fields'] if x['name']==field)
        bh,br=table(samples[(pair,baseline)]/'results/step-1/native'/be['file']);ch,cr=table(samples[(pair,candidate)]/'results/step-1/native'/ce['file'])
        aligned=bh==ch and [x[0] for x in br]==[x[0] for x in cr] and all(len(a[1])==len(b[1]) for a,b in zip(br,cr))
        bv=[x for _,row in br for x in row];cv=[x for _,row in cr for x in row]
        bn=math.sqrt(sum(x*x for x in bv));dn=math.sqrt(sum((a-b)*(a-b) for a,b in zip(bv,cv)));limit=256*sys.float_info.epsilon*scales[field]+rtol*max(bn,scales[field])
        checks[f'pair_{pair}_{field}_within_frozen_tolerance']=aligned and math.isfinite(dn) and dn<=limit
    base_calls=profiles[(pair,baseline)]['phases']['solver_setup']['calls'];c=reuse[(pair,candidate)]
    checks[f'pair_{pair}_reuse_visible_and_fewer_builds']=c['attempts']>0 and c['accepts']>0 and c['fresh']<base_calls and c['fresh']+c['accepts']==base_calls
    checks[f'pair_{pair}_setup_reduced']=profiles[(pair,candidate)]['phases']['solver_setup']['exclusive_s']<profiles[(pair,baseline)]['phases']['solver_setup']['exclusive_s']
reductions=[(elapsed[(p,baseline)]-elapsed[(p,candidate)])/elapsed[(p,baseline)] for p in (1,2,3)]
checks['candidate_faster_every_pair']=all(x>0 for x in reductions)
checks['median_elapsed_reduction_at_least_5_percent']=statistics.median(reductions)>=.05
status='passed' if all(checks.values()) else 'failed'
timing={'schema_version':1,'baseline':baseline,'candidate':candidate,'elapsed_s':{f'pair-{p}':{'baseline':elapsed[(p,baseline)],'candidate':elapsed[(p,candidate)],'reduction_fraction':reductions[p-1]} for p in (1,2,3)},'median_reduction_fraction':statistics.median(reductions),'setup_s':{f'pair-{p}':{'baseline':profiles[(p,baseline)]['phases']['solver_setup']['exclusive_s'],'candidate':profiles[(p,candidate)]['phases']['solver_setup']['exclusive_s']} for p in (1,2,3)},'reuse':{f'pair-{p}':reuse[(p,candidate)] for p in (1,2,3)}}
(root/'timing.json').write_text(json.dumps(timing,indent=2)+'\n');(root/'correctness.json').write_text(json.dumps({'status':status,'checks':checks,'formal_physics_acceptance':False,'scope':'three independent original-Y one-step matched process pairs; frozen per-quantity rtol 1e-8'},indent=2)+'\n')
(root/'summary.md').write_text(f'# P5 matched W1: {status}\n\nPair reductions: {reductions}; median {statistics.median(reductions):.3%}. Failed: {[k for k,v in checks.items() if not v]}.\n')
sys.exit(0 if status=='passed' else 1)
PY
stage=complete
