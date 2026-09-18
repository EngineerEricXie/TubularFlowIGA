#!/bin/bash
set -euo pipefail
while test $# -gt 0; do case "$1" in --root) root=$2;shift 2;;--variants) variant=$2;shift 2;;--case-list) cases=$2;shift 2;;*)exit 2;;esac;done
: "${LOCAL:?}" "${SLURM_JOB_ID:?}" "${root:?}" "${variant:?}" "${cases:?}"
work="$LOCAL/fsi-p4-w1-$SLURM_JOB_ID";out="$root/job-$SLURM_JOB_ID";mkdir -p "$work/evidence" "$out" "$work/$variant";stage=staging
finish(){ result=$?;trap - EXIT TERM;printf '{"job_id":"%s","stage":"%s","status":"%s","exit_code":%s}\n' "$SLURM_JOB_ID" "$stage" "$([ "$result" -eq 0 ]&&echo passed||echo failed)" "$result">"$work/evidence/status.json";cp -au "$work/evidence/." "$out/"||{ test "$result" -ne 0||result=1;};exit "$result";};trap finish EXIT;trap 'exit 124' TERM
sha256sum -c "$root/$variant-source.sha256";tar -xf "$root/$variant-source.tar" -C "$work/$variant";cp "$root/$variant-source.sha256" "$cases" "$work/evidence/"
src="$work/$variant";tools="$src/benchmarks/fsi-time-opt";stage=build
timeout -k 30s 900s /usr/bin/time -v -o "$work/evidence/build.time" make -f "$tools/fsi_case.mk" SOURCE="$src" BUILD="$src" CASE_SOURCE="$src/solvers/cpu/tests/test_immersed_transient_flow.cpp" CASE_NAME=p4-w1 fsi_case >"$work/evidence/build.stdout" 2>"$work/evidence/build.stderr"
cp "$src/p4-w1" "$work/evidence/";sha256sum "$src/p4-w1">"$work/evidence/binary.sha256"
stage=w1;mkdir -p "$work/evidence/w1";cd "$work/evidence/w1"
orders=(off on on off off on)
for i in 0 1 2 3 4 5; do mode=${orders[$i]}; IGA_PROFILE_DETAIL=1 timeout -k 30s 2400s /usr/bin/time -v -o "$mode-$i.time" mpiexec -np 1 --map-by slot:PE="$OMP_NUM_THREADS" --bind-to core "$src/p4-w1" "--p4-cache-w1-$mode" >"$mode-$i.stdout" 2>"$mode-$i.stderr"; done
stage=verification
python3 - "$work/evidence" <<'PY'
import json,os,pathlib,statistics,sys
r=pathlib.Path(sys.argv[1]);w=r/'w1';records=[]
for i,mode in enumerate(['off','on','on','off','off','on']):
 lines=[x for x in (w/f'{mode}-{i}.stdout').read_text().splitlines() if x.startswith('volume_basis_cache_w1 ')]
 if len(lines)!=1: raise SystemExit(f'missing W1 marker {mode}-{i}')
 fields=dict(x.split('=',1) for x in lines[0].split()[1:]); rec={'mode':mode,'process':i,'construct_s':float(fields['construct_s']),'samples_s':[float(fields[f'sample{j}_s']) for j in range(3)],'residual_sha':fields['residual_sha'],'jacobian_action_sha':fields['jacobian_action_sha'],'points':int(fields['points']),'hits':int(fields['hits']),'misses':int(fields['misses']),'bytes':int(fields['bytes'])};rec['median_s']=statistics.median(rec['samples_s']);records.append(rec)
checks={}
checks['all_processes_exit_zero']=all('Exit status: 0' in (w/f"{x['mode']}-{x['process']}.time").read_text() for x in records)
checks['exact_residual_hash']=len({x['residual_sha'] for x in records})==1
checks['exact_jacobian_action_hash']=len({x['jacobian_action_sha'] for x in records})==1
checks['point_count_exact']=len({x['points'] for x in records})==1 and records[0]['points']>0
checks['enabled_hits_only']=all(x['hits']==4*x['points'] and x['misses']==0 and x['bytes']>0 for x in records if x['mode']=='on')
checks['disabled_misses_only']=all(x['hits']==0 and x['misses']==4*x['points'] and x['bytes']==0 for x in records if x['mode']=='off')
pairs=[(records[0],records[1]),(records[3],records[2]),(records[4],records[5])]
pair_reductions=[(off['median_s']-on['median_s'])/off['median_s'] for off,on in pairs]
checks['cached_faster_every_pair']=all(x>0 for x in pair_reductions)
checks['cached_median_reduction_at_least_5pct']=statistics.median(pair_reductions)>=.05
status='passed' if all(checks.values()) else 'failed'
timing={'schema_version':1,'records':records,'pair_reduction_fractions':pair_reductions,'median_pair_reduction_fraction':statistics.median(pair_reductions),'scope':'three independent process pairs; each median excludes one warmup assembly; construction reported separately'}
(r/'timing.json').write_text(json.dumps(timing,indent=2)+'\n');(r/'correctness.json').write_text(json.dumps({'status':status,'checks':checks,'job_id':os.environ['SLURM_JOB_ID'],'formal_physics_acceptance':False},indent=2)+'\n');(r/'summary.md').write_text(f"# P4 volume basis cache W1: {status}\n\nMedian paired assembly reduction: {100*timing['median_pair_reduction_fraction']:.3f}%. Failed: {[k for k,v in checks.items() if not v]}. W2 pending.\n");sys.exit(0 if status=='passed' else 1)
PY
stage=complete
