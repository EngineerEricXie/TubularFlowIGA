#!/usr/bin/env python3
"""P1B gates: original full-assembly inputs, local readiness, exact physics."""
import hashlib,json,os,sys,tarfile
from pathlib import Path

def rows(path):
    return [json.loads(line) for line in path.read_text().splitlines()
            if json.loads(line)['kind'] in ('assembly','assembly-reuse')]

inputs=['step','coupling_iteration','active_dofs','geometry_epoch','input_identity',
        'state_identity','history_identity','moving_map_identity','state_generation',
        'residual_generation','jacobian_generation']

def prove(full):
    pairs=[]
    for previous,current in zip(full,full[1:]):
        if current['reason']=='next-Newton':
            assert previous['reason']=='line-search' and previous['status']==current['status']=='completed'
            assert current['newton_iteration']==previous['newton_iteration']+1
            assert all(previous[k]==current[k] for k in inputs)
            pairs.append({'step':current['step'],'coupling':current['coupling_iteration'],
                          'newton':current['newton_iteration'],'unused_wall_s':current['wall_s']})
    assert pairs
    return pairs

if sys.argv[1]=='--proof':
    old,out=map(Path,sys.argv[2:4]); path=old/'W2-P1A/assembly-detail.jsonl'
    pairs=prove(rows(path)); result={'status':'passed','source_job':46071831,
        'input_log_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'identical_pairs':pairs,'unused_next_Newton_wall_s':sum(p['unused_wall_s'] for p in pairs)}
    (out/'p1b-before-proof.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','identical_pairs':len(pairs)})); sys.exit(0)

if sys.argv[1]=='--prefix':
    root,out=map(Path,sys.argv[2:4]); case=json.loads(Path(sys.argv[4]).read_text()); common=Path(sys.argv[5])
    if not case.get('passed_prefix_job'): sys.exit(0)
    old=root/f'job-{case["passed_prefix_job"]}'
    assert json.loads((old/'status.json').read_text())=={'job_id':'46083059','stage':'targeted-regression','exit_code':124,'status':'failed'}
    assert (old/'P1B-source.sha256').read_text().split()[0]==(out/'P1B-source.sha256').read_text().split()[0]
    digest=(old/'binary.sha256').read_text().split()[0]
    assert hashlib.sha256((old/'binaries/accepted-full').read_bytes()).hexdigest()==digest
    name='solvers/cpu/tests/test_accepted_line_search_assembly.cpp'
    with tarfile.open(root/'P1B-source.tar') as archive: original=archive.extractfile(name).read().decode()
    current=(common/name).read_text()
    # Only a CLI selector in main changes; ALL case assertions and solver fixtures
    # are identical to those compiled into the saved completed prefix.
    additions='''    std::string selected;
    for(int i=1;i<argc;++i) if(std::string(argv[i])=="--reuse-case") {
        if(++i==argc) return 2;
        selected=argv[i];
    }
    if(!selected.empty() && selected!="zero" && selected!="one-step" && selected!="normal"
        && selected!="half" && selected!="nan-first" && selected!="failure"
        && selected!="nan-all" && selected!="throw-ready") return 2;
'''
    assert current.replace(additions,'',1).replace('            if(!selected.empty() && selected!=name) continue;\n','',1)==original
    text=(old/'accepted-full/stdout.txt').read_text()
    cases=case['passed_prefix_cases']
    assert cases==['zero','one-step','normal','half','nan-first','failure','nan-all']
    for name in cases:
        assert sum(line.startswith(f'accepted_full_pair {name} ') and 'bitwise=1 passed' in line for line in text.splitlines())==1
        for setting in (0,1):
            assert sum(line.startswith(f'accepted_full_case {name} reuse={setting} ') and 'exact_retry_operator=1 passed' in line for line in text.splitlines())==1
    (out/'passed-prefix-stdout.txt').write_text(text)
    (out/'passed-prefix.json').write_text(json.dumps({'status':'passed','scope':'seven completed assertion-bearing case pairs ONLY; old overall process timed out',
        'job_id':case['passed_prefix_job'],'cases':cases,'binary_sha256':digest,
        'stdout_sha256':hashlib.sha256(text.encode()).hexdigest(),'exact_test_selector_only_delta':True},indent=2)+'\n')
    print('seven_completed_prefix_pairs_source_binary_and_assertions_verified'); sys.exit(0)

root=Path(sys.argv[1]); before,after=sys.argv[2:4]; checks={}; profiles={}; raw={}; native={}
test=root/'accepted-full/stdout.txt'; text=test.read_text()
if (root/'passed-prefix.json').is_file():
    prefix=json.loads((root/'passed-prefix.json').read_text()); saved=(root/'passed-prefix-stdout.txt').read_text()
    checks['prefix_saved_stdout_sha']=hashlib.sha256(saved.encode()).hexdigest()==prefix['stdout_sha256']
    text=saved+'\n'+text
checks['eight_exact_branch_pairs']=sum(line.startswith('accepted_full_pair ') and 'bitwise=1 passed' in line for line in text.splitlines())==8
checks['targeted_regression_exit_zero']='Exit status: 0' in (test.parent/'timing.txt').read_text()
for variant in (before,after):
    sample=root/f'{variant}-fsi'; stdout=(sample/'stdout.txt').read_text()
    profiles[variant]=json.loads(stdout.split('hpc_profile ')[-1].splitlines()[0])
    checks[variant+'_full_FSI_passed']='compliant_channel_fsi converged=true' in stdout and 'Exit status: 0' in (sample/'timing.txt').read_text()
    manifest=json.loads((sample/'reference/manifest.json').read_text()); native[variant]=manifest
    checks[variant+'_native_gates']=manifest['native_gates_passed'] is True
    for field in manifest['fields']:
        checks[variant+'_'+field['name']+'_saved_sha']=hashlib.sha256((sample/'reference'/field['file']).read_bytes()).hexdigest()==field['sha256']
    raw[variant]=rows(sample/'assembly-detail.jsonl')
checks['all_native_fields_exact']=native[before]==native[after]
baseline=raw[before]; candidate=raw[after]; pairs=prove(baseline)
reuse=[r for r in candidate if r['kind']=='assembly-reuse']; expanded=[]; previous=None
for row in candidate:
    if row['kind']=='assembly-reuse':
        checks[f'reuse_{len(expanded)}_exact_accepted_full_inputs']=previous is not None and previous['reason']=='line-search' and row['newton_iteration']==previous['newton_iteration']+1 and all(previous[k]==row[k] for k in inputs)
        clone=dict(previous); clone.update(row); clone.update(kind='assembly',reason='next-Newton',damping=0)
        expanded.append(clone)
    else: expanded.append(row)
    previous=row
keys=inputs+['reason','status','newton_iteration','damping','inside_cells','cut_cells',
             'volume_quadrature_points','surface_quadrature_points']
checks['remaining_full_and_reused_inputs_exact']=[{k:r.get(k) for k in keys} for r in baseline]==[{k:r.get(k) for k in keys} for r in expanded]
checks['exact_real_full_assembly_reduction']=len(reuse)==len(pairs)>0 and len(baseline)-sum(r['kind']=='assembly' for r in candidate)==len(reuse)
def physical_line(v):
    return next(line for line in (root/f'{v}-fsi/stdout.txt').read_text().splitlines() if line.startswith('compliant_channel_fsi converged=true'))
checks['coupling_conservation_summary_exact']=physical_line(before)==physical_line(after)
status='passed' if all(checks.values()) else 'failed_numerical'
result={'status':status,'checks':checks,'job_id':os.environ['SLURM_JOB_ID'],
        'formal_physics_acceptance':False,'promoted':False,'W2':'pending',
        'before_elapsed_s':profiles[before]['elapsed_s'],'after_elapsed_s':profiles[after]['elapsed_s'],
        'removed_full_calls':len(reuse),'timing_scope':'ONE screening pair only'}
(root/'correctness.json').write_text(json.dumps(result,indent=2)+'\n')
failed=[k for k,v in checks.items() if not v]
(root/'summary.md').write_text(f'# P1-B targeted + small FSI: {status}\n\nRemoved full calls {len(reuse)}. Failed checks: {failed}. W2 pending; not promoted.\n')
print(json.dumps({'status':status,'failed':failed})); sys.exit(0 if status=='passed' else 1)
