#!/usr/bin/env python3
"""P1-A gates: existing assertions plus exact accepted FSI reference and call removal."""
import hashlib
import json
import os
from pathlib import Path
import re
import statistics
import sys

root=Path(sys.argv[1]); before,after=sys.argv[2:4]
checks={}; profiles={}; assembly={}; references={}
case=json.loads((root/'case.json').read_text())
selected=case.get('regressions',['immersed','adapter','volume'])
if case.get('passed_paired_adapter_job'):
    evidence_store=Path(sys.argv[4]) if len(sys.argv)>4 else root.parent
    previous=evidence_store/f'job-{case["passed_paired_adapter_job"]}'
    current_test=Path(__file__).resolve().parents[2]/'solvers/cpu/tests/test_moving_immersed_transient_flow_fsi_runtime.cpp'
    checks['prior_adapter_exact_common_test']=hashlib.sha256(current_test.read_bytes()).hexdigest()==(previous/'workload-test.sha256').read_text().split()[0]
    for variant,label in ((before,'adapter-before'),(after,'adapter')):
        text=(previous/label/'stdout.txt').read_text()
        checks[variant+'_prior_adapter_zero_and_first_fault']='fsi_adapter_zero_initial_residual passed' in text and 'fsi_adapter_first_assembly_failure_retry passed' in text and 'channel patch_nodes=9' in text and 'converged=true' in text
        checks[variant+'_prior_adapter_exit_zero']='Exit status: 0' in (previous/label/'timing.txt').read_text()
        checks[variant+'_prior_numerical_source_same']=(previous/(variant+'-source.sha256')).read_text().split()[0]==(root/(variant+'-source.sha256')).read_text().split()[0]
        binary_hash=next(line.split()[0] for line in (previous/(variant+'-binary.sha256')).read_text().splitlines() if Path(line.split()[1]).name=='adapter')
        checks[variant+'_prior_saved_adapter_binary_sha']=hashlib.sha256((previous/(variant+'-adapter')).read_bytes()).hexdigest()==binary_hash
operator_records=[]
if case.get('w1_operator'):
    expected={0:'d67dfc4ff72363aa99b0b7228a3328fdfb162f687f78d74189eb8816a048cf8f',1:'982aa0be0e6fc99bfc949526a1b227a963bed900c01bb459e469c9366790c373'}
    pattern=re.compile(r'benchmark_sample trial=(-?\d+) state=(\d+) assembly_s=(\S+) sha256=(\w+) mallocs=(\S+)')
    for pair in range(3):
        for variant in (before,after):
            sample=root/f'W1-pair-{pair}-{variant}'; text=(sample/'stdout.txt').read_text()
            records=list(pattern.finditer(text)); prefix=f'W1-{pair}-{variant}'
            checks[prefix+'_seven_including_warmup']=len(records)==7
            for record in records:
                checks[prefix+'_'+record[1]+'_exact_archived_R_J']=record[4]==expected[int(record[2])]
                checks[prefix+'_'+record[1]+'_mallocs_zero']=float(record[5])==0.
            detail=[json.loads(line) for line in (sample/'assembly-detail.jsonl').read_text().splitlines() if json.loads(line)['kind']=='assembly']
            checks[prefix+'_team_and_generation']=len(detail)==7 and all(r['status']=='completed' and r['team_size']==int(os.environ['OMP_NUM_THREADS']) and r['state_generation']==r['residual_generation']==r['jacobian_generation'] for r in detail)
            profile=json.loads(text.split('hpc_profile ')[-1].splitlines()[0])
            checks[prefix+'_status_zero']=profile['status']==0
            seconds=[float(r[3]) for r in records if int(r[1])>=0]
            operator_records.append({'pair':pair,'variant':variant,'median_assembly_s':statistics.median(seconds),'samples_s':seconds,'profile':profile})
for name in selected:
    sample=root/name
    text=(sample/'stdout.txt').read_text()
    if name=='immersed':
        checks[name+'_profile_passed']=json.loads(text.split('hpc_profile ')[-1].splitlines()[0])['status']==0
    elif name=='adapter':
        checks['adapter_complete']='channel patch_nodes=9' in text and 'converged=true' in text
    else:
        checks['volume_six_exact_fault_retry_cases']=sum('exact_fields=1 rollback_retry=1 passed' in line for line in text.splitlines())==6
if 'adapter' in selected:
    adapter=(root/'adapter/stdout.txt').read_text()
    for marker in ('fsi_adapter_zero_initial_residual passed','fsi_adapter_first_assembly_failure_retry passed'):
        checks[marker]=marker in adapter
    if case.get('paired_regressions'):
        baseline_adapter=(root/'adapter-before/stdout.txt').read_text()
        checks['baseline_adapter_complete']='channel patch_nodes=9' in baseline_adapter and 'converged=true' in baseline_adapter
        for marker in ('fsi_adapter_zero_initial_residual passed','fsi_adapter_first_assembly_failure_retry passed'):
            checks['baseline_'+marker]=marker in baseline_adapter
if case.get('matched_complete_fsi') is False:
    status='passed' if checks and all(checks.values()) else 'failed_numerical'
    result={'schema_version':1,'status':status,'checks':checks,'phase':case['phase'],'job_id':os.environ['SLURM_JOB_ID'],'variant':after,'formal_physics_acceptance':False}
    (root/'correctness.json').write_text(json.dumps(result,indent=2)+'\n')
    (root/'summary.md').write_text(f'# {case["phase"]}: {status}\n\nSelected complete original regression suites: {selected}.\nAll mathematical gates unchanged. Other P1-A gates remain separate.\n')
    print(json.dumps(result))
    sys.exit(0 if status=='passed' else 1)
for variant in (before,after):
    sample=root/(variant+'-fsi'); text=(sample/'stdout.txt').read_text()
    profiles[variant]=json.loads(text.split('hpc_profile ')[-1].splitlines()[0])
    checks[variant+'_complete_fsi']='compliant_channel_fsi converged=true' in text
    references[variant]=json.loads((sample/'reference/manifest.json').read_text())
    checks[variant+'_native_gates']=references[variant]['native_gates_passed'] is True
    for field in references[variant]['fields']:
        checks[variant+'_'+field['name']+'_saved_sha']=hashlib.sha256((sample/'reference'/field['file']).read_bytes()).hexdigest()==field['sha256']
    assembly[variant]=[json.loads(line) for line in (sample/'assembly-detail.jsonl').read_text().splitlines() if json.loads(line)['kind']=='assembly']
checks['accepted_reference_bitwise_equal']=references[before]==references[after]
def completed(rows,reason):
    return [r for r in rows if r['status']=='completed' and r['reason']==reason]
removed=completed(assembly[before],'adapter-initial')
checks['baseline_preassembly_proven']=len(removed)>0
checks['candidate_no_preassembly']=not completed(assembly[after],'adapter-initial')
checks['one_removed_per_fluid_trial']=len(removed)==len(completed(assembly[before],'initial'))==len(completed(assembly[after],'initial'))
checks['exact_total_call_reduction']=len(assembly[before])-len(assembly[after])==len(removed)
# Removing an unused scratch rebuild must not alter Newton/line-search trajectories.
identity=['reason','status','step','coupling_iteration','newton_iteration','damping','active_dofs','geometry_epoch','input_identity','state_identity','history_identity','moving_map_identity','state_generation','inside_cells','cut_cells','volume_quadrature_points','surface_quadrature_points']
remaining=[r for r in assembly[before] if r not in removed]
checks['remaining_assembly_inputs_identical']=[{k:r.get(k) for k in identity} for r in remaining]==[{k:r.get(k) for k in identity} for r in assembly[after]]
def physics_line(variant):
    return next(x for x in (root/(variant+'-fsi')/'stdout.txt').read_text().splitlines() if x.startswith('compliant_channel_fsi converged=true'))
checks['coupling_conservation_summary_identical']=physics_line(before)==physics_line(after)
status='passed' if all(checks.values()) else 'failed_numerical'
result={'schema_version':1,'status':status,'checks':checks,'removed_calls':len(removed),'baseline_calls':len(assembly[before]),'candidate_calls':len(assembly[after]),'formal_physics_acceptance':False,'note':'Small complete FSI and V0/V1 gate only; W2 original Y and final full cycle still pending'}
timing={'profiles':profiles,'W1_operator_records':operator_records,'baseline_assembly_wall_s':sum(r['wall_s'] for r in assembly[before]),'candidate_assembly_wall_s':sum(r['wall_s'] for r in assembly[after]),'baseline_unused_preassembly_wall_s':sum(r['wall_s'] for r in removed),'paired_processes':1,'speedup_claim':'screening only, not final repeated end-to-end acceptance'}
manifest={'schema_version':1,'phase':'P1-A-V0-V1-V2','variants':[before,after],'job_id':os.environ.get('IGA_COMPUTE_JOB_ID',os.environ['SLURM_JOB_ID']),'verification_job_id':os.environ['SLURM_JOB_ID'],'verification_only':os.environ.get('IGA_VERIFICATION_ONLY')=='1','node':os.environ.get('IGA_COMPUTE_NODE',os.environ.get('SLURM_JOB_NODELIST')),'threads':int(os.environ['OMP_NUM_THREADS']),'batch':int(os.environ['IGA_ASSEMBLY_BATCH_SIZE']),'source_hashes':[before+'-source.sha256',after+'-source.sha256'],'assertions':True,'same_allocation':True,'all_compute_node_local':True}
for name,obj in [('correctness.json',result),('timing.json',timing),('manifest.json',manifest)]:
    (root/name).write_text(json.dumps(obj,indent=2)+'\n')
failed=[k for k,v in checks.items() if not v]
(root/'summary.md').write_text(f'# P1-A V0/V1/V2: {status}\n\nRemoved calls: {len(removed)} ({len(assembly[before])} → {len(assembly[after])}).\nExact accepted reference equality: {checks["accepted_reference_bitwise_equal"]}.\nFailed checks: {failed}.\nW2 and final full-cycle gates remain pending.\n')
print(json.dumps({'status':status,'failed_checks':failed}))
sys.exit(0 if status=='passed' else 1)
