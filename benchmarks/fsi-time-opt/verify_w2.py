#!/usr/bin/env python3
"""Exact original Y first-two-step regression, unchanged formal gates, no shrink gate."""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys

root=Path(sys.argv[1]); before,after=sys.argv[2:4]
case=json.loads((root/'case.json').read_text()); card=case.get('optimization_card','P1A')
verification_threads=int(os.environ['OMP_NUM_THREADS'])
evidence_threads=int(case.get('threads',verification_threads))
checks={}; histories={}; profiles={}; assembly={}; references={}
historical=list(csv.DictReader((root/'historical-history.csv').open()))[:2]
required={'fluid_velocity','fluid_pressure','membrane_displacement','membrane_velocity','surface_traction','surface_force','material_position','material_displacement','material_velocity','port_flow','port_mean_pressure','port_mean_normal_traction','port_area'}
physical_scales={'fluid_velocity':1.0,'fluid_pressure':100.0,'membrane_displacement':1e-4,'membrane_velocity':1e-3,'surface_traction':100.0,'surface_force':1e-2,'material_position':1e-1,'material_displacement':1e-4,'material_velocity':1e-3,'port_flow':1e-6,'port_mean_pressure':100.0,'port_mean_normal_traction':100.0,'port_area':1e-4}
frozen_rtol=1e-8
for variant in (before,after):
    sample=root/f'W2-{variant}'; result=sample/'results'
    histories[variant]=list(csv.DictReader((result/'history.csv').open()))
    if card=='P5' and variant==after:
        checks[variant+'_original_step_time_sequence']=[(r['step'],r['time_s']) for r in histories[variant]]==[(r['step'],r['time_s']) for r in historical]
    else:
        checks[variant+'_exact_original_first_two_history']=histories[variant]==historical
    checks[variant+'_two_accepted_steps']=len(histories[variant])==2
    run=json.loads((result/'run.json').read_text())
    checks[variant+'_original_parameters']=all(run.get(k)==v for k,v in {'status':'passed','conservation_failed_steps':0,'steps':2,'dt_s':.05,'grid_xy':16,'cut_depth':1,'pulse_amplitude_pa':30,'period_s':1,'inlet_pressure_pa':5,'density_kg_m3':1060,'viscosity_pa_s':.0035,'wall_inertial_gamma0':1}.items())
    for row in histories[variant]:
        step=int(row['step']); prefix=f'{variant}-{step}'
        checks[prefix+'_formal_gates']=int(row['conservation_passed'])==1 and float(row['mass_defect'])<.03 and float(row['wall_leakage'])<.03 and float(row['continuity'])<1e-8 and float(row['inlet_m3_s'])<0. and float(row['lower_m3_s'])>0. and float(row['upper_m3_s'])>0.
        reference=result/f'step-{step}'/'native'
        manifest=json.loads((reference/'manifest.json').read_text()); references[(variant,step)]=manifest
        checks[prefix+'_native_schema_gate']=manifest['kind']=='native_fsi_reference_state' and manifest['checkpoint'] is False and manifest['native_gates_passed'] is True and manifest['step']==step and manifest['time_s']==float(row['time_s'])
        checks[prefix+'_native_field_coverage']=required.issubset({f['name'] for f in manifest['fields']})
        for field in manifest['fields']:
            path=reference/field['file']; table=list(csv.reader(path.open()))
            checks[prefix+'_'+field['name']+'_hash']=hashlib.sha256(path.read_bytes()).hexdigest()==field['sha256']
            ids=[int(r[0]) for r in table[1:]]
            checks[prefix+'_'+field['name']+'_ids_dimensions_finite']=len(table)==field['rows']+1 and len(ids)==len(set(ids)) and ids==sorted(ids) and all(i>=0 for i in ids) and all(len(r)==field['columns']+1 for r in table) and all(math.isfinite(float(v)) for r in table[1:] for v in r[1:])
    text=(sample/'stdout.txt').read_text()
    profiles[variant]=json.loads(text.split('hpc_profile ')[-1].splitlines()[0])
    checks[variant+'_complete_exit_zero']='bifurcation_fsi completed steps=2' in text and 'Exit status: 0' in (sample/'timing.txt').read_text()
    assembly[variant]=[json.loads(line) for line in (sample/'assembly-detail.jsonl').read_text().splitlines() if json.loads(line)['kind']=='assembly']
    checks[variant+'_all_assembly_completed']=all(r['status']=='completed' for r in assembly[variant])
    expected_batch=case.get('batch_by_variant',{}).get(variant,case.get('batch',int(os.environ['IGA_ASSEMBLY_BATCH_SIZE'])))
    checks[variant+'_saved_batch_configuration']=int((sample/'batch-size.txt').read_text())==expected_batch
    checks[variant+'_assembly_configuration']=all(r['team_size']==evidence_threads and r['batch_capacity']==expected_batch for r in assembly[variant])
def csv_columns_within_tolerance(left,right,scales,exact,roundoff_scales={}):
    a=list(csv.DictReader(left.open()));b=list(csv.DictReader(right.open()))
    if len(a)!=len(b) or (a and b and a[0].keys()!=b[0].keys()): return False
    for key in exact:
        if [r[key] for r in a]!=[r[key] for r in b]: return False
    for key,scale in scales.items():
        av=[float(r[key]) for r in a];bv=[float(r[key]) for r in b]
        if not all(math.isfinite(x) for x in av+bv): return False
        an=math.sqrt(sum(x*x for x in av));dn=math.sqrt(sum((x-y)*(x-y) for x,y in zip(av,bv)))
        if dn>256*sys.float_info.epsilon*scale+frozen_rtol*max(an,scale): return False
    for key,scale in roundoff_scales.items():
        av=[float(r[key]) for r in a];bv=[float(r[key]) for r in b]
        if not all(math.isfinite(x) for x in av+bv): return False
        dn=math.sqrt(sum((x-y)*(x-y) for x,y in zip(av,bv)))
        if dn>256*sys.float_info.epsilon*scale*math.sqrt(max(1,len(av))): return False
    return True
def native_field_within_tolerance(step,name):
    left=root/f'W2-{before}/results'/f'step-{step}'/'native'
    right=root/f'W2-{after}/results'/f'step-{step}'/'native'
    lm=references[(before,step)];rm=references[(after,step)]
    le=next(x for x in lm['fields'] if x['name']==name);re=next(x for x in rm['fields'] if x['name']==name)
    la=list(csv.reader((left/le['file']).open()));ra=list(csv.reader((right/re['file']).open()))
    if not la or not ra or la[0]!=ra[0] or len(la)!=len(ra): return False
    if [row[0] for row in la[1:]]!=[row[0] for row in ra[1:]]: return False
    lv=[float(x) for row in la[1:] for x in row[1:]];rv=[float(x) for row in ra[1:] for x in row[1:]]
    if len(lv)!=len(rv) or not all(math.isfinite(x) for x in lv+rv): return False
    ln=math.sqrt(sum(x*x for x in lv));dn=math.sqrt(sum((x-y)*(x-y) for x,y in zip(lv,rv)));scale=physical_scales[name]
    return dn<=256*sys.float_info.epsilon*scale+frozen_rtol*max(ln,scale)
if card=='P5':
    history_scales={'time_s':1.0,'max_displacement_m':1e-4,'mass_defect':.03,'wall_leakage':.03,'inlet_m3_s':1e-6,'lower_m3_s':1e-6,'upper_m3_s':1e-6,'normalization_scale_m3_s':1e-6,'moving_mass_defect_m3_s':1e-6,'wall_relative_leakage_m3_s':1e-6,'backward_euler_volume_rate_m3_s':1e-6,'total_fluid_surface_outward_flow_m3_s':1e-6,'total_material_surface_outward_flow_m3_s':1e-6,'inlet_pressure_pa':100.0}
    checks['history_within_frozen_per_quantity_tolerance']=csv_columns_within_tolerance(root/f'W2-{before}/results/history.csv',root/f'W2-{after}/results/history.csv',history_scales,{'step','iterations','conservation_passed'},{'continuity':1.0})
    checks['coupling_within_frozen_per_quantity_tolerance']=csv_columns_within_tolerance(root/f'W2-{before}/results/coupling.csv',root/f'W2-{after}/results/coupling.csv',{'rms_m':1e-4,'threshold_m':1e-9},{'step','iteration'})
    for step in (1,2):
        for name in required: checks[f'native_step_{step}_{name}_within_frozen_tolerance']=native_field_within_tolerance(step,name)
else:
    checks['coupling_history_exact']=(root/f'W2-{before}/results/coupling.csv').read_bytes()==(root/f'W2-{after}/results/coupling.csv').read_bytes()
    checks['all_native_fields_exact']=all(references[(before,step)]==references[(after,step)] for step in (1,2))
keys=['reason','status','step','coupling_iteration','newton_iteration','damping','active_dofs','geometry_epoch','input_identity','state_identity','history_identity','moving_map_identity','state_generation','inside_cells','cut_cells','volume_quadrature_points','surface_quadrature_points']
candidate_full_count=len(assembly[after]); candidate_residual_count=0
candidate_cache_bytes=0; candidate_cache_hits=0; candidate_cache_build_misses=0
candidate_preconditioner_builds=0; candidate_preconditioner_reuse_attempts=0
candidate_preconditioner_reuse_accepts=0; candidate_preconditioner_rebuilds=0
if card=='P1B':
    removed=[r for r in assembly[before] if r['reason']=='next-Newton']
    proof_keys=['step','coupling_iteration','active_dofs','geometry_epoch','input_identity','state_identity','history_identity','moving_map_identity','state_generation','residual_generation','jacobian_generation']
    for index,(previous,current) in enumerate(zip(assembly[before],assembly[before][1:])):
        if current['reason']=='next-Newton':
            checks[f'before_accepted_full_proof_{index}']=previous['reason']=='line-search' and previous['status']==current['status']=='completed' and current['newton_iteration']==previous['newton_iteration']+1 and all(previous[k]==current[k] for k in proof_keys)
    raw=[json.loads(line) for line in (root/f'W2-{after}/assembly-detail.jsonl').read_text().splitlines() if json.loads(line)['kind'] in ('assembly','assembly-reuse')]
    expanded=[]; previous=None; reused=0
    for row in raw:
        if row['kind']=='assembly-reuse':
            reused+=1
            checks[f'candidate_reuse_inputs_{reused}']=previous is not None and previous['reason']=='line-search' and row['status']==previous['status']=='completed' and row['newton_iteration']==previous['newton_iteration']+1 and all(previous[k]==row[k] for k in proof_keys)
            clone=dict(previous); clone.update(row); clone.update(kind='assembly',reason='next-Newton',damping=0)
            expanded.append(clone)
        else: expanded.append(row)
        previous=row
    checks['remaining_full_and_reused_inputs_exact']=[{k:r.get(k) for k in keys} for r in assembly[before]]==[{k:r.get(k) for k in keys} for r in expanded]
    checks['exact_accepted_full_call_reduction']=len(removed)==reused>0 and len(assembly[before])-len(assembly[after])==reused and not any(r['reason']=='next-Newton' for r in assembly[after])
elif card=='P2':
    removed=[]
    def all_calls(variant):
        return [json.loads(line) for line in (root/f'W2-{variant}/assembly-detail.jsonl').read_text().splitlines() if json.loads(line)['kind'] in ('assembly','assembly-reuse')]
    sequence_keys=['kind','reason','status','step','coupling_iteration','newton_iteration','damping','active_dofs','geometry_epoch','input_identity','state_identity','history_identity','moving_map_identity','state_generation','residual_generation','jacobian_generation','inside_cells','cut_cells','volume_quadrature_points','surface_quadrature_points']
    checks['exact_full_and_reuse_sequence_inputs']=[{k:r.get(k) for k in sequence_keys} for r in all_calls(before)]==[{k:r.get(k) for k in sequence_keys} for r in all_calls(after)]
    checks['batch_only_no_call_count_change']=len(all_calls(before))==len(all_calls(after))
elif card=='P3':
    def all_calls(variant):
        return [json.loads(line) for line in (root/f'W2-{variant}/assembly-detail.jsonl').read_text().splitlines() if json.loads(line)['kind'] in ('assembly','assembly-reuse')]
    baseline_raw=all_calls(before); candidate_raw=all_calls(after)
    expanded=[]; previous=None
    for row in baseline_raw:
        if row['kind']=='assembly-reuse':
            checks[f'baseline_reuse_{len(expanded)}_follows_line_search']=(previous is not None and previous['kind']=='assembly'
                and previous['reason']=='line-search' and previous['status']==row['status']=='completed'
                and row['newton_iteration']==previous['newton_iteration']+1)
            clone=dict(previous); clone.update(row); clone.update(kind='assembly',reason='next-Newton',damping=0)
            expanded.append(clone)
        else: expanded.append(row)
        previous=row
    candidate_full=[r for r in candidate_raw if r.get('assembly_request')=='ResidualAndJacobian']
    candidate_residual=[r for r in candidate_raw if r.get('assembly_request')=='ResidualOnly']
    candidate_full_count=len(candidate_full); candidate_residual_count=len(candidate_residual)
    removed=[None]*(len(assembly[before])-candidate_full_count)
    checks['candidate_contains_no_assembly_reuse']=not any(r['kind']=='assembly-reuse' for r in candidate_raw)
    checks['candidate_residual_only_is_line_search']=bool(candidate_residual) and all(r['reason']=='line-search' for r in candidate_residual)
    checks['candidate_full_is_not_line_search']=bool(candidate_full) and all(r['reason']!='line-search' for r in candidate_full)
    checks['candidate_residual_generation_current']=all(r['residual_generation']==r['state_generation'] for r in candidate_raw)
    checks['candidate_line_search_jacobian_stale']=all(r['jacobian_generation']!=r['state_generation'] for r in candidate_residual)
    checks['candidate_full_jacobian_current']=all(r['jacobian_generation']==r['state_generation'] for r in candidate_full)
    checks['line_search_count_preserved']=sum(r['reason']=='line-search' for r in assembly[before])==candidate_residual_count
    checks['next_newton_rebuild_matches_baseline_reuse']=sum(r['kind']=='assembly-reuse' for r in baseline_raw)==sum(r['reason']=='next-Newton' for r in candidate_raw)
    checks['logical_work_event_count_preserved']=len(expanded)==len(candidate_raw)
    checks['real_full_assembly_count_reduced']=0<candidate_full_count<len(assembly[before])
    p3_keys=keys+['residual_generation']
    checks['exact_state_aligned_assembly_sequence']=[{k:r.get(k) for k in p3_keys} for r in expanded]==[{k:r.get(k) for k in p3_keys} for r in candidate_raw]
    checks['exact_current_j_generation_on_full_events']=all(b.get('jacobian_generation')==c.get('jacobian_generation') for b,c in zip(expanded,candidate_raw) if c['reason']!='line-search')
    residual_payload=[r.get('maximum_resident_result_payload_bytes') for r in candidate_residual]
    full_payload=[r.get('maximum_resident_result_payload_bytes') for r in candidate_full]
    checks['residual_payload_recorded_and_smaller']=(all(v is not None for v in residual_payload+full_payload) and max(residual_payload)<min(full_payload))
elif card=='P4':
    removed=[]
    baseline_raw=assembly[before]; candidate_raw=assembly[after]
    candidate_full=[r for r in candidate_raw if r.get('assembly_request')=='ResidualAndJacobian']
    candidate_residual=[r for r in candidate_raw if r.get('assembly_request')=='ResidualOnly']
    candidate_full_count=len(candidate_full); candidate_residual_count=len(candidate_residual)
    sequence_keys=keys+['residual_generation','jacobian_generation','assembly_request']
    checks['exact_state_aligned_assembly_sequence']=[{k:r.get(k) for k in sequence_keys} for r in baseline_raw]==[{k:r.get(k) for k in sequence_keys} for r in candidate_raw]
    checks['cache_only_no_call_count_change']=len(baseline_raw)==len(candidate_raw)>0
    checks['baseline_has_no_volume_basis_cache_fields']=all('volume_basis_cache_enabled' not in r for r in baseline_raw)
    cache_keys={r.get('volume_basis_cache_key') for r in candidate_raw}
    cache_bytes={r.get('volume_basis_cache_bytes') for r in candidate_raw}
    build_misses_by_key={}
    for row in candidate_raw:
        key=row.get('volume_basis_cache_key'); misses=row.get('volume_basis_cache_build_misses')
        if key in build_misses_by_key: checks[f'candidate_cache_build_misses_stable_{len(build_misses_by_key)}']=build_misses_by_key[key]==misses
        else: build_misses_by_key[key]=misses
    checks['candidate_cache_enabled']=all(r.get('volume_basis_cache_enabled')==1 for r in candidate_raw)
    checks['candidate_cache_key_epoch_bound']=len(cache_keys)==len({r['geometry_epoch'] for r in candidate_raw}) and all(isinstance(k,str) and len(k)==64 for k in cache_keys)
    checks['candidate_cache_build_misses_match_points']=all(r.get('volume_basis_cache_build_misses')==r.get('volume_quadrature_points') for r in candidate_raw)
    checks['candidate_cache_hits_match_points']=all(r.get('volume_basis_cache_hits')==r.get('volume_quadrature_points') for r in candidate_raw)
    checks['candidate_cache_runtime_misses_zero']=all(r.get('volume_basis_cache_misses')==0 for r in candidate_raw)
    checks['candidate_cache_bytes_positive']=all(isinstance(v,int) and v>0 for v in cache_bytes)
    candidate_cache_bytes=max(cache_bytes) if cache_bytes else 0
    candidate_cache_hits=sum(r.get('volume_basis_cache_hits',0) for r in candidate_raw)
    candidate_cache_build_misses=sum(build_misses_by_key.values())
    def maximum_rss_kib(sample):
        for line in (sample/'timing.txt').read_text().splitlines():
            if 'Maximum resident set size (kbytes):' in line: return int(line.rsplit(':',1)[1])
        raise AssertionError('missing maximum RSS')
    checks['candidate_cache_rss_below_fraction']=maximum_rss_kib(root/f'W2-{after}') < case.get('memory_mb',15200)*1024*case.get('maximum_rss_fraction',.7)
elif card=='P5':
    removed=[]
    baseline_raw=assembly[before];candidate_raw=assembly[after]
    candidate_full=[r for r in candidate_raw if r.get('assembly_request')=='ResidualAndJacobian']
    candidate_residual=[r for r in candidate_raw if r.get('assembly_request')=='ResidualOnly']
    candidate_full_count=len(candidate_full);candidate_residual_count=len(candidate_residual)
    logical_keys=['reason','status','step','coupling_iteration','newton_iteration','damping','active_dofs','inside_cells','cut_cells','volume_quadrature_points','surface_quadrature_points','assembly_request']
    checks['same_logical_assembly_path']=[{k:r.get(k) for k in logical_keys} for r in baseline_raw]==[{k:r.get(k) for k in logical_keys} for r in candidate_raw]
    checks['solver_path_identity_is_distinct']=set(r['input_identity'] for r in baseline_raw).isdisjoint(set(r['input_identity'] for r in candidate_raw))
    raw=[json.loads(line) for line in (root/f'W2-{after}/assembly-detail.jsonl').read_text().splitlines()]
    linear=[row for row in raw if row['kind']=='linear-solve-attempt']
    candidate_preconditioner_builds=sum(row.get('reuse_preconditioner')==0 for row in linear)
    candidate_preconditioner_reuse_attempts=sum(row.get('reuse_preconditioner')==1 for row in linear)
    candidate_preconditioner_rebuilds=sum(row['kind']=='linear-solve-attempt' and row.get('reuse_preconditioner')==1 and i+1<len(raw) and raw[i+1]['kind']=='linear-solve-attempt' and raw[i+1].get('reuse_preconditioner')==0 for i,row in enumerate(raw))
    candidate_preconditioner_reuse_accepts=candidate_preconditioner_reuse_attempts-candidate_preconditioner_rebuilds
    baseline_linear_systems=profiles[before]['phases']['solver_setup']['calls']
    checks['candidate_linear_attempt_records_complete']=bool(linear) and all(row['status']=='completed' for row in linear)
    checks['candidate_true_linear_residual_gate']=all(row['ksp_reason']>0 and math.isfinite(row['true_linear_relative_residual']) and row['true_linear_relative_residual']<=case.get('reused_preconditioner_true_linear_relative_tolerance',1e-10) for row in linear)
    checks['candidate_reuse_visible_and_fewer_factor_builds']=candidate_preconditioner_reuse_attempts>0 and candidate_preconditioner_reuse_accepts>0 and candidate_preconditioner_builds<baseline_linear_systems
    checks['candidate_linear_system_accounting']=candidate_preconditioner_builds+candidate_preconditioner_reuse_accepts==baseline_linear_systems
    checks['candidate_setup_time_reduced']=profiles[after]['phases']['solver_setup']['exclusive_s']<profiles[before]['phases']['solver_setup']['exclusive_s']
else:
    assert card=='P1A'
    removed=[r for r in assembly[before] if r['reason']=='adapter-initial']
    remaining=[r for r in assembly[before] if r['reason']!='adapter-initial']
    checks['remaining_assembly_inputs_exact']=[{k:r.get(k) for k in keys} for r in remaining]==[{k:r.get(k) for k in keys} for r in assembly[after]]
    checks['exact_one_preassembly_removed_per_trial']=len(removed)>0 and len(removed)==sum(r['reason']=='initial' for r in assembly[before])==sum(r['reason']=='initial' for r in assembly[after]) and len(assembly[before])-len(assembly[after])==len(removed) and not any(r['reason']=='adapter-initial' for r in assembly[after])
elapsed_before=profiles[before]['elapsed_s']; elapsed_after=profiles[after]['elapsed_s']
checks['W2_no_timing_regression']=elapsed_after<=elapsed_before*1.03
numerical=all(v for k,v in checks.items() if k!='W2_no_timing_regression')
status='passed' if numerical and checks['W2_no_timing_regression'] else 'inconclusive' if numerical else 'failed_numerical'
correctness={'schema_version':1,'status':status,'checks':checks,'formal_physics_acceptance':False,'scope':'original Y two-step only; original full-cycle conservation failures not repaired'}
timing={'schema_version':1,'profiles':profiles,'before_elapsed_s':elapsed_before,'after_elapsed_s':elapsed_after,'elapsed_reduction_fraction':1-elapsed_after/elapsed_before,'baseline_assembly_calls':len(assembly[before]),'candidate_assembly_calls':len(assembly[after]),'candidate_full_assembly_calls':candidate_full_count,'candidate_residual_only_calls':candidate_residual_count,'removed_full_calls':len(removed),'candidate_cache_bytes':candidate_cache_bytes,'candidate_cache_hits':candidate_cache_hits,'candidate_cache_build_misses':candidate_cache_build_misses,'candidate_preconditioner_builds':candidate_preconditioner_builds,'candidate_preconditioner_reuse_attempts':candidate_preconditioner_reuse_attempts,'candidate_preconditioner_reuse_accepts':candidate_preconditioner_reuse_accepts,'candidate_preconditioner_rebuilds':candidate_preconditioner_rebuilds,'optimization_card':card,'paired_processes':1,'repeat_gate':'screening pair; winning final candidate requires independent repeated pairs'}
manifest={'schema_version':1,'phase':card+'-W2','job_id':os.environ['SLURM_JOB_ID'],'variants':[before,after],'node':os.environ['SLURM_JOB_NODELIST'],'threads':evidence_threads,'verification_threads':verification_threads,'batch':int(os.environ['IGA_ASSEMBLY_BATCH_SIZE']),'binary_build_job':(root/'binary-build-job.txt').read_text().strip(),'same_allocation':True,'native_reference':True,'VTU_restart':False,'full_cycle_shrink_gate_applied':False}
for name,data in [('correctness.json',correctness),('timing.json',timing),('manifest.json',manifest)]:
    (root/name).write_text(json.dumps(data,indent=2)+'\n')
failed=[k for k,v in checks.items() if not v]
(root/'summary.md').write_text(f'# {card} original Y W2: {status}\n\nBefore {elapsed_before:.3f} s; after {elapsed_after:.3f} s; reduction {1-elapsed_after/elapsed_before:.3%}.\nFull assembly calls {len(assembly[before])} → {candidate_full_count}; residual-only calls {candidate_residual_count}; removed {len(removed)}.\nFailed checks: {failed}.\nFormal acceptance remains false for the full original cycle.\n')
print(json.dumps({'status':status,'failed_checks':failed}))
sys.exit(0 if status=='passed' else 1)
