#!/usr/bin/env python3
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys

root=Path(sys.argv[1]); case=json.loads((root/'case.json').read_text())
threads=case['expected_threads']; checks={}; histories={}; references={}; assemblies={}; profiles={}; linear={}; rss={}
required={'fluid_velocity','fluid_pressure','membrane_displacement','membrane_velocity','surface_traction','surface_force','material_position','material_displacement','material_velocity','port_flow','port_mean_pressure','port_mean_normal_traction','port_area'}
for count in threads:
    sample=root/f'W2-threads-{count}'; result=sample/'results'; prefix=f't{count}'
    histories[count]=(result/'history.csv').read_bytes()
    rows=list(csv.DictReader((result/'history.csv').open()))
    checks[prefix+'_two_accepted_steps']=len(rows)==2 and [(r['step'],r['time_s']) for r in rows]==[('1','0.050000000000000003'),('2','0.10000000000000001')]
    checks[prefix+'_formal_gates']=all(int(r['conservation_passed'])==1 and float(r['mass_defect'])<.03 and float(r['wall_leakage'])<.03 and float(r['continuity'])<1e-8 and float(r['inlet_m3_s'])<0 and float(r['lower_m3_s'])>0 and float(r['upper_m3_s'])>0 for r in rows)
    run=json.loads((result/'run.json').read_text())
    checks[prefix+'_original_parameters']=all(run.get(k)==v for k,v in {'status':'passed','conservation_failed_steps':0,'steps':2,'dt_s':.05,'grid_xy':16,'cut_depth':1,'pulse_amplitude_pa':30,'period_s':1,'inlet_pressure_pa':5,'density_kg_m3':1060,'viscosity_pa_s':.0035,'wall_inertial_gamma0':1}.items())
    native={}
    for step in (1,2):
        directory=result/f'step-{step}'/'native'; manifest=json.loads((directory/'manifest.json').read_text())
        checks[f'{prefix}_step_{step}_native_schema']=manifest['kind']=='native_fsi_reference_state' and manifest['checkpoint'] is False and manifest['native_gates_passed'] is True and manifest['step']==step
        checks[f'{prefix}_step_{step}_native_coverage']=required.issubset({x['name'] for x in manifest['fields']})
        for field in manifest['fields']:
            path=directory/field['file']; table=list(csv.reader(path.open())); ids=[int(r[0]) for r in table[1:]]
            checks[f'{prefix}_step_{step}_{field["name"]}_integrity']=hashlib.sha256(path.read_bytes()).hexdigest()==field['sha256'] and len(table)==field['rows']+1 and ids==sorted(set(ids)) and all(len(r)==field['columns']+1 for r in table) and all(math.isfinite(float(v)) for r in table[1:] for v in r[1:])
            native[(step,field['name'])]=path.read_bytes()
    references[count]=native
    detail=[json.loads(line) for line in (sample/'assembly-detail.jsonl').read_text().splitlines()]
    assemblies[count]=[x for x in detail if x['kind']=='assembly']
    linear[count]=[x for x in detail if x['kind']=='linear-solve-attempt']
    checks[prefix+'_assembly_configuration']=len(assemblies[count])==60 and all(x['status']=='completed' and x['team_size']==count and x['batch_capacity']==case['batch'] for x in assemblies[count])
    full=[x for x in assemblies[count] if x.get('assembly_request')=='ResidualAndJacobian']; residual=[x for x in assemblies[count] if x.get('assembly_request')=='ResidualOnly']
    checks[prefix+'_residual_full_path']=len(full)==30 and len(residual)==30 and all(x['reason']=='line-search' for x in residual) and all(x['reason']!='line-search' for x in full)
    checks[prefix+'_cache_contract']=all(x.get('volume_basis_cache_enabled')==1 and x.get('volume_basis_cache_hits')==x.get('volume_quadrature_points') and x.get('volume_basis_cache_misses')==0 for x in assemblies[count])
    builds=sum(x.get('reuse_preconditioner')==0 for x in linear[count]); attempts=sum(x.get('reuse_preconditioner')==1 for x in linear[count])
    rebuilds=sum(x.get('reuse_preconditioner')==1 and i+1<len(detail) and detail[i+1]['kind']=='linear-solve-attempt' and detail[i+1].get('reuse_preconditioner')==0 for i,x in enumerate(detail))
    checks[prefix+'_preconditioner_reuse_accounting']=len(linear[count])==30 and builds==10 and attempts==20 and rebuilds==0
    checks[prefix+'_true_linear_residual']=all(x['status']=='completed' and x['ksp_reason']>0 and math.isfinite(x['true_linear_relative_residual']) and x['true_linear_relative_residual']<=case['reused_preconditioner_true_linear_relative_tolerance'] for x in linear[count])
    text=(sample/'stdout.txt').read_text(); found=re.findall(r'hpc_profile (\{.*\})',text)
    checks[prefix+'_complete_exit_profile']=bool(found) and 'bifurcation_fsi completed steps=2' in text and 'Exit status: 0' in (sample/'timing.txt').read_text()
    profiles[count]=json.loads(found[-1])
    rss[count]=next(int(line.rsplit(':',1)[1]) for line in (sample/'timing.txt').read_text().splitlines() if 'Maximum resident set size (kbytes):' in line)

first,second=threads
checks['history_bitwise_exact']=histories[first]==histories[second]
checks['coupling_bitwise_exact']=(root/f'W2-threads-{first}/results/coupling.csv').read_bytes()==(root/f'W2-threads-{second}/results/coupling.csv').read_bytes()
checks['all_native_fields_bitwise_exact']=references[first].keys()==references[second].keys() and all(references[first][k]==references[second][k] for k in references[first])
logical_keys=['reason','status','step','coupling_iteration','newton_iteration','damping','active_dofs','inside_cells','cut_cells','volume_quadrature_points','surface_quadrature_points','assembly_request','geometry_epoch','input_identity','state_identity','history_identity','moving_map_identity','state_generation','residual_generation','jacobian_generation']
checks['assembly_path_and_inputs_bitwise_exact']=[{k:x.get(k) for k in logical_keys} for x in assemblies[first]]==[{k:x.get(k) for k in logical_keys} for x in assemblies[second]]
linear_keys=['reason','status','coupling_iteration','ksp_iterations','ksp_reason','reuse_preconditioner','true_linear_relative_residual']
checks['linear_path_bitwise_exact']=[{k:x.get(k) for k in linear_keys} for x in linear[first]]==[{k:x.get(k) for k in linear_keys} for x in linear[second]]
elapsed={t:profiles[t]['elapsed_s'] for t in threads}; fastest=min(threads,key=lambda t:(elapsed[t],t)); limit=elapsed[fastest]*(1+case['closest_to_fastest_fraction'])
preferred=min(t for t in threads if elapsed[t]<=limit)
checks['selected_configuration_is_measured']=preferred in threads and fastest in threads
checks['rss_below_70_percent_allocation']=max(rss.values())<15200*1024*.7
status='passed' if all(checks.values()) else 'failed_numerical'
timing={'schema_version':1,'elapsed_s':{str(k):v for k,v in elapsed.items()},'phase_profiles':{str(k):v for k,v in profiles.items()},'max_rss_kib':{str(k):v for k,v in rss.items()},'fastest_thread':fastest,'preferred_within_5_percent_thread':preferred,'closest_to_fastest_fraction':case['closest_to_fastest_fraction'],'batch':case['batch'],'same_allocation':True,'paired_processes':1}
correctness={'schema_version':1,'status':status,'checks':checks,'job_id':os.environ['SLURM_JOB_ID'],'formal_physics_acceptance':False,'scope':'same accepted P5 code, original-Y two-step W2, top-two production thread configurations'}
(root/'selected-production-thread.txt').write_text(str(preferred)+'\n')
(root/'timing.json').write_text(json.dumps(timing,indent=2)+'\n')
(root/'correctness.json').write_text(json.dumps(correctness,indent=2)+'\n')
(root/'summary.md').write_text(f'# R6 W2 top-two thread scaling: {status}\n\nElapsed: {elapsed}; fastest: {fastest}; preferred within 5%: {preferred}. Failed: {[k for k,v in checks.items() if not v]}.\n')
print(json.dumps({'status':status,'elapsed_s':elapsed,'fastest':fastest,'preferred':preferred,'failed_checks':[k for k,v in checks.items() if not v]}))
sys.exit(0 if status=='passed' else 1)
