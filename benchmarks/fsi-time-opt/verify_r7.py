#!/usr/bin/env python3
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys

root=Path(sys.argv[1]); sample=root/'final'; result=sample/'results'
case=json.loads((root/'case.json').read_text()); historical=root/'historical'
checks={}; rtol=case['frozen_per_quantity_rtol']; eps=sys.float_info.epsilon

def csv_within(left,right,scales,exact,roundoff={}):
    a=list(csv.DictReader(left.open())); b=list(csv.DictReader(right.open()))
    if len(a)!=len(b) or (a and b and a[0].keys()!=b[0].keys()): return False
    for key in exact:
        if [x[key] for x in a]!=[x[key] for x in b]: return False
    for key,scale in scales.items():
        av=[float(x[key]) for x in a]; bv=[float(x[key]) for x in b]
        if not all(math.isfinite(x) for x in av+bv): return False
        an=math.sqrt(sum(x*x for x in av)); dn=math.sqrt(sum((x-y)**2 for x,y in zip(av,bv)))
        if dn>256*eps*scale+rtol*max(an,scale): return False
    for key,scale in roundoff.items():
        av=[float(x[key]) for x in a]; bv=[float(x[key]) for x in b]
        if not all(math.isfinite(x) for x in av+bv): return False
        if math.sqrt(sum((x-y)**2 for x,y in zip(av,bv)))>256*eps*scale*math.sqrt(max(1,len(av))): return False
    return True

run=json.loads((result/'run.json').read_text()); old_run=json.loads((historical/'run.json').read_text())
expected={'status':'visualization_only','conservation_failed_steps':case['expected_historical_conservation_failed_steps'],'pulse_amplitude_pa':30,'period_s':1,'display_scale':100,'cut_depth':1,'steps':20,'dt_s':.05,'grid_xy':16,'inlet_pressure_pa':5,'density_kg_m3':1060,'viscosity_pa_s':.0035,'wall_gamma0':2,'wall_inertial_gamma0':1,'model':'small-displacement pretensioned membrane'}
checks['run_parameters_and_status']=all(run.get(k)==v for k,v in expected.items())
checks['historical_run_identity']=old_run==expected
history_scales={'time_s':1.0,'max_displacement_m':1e-4,'mass_defect':.03,'wall_leakage':.03,'inlet_m3_s':1e-6,'lower_m3_s':1e-6,'upper_m3_s':1e-6,'normalization_scale_m3_s':1e-6,'moving_mass_defect_m3_s':1e-6,'wall_relative_leakage_m3_s':1e-6,'backward_euler_volume_rate_m3_s':1e-6,'total_fluid_surface_outward_flow_m3_s':1e-6,'total_material_surface_outward_flow_m3_s':1e-6,'inlet_pressure_pa':100.0}
checks['full_history_within_frozen_per_quantity_tolerance']=csv_within(historical/'history.csv',result/'history.csv',history_scales,{'step','iterations','conservation_passed'},{'continuity':1.0})
checks['full_coupling_within_frozen_per_quantity_tolerance']=csv_within(historical/'coupling.csv',result/'coupling.csv',{'rms_m':1e-4,'threshold_m':1e-9},{'step','iteration'})
history=list(csv.DictReader((result/'history.csv').open()))
checks['twenty_finite_steps']=len(history)==20 and [int(x['step']) for x in history]==list(range(1,21)) and all(math.isfinite(float(v)) for x in history for k,v in x.items() if k not in ('step','iterations','conservation_passed'))
checks['conservation_failure_count_preserved']=sum(int(x['conservation_passed'])==0 for x in history)==case['expected_historical_conservation_failed_steps']

heartbeat=json.loads((result/'heartbeat-verification.json').read_text()); old_heartbeat=json.loads((historical/'heartbeat-verification.json').read_text())
checks['expansion_contraction_verifier']=heartbeat['status']=='visualization_verified' and heartbeat['formal_acceptance'] is False and heartbeat['frames']==20 and heartbeat['expansion_observed'] is True and heartbeat['contraction_observed'] is True and heartbeat['conservation_failed_steps']==case['expected_historical_conservation_failed_steps']
av=old_heartbeat['mean_normal_displacement_m']; bv=heartbeat['mean_normal_displacement_m']
checks['mean_normal_displacement_within_frozen_tolerance']=len(av)==len(bv)==20 and math.sqrt(sum((x-y)**2 for x,y in zip(av,bv)))<=256*eps*1e-4+rtol*max(math.sqrt(sum(x*x for x in av)),1e-4)
vtk=json.loads((result/'heartbeat-vtk-verification.json').read_text()); old_vtk=json.loads((historical/'heartbeat-vtk-verification.json').read_text())
checks['native_vtk_verifier']=vtk['status']=='visualization_vtk_verified' and vtk['formal_acceptance'] is False and len(vtk['frames'])==20
checks['native_vtk_frame_topology_exact']=[x['points'] for x in vtk['frames']]==[x['points'] for x in old_vtk['frames']]
checks['native_vtk_frame_times_exact']=[x['time_s'] for x in vtk['frames']]==[x['time_s'] for x in old_vtk['frames']]

detail=[json.loads(line) for line in (sample/'assembly-detail.jsonl').read_text().splitlines()]
assembly=[x for x in detail if x['kind']=='assembly']; linear=[x for x in detail if x['kind']=='linear-solve-attempt']
checks['all_assembly_records_complete']=bool(assembly) and all(x['status']=='completed' and x['team_size']==case['production_threads'] and x['batch_capacity']==case['batch'] for x in assembly)
full=[x for x in assembly if x.get('assembly_request')=='ResidualAndJacobian']; residual=[x for x in assembly if x.get('assembly_request')=='ResidualOnly']
checks['residual_full_paths_visible']=bool(full) and bool(residual) and all(x['reason']=='line-search' for x in residual) and all(x['reason']!='line-search' for x in full)
checks['volume_basis_cache_contract']=all(x.get('volume_basis_cache_enabled')==1 and x.get('volume_basis_cache_hits')==x.get('volume_quadrature_points') and x.get('volume_basis_cache_misses')==0 for x in assembly)
builds=sum(x.get('reuse_preconditioner')==0 for x in linear); attempts=sum(x.get('reuse_preconditioner')==1 for x in linear)
rebuilds=sum(x.get('reuse_preconditioner')==1 and i+1<len(detail) and detail[i+1]['kind']=='linear-solve-attempt' and detail[i+1].get('reuse_preconditioner')==0 for i,x in enumerate(detail))
checks['preconditioner_reuse_visible']=bool(linear) and builds>0 and attempts>0 and attempts>=rebuilds and builds+attempts==len(linear)
checks['true_linear_residual_gate']=all(x['status']=='completed' and x['ksp_reason']>0 and math.isfinite(x['true_linear_relative_residual']) and x['true_linear_relative_residual']<=case['reused_preconditioner_true_linear_relative_tolerance'] for x in linear)
text=(sample/'stdout.txt').read_text(); profiles=re.findall(r'hpc_profile (\{.*\})',text); profile=json.loads(profiles[-1]) if profiles else None
checks['solver_complete_and_profiled']=bool(profile) and 'bifurcation_fsi completed steps=20' in text and 'Exit status: 0' in (sample/'timing.txt').read_text()
checks['archive_integrity']=Path(root/'heartbeat-paraview.tar.gz').is_file() and Path(root/'archive.sha256').is_file() and 'results/run.json' in (root/'archive-files.txt').read_text() and 'results/step-20/fluid.vtu' in (root/'archive-files.txt').read_text() and 'results/wall-display.pvd' in (root/'archive-files.txt').read_text()
archive_digest=hashlib.sha256()
with (root/'heartbeat-paraview.tar.gz').open('rb') as stream:
    for chunk in iter(lambda:stream.read(1024*1024),b''): archive_digest.update(chunk)
archive_hash=archive_digest.hexdigest()
checks['archive_sha256_matches']=archive_hash==(root/'archive.sha256').read_text().split()[0]

elapsed=profile['elapsed_s'] if profile else math.nan; historical_elapsed=case['historical_solver_elapsed_s']
status='passed' if all(checks.values()) else 'failed_numerical'
correctness={'schema_version':1,'status':status,'checks':checks,'job_id':os.environ['SLURM_JOB_ID'],'formal_physics_acceptance':False,'conservation_failed_steps':run.get('conservation_failed_steps'),'scope':'full 20-step visualization-only heartbeat with historical per-quantity comparison, native VTK, and package'}
timing={'schema_version':1,'solver_elapsed_s':elapsed,'historical_solver_elapsed_s':historical_elapsed,'historical_elapsed_reduction_fraction':1-elapsed/historical_elapsed,'phase_profile':profile,'assembly_records':len(assembly),'full_assembly_records':len(full),'residual_only_records':len(residual),'preconditioner_builds':builds,'preconditioner_reuse_attempts':attempts,'preconditioner_rebuilds':rebuilds,'threads':case['production_threads'],'batch':case['batch'],'historical_comparison':'cross-job/cross-node reference, not paired full-cycle timing','archive_sha256':archive_hash}
(root/'correctness.json').write_text(json.dumps(correctness,indent=2)+'\n')
(root/'timing.json').write_text(json.dumps(timing,indent=2)+'\n')
(root/'summary.md').write_text(f'# R7 full heartbeat: {status}\n\nSolver {elapsed:.3f} s versus historical {historical_elapsed:.3f} s ({1-elapsed/historical_elapsed:.3%} reduction; cross-job reference). Conservation failures: {run.get("conservation_failed_steps")}/20; formal acceptance remains false. Archive SHA256: {archive_hash}. Failed: {[k for k,v in checks.items() if not v]}.\n')
print(json.dumps({'status':status,'solver_elapsed_s':elapsed,'historical_reduction_fraction':1-elapsed/historical_elapsed,'conservation_failed_steps':run.get('conservation_failed_steps'),'archive_sha256':archive_hash,'failed_checks':[k for k,v in checks.items() if not v]}))
sys.exit(0 if status=='passed' else 1)
