#!/usr/bin/env python3
"""P3 matched-FSI gate: exact physics and state-aligned assembly semantics."""
import hashlib
import json
import os
from pathlib import Path
import statistics
import sys

root=Path(sys.argv[1]); before,after=sys.argv[2:4]
case=json.loads((root/'case.json').read_text()); checks={}; manifests={}; profiles={}; raw={}

def records(path):
    return [json.loads(line) for line in path.read_text().splitlines()
            if json.loads(line)['kind'] in ('assembly','assembly-reuse')]

def normalized(rows):
    expanded=[]; previous=None
    for row in rows:
        if row['kind']=='assembly-reuse':
            checks[f'baseline_reuse_{len(expanded)}_follows_accepted_line_search']=(
                previous is not None and previous['kind']=='assembly' and previous['reason']=='line-search'
                and previous['status']==row['status']=='completed'
                and row['newton_iteration']==previous['newton_iteration']+1)
            clone=dict(previous); clone.update(row)
            clone.update(kind='assembly',reason='next-Newton',damping=0)
            expanded.append(clone)
        else:
            expanded.append(row)
        previous=row
    return expanded

for variant in (before,after):
    sample=root/f'{variant}-fsi'; stdout=(sample/'stdout.txt').read_text()
    checks[variant+'_process_passed']='compliant_channel_fsi converged=true' in stdout and 'Exit status: 0' in (sample/'timing.txt').read_text()
    profiles[variant]=json.loads(stdout.split('hpc_profile ')[-1].splitlines()[0])
    checks[variant+'_profile_status_zero']=profiles[variant]['status']==0
    manifest=json.loads((sample/'reference/manifest.json').read_text()); manifests[variant]=manifest
    checks[variant+'_native_gates']=manifest['native_gates_passed'] is True
    for field in manifest['fields']:
        checks[f'{variant}_{field["name"]}_saved_sha']=hashlib.sha256((sample/'reference'/field['file']).read_bytes()).hexdigest()==field['sha256']
    raw[variant]=records(sample/'assembly-detail.jsonl')
    checks[variant+'_all_records_completed']=all(row['status']=='completed' for row in raw[variant])
    expected=case.get('batch_by_variant',{}).get(variant,case['batch'])
    checks[variant+'_saved_batch']=int((sample/'batch-size.txt').read_text())==expected
    checks[variant+'_assembly_team_batch']=all(row['kind']!='assembly' or (row['team_size']==case['threads'] and row['batch_capacity']==expected) for row in raw[variant])

checks['all_native_fields_bitwise_exact']=manifests[before]==manifests[after]
def physics_line(variant):
    return next(line for line in (root/f'{variant}-fsi/stdout.txt').read_text().splitlines() if line.startswith('compliant_channel_fsi converged=true'))
checks['coupling_conservation_summary_bitwise_exact']=physics_line(before)==physics_line(after)

baseline=normalized(raw[before]); candidate=raw[after]
checks['baseline_contains_reuse']=any(row['kind']=='assembly-reuse' for row in raw[before])
checks['candidate_contains_no_reuse']=not any(row['kind']=='assembly-reuse' for row in candidate)
candidate_full=[row for row in candidate if row.get('assembly_request')=='ResidualAndJacobian']
candidate_residual=[row for row in candidate if row.get('assembly_request')=='ResidualOnly']
baseline_full=[row for row in raw[before] if row['kind']=='assembly']
checks['candidate_residual_only_is_line_search']=bool(candidate_residual) and all(row['reason']=='line-search' for row in candidate_residual)
checks['candidate_full_is_not_line_search']=bool(candidate_full) and all(row['reason']!='line-search' for row in candidate_full)
checks['candidate_residual_generation_current']=all(row['residual_generation']==row['state_generation'] for row in candidate)
checks['candidate_line_search_jacobian_stale']=all(row['jacobian_generation']!=row['state_generation'] for row in candidate_residual)
checks['candidate_full_jacobian_current']=all(row['jacobian_generation']==row['state_generation'] for row in candidate_full)
checks['line_search_call_count_preserved']=sum(row['reason']=='line-search' for row in baseline_full)==len(candidate_residual)
checks['next_newton_rebuild_count_matches_baseline_reuse']=sum(row['kind']=='assembly-reuse' for row in raw[before])==sum(row['reason']=='next-Newton' for row in candidate)
checks['logical_work_event_count_preserved']=len(baseline)==len(candidate)
checks['real_full_assembly_count_reduced']=len(candidate_full)<len(baseline_full)

keys=['reason','status','step','coupling_iteration','newton_iteration','damping','active_dofs',
      'geometry_epoch','input_identity','state_identity','history_identity','moving_map_identity',
      'state_generation','residual_generation','inside_cells','cut_cells',
      'volume_quadrature_points','surface_quadrature_points']
checks['exact_state_aligned_assembly_sequence']=[{k:r.get(k) for k in keys} for r in baseline]==[{k:r.get(k) for k in keys} for r in candidate]
non_line=[(b,c) for b,c in zip(baseline,candidate) if c['reason']!='line-search']
checks['exact_current_j_generation_on_full_events']=all(b.get('jacobian_generation')==c.get('jacobian_generation') for b,c in non_line)
residual_payload=[row.get('maximum_resident_result_payload_bytes') for row in candidate_residual]
full_payload=[row.get('maximum_resident_result_payload_bytes') for row in candidate_full]
checks['residual_payload_recorded_and_smaller']=(all(v is not None for v in residual_payload+full_payload)
    and max(residual_payload)<min(full_payload))

status='passed' if all(checks.values()) else 'failed_numerical'
timing={'schema_version':1,'profiles':profiles,'before_elapsed_s':profiles[before]['elapsed_s'],
        'after_elapsed_s':profiles[after]['elapsed_s'],'baseline_real_full_calls':len(baseline_full),
        'candidate_real_full_calls':len(candidate_full),'candidate_residual_only_calls':len(candidate_residual),
        'removed_real_full_calls':len(baseline_full)-len(candidate_full),
        'candidate_full_wall_median_s':statistics.median(row['wall_s'] for row in candidate_full),
        'candidate_residual_wall_median_s':statistics.median(row['wall_s'] for row in candidate_residual),
        'timing_scope':'ONE matched small-FSI screening pair only'}
result={'schema_version':1,'status':status,'checks':checks,'phase':case['phase'],
        'job_id':os.environ['SLURM_JOB_ID'],'formal_physics_acceptance':False,
        'promoted':False,'W2':'pending'}
manifest={'schema_version':1,'phase':'P3-matched-small-FSI','variants':[before,after],
          'job_id':os.environ['SLURM_JOB_ID'],'node':os.environ['SLURM_JOB_NODELIST'],
          'threads':case['threads'],'batch_by_variant':case.get('batch_by_variant',{}),
          'same_allocation':True,'native_reference':True}
for name,data in [('correctness.json',result),('timing.json',timing),('manifest.json',manifest)]:
    (root/name).write_text(json.dumps(data,indent=2)+'\n')
failed=[key for key,value in checks.items() if not value]
(root/'summary.md').write_text(
    f'# P3 matched small FSI: {status}\n\nNative fields and physical summary exact: {checks["all_native_fields_bitwise_exact"] and checks["coupling_conservation_summary_bitwise_exact"]}.\n'
    f'Full assemblies {len(baseline_full)} -> {len(candidate_full)}; residual-only calls {len(candidate_residual)}.\nFailed checks: {failed}. W2 pending; P3 not promoted.\n')
print(json.dumps({'status':status,'failed':failed,'full_before':len(baseline_full),'full_after':len(candidate_full),'residual_only':len(candidate_residual)}))
sys.exit(0 if status=='passed' else 1)
