#!/usr/bin/env python3
"""Verify actual preassembly/initial input equality before allowing P1-A edits."""
import json
import os
from pathlib import Path
import sys

root=Path(sys.argv[1])
rows=[json.loads(x) for x in (root/'assembly-detail.jsonl').read_text().splitlines()]
assembly=[r for r in rows if r['kind']=='assembly']
keys=['step','coupling_iteration','active_dofs','geometry_epoch','input_identity','state_identity',
      'history_identity','moving_map_identity','state_generation','inside_cells','cut_cells',
      'volume_quadrature_points','surface_quadrature_points']
pairs=[]
for i,a in enumerate(assembly):
    if a['reason']!='adapter-initial' or a['status']!='completed': continue
    b=assembly[i+1] if i+1<len(assembly) else {}
    pairs.append({'index':i,'step':a['step'],'coupling_iteration':a['coupling_iteration'],
        'next_is_initial':b.get('reason')=='initial' and b.get('status')=='completed',
        'identical':all(k in b and a[k]==b[k] for k in keys),
        'identity_fields':{k:a[k] for k in keys},'preassembly_seconds':a['wall_s']})
stdout=(root/'stdout.txt').read_text()
checks={'at_least_one_pair':bool(pairs),'all_preassemblies_have_identical_initial':all(p['next_is_initial'] and p['identical'] for p in pairs),
        'existing_complete_fsi_regression':'compliant_channel_fsi converged=true' in stdout,
        'accepted_state_reference':(root/'reference/manifest.json').is_file()}
status='passed' if all(checks.values()) else 'failed_numerical'
profile=json.loads(stdout.split('hpc_profile ')[-1].splitlines()[0])
result={'schema_version':1,'status':status,'job_id':os.environ['SLURM_JOB_ID'],'checks':checks,'pairs':pairs,
        'preassembly_calls':len(pairs),'assembly_calls':len(assembly),
        'avoidable_preassembly_wall_s':sum(p['preassembly_seconds'] for p in pairs),
        'note':'Before-change proof only; no speedup claim; adapter publishes only after SolveTrial',
        'formal_physics_acceptance':False}
manifest={'schema_version':1,'phase':'P1-A-proof-before-change','variant':sys.argv[2],'job_id':os.environ['SLURM_JOB_ID'],
          'node':os.environ.get('SLURM_JOB_NODELIST'),'threads':int(os.environ['OMP_NUM_THREADS']),
          'batch':int(os.environ['IGA_ASSEMBLY_BATCH_SIZE']),'source_hashes':'source.sha256',
          'binary_hashes':'binary.sha256','runner_hashes':'runner-input.sha256','assertions':True}
for name,obj in [('correctness.json',result),('timing.json',profile),('manifest.json',manifest)]:
    (root/name).write_text(json.dumps(obj,indent=2)+'\n')
(root/'summary.md').write_text(f'# P1-A before-change proof: {status}\n\nIdentical preassembly/initial pairs: {len(pairs)}; total assemblies: {len(assembly)}.\nPreassembly wall seconds: {result["avoidable_preassembly_wall_s"]:.6f}.\nNo production algorithm change yet.\n')
print(json.dumps({'status':status,'pairs':len(pairs),'failed_checks':[k for k,v in checks.items() if not v]}))
sys.exit(0 if status=='passed' else 1)
