#!/usr/bin/env python3
"""Compute-node analysis of independent process pairs; no physics tolerance relaxation."""
import json
import os
from pathlib import Path
import re
import statistics
import sys

root=Path(sys.argv[1])
baseline,candidate=sys.argv[2].split(',')
expected={0:'d67dfc4ff72363aa99b0b7228a3328fdfb162f687f78d74189eb8816a048cf8f',
          1:'982aa0be0e6fc99bfc949526a1b227a963bed900c01bb459e469c9366790c373'}
checks={}
records=[]
details=[]
pattern=re.compile(r'benchmark_sample trial=(-?\d+) state=(\d+) assembly_s=(\S+) sha256=(\w+) mallocs=(\S+)')
labels=[baseline,candidate+'-off',candidate+'-on']
for pair in range(3):
    for label in labels:
        folder=root/f'W1-pair-{pair}-{label}'
        stdout=(folder/'stdout.txt').read_text()
        samples=[]
        for m in pattern.finditer(stdout):
            sample,state=int(m[1]),int(m[2])
            checks[f'{pair}-{label}-{sample}-exact_operator']=m[4]==expected[state]
            checks[f'{pair}-{label}-{sample}-mallocs_zero']=float(m[5])==0.
            if sample>=0: samples.append({'sample':sample,'state':state,'seconds':float(m[3]),'sha256':m[4]})
        checks[f'{pair}-{label}-sample_count']=len(samples)==6
        checks[f'{pair}-{label}-actual_team']=f'team_size={os.environ["OMP_NUM_THREADS"]}' in stdout
        profile=json.loads(stdout.split('hpc_profile ')[-1].splitlines()[0])
        entry={'pair':pair,'label':label,'samples':samples,'median_s':statistics.median(s['seconds'] for s in samples),'coarse_profile':profile}
        records.append(entry)
        if label.endswith('-on'):
            raw=[json.loads(line) for line in (folder/'assembly-detail.jsonl').read_text().splitlines()]
            assembly=[r for r in raw if r['kind']=='assembly' and r['status']=='completed']
            checks[f'{pair}-detail_count']=len(assembly)==7
            detailed=sum(r['wall_s'] for r in assembly)
            coarse=profile['phases']['assembly']['exclusive_s']
            checks[f'{pair}-detail_vs_coarse']=abs(detailed-coarse)<=max(.02,.003*coarse)
            for r in assembly:
                # Nested scatter/ports are already included in consume; worker work is NOT wall time.
                accounted=sum(r[k] for k in ['prepare_wall_s','compute_wall_s','consume_wall_s'])+sum(r['timers_s'].get(k,0.) for k in ['ghost','mat_vec_assembly','port_measurement_gauge','hash_diagnostic'])
                checks[f'{pair}-accounting-{len(details)}']=accounted<=r['wall_s']*1.003 and accounted>=r['wall_s']*.95
                checks[f'{pair}-worker-total-{len(details)}']=r['worker_work_s']>=r['max_worker_s']
                checks[f'{pair}-generation-{len(details)}']=r['state_generation']==r['residual_generation']==r['jacobian_generation']
                details.append(r)
        else:
            checks[f'{pair}-{label}-disabled_no_detail']=not (folder/'assembly-detail.jsonl').exists()
checks['W0_batch_serial']='passed' in (root/'W0-parallel_element_batch_test.stdout').read_text()
checks['W0_batch_openmp']='passed' in (root/'W0-parallel_element_batch_openmp_test.stdout').read_text()
checks['W0_caller_funneled']='passed' in (root/'W0-element_assembly_execution_openmp_test.stdout').read_text()
checks['W0_thread_single_negative']='passed' in (root/'W0-mpi-single.stdout').read_text()
volume_log=(root/'W0-volume/stdout.txt').read_text()
checks['W0_volume_expanded_compact_exact_retry_abort']='passed' in volume_log
raw=[json.loads(line) for line in (root/'W0-volume/assembly-detail.jsonl').read_text().splitlines()]
checks['W0_fault_log']=any(r['kind']=='assembly' and r['status']=='failed' for r in raw)
off=[r['median_s'] for r in records if r['label']==candidate+'-off']
on=[r['median_s'] for r in records if r['label']==candidate+'-on']
base=[r['median_s'] for r in records if r['label']==baseline]
overheads=[b/a-1 for a,b in zip(off,on)]
overhead=statistics.median(overheads)
numerical=all(checks.values())
status='passed' if numerical and overhead<.03 else 'inconclusive' if numerical else 'failed_numerical'
hot_worker={k:sum(r[k] for r in details) for k in ['volume_worker_work_s','trace_worker_work_s','wall_worker_work_s']}
hot_caller={k:sum(r['timers_s'].get(k,0.) for r in details) for k in ['ghost','scatter','ports','mat_vec_assembly','port_measurement_gauge','hash_diagnostic']}
correctness={'schema_version':1,'status':status,'checks':checks,'profiling_overhead_fraction':overhead,'profiling_overhead_gate':overhead<.03,'formal_physics_acceptance':False}
timing={'schema_version':1,'records':records,'same_binary_overhead_pairs':overheads,'same_binary_median_overhead':overhead,'baseline_medians_s':base,'candidate_disabled_medians_s':off,'top_worker_hot_paths_work_seconds':sorted(hot_worker.items(),key=lambda x:-x[1]),'caller_hot_paths_wall_seconds':sorted(hot_caller.items(),key=lambda x:-x[1]),'accounting_note':'worker-work seconds overlap; scatter/ports nested inside consume, not additive to it; six samples in two states per process'}
manifest={'schema_version':1,'variant':candidate,'baseline':baseline,'test_id':'P0-W0-W1','job_id':os.environ['SLURM_JOB_ID'],'node':os.environ.get('SLURM_JOB_NODELIST'),'threads':int(os.environ['OMP_NUM_THREADS']),'batch':int(os.environ['IGA_ASSEMBLY_BATCH_SIZE']),'source_manifests':[f'{baseline}-source.sha256',f'{candidate}-source.sha256'],'binary_manifests':[f'{baseline}-binary.sha256',f'{candidate}-binary.sha256'],'input_manifest':'input.sha256','workload_manifest':'workload-source.sha256','compile_logs':[f'{baseline}-build.stdout',f'{candidate}-build.stdout'],'assertions':'-UNDEBUG in benchmark/volume; default Makefile has no NDEBUG','physics_unchanged':True}
for name,obj in [('correctness.json',correctness),('timing.json',timing),('manifest.json',manifest)]:
    (root/name).write_text(json.dumps(obj,indent=2)+'\n')
(root/'summary.md').write_text(f'# P0 W0/W1: {status}\n\nJob {manifest["job_id"]}; threads={manifest["threads"]}, batch={manifest["batch"]}.\nSame-binary profiling overhead median: {overhead:.3%} (gate <3%).\nW1 uses three independent process triples, two states, three measured samples/state after one warmup.\nNo full-FSI performance or formal-physics acceptance claimed.\n\nWorker hot paths (work seconds, not wall):\n'+''.join(f'- {k}: {v:.6f}\n' for k,v in timing['top_worker_hot_paths_work_seconds'])+'\nFailed checks: '+str([k for k,v in checks.items() if not v])+'\n')
print(json.dumps({'status':status,'job_id':manifest['job_id'],'overhead':overhead,'failed_checks':[k for k,v in checks.items() if not v]}))
sys.exit(0 if status=='passed' else 3 if status=='inconclusive' else 1)
