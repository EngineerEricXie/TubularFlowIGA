#!/usr/bin/env python3
import json,os,re,statistics,sys
from pathlib import Path
root=Path(sys.argv[-1]); expected={0:'d67dfc4ff72363aa99b0b7228a3328fdfb162f687f78d74189eb8816a048cf8f',1:'982aa0be0e6fc99bfc949526a1b227a963bed900c01bb459e469c9366790c373'}
pattern=re.compile(r'benchmark_sample trial=(-?\d+) state=(\d+) assembly_s=(\S+) sha256=(\w+) mallocs=(\S+)')
checks={}; records={}; rss={}; details={}
for batch in (8,16,32):
    records[batch]=[]; details[batch]=[]; rss[batch]=[]
    for process in 'AB':
        folder=root/f'W1-batch-{batch}-{process}'; text=(folder/'stdout.txt').read_text(); matches=list(pattern.finditer(text))
        checks[f'b{batch}{process}_seven_records']=len(matches)==7
        checks[f'b{batch}{process}_operator_exact']=all(m[4]==expected[int(m[2])] and float(m[5])==0 for m in matches)
        records[batch].append(statistics.median(float(m[3]) for m in matches if int(m[1])>=0))
        rows=[json.loads(line) for line in (folder/'assembly-detail.jsonl').read_text().splitlines() if json.loads(line)['kind']=='assembly']
        checks[f'b{batch}{process}_details']=len(rows)==7 and all(r['status']=='completed' and r['team_size']==8 and r['batch_capacity']==batch and r['maximum_resident_result_payload_bytes']>0 and r['maximum_worker_completion_spread_s']>=0 for r in rows)
        details[batch]+=rows
        value=next(int(line.split(':',1)[1]) for line in (folder/'timing.txt').read_text().splitlines() if 'Maximum resident set size' in line)
        rss[batch].append(value)
fastest=min((statistics.median(v),b) for b,v in records.items())[1]
limit=statistics.median(records[fastest])*1.02
selected=min(b for b,v in records.items() if statistics.median(v)<=limit)
checks['selected_reproducibly_faster_than_8']=selected!=8 and all(x<y for x,y in zip(records[selected],records[8]))
checks['rss_below_70_percent_allocation']=max(rss[selected])<10640*1024
summary={b:{'process_medians_s':records[b],'median_s':statistics.median(records[b]),'max_rss_kib':max(rss[b]),
    'prepare_fraction':sum(r['prepare_wall_s'] for r in details[b])/sum(r['wall_s'] for r in details[b]),
    'consume_fraction':sum(r['consume_wall_s'] for r in details[b])/sum(r['wall_s'] for r in details[b]),
    'maximum_worker_completion_spread_s':max(r['maximum_worker_completion_spread_s'] for r in details[b]),
    'maximum_resident_result_payload_bytes':max(r['maximum_resident_result_payload_bytes'] for r in details[b])} for b in records}
if sys.argv[1]=='--select':
    (root/'selected-batch.txt').write_text(str(selected)+'\n'); (root/'selection.json').write_text(json.dumps({'selected':selected,'fastest':fastest,'summary':summary,'prechecks':checks},indent=2)+'\n')
    print(json.dumps({'selected':selected,'fastest':fastest,'reproducible':checks['selected_reproducibly_faster_than_8']})); sys.exit(0 if all(checks.values()) else 3)
for name in ('parallel_element_batch_test','parallel_element_batch_openmp_test','element_assembly_execution_openmp_test'):
    checks[name]='passed' in (root/f'{name}.stdout').read_text()
checks['mpi_single_negative']='passed' in (root/'mpi-single.stdout').read_text()
volume=(root/'W0-volume/stdout.txt').read_text(); checks['expanded_compact_ordered_fault_retry']='passed' in volume
rows=[json.loads(line) for line in (root/'W0-volume/assembly-detail.jsonl').read_text().splitlines()]
checks['fault_record']=any(r['kind']=='assembly' and r['status']=='failed' for r in rows)
status='passed' if all(checks.values()) else 'failed_numerical'
result={'status':status,'checks':checks,'selected_batch':selected,'summary':summary,'job_id':os.environ['SLURM_JOB_ID'],'formal_physics_acceptance':False,'W2':'pending'}
(root/'correctness.json').write_text(json.dumps(result,indent=2)+'\n'); (root/'summary.md').write_text(f'# P2 W0/W1: {status}\n\nSelected batch {selected}; fastest {fastest}. Failed: {[k for k,v in checks.items() if not v]}. W2 pending.\n')
print(json.dumps({'status':status,'selected':selected})); sys.exit(0 if status=='passed' else 1)
