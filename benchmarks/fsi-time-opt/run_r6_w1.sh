#!/bin/bash
set -euo pipefail
while test $# -gt 0; do
    case "$1" in
        --root) root=$2; shift 2;;
        --variants) variant=$2; shift 2;;
        --case-list) cases=$2; shift 2;;
        *) exit 2;;
    esac
done
: "${root:?}" "${variant:?}" "${cases:?}" "${LOCAL:?}" "${SLURM_JOB_ID:?}"
work="$LOCAL/fsi-r6-w1-$SLURM_JOB_ID"
out="$root/job-$SLURM_JOB_ID"
mkdir -p "$work/evidence" "$out" "$work/$variant"
stage=staging
failure_status=failed_environment
simulation_executed=false
started=$(date -u +%FT%TZ)
finish() {
    result=$?
    trap - EXIT TERM USR1
    test "$result" -ne 124 || failure_status=timeout
    printf '{"job_id":"%s","stage":"%s","status":"%s","exit_code":%s,"simulation_executed":%s,"started":"%s","ended":"%s"}\n' "$SLURM_JOB_ID" "$stage" "$([ "$result" -eq 0 ] && echo passed || echo "$failure_status")" "$result" "$simulation_executed" "$started" "$(date -u +%FT%TZ)" > "$work/evidence/status.json"
    cp -au "$work/evidence/." "$out/"
    exit "$result"
}
trap finish EXIT
trap 'failure_status=interrupted; exit 124' TERM USR1

read -r expected_variant < <(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["variant"])' "$cases")
test "$variant" = "$expected_variant"
sha256sum -c "$root/$variant-source.sha256"
tar -xf "$root/$variant-source.tar" -C "$work/$variant"
sha256sum -c "$root/surface.sha256"
sha256sum -c "$root/archived-workload.sha256"
cp "$root/surface.txt" "$work/surface.txt"
src="$work/$variant"
tools="$src/benchmarks/fsi-time-opt"
tar -xOf "$root/archived-workload.tar.gz" ./driver.cpp > "$tools/driver.cpp"
cp "$root/$variant-source.sha256" "$cases" "$0" "$work/evidence/"
sha256sum "$tools/driver.cpp" "$tools/benchmark.cpp" "$work/surface.txt" "$cases" > "$work/evidence/input.sha256"
{
    date -u +%FT%TZ
    module list
    mpicxx --version
    lscpu
    scontrol show job "$SLURM_JOB_ID"
    grep Cpus_allowed_list /proc/self/status
} > "$work/evidence/environment.txt" 2>&1

stage=build
timeout -k 30s 900s /usr/bin/time -v -o "$work/evidence/build.time" make -f "$tools/runner.mk" SOURCE="$src" WORKLOAD_SOURCE="$src" BUILD="$src" benchmark > "$work/evidence/build.stdout" 2> "$work/evidence/build.stderr"
cp "$src/benchmark" "$work/evidence/benchmark"
sha256sum "$src/benchmark" > "$work/evidence/binary.sha256"

stage=W1
failure_status=failed_numerical
simulation_executed=true
for order in 1,2,4,8 8,4,2,1; do
    IFS=, read -r -a thread_counts <<< "$order"
    suffix=$(test "$order" = 1,2,4,8 && echo A || echo B)
    for threads in "${thread_counts[@]}"; do
        sample="$work/evidence/W1-threads-$threads-$suffix"
        mkdir -p "$sample"
        cd "$sample"
        OMP_NUM_THREADS=$threads OMP_THREAD_LIMIT=$threads IGA_ASSEMBLY_THREADS=$threads IGA_ASSEMBLY_BATCH_SIZE=32 IGA_PROFILE_DETAIL=1 timeout -k 30s 900s /usr/bin/time -v -o timing.txt mpiexec -np 1 --map-by slot:PE="$threads" --bind-to core "$src/benchmark" "$work/surface.txt" > stdout.txt 2> stderr.txt
    done
done

stage=verification
python3 - "$work/evidence" "$cases" <<'PY'
import json,os,pathlib,re,statistics,sys
root=pathlib.Path(sys.argv[1]); case=json.load(open(sys.argv[2]))
expected={0:case['operator_hashes']['zero'],1:case['operator_hashes']['nonzero']}
pattern=re.compile(r'benchmark_sample trial=(-?\d+) state=(\d+) assembly_s=(\S+) sha256=(\w+) mallocs=(\S+)')
checks={}; summaries={}; all_hashes=set()
for threads in case['threads']:
    medians=[]; rss=[]; rows=[]
    for suffix in 'AB':
        folder=root/f'W1-threads-{threads}-{suffix}'
        matches=list(pattern.finditer((folder/'stdout.txt').read_text()))
        checks[f't{threads}{suffix}_seven_records']=len(matches)==7
        checks[f't{threads}{suffix}_operator_exact']=all(m[4]==expected[int(m[2])] and float(m[5])==0 for m in matches)
        all_hashes.update(m[4] for m in matches)
        medians.append(statistics.median(float(m[3]) for m in matches if int(m[1])>=0))
        detail=[json.loads(line) for line in (folder/'assembly-detail.jsonl').read_text().splitlines() if json.loads(line)['kind']=='assembly']
        checks[f't{threads}{suffix}_details']=len(detail)==7 and all(x['status']=='completed' and x['team_size']==threads and x['batch_capacity']==case['batch'] and x['maximum_resident_result_payload_bytes']>0 and x['maximum_worker_completion_spread_s']>=0 for x in detail)
        rows.extend(detail)
        rss.append(next(int(line.split(':',1)[1]) for line in (folder/'timing.txt').read_text().splitlines() if 'Maximum resident set size' in line))
    summaries[threads]={
        'process_medians_s':medians,
        'median_s':statistics.median(medians),
        'min_process_median_s':min(medians),
        'max_process_median_s':max(medians),
        'max_rss_kib':max(rss),
        'maximum_worker_completion_spread_s':max(x['maximum_worker_completion_spread_s'] for x in rows),
        'maximum_resident_result_payload_bytes':max(x['maximum_resident_result_payload_bytes'] for x in rows)
    }
checks['exact_two_expected_operator_hashes']=all_hashes==set(expected.values())
ranked=sorted(case['threads'],key=lambda t:(summaries[t]['median_s'],t))
top_two=ranked[:2]
fastest=ranked[0]; limit=summaries[fastest]['median_s']*(1+case['closest_to_fastest_fraction'])
preferred=min(t for t in case['threads'] if summaries[t]['median_s']<=limit)
base=summaries[1]['median_s']
for threads in case['threads']:
    summaries[threads]['speedup_vs_1']=base/summaries[threads]['median_s']
    summaries[threads]['parallel_efficiency_vs_1']=(base/summaries[threads]['median_s'])/threads
checks['top_two_distinct']=len(top_two)==2 and top_two[0]!=top_two[1]
checks['rss_below_70_percent_allocation']=max(summaries[t]['max_rss_kib'] for t in case['threads'])<10640*1024
status='passed' if all(checks.values()) else 'failed_numerical'
(root/'selected-w2-threads.txt').write_text(','.join(map(str,top_two))+'\n')
(root/'preferred-w1-thread.txt').write_text(str(preferred)+'\n')
timing={'schema_version':1,'summaries':{str(k):v for k,v in summaries.items()},'ranked_threads':ranked,'top_two_threads':top_two,'fastest_thread':fastest,'preferred_within_5_percent_thread':preferred,'batch':case['batch'],'paired_processes_per_configuration':2,'allocation_cpus':case['allocation_cpus']}
correctness={'schema_version':1,'status':status,'checks':checks,'job_id':os.environ['SLURM_JOB_ID'],'W2':'pending_for_top_two','formal_physics_acceptance':False}
(root/'timing.json').write_text(json.dumps(timing,indent=2)+'\n')
(root/'correctness.json').write_text(json.dumps(correctness,indent=2)+'\n')
(root/'summary.md').write_text(f'# R6 W1 thread scaling: {status}\n\nRanked: {ranked}; W2 top two: {top_two}; preferred within 5%: {preferred}. Failed: {[k for k,v in checks.items() if not v]}.\n')
print(json.dumps({'status':status,'ranked':ranked,'top_two':top_two,'preferred':preferred}))
sys.exit(0 if status=='passed' else 1)
PY
stage=complete
