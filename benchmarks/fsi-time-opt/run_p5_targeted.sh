#!/bin/bash
set -euo pipefail
while test $# -gt 0; do case "$1" in --root) root=$2;shift 2;;--variants) variant=$2;shift 2;;--case-list) cases=$2;shift 2;;*)exit 2;;esac;done
: "${LOCAL:?}" "${SLURM_JOB_ID:?}" "${root:?}" "${variant:?}" "${cases:?}"
work="$LOCAL/fsi-p5-$SLURM_JOB_ID";out="$root/job-$SLURM_JOB_ID";mkdir -p "$work/evidence" "$out" "$work/$variant";stage=staging
finish(){ result=$?;trap - EXIT TERM;printf '{"job_id":"%s","stage":"%s","status":"%s","exit_code":%s}\n' "$SLURM_JOB_ID" "$stage" "$([ "$result" -eq 0 ]&&echo passed||echo failed)" "$result">"$work/evidence/status.json";cp -au "$work/evidence/." "$out/"||{ test "$result" -ne 0||result=1;};exit "$result";};trap finish EXIT;trap 'exit 124' TERM
sha256sum -c "$root/$variant-source.sha256";tar -xf "$root/$variant-source.tar" -C "$work/$variant";cp "$root/$variant-source.sha256" "$cases" "$work/evidence/"
src="$work/$variant";tools="$src/benchmarks/fsi-time-opt";stage=build
timeout -k 30s 900s /usr/bin/time -v -o "$work/evidence/build.time" make -f "$tools/fsi_case.mk" SOURCE="$src" BUILD="$src" CASE_SOURCE="$src/solvers/cpu/tests/test_immersed_transient_flow.cpp" CASE_NAME=p5-targeted fsi_case >"$work/evidence/build.stdout" 2>"$work/evidence/build.stderr"
cp "$src/p5-targeted" "$work/evidence/";sha256sum "$src/p5-targeted">"$work/evidence/binary.sha256"
stage=targeted;mkdir -p "$work/evidence/targeted";cd "$work/evidence/targeted"
IGA_PROFILE_DETAIL=1 timeout -k 30s 2400s /usr/bin/time -v -o timing.txt mpiexec -np 1 --map-by slot:PE="$OMP_NUM_THREADS" --bind-to core "$src/p5-targeted" --p5-preconditioner-reuse >stdout.txt 2>stderr.txt
stage=verification
python3 - "$work/evidence" <<'PY'
import json,os,pathlib,sys
r=pathlib.Path(sys.argv[1]);text=(r/'targeted/stdout.txt').read_text();timing=(r/'targeted/timing.txt').read_text()
records=[json.loads(line) for line in (r/'targeted/assembly-detail.jsonl').read_text().splitlines() if json.loads(line)['kind']=='linear-solve-attempt']
checks={
    'reuse_contract_marker':sum(x.startswith('preconditioner_reuse_contract ') for x in text.splitlines())==1,
    'fresh_and_reused_attempts_visible':any(x['reuse_preconditioner']==0 for x in records) and any(x['reuse_preconditioner']==1 for x in records),
    'all_solve_attempt_records_completed':bool(records) and all(x['status']=='completed' for x in records),
    'exit_zero':'Exit status: 0' in timing
}
status='passed' if all(checks.values()) else 'failed'
(r/'correctness.json').write_text(json.dumps({'status':status,'checks':checks,'job_id':os.environ['SLURM_JOB_ID'],'scope':'P5 same-epoch preconditioner reuse disabled/enabled/forced-fallback contract; W1/W2 pending','formal_physics_acceptance':False},indent=2)+'\n')
(r/'summary.md').write_text(f'# P5 preconditioner reuse targeted: {status}\n\nFailed: {[k for k,v in checks.items() if not v]}. W1/W2 pending.\n')
sys.exit(0 if status=='passed' else 1)
PY
stage=complete
