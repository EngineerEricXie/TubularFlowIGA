#!/bin/bash
set -euo pipefail
while test $# -gt 0; do case "$1" in --root) root=$2;shift 2;;--variants) variant=$2;shift 2;;--case-list) cases=$2;shift 2;;*)exit 2;;esac;done
: "${LOCAL:?}" "${SLURM_JOB_ID:?}" "${root:?}" "${variant:?}" "${cases:?}"
work="$LOCAL/fsi-p4-$SLURM_JOB_ID";out="$root/job-$SLURM_JOB_ID";mkdir -p "$work/evidence" "$out" "$work/$variant";stage=staging
finish(){ result=$?;trap - EXIT TERM;printf '{"job_id":"%s","stage":"%s","status":"%s","exit_code":%s}\n' "$SLURM_JOB_ID" "$stage" "$([ "$result" -eq 0 ]&&echo passed||echo failed)" "$result">"$work/evidence/status.json";cp -au "$work/evidence/." "$out/"||{ test "$result" -ne 0||result=1;};exit "$result";};trap finish EXIT;trap 'exit 124' TERM
sha256sum -c "$root/$variant-source.sha256";tar -xf "$root/$variant-source.tar" -C "$work/$variant";cp "$root/$variant-source.sha256" "$cases" "$work/evidence/"
src="$work/$variant";tools="$src/benchmarks/fsi-time-opt";stage=build
timeout -k 30s 900s /usr/bin/time -v -o "$work/evidence/build.time" make -f "$tools/fsi_case.mk" SOURCE="$src" BUILD="$src" CASE_SOURCE="$src/solvers/cpu/tests/test_immersed_transient_flow.cpp" CASE_NAME=p4-targeted fsi_case >"$work/evidence/build.stdout" 2>"$work/evidence/build.stderr"
cp "$src/p4-targeted" "$work/evidence/";sha256sum "$src/p4-targeted">"$work/evidence/binary.sha256"
stage=targeted;mkdir -p "$work/evidence/targeted";cd "$work/evidence/targeted"
IGA_PROFILE_DETAIL=1 timeout -k 30s 2400s /usr/bin/time -v -o timing.txt mpiexec -np 1 --map-by slot:PE="$OMP_NUM_THREADS" --bind-to core "$src/p4-targeted" --p4-cache >stdout.txt 2>stderr.txt
stage=verification
python3 - "$work/evidence" <<'PY'
import json,os,pathlib,sys
r=pathlib.Path(sys.argv[1]);text=(r/'targeted/stdout.txt').read_text();timing=(r/'targeted/timing.txt').read_text();checks={'cache_contract_marker':sum(x.startswith('volume_basis_cache_contract passed ') for x in text.splitlines())==1,'exit_zero':'Exit status: 0' in timing};status='passed' if all(checks.values()) else 'failed';(r/'correctness.json').write_text(json.dumps({'status':status,'checks':checks,'job_id':os.environ['SLURM_JOB_ID'],'scope':'P4 epoch-local volume basis cache targeted contract; W1/W2 pending','formal_physics_acceptance':False},indent=2)+'\n');(r/'summary.md').write_text(f'# P4 volume basis cache targeted: {status}\n\nFailed: {[k for k,v in checks.items() if not v]}. W1/W2 pending.\n');sys.exit(0 if status=='passed' else 1)
PY
stage=complete
