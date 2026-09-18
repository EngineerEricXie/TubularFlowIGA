#!/bin/bash
set -euo pipefail
while test $# -gt 0;do case "$1" in --root)root=$2;shift 2;;--variants)variants=$2;shift 2;;--case-list)cases=$2;shift 2;;*)exit 2;;esac;done
: "${LOCAL:?}" "${SLURM_JOB_ID:?}" "${root:?}" "${variants:?}" "${cases:?}"
IFS=, read -r baseline candidate <<< "$variants"
work="$LOCAL/fsi-w2-reverify-$SLURM_JOB_ID";out="$root/job-$SLURM_JOB_ID";mkdir -p "$work/evidence" "$out";stage=staging
finish(){ result=$?;trap - EXIT TERM;printf '{"job_id":"%s","stage":"%s","status":"%s","exit_code":%s,"simulation_executed":false}\n' "$SLURM_JOB_ID" "$stage" "$([ "$result" -eq 0 ]&&echo passed||echo failed)" "$result">"$work/evidence/status.json";cp -au "$work/evidence/." "$out/"||{ test "$result" -ne 0||result=1;};exit "$result";};trap finish EXIT;trap 'exit 124' TERM
read -r source_job verifier expected_baseline expected_candidate < <(python3 -c 'import json,sys;c=json.load(open(sys.argv[1]));print(c["source_job"],c["verifier_variant"],c["baseline"],c["candidate"])' "$cases")
test "$baseline" = "$expected_baseline";test "$candidate" = "$expected_candidate"
sha256sum -c "$root/$verifier-source.sha256";mkdir -p "$work/verifier";tar -xf "$root/$verifier-source.tar" -C "$work/verifier"
python3 - "$root/job-$source_job/status.json" <<'PY'
import json,sys
s=json.load(open(sys.argv[1]));assert s['simulation_executed'] is True and s['stage']=='verification' and s['exit_code']==1
PY
cp -a "$root/job-$source_job/." "$work/evidence/"
mv "$work/evidence/correctness.json" "$work/evidence/initial-verifier-correctness.json"
mv "$work/evidence/timing.json" "$work/evidence/initial-verifier-timing.json"
mv "$work/evidence/summary.md" "$work/evidence/initial-verifier-summary.md"
cp "$cases" "$work/evidence/reverification-case.json";cp "$root/$verifier-source.sha256" "$work/evidence/$verifier-source.sha256";cp "$0" "$work/evidence/reverification-runner.sh"
stage=verification
python3 "$work/verifier/benchmarks/fsi-time-opt/verify_w2.py" "$work/evidence" "$baseline" "$candidate"
stage=complete
