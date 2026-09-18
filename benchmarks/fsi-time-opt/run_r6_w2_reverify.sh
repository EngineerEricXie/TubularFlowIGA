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
: "${LOCAL:?Slurm LOCAL required}" "${SLURM_JOB_ID:?Slurm required}" "${root:?}" "${variant:?}" "${cases:?}"
work="$LOCAL/fsi-r6-w2-reverify-$SLURM_JOB_ID"
out="$root/job-$SLURM_JOB_ID"
mkdir -p "$work/evidence" "$out" "$work/$variant"
stage=staging
failure_status=failed_environment
started=$(date -u +%FT%TZ)
finish() {
    result=$?
    trap - EXIT TERM USR1
    test "$result" -ne 124 || failure_status=timeout
    printf '{"job_id":"%s","stage":"%s","status":"%s","exit_code":%s,"simulation_executed":false,"started":"%s","ended":"%s"}\n' "$SLURM_JOB_ID" "$stage" "$([ "$result" -eq 0 ] && echo passed || echo "$failure_status")" "$result" "$started" "$(date -u +%FT%TZ)" > "$work/evidence/status.json"
    cp -au "$work/evidence/." "$out/"
    exit "$result"
}
trap finish EXIT
trap 'failure_status=interrupted; exit 124' TERM USR1

read -r source_job expected_variant < <(python3 -c 'import json,sys; c=json.load(open(sys.argv[1])); print(c["source_job"],c["verifier_variant"])' "$cases")
test "$variant" = "$expected_variant"
sha256sum -c "$root/$variant-source.sha256"
tar -xf "$root/$variant-source.tar" -C "$work/$variant"
python3 - "$root/job-$source_job" <<'PY'
import json,pathlib,sys
r=pathlib.Path(sys.argv[1]); status=json.loads((r/'status.json').read_text())
assert status['stage']=='verification' and status['exit_code']==1 and status['simulation_executed'] is True
for count in (8,4):
    timing=(r/f'W2-threads-{count}/timing.txt').read_text()
    assert 'Exit status: 0' in timing and (r/f'W2-threads-{count}/results/history.csv').is_file()
PY
cp -a "$root/job-$source_job/." "$work/evidence/"
mv "$work/evidence/status.json" "$work/evidence/initial-status.json"
cp "$work/evidence/cases-r6-w2.json" "$work/evidence/case.json"
cp "$cases" "$work/evidence/reverification-case.json"
cp "$root/$variant-source.sha256" "$work/evidence/$variant-source.sha256"
cp "$0" "$work/evidence/reverification-runner.sh"
stage=verification
failure_status=failed_numerical
python3 "$work/$variant/benchmarks/fsi-time-opt/verify_r6_w2.py" "$work/evidence"
stage=complete
