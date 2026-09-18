#!/bin/bash
# Login-node coordination only. Silent until a terminal scheduler state or deadline.
set -euo pipefail
job=${1:?JobID required}
root=${2:?evidence root required}
limit=${3:-86400}
deadline=$((SECONDS+limit))
while true; do
    row=$(sacct -X -n -P -j "$job" --format=JobID,State,ExitCode,ElapsedRaw,AllocCPUS 2>/dev/null | awk -F'|' -v id="$job" '$1==id {print; exit}') || row=
    IFS='|' read -r actual state code elapsed cpus <<< "$row"
    state=${state%% *}
    case "$state" in
        COMPLETED|FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|PREEMPTED|BOOT_FAIL|DEADLINE|REVOKED)
            printf '{"job_id":"%s","state":"%s","exit_code":"%s","elapsed_s":%s,"allocated_cpus":%s,"observed_utc":"%s"}\n' "$job" "$state" "$code" "${elapsed:-0}" "${cpus:-0}" "$(date -u +%FT%TZ)" > "$root/job-$job-terminal.json"
            sacct -j "$job" --format=JobID,State,ElapsedRaw,AllocCPUS,AllocTRES,TotalCPU,MaxRSS,ExitCode -P > "$root/job-$job-accounting.txt"
            printf 'Job %s terminal: %s (%s). Evidence: %s/job-%s\n' "$job" "$state" "$code" "$root" "$job"
            exit 0;;
    esac
    if test "$SECONDS" -ge "$deadline"; then
        printf '{"job_id":"%s","monitor_status":"monitor_timeout","last_scheduler_state":"%s","observed_utc":"%s"}\n' "$job" "$state" "$(date -u +%FT%TZ)" > "$root/job-$job-terminal.json"
        printf 'Job %s monitor timeout; last state=%s; do not assume solver timed out.\n' "$job" "$state"
        exit 3
    fi
    sleep 300
done
