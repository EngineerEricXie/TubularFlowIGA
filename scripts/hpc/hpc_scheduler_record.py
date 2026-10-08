#!/usr/bin/env python3
"""Write a machine-readable record for one scheduler attempt."""

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess


def probe(argv):
    try:
        result = subprocess.run(argv, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, timeout=20)
        return {"argv": argv, "returncode": result.returncode,
                "output": result.stdout.rstrip()}
    except (OSError, subprocess.TimeoutExpired) as error:
        return {"argv": argv, "returncode": 127, "output": str(error)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--status", choices=("starting", "passed", "failed", "checkpointed"),
                        required=True)
    parser.add_argument("--returncode", type=int)
    parser.add_argument("--checkpoint-requested", action="store_true")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.command[:1] == ["--"]:
        args.command = args.command[1:]
    job = os.environ.get("SLURM_JOB_ID")
    report = {
        "schema_version": 1, "kind": "hpc_scheduler_attempt",
        "updated_utc": datetime.now(timezone.utc).isoformat(),
        "status": args.status, "returncode": args.returncode,
        "checkpoint_requested": args.checkpoint_requested,
        "command": args.command,
        "slurm": {name: os.environ.get(name) for name in (
            "SLURM_JOB_ID", "SLURM_RESTART_COUNT", "SLURM_JOB_NUM_NODES",
            "SLURM_NODELIST", "SLURM_NTASKS", "SLURM_NTASKS_PER_NODE",
            "SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE", "SLURM_MEM_PER_CPU")},
        "job": probe(["scontrol", "show", "job", "-dd", job]) if job else None,
        "step_accounting": probe(["sstat", "-j", f"{job}.batch", "--parsable2",
                                  "--noheader", "--format=JobID,NTasks,MaxRSS,AveRSS,AveCPU"])
                           if job else None,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(args.output.name + ".tmp")
    temporary.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    temporary.replace(args.output)


if __name__ == "__main__":
    main()
