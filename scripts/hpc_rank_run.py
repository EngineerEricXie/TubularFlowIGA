#!/usr/bin/env python3
"""Measure one executable per MPI rank; never label a failed run as a baseline."""

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

from hpc_inventory import RESOURCE_ENV, digest


def rank_identity(environment):
    for rank_key, size_key in (("OMPI_COMM_WORLD_RANK", "OMPI_COMM_WORLD_SIZE"),
                               ("PMI_RANK", "PMI_SIZE"), ("SLURM_PROCID", "SLURM_NTASKS")):
        if rank_key in environment:
            rank, size = int(environment[rank_key]), int(environment[size_key])
            if size < 1 or rank < 0 or rank >= size:
                raise ValueError("invalid rank environment")
            return rank, size
    return 0, 1


def run(command, directory, expected_ranks, timeout):
    rank, size = rank_identity(os.environ)
    if size != expected_ranks:
        raise ValueError(f"launcher reports {size} ranks, expected {expected_ranks}")
    # The controller creates the parent. Exclusive rank directories prevent reruns
    # from mixing old and new evidence. Each fresh-process repetition uses a new parent.
    directory = directory.resolve()
    destination = directory / f"rank-{rank}"
    destination.mkdir()
    timestamp = datetime.now(timezone.utc).isoformat()
    wall_start = time.perf_counter()
    timed_out = False
    with (destination / "stdout.log").open("x") as stdout, (destination / "stderr.log").open("x") as stderr:
        process = subprocess.Popen(
            ["/usr/bin/time", "-f", "%M %U %S", "-o", str(destination / "resource.txt"), "--"] + command,
            stdout=stdout, stderr=stderr, start_new_session=True)
        try:
            returncode = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=2)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            returncode = 124
    elapsed = time.perf_counter() - wall_start
    resource = destination / "resource.txt"
    usage = None
    if resource.is_file():
        # GNU time emits an explanatory line before the last line on nonzero exit.
        lines = resource.read_text().splitlines()
        if lines:
            fields = lines[-1].split()
            if len(fields) == 3:
                try:
                    usage = {"peak_rss_bytes": int(fields[0]) * 1024,
                             "user_cpu_s": float(fields[1]), "system_cpu_s": float(fields[2])}
                except ValueError:
                    pass
    report = {"schema_version": 1, "kind": "hpc_rank_run", "created_utc": timestamp,
              "command_argv": command, "cwd": str(Path.cwd()), "rank": rank, "ranks": size,
              "hostname": os.uname().nodename, "pid": process.pid,
              "affinity_cpus": sorted(os.sched_getaffinity(0)),
              "environment": {name: os.environ.get(name) for name in RESOURCE_ENV},
              "returncode": returncode, "timed_out": timed_out, "wall_s": elapsed,
              "resource": usage,
              "status": "process_passed" if returncode == 0 and usage is not None else "failed",
              "notes": ["Wall time includes process startup, numerical work, and output; this is not assembly time.",
                        "A zero exit code alone does not certify fields, convergence, or baseline acceptance.",
                        "Per-rank maxima summed across ranks are not a simultaneous aggregate memory peak."],
              "logs": {name: digest(destination / name) for name in ("stdout.log", "stderr.log")}}
    (destination / "run.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return returncode if returncode else (0 if usage is not None else 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-ranks", type=int, default=1)
    parser.add_argument("--timeout", type=float, default=600)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    options = parser.parse_args()
    command = options.command[1:] if options.command[:1] == ["--"] else options.command
    if not command or options.expected_ranks < 1 or not 0 < options.timeout < float("inf"):
        parser.error("require command, positive rank count, and finite positive timeout")
    return run(command, options.output_dir, options.expected_ranks, options.timeout)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError) as error:
        print(f"hpc_rank_run: {error}", file=sys.stderr)
        sys.exit(1)
