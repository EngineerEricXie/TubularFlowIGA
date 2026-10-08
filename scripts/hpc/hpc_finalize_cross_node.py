#!/usr/bin/env python3
"""Consolidate cross-node graph, FSI, test-tier, and scaling evidence."""

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path


def require(value, message):
    if not value:
        raise ValueError(message)


def read(path, kind=None):
    require(path.is_file(), f"missing evidence: {path}")
    value = json.loads(path.read_text())
    if kind is not None:
        require(value.get("kind") == kind, f"unexpected evidence kind: {path}")
    return value


def integer(value, name):
    try:
        result = int(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"invalid {name}") from error
    require(result > 0, f"invalid {name}")
    return result


def build_revision(path):
    value = read(path, "hpc_build_manifest")
    require(value.get("dependency_status") == "passed", f"dependencies did not pass: {path}")
    commit = value["source"]["commit"]
    status = value["source"]["status_porcelain"]
    require(commit.get("returncode") == 0 and commit.get("output"), f"missing source revision: {path}")
    require(status.get("returncode") == 0 and not status.get("output"),
        f"source tree was not clean: {path}")
    require(value.get("binaries") and all(item.get("exists") for item in value["binaries"]),
        f"missing recorded binary: {path}")
    return commit["output"]


def hosts(directory):
    values = set()
    for path in directory.glob("rank-*/run.json"):
        record = read(path, "hpc_rank_run")
        require(record.get("status") == "process_passed" and record.get("returncode") == 0
            and not record.get("timed_out"), f"invalid rank run: {path}")
        values.add(record.get("hostname"))
    values.discard(None)
    return values


def graph_evidence(root):
    attempts = []
    for path in sorted(root.glob("attempt-*/scheduler.json")):
        scheduler = read(path, "hpc_scheduler_attempt")
        nodes = integer(scheduler.get("slurm", {}).get("SLURM_JOB_NUM_NODES"), "graph node count")
        require(nodes >= 2 and scheduler.get("slurm", {}).get("SLURM_JOB_ID"),
            "graph attempt is not from a cross-node Slurm job")
        attempt = integer(int(path.parent.name.removeprefix("attempt-")) + 1,
            "graph restart count") - 1
        require(int(scheduler["slurm"].get("SLURM_RESTART_COUNT") or 0) == attempt,
            "graph attempt directory differs from Slurm restart count")
        attempts.append((attempt, path.parent, scheduler))
    require(attempts, "missing graph attempts")
    require(all(item[2].get("status") in ("checkpointed", "passed") for item in attempts),
        "graph evidence contains an unsuccessful attempt")
    for _, directory, _ in attempts:
        node_file = directory / "nodes.txt"
        require(node_file.is_file() and len(set(node_file.read_text().splitlines())) >= 2,
            "graph attempt lacks two distinct allocated hosts")
    checkpointed = [item for item in attempts if item[2].get("status") == "checkpointed"]
    passed = [item for item in attempts if item[2].get("status") == "passed"]
    require(checkpointed and passed, "graph evidence requires checkpointed and passed attempts")
    first = min(checkpointed, key=lambda item: item[0])
    last = max(passed, key=lambda item: item[0])
    require(first[0] < last[0] and first[2].get("checkpoint_requested"),
        "graph did not checkpoint and resume in a later allocation attempt")
    require("--restart-dir" in last[2].get("command", []),
        "completed graph attempt did not load a checkpoint")
    require((last[1] / "results/graph_binding_manifest.json").is_file(),
        "resumed graph lacks its completion manifest")
    revisions = {build_revision(item[1] / "build.json") for item in attempts}
    require(len(revisions) == 1, "graph attempts used different revisions")
    return {"revision": revisions.pop(), "attempts": len(attempts),
        "checkpoint_attempt": first[0], "completion_attempt": last[0],
        "job_ids": [item[2]["slurm"]["SLURM_JOB_ID"] for item in attempts]}


def fsi_evidence(root):
    scheduler = read(root / "scheduler.json", "hpc_scheduler_attempt")
    require(scheduler.get("status") == "passed" and scheduler.get("returncode") == 0,
        "cross-node FSI scheduler attempt did not pass")
    require(integer(scheduler.get("slurm", {}).get("SLURM_JOB_NUM_NODES"), "FSI node count") >= 2,
        "FSI evidence is not cross-node")
    tier = read(root / "scheduled/result.json", "hpc_test_tier")
    require(tier.get("status") == "passed" and tier.get("tier") == "scheduled"
        and tier.get("commands"), "scheduled tier was skipped or failed")
    strong = read(root / "strong-fsi.json")
    restart = read(root / "paired-restart.json")
    require(strong.get("status") == "passed", "strong FSI comparison did not pass")
    require(restart.get("status") == "passed" and restart.get("source_ranks") == 4
        and restart.get("target_ranks") == 2, "paired 4-to-2 FSI restart did not pass")
    require(len(hosts(root / "writer-4")) >= 2 and len(hosts(root / "reader-2")) >= 2,
        "FSI writer and restart reader must each span two hosts")
    revision = build_revision(root / "build.json")
    require(tier.get("source_commit") == revision, "scheduled tier used a different revision")
    return {"revision": revision, "job_id": scheduler["slurm"]["SLURM_JOB_ID"],
        "writer_hosts": sorted(hosts(root / "writer-4")),
        "reader_hosts": sorted(hosts(root / "reader-2"))}


def scaling_evidence(root):
    report = read(root / "summary.json", "hpc_cross_node_scaling")
    require(report.get("status") == "passed", "scaling collector did not pass")
    require(not report.get("source_status"), "scaling used a dirty source tree")
    require(integer(report.get("allocation", {}).get("SLURM_JOB_NUM_NODES"),
        "scaling node count") >= 2, "scaling evidence is not cross-node")
    require(report.get("allocation", {}).get("SLURM_JOB_ID"),
        "scaling evidence lacks a Slurm job ID")
    ranks = report.get("ranks", [])
    require(ranks and ranks[0] == 1 and max(ranks) > 128 and report.get("repetitions", 0) >= 3,
        "scaling must include 1 and more than 128 ranks with three repetitions")
    strong_elements = {report["inputs"][f"strong-{count}"]["elements"] for count in ranks}
    require(len(strong_elements) == 1 and next(iter(strong_elements)) >= 16384,
        "strong scaling problem must be fixed and contain at least 16,384 elements")
    weak_work = [report["inputs"][f"weak-{count}"]["elements_per_rank"] for count in ranks]
    require(min(weak_work) >= 256 and max(weak_work) / min(weak_work) - 1 <= 0.15,
        "weak scaling must retain at least 256 elements per rank within 15 percent")
    expected_runs = 2 * len(ranks) * report["repetitions"]
    require(len(report.get("runs", [])) == expected_runs,
        "scaling repetition matrix is incomplete")
    for run in report["runs"]:
        require(run.get("status") == "passed" and
            run.get("physical_validation", {}).get("returncode") == 0,
            "a scaling run or physical validation did not pass")
        if run.get("mode") == "strong":
            require(all(value.get("passed") for value in run.get("field_comparison", {}).values()),
                "strong-scaling field comparison did not pass")
    for mode in ("strong", "weak"):
        for count in ranks:
            item = report.get("summary", {}).get(mode, {}).get(str(count), {})
            require(all(name in item for name in ("elements_per_rank", "max_rank_process_wall_s",
                "max_rank_peak_rss_bytes", "solver_iterations", "speedup_vs_one_rank",
                "parallel_efficiency")), "scaling summary lacks a required metric")
            require("communication" in item.get("phases_max_rank_exclusive_s", {}),
                "scaling summary lacks communication timing")
    large_runs = [Path(run["directory"]) for run in report["runs"] if run["ranks"] == max(ranks)]
    large_hosts = set().union(*(hosts(path) for path in large_runs))
    require(len(large_hosts) >= 2, "largest scaling runs do not contain two hostnames")
    revision = build_revision(root.with_name(root.name + ".build.json"))
    require(report.get("source_commit") == revision, "scaling report and build revision differ")
    return {"revision": revision, "job_id": report["allocation"]["SLURM_JOB_ID"],
        "ranks": ranks, "repetitions": report["repetitions"],
        "largest_run_hosts": sorted(large_hosts), "summary": report["summary"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph-output", type=Path, required=True)
    parser.add_argument("--fsi-output", type=Path, required=True)
    parser.add_argument("--scaling-output", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = {"schema_version": 1, "kind": "hpc_cross_node_acceptance",
        "created_utc": datetime.now(timezone.utc).isoformat(), "status": "running"}
    try:
        report["graph"] = graph_evidence(args.graph_output.resolve())
        report["fsi"] = fsi_evidence(args.fsi_output.resolve())
        report["scaling"] = scaling_evidence(args.scaling_output.resolve())
        revisions = {report[name]["revision"] for name in ("graph", "fsi", "scaling")}
        require(len(revisions) == 1, "cross-node workflows used different source revisions")
        report["source_commit"] = revisions.pop()
        report["completed_items"] = ["HPC-05C", "HPC-07C", "HPC-07D",
            "HPC-09A", "HPC-09B", "HPC-09C", "HPC-09D"]
        report["status"] = "passed"
    except (KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
        report["status"] = "failed"
        report["failure"] = str(error)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    if report["status"] != "passed":
        raise SystemExit(report["failure"])
    print(f"cross-node acceptance passed for {report['source_commit']}")


if __name__ == "__main__":
    main()
