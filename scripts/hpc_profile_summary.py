#!/usr/bin/env python3
"""Aggregate one completed run's per-rank application phases and process RSS."""

import argparse
import json
import math
from pathlib import Path
import statistics
import sys

from hpc_inventory import digest


def distribution(values):
    if not values or any(not isinstance(value, (int, float)) or
                         not math.isfinite(value) or value < 0 for value in values):
        raise ValueError("metrics must be finite and nonnegative")
    mean = statistics.mean(values)
    return {"per_rank": values, "min": min(values), "mean": mean, "max": max(values),
            "max_over_mean": max(values) / mean if mean else None}


def summarize(directory):
    records = []
    for path in directory.glob("rank-*/run.json"):
        record = json.loads(path.read_text())
        if record.get("kind") != "hpc_rank_run" or record.get("schema_version") != 1:
            raise ValueError("unsupported rank record")
        if record["returncode"] != 0 or record["timed_out"] or record["status"] != "process_passed":
            raise ValueError("cannot summarize a failed or incomplete rank as successful")
        if path.parent.name != f"rank-{record['rank']}":
            raise ValueError("rank directory does not match its record")
        log = path.parent / "stdout.log"
        if digest(log) != record["logs"]["stdout.log"]:
            raise ValueError("stdout hash differs from the completed rank record")
        profiles = [json.loads(line[len("hpc_profile "):]) for line in log.read_text().splitlines()
                    if line.startswith("hpc_profile ")]
        if len(profiles) != 1:
            raise ValueError("expected exactly one hpc_profile record per rank; enable IGA_PROFILE=1")
        profile = profiles[0]
        if (profile["schema_version"] != 1 or profile["status"] != 0 or
                profile["rank"] != record["rank"] or profile["ranks"] != record["ranks"]):
            raise ValueError("profile identity or status differs from launcher record")
        distribution([profile["elapsed_s"], profile["unscoped_s"]])
        total = profile["unscoped_s"]
        for phase in profile["phases"].values():
            distribution([phase["inclusive_s"], phase["exclusive_s"]])
            if (not isinstance(phase["calls"], int) or phase["calls"] < 0 or
                    phase["exclusive_s"] > phase["inclusive_s"] + 1e-9):
                raise ValueError("invalid phase count or exclusive duration")
            total += phase["exclusive_s"]
        if not math.isclose(total, profile["elapsed_s"], rel_tol=1e-9, abs_tol=1e-7):
            raise ValueError("exclusive phases and unscoped time do not partition elapsed time")
        records.append((record, profile))
    records.sort(key=lambda pair: pair[0]["rank"])
    if not records:
        raise ValueError("no completed rank records")
    ranks = records[0][0]["ranks"]
    if ([record["rank"] for record, _ in records] != list(range(ranks)) or
            any(record["ranks"] != ranks for record, _ in records)):
        raise ValueError("missing or inconsistent rank records")
    names = set(records[0][1]["phases"])
    if any(set(profile["phases"]) != names for _, profile in records):
        raise ValueError("phase schemas differ between ranks")
    rss = [record["resource"]["peak_rss_bytes"] for record, _ in records]
    return {"schema_version": 1, "kind": "hpc_profile_summary", "ranks": ranks,
            "source_directory": str(directory.resolve()),
            "process_wall_s": distribution([record["wall_s"] for record, _ in records]),
            "application_elapsed_s": distribution([profile["elapsed_s"] for _, profile in records]),
            "unscoped_s": distribution([profile["unscoped_s"] for _, profile in records]),
            "peak_rss_bytes": distribution(rss), "sum_individual_peak_rss_bytes": sum(rss),
            "phases": {name: {
                metric: distribution([profile["phases"][name][metric] for _, profile in records])
                for metric in ("inclusive_s", "exclusive_s", "calls")} for name in sorted(names)},
            "notes": ["This is one run, not a repeated scaling study or numerical acceptance.",
                      "Sum of individual RSS peaks is not a simultaneous aggregate peak.",
                      "Only exclusive phases plus unscoped time partition application elapsed time.",
                      "Library-internal communication stays in the caller's phase.",
                      "Zero calls means uninstrumented or unexecuted, not proof of zero cost."]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    options = parser.parse_args()
    print(json.dumps(summarize(options.directory), indent=2, allow_nan=False))


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, TypeError) as error:
        print(f"hpc_profile_summary: {error}", file=sys.stderr)
        sys.exit(1)
