#!/usr/bin/env python3
"""Run sequential fresh-process FSI OpenMP repetitions with complete field gates."""

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import sys
from types import SimpleNamespace

from hpc_cpu_matrix import ROOT, write_json
from hpc_reference_states import POLICY, accepted_collection, collect, compare_states, execution_configuration
from hpc_serial_matrix import aggregate_serial, file_hashes


def publish(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    write_json(temporary, value)
    temporary.replace(path)


def aggregate_openmp(records, threads, repetitions, memory_limit):
    if (not threads or threads[0] != 1 or len(set(threads)) != len(threads)
            or len(records) != len(threads) * (repetitions + 1)
            or {r["threads"] for r in records} != set(threads)
            or any(not r["field_comparison"]["passed"] for r in records)):
        raise ValueError("incomplete OpenMP matrix or failed field comparison")
    modes = {str(count): aggregate_serial([r for r in records if r["threads"] == count], repetitions, False)
             for count in threads}
    reference = modes["1"]
    if reference["launcher_wall_s"]["median"] <= 0 or reference["peak_rss_bytes"]["max"] <= 0:
        raise ValueError("reference time and memory must be positive")
    for count in threads:
        mode = modes[str(count)]
        if mode["launcher_wall_s"]["median"] <= 0:
            raise ValueError("candidate time must be positive")
        ratio = mode["launcher_wall_s"]["median"] / reference["launcher_wall_s"]["median"]
        memory_ratio = mode["peak_rss_bytes"]["max"] / reference["peak_rss_bytes"]["max"]
        mode.update(observed_speedup=1 / ratio, parallel_efficiency=1 / ratio / count,
                    peak_rss_ratio_to_one_thread=memory_ratio,
                    memory_gate_passed=memory_ratio <= memory_limit,
                    # Matches the >=10% pre-optimization hpc_baselines policy.
                    speedup_claim_allowed=count > 1 and ratio <= .90 and memory_ratio <= memory_limit)
    return modes


def run(options):
    if (options.repetitions < 3 or not options.threads or options.threads[0] != 1
            or len(set(options.threads)) != len(options.threads)
            or not math.isfinite(options.timeout) or options.timeout <= 0):
        raise ValueError("require unique thread configurations starting at 1, at least three repetitions and positive timeout")
    policy = json.loads(POLICY.read_text())
    catalog = ROOT / "benchmarks/hpc_baselines.json"
    comparison_policy = json.loads(catalog.read_text())["comparison_policy"]
    reference = options.reference.resolve(strict=True)
    baseline = accepted_collection(reference, policy)
    if baseline["fixture"] != "fsi":
        raise ValueError("this OpenMP matrix requires the complete FSI reference")
    binary = options.binary.resolve(strict=True)
    output = options.output_dir.resolve()
    configurations = []
    for count in options.threads:
        settings = SimpleNamespace(cpu=options.cpus[0], cpus=options.cpus[:count],
                                   threads=count, timeout=options.timeout)
        cpus, environment = execution_configuration(settings)
        configurations.append(dict(threads=count, cpus=cpus, environment=environment))
    paths = [Path(__file__).resolve(), catalog, POLICY, binary,
             reference / "collection.json", reference / "state/manifest.json"]
    before = file_hashes(paths)
    output.mkdir(parents=True, exist_ok=False)
    metadata = dict(schema_version=1, kind="fsi_openmp_matrix", status="running",
                    created_utc=datetime.now(timezone.utc).isoformat(),
                    reference=str(reference), repetitions=options.repetitions,
                    configurations=configurations, input_binary_tool_sha256=before,
                    comparison_policy=comparison_policy, runs=[],
                    note="Runs are serialized by this controller. Use otherwise idle allocated CPUs; external workload isolation is the operator's responsibility. Repetition 0 is excluded; it is not a controlled cold-cache experiment. Full state output is enabled in every run.")
    publish(output / "matrix.json", metadata)
    try:
        # Alternate configurations between repetitions to reduce order bias.
        for repetition in range(options.repetitions + 1):
            for config in configurations:
                directory = output / ("threads-" + str(config["threads"])) / ("repetition-" + str(repetition))
                settings = SimpleNamespace(case="fsi", case_dir=None, output_dir=directory,
                                           cpu=config["cpus"][0], cpus=config["cpus"], threads=config["threads"],
                                           binary=binary, timeout=options.timeout)
                if collect(settings):
                    raise ValueError("native collection failed: " + str(directory))
                record = accepted_collection(directory, policy)
                if "collection_source_sha256" not in metadata:
                    metadata["collection_source_sha256"] = record["source_binary_sha256"]
                if record["source_binary_sha256"] != metadata["collection_source_sha256"]:
                    raise ValueError("collection source snapshot changed between repetitions")
                if record["input_identity"] != baseline["input_identity"]:
                    raise ValueError("FSI physical inputs differ from reference")
                comparison = compare_states(reference / "state", directory / "state", policy)
                write_json(directory / "field-comparison.json", comparison)
                record.update(directory=str(directory), repetition=repetition,
                              threads=config["threads"], field_comparison=comparison)
                metadata["runs"].append(record)
                publish(output / "matrix.json", metadata)
                if not comparison["passed"]:
                    raise ValueError("complete FSI field comparison failed")
                if file_hashes(paths) != before:
                    raise ValueError("matrix driver, policy, reference or executable changed")
        metadata["summary"] = aggregate_openmp(metadata["runs"], options.threads, options.repetitions,
                                               comparison_policy["memory_regression_limit_ratio"])
        metadata["status"] = "accepted"
    except (OSError, ValueError, KeyError, TypeError) as error:
        metadata["status"] = "failed"
        metadata["error"] = str(error)
    publish(output / "matrix.json", metadata)
    print(json.dumps(dict(directory=str(output), status=metadata["status"], error=metadata.get("error"))), flush=True)
    return 0 if metadata["status"] == "accepted" else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--binary", type=Path, default=ROOT / "solvers/cpu/compliant_channel_fsi_openmp_test")
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 4])
    parser.add_argument("--cpus", type=int, nargs="+", required=True)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=1800)
    return run(parser.parse_args())


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError, TypeError) as error:
        print("hpc_openmp_matrix: " + str(error), file=sys.stderr)
        sys.exit(1)
