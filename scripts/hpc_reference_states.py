#!/usr/bin/env python3
"""Collect and compare complete accepted immersed/FSI fixture states."""

import argparse
import json
import math
import os
import re
from pathlib import Path
import sys

from hpc_compare_fields import compare
from hpc_cpu_matrix import ROOT, launch, write_json
from hpc_inventory import digest
from hpc_profile_summary import summarize
from hpc_serial_matrix import file_hashes, native_diagnostics


POLICY = ROOT / "benchmarks/hpc_reference_states.json"


def state_manifest(directory, policy):
    data = json.loads((directory / "manifest.json").read_text())
    if (data.get("schema_version") != 1 or data.get("kind") != "fixture_reference_state"
            or data.get("native_gates_passed") is not True or data.get("fixture") not in policy["fixtures"]):
        raise ValueError("invalid or incomplete fixture state manifest")
    expected = policy["fixtures"][data["fixture"]]
    if data["time_s"] != expected["time_s"] or data["step"] != expected["step"]:
        raise ValueError("reference state epoch differs from registered fixture")
    fields = data["fields"]
    if len(fields) != len(expected["fields"]) or {f["name"] for f in fields} != set(expected["fields"]):
        raise ValueError("reference state has missing, duplicate, or unexpected fields")
    for field in fields:
        schema = expected["fields"][field["name"]]
        if (field["file"] != field["name"] + ".txt" or field["units"] != schema["units"]
                or field["columns"] != schema["columns"] or type(field["rows"]) is not int or field["rows"] <= 0):
            raise ValueError("invalid reference field schema")
        path = directory / field["file"]
        if not path.resolve().is_relative_to(directory.resolve()) or digest(path) != field["sha256"]:
            raise ValueError("reference field hash or path differs")
        shape = compare(path, path, policy["relative_l2_max"], policy["zero_reference_absolute_l2_max"], True)
        if shape["rows"] != field["rows"] or shape["columns"] != field["columns"]:
            raise ValueError("reference field dimensions differ from manifest")
    return data


def compare_states(reference, candidate, policy):
    first = state_manifest(reference, policy)
    second = state_manifest(candidate, policy)
    if first["fixture"] != second["fixture"]:
        raise ValueError("cannot compare different fixture models")
    results = {}
    for field in first["fields"]:
        name = field["name"]
        results[name] = compare(reference / field["file"], candidate / field["file"],
                                policy["relative_l2_max"], policy["zero_reference_absolute_l2_max"], True)
    return {"fixture": first["fixture"], "passed": all(r["passed"] for r in results.values()), "fields": results}


def execution_configuration(options):
    cpus = [options.cpu] if options.cpus is None else options.cpus
    if (not math.isfinite(options.timeout) or options.timeout <= 0
            or options.threads < 1 or not cpus or len(set(cpus)) != len(cpus)
            or any(cpu not in os.sched_getaffinity(0) for cpu in cpus)
            or options.threads > len(cpus)):
        raise ValueError("require positive timeout, available unique CPUs, and at least one CPU per thread")
    environment = {"OMP_NUM_THREADS": str(options.threads), "OMP_THREAD_LIMIT": str(options.threads),
                   "OMP_DYNAMIC": "FALSE", "OMP_PROC_BIND": "close", "OMP_PLACES": "threads",
                   "IGA_ASSEMBLY_THREADS": str(options.threads), "IGA_ASSEMBLY_BATCH_SIZE": str(min(options.threads, 8)),
                   "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "BLIS_NUM_THREADS": "1",
                   "PETSC_OPTIONS": "", "IGA_PROFILE": "1"}
    return cpus, environment



def assembly_diagnostics(log, threads):
    lines = [line for line in log.splitlines() if line.startswith("element_assembly ")]
    rows = []
    pattern = (r"element_assembly threads_requested=(\d+) team_size=(\d+) cells=(\d+) "
               r"batches=(\d+) maximum_resident_items=(\d+)")
    for line in lines:
        match = re.fullmatch(pattern, line)
        if match is None:
            raise ValueError("invalid element assembly diagnostic")
        requested, team, cells, batches, resident = map(int, match.groups())
        if (requested != threads or team != threads or cells < 1 or batches < 1
                or resident < 1 or resident > min(threads, 8) or resident > cells):
            raise ValueError("actual assembly team or batch capacity differs from requested execution")
        rows.append(dict(threads_requested=requested, team_size=team, cells=cells,
                         batches=batches, maximum_resident_items=resident))
    if threads > 1 and not rows:
        raise ValueError("multithread collection lacks actual element assembly evidence")
    return rows


def collect(options):
    cpus, settings = execution_configuration(options)
    policy = json.loads(POLICY.read_text())
    output = options.output_dir.resolve()
    test = ROOT / "solvers/cpu/tests" / ("test_compliant_channel_fsi.cpp" if options.case == "fsi"
                                        else "test_immersed_aneurysm_jacobian.cpp")
    binary = ROOT / "solvers/cpu" / ("compliant_channel_fsi_test" if options.case == "fsi"
                                    else "immersed_aneurysm_jacobian_test")
    if options.binary is not None:
        binary = options.binary.resolve(strict=True)
    inputs = [test]
    command = [str(binary)]
    if options.case == "immersed":
        if options.case_dir is None:
            raise ValueError("immersed reference requires an explicit depth-2 case directory")
        case = options.case_dir.resolve(strict=True)
        if output.is_relative_to(case) or case.is_relative_to(output):
            raise ValueError("input and output directories must be separate")
        geometry = json.loads((case / "immersed_geometry.json").read_text())
        if geometry["volume_quadrature"]["max_depth"] != 2:
            raise ValueError("selected immersed reference requires depth 2")
        inputs += [p for p in case.rglob("*") if p.is_file()]
        command.append(str(case))
    else:
        if options.case_dir is not None:
            raise ValueError("FSI reference uses the in-memory fixture")
        inputs.append(ROOT / "solvers/cpu/tests/CompliantChannelFsiFixture.hpp")
    # Input paths are retained for provenance, while logical names allow an
    # identical immutable case to be staged at a different filesystem path.
    input_identity = {"fixture_source/" + test.name: digest(test)}
    for path in inputs[1:]:
        name = str(path.relative_to(case)) if options.case == "immersed" else path.name
        input_identity[name] = digest(path)
    sources = [POLICY, binary, ROOT / "solvers/cpu/tests/ReferenceStateOutput.hpp"] + inputs
    sources += list((ROOT / "include").glob("*.hpp")) + list((ROOT / "solvers/cpu/include").glob("*.hpp"))
    sources += [ROOT / "scripts" / name for name in ("hpc_reference_states.py", "hpc_compare_fields.py",
                "hpc_cpu_matrix.py", "hpc_inventory.py", "hpc_profile_summary.py", "hpc_serial_matrix.py", "hpc_rank_run.py")]
    before = file_hashes(sources)
    environment = os.environ.copy()
    environment.update(settings)
    output.mkdir(parents=True, exist_ok=False)
    record = {"schema_version": 1, "kind": "fixture_reference_collection", "fixture": options.case,
              "input_identity": input_identity, "source_binary_sha256": before,
              "execution": {"cpus": cpus, "threads": options.threads, "environment": settings},
              "policy_sha256": digest(POLICY), "status": "failed",
              "note": "Correctness collection with optional state output; not an isolated timing matrix."}
    write_json(output / "inputs.json", record)
    command += ["--reference-output", str(output / "state")]
    argv = ["taskset", "--cpu-list", ",".join(map(str, cpus)), sys.executable, str(ROOT / "scripts/hpc_rank_run.py"),
            "--output-dir", str(output), "--expected-ranks", "1", "--timeout", str(options.timeout), "--"] + command
    try:
        record["launch"] = launch(argv, output, environment, options.timeout + 30)
        if record["launch"]["returncode"] != 0:
            raise ValueError("native reference fixture failed")
        record["profile"] = summarize(output)
        native_log = (output / "rank-0/stdout.log").read_text()
        record["native_diagnostics"] = native_diagnostics(options.case, native_log)
        record["element_assembly"] = assembly_diagnostics(native_log, options.threads)
        record["state"] = state_manifest(output / "state", policy)
        record["manifest_sha256"] = digest(output / "state/manifest.json")
        if file_hashes(sources) != before:
            raise ValueError("reference inputs, code or binary changed during collection")
        record["status"] = "accepted"
    except (OSError, ValueError, KeyError, TypeError) as error:
        record["error"] = str(error)
    write_json(output / "collection.json", record)
    print(json.dumps({"directory": str(output), "status": record["status"], "error": record.get("error")}), flush=True)
    return 0 if record["status"] == "accepted" else 1


def accepted_collection(directory, policy):
    record = json.loads((directory / "collection.json").read_text())
    if (record.get("kind") != "fixture_reference_collection" or record.get("status") != "accepted"
            or record["policy_sha256"] != digest(POLICY)):
        raise ValueError("require accepted collection under the current reference policy")
    if summarize(directory) != record["profile"]:
        raise ValueError("recorded native profile differs")
    diagnostics = native_diagnostics(record["fixture"], (directory / "rank-0/stdout.log").read_text())
    if "element_assembly" in record:
        observed = assembly_diagnostics((directory / "rank-0/stdout.log").read_text(), record["execution"]["threads"])
        if observed != record["element_assembly"]:
            raise ValueError("recorded assembly execution differs")
    if diagnostics != record["native_diagnostics"]:
        raise ValueError("native acceptance diagnostics differ")
    if digest(directory / "state/manifest.json") != record["manifest_sha256"]:
        raise ValueError("accepted state manifest changed")
    if state_manifest(directory / "state", policy) != record["state"]:
        raise ValueError("accepted state metadata changed")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("collect")
    run.add_argument("--case", choices=("immersed", "fsi"), required=True)
    run.add_argument("--case-dir", type=Path)
    run.add_argument("--output-dir", type=Path, required=True)
    run.add_argument("--cpu", type=int, default=0)
    run.add_argument("--timeout", type=float, default=1500)
    run.add_argument("--binary", type=Path, help="fixture-compatible executable with identical native gates and state export")
    run.add_argument("--threads", type=int, default=1)
    run.add_argument("--cpus", type=int, nargs="+", help="explicit CPU IDs; overrides --cpu")
    check = commands.add_parser("compare")
    check.add_argument("reference", type=Path)
    check.add_argument("candidate", type=Path)
    options = parser.parse_args()
    if options.command == "collect":
        return collect(options)
    policy = json.loads(POLICY.read_text())
    reference = options.reference.resolve(strict=True)
    candidate = options.candidate.resolve(strict=True)
    first = accepted_collection(reference, policy)
    second = accepted_collection(candidate, policy)
    if first["input_identity"] != second["input_identity"]:
        raise ValueError("fixture definitions or physical inputs differ")
    result = compare_states(reference / "state", candidate / "state", policy)
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0 if result["passed"] else 2


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError, TypeError) as error:
        print(f"hpc_reference_states: {error}", file=sys.stderr)
        sys.exit(1)
