#!/usr/bin/env python3
"""Validate complete baseline transport histories with independent species budgets."""

import argparse
import json
import os
from pathlib import Path
import sys

from hpc_compare_fields import compare
from hpc_cpu_matrix import ROOT, launch, snapshot, write_json
from hpc_profile_summary import summarize
from hpc_serial_matrix import cpu_reference


def run(options):
    metadata, prior_reference, _ = cpu_reference(options.cpu_matrix.resolve(strict=True), "neuron_transport")
    case = Path(metadata["case_dir"])
    configuration = json.loads((case / "simulation_config.json").read_text())
    if configuration["time"]["steps"] != 2:
        raise ValueError("this regression requires the catalog's two-step transport fixture")
    systems = [s for s in configuration["equation_systems"] if s["name"] == "neuron_transport"]
    if len(systems) != 1 or systems[0]["unknowns"] != ["N0", "Nplus"]:
        raise ValueError("this regression requires the catalog's ordered N0/Nplus system")
    ranks = sorted(set(options.ranks))
    if not ranks or ranks[0] != 1 or len(ranks) != len(options.ranks):
        raise ValueError("require unique positive ranks including 1")
    if any(n not in metadata["ranks"] for n in ranks):
        raise ValueError("requested ranks were not covered by the accepted CPU reference matrix")
    policy_path = ROOT / "benchmarks/hpc_transport_budget.json"
    policy = json.loads(policy_path.read_text())
    comparison_policy = metadata["comparison_policy"]
    validator = ROOT / "solvers/cpu/iga_transport_validate"
    cpu = ROOT / "solvers/cpu/iga_solve"
    gpu = ROOT / "solvers/cuda/iga_cuda"
    output = options.output_dir.resolve()
    if output.is_relative_to(case) or case.is_relative_to(output):
        raise ValueError("case and output directories must be separate")
    provenance = [validator, cpu, policy_path, Path(__file__).resolve(),
                  ROOT / "solvers/cpu/include/TransportBudget.hpp", ROOT / "solvers/cpu/src/iga_transport_validate.cpp",
                  ROOT / "solvers/cpu/tests/test_transport_budget.cpp", ROOT / "solvers/cpu/Makefile",
                  options.cpu_matrix.resolve() / "matrix.json", options.cpu_matrix.resolve() / "summary.json"]
    # The independent integrator still shares geometry and configuration
    # primitives; record those sources as well as its own implementation.
    provenance += list((ROOT / "solvers/cpu/include").glob("*.hpp"))
    provenance += list((ROOT / "include").glob("*.hpp"))
    provenance += [ROOT / "scripts" / name for name in
                   ("hpc_compare_fields.py", "hpc_cpu_matrix.py", "hpc_serial_matrix.py",
                    "hpc_profile_summary.py", "hpc_rank_run.py", "hpc_inventory.py")]
    if options.include_cuda:
        provenance.append(gpu)
    before = snapshot(case, provenance)
    environment = os.environ.copy()
    environment.update({"OMP_NUM_THREADS": "1", "OMP_THREAD_LIMIT": "1", "OMP_DYNAMIC": "FALSE",
                        "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "BLIS_NUM_THREADS": "1",
                        "IGA_PROFILE": "1", "PETSC_OPTIONS": ""})
    output.mkdir(parents=True, exist_ok=False)
    write_json(output / "inputs.json", {"schema_version": 1, "policy": policy,
                                       "case": str(case), "ranks": ranks, "include_cuda": options.include_cuda,
                                       "input_binary_and_source_sha256": before,
                                       "note": "Validation runs include per-step output; these are not timing comparisons with the earlier matrices."})
    observations = []
    negative = []
    error = None
    reference = output / "cpu1"
    try:
        modes = [(f"cpu{n}", n, False) for n in ranks]
        if options.include_cuda:
            modes.append(("cuda", 1, True))
        for name, n, is_cuda in modes:
            directory = output / name
            directory.mkdir()
            database = case / f"straight_neurite-{n}.ntiga"
            command = ([str(gpu), "solve"] if is_cuda else [str(cpu)]) + [
                str(database), str(case), "--system", "neuron_transport", "--output", str(directory / "field.txt"),
                "--output-every", "1"]
            prefix = [] if is_cuda else ["mpiexec", "--bind-to", "core", "--map-by", "core", "--nooversubscribe", "-np", str(n)]
            argv = prefix + [sys.executable, str(ROOT / "scripts/hpc_rank_run.py"),
                             "--output-dir", str(directory), "--expected-ranks", str(n), "--timeout", "120", "--"] + command
            record = {"mode": name, "status": "failed", "budgets": [], "comparisons": []}
            observations.append(record)
            print("running transport history " + name, flush=True)
            record["launch"] = launch(argv, directory, environment, 150)
            if record["launch"]["returncode"] != 0:
                raise ValueError("history generation failed: " + name)
            record["profile"] = summarize(directory)
            rtol = comparison_policy["cpu_cuda_field_relative_l2_max" if is_cuda else "cpu_rank_field_relative_l2_max"]
            for step in range(3):
                field = f"field.step{step:06d}.txt"
                record["comparisons"].append(compare(reference / field, directory / field, rtol,
                                                       comparison_policy["zero_reference_absolute_l2_max"], True))
                if step == 0:
                    continue
                validation_dir = directory / f"budget-step{step}"
                validation_dir.mkdir()
                command = [str(validator), str(database), str(case), "neuron_transport",
                           str(directory / f"field.step{step-1:06d}.txt"), str(directory / field),
                           "--rtol", str(policy["relative_tolerance"]), "--zero-atol", str(policy["zero_reference_absolute_tolerance"])]
                launched = launch(command, validation_dir, environment, 120)
                result = json.loads((validation_dir / "launcher.stdout").read_text()) if launched["returncode"] in (0, 2) else None
                passed = (launched["returncode"] == 0 and result is not None and result["accepted"]
                          and abs(result["net_linear_reaction"]) <= policy["zero_reference_absolute_tolerance"])
                record["budgets"].append({"step": step, "launch": launched, "result": result, "passed": bool(passed)})
                print(f"{name} step {step}: budget {'passed' if passed else 'FAILED'}", flush=True)
            record["prior_reference_comparison"] = compare(prior_reference / "field.txt", directory / "field.txt", rtol,
                                                             comparison_policy["zero_reference_absolute_l2_max"], True)
            record["status"] = "accepted" if (all(b["passed"] for b in record["budgets"])
                                                and all(c["passed"] for c in record["comparisons"])
                                                and record["prior_reference_comparison"]["passed"]) else "failed"
        # Exercise actual validator error paths after numerical observations;
        # never replace the original fields or relax failed budget thresholds.
        original = (reference / "field.step000001.txt").read_text().splitlines()
        tokens = (case / "controlmesh.vtk").read_text().split()
        start = tokens.index("POINT_DATA")
        nodes = int(tokens[start + 1])
        labels_start = tokens.index("LOOKUP_TABLE", start) + 2
        labels = [int(float(x)) for x in tokens[labels_start:labels_start+nodes]]
        constrained_labels = {b["label"] for b in configuration["boundaries"]
                              if any(c["field"] == "N0" and c["type"] == "dirichlet" for c in b["conditions"])}
        free_node = next((i for i, label in enumerate(labels) if label not in constrained_labels), None)
        if free_node is None:
            raise ValueError("baseline requires a free N0 coefficient for the corruption check")
        corrupt = original.copy()
        values = corrupt[free_node].split()
        values[1] = format(float(values[1]) + 0.1, ".17g")
        corrupt[free_node] = " ".join(values)
        bad_id = original.copy()
        bad_id[0] = "0.0 " + " ".join(bad_id[0].split()[1:])
        nonfinite = original.copy()
        nonfinite[0] = "0 nan " + " ".join(nonfinite[0].split()[2:])
        extra = original.copy()
        extra[0] += " 99"
        cases = [("free-field-corruption", corrupt, 2), ("fractional-node-id", bad_id, 1),
                 ("truncated-field", original[:-1], 1), ("nonfinite-field", nonfinite, 1), ("extra-column", extra, 1)]
        for name, lines, expected in cases:
            directory = output / name
            directory.mkdir()
            field = directory / "field.txt"
            field.write_text("\n".join(lines) + "\n")
            command = [str(validator), str(case / "straight_neurite-1.ntiga"), str(case), "neuron_transport",
                       str(reference / "field.step000000.txt"), str(field), "--rtol", str(policy["relative_tolerance"]),
                       "--zero-atol", str(policy["zero_reference_absolute_tolerance"])]
            launched = launch(command, directory, environment, 120)
            check = {"case": name, "expected_exit": expected, "launch": launched,
                     "passed": launched["returncode"] == expected}
            if name == "free-field-corruption" and launched["returncode"] == 2:
                result = json.loads((directory / "launcher.stdout").read_text())
                check["free_residual_rejected"] = any(not f["gates"]["free_residual"] for f in result["fields"])
                check["passed"] = check["passed"] and check["free_residual_rejected"]
            negative.append(check)
        if snapshot(case, provenance) != before:
            raise ValueError("budget inputs or binaries changed during validation")
    except (OSError, ValueError, KeyError, TypeError) as exception:
        error = str(exception)
    accepted = (error is None and bool(observations) and len(negative) == 5
                and all(r["status"] == "accepted" for r in observations) and all(r["passed"] for r in negative))
    write_json(output / "summary.json", {"schema_version": 1, "kind": "hpc_transport_budget_regression",
                                        "status": "requested_checks_passed" if accepted else "failed", "error": error,
                                        "include_cuda": options.include_cuda, "observations": observations, "negative_tests": negative})
    return 0 if accepted else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cpu-matrix", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--ranks", type=int, nargs="+", default=[1, 2, 4, 8])
    parser.add_argument("--include-cuda", action="store_true")
    return run(parser.parse_args())


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError, TypeError) as error:
        print(f"hpc_transport_budget_regression: {error}", file=sys.stderr)
        sys.exit(1)
