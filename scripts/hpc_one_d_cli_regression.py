#!/usr/bin/env python3
"""Validate native 1D CLI rank-local faults and serial/MPI field agreement."""

import argparse
import csv
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

from hpc_inventory import digest


def compare(reference, current, start_time=None, include_species=False):
    errors = {}
    filenames = ["profile_1d.csv", "branch_timeseries.csv"]
    if include_species:
        filenames += ["species_profile_1d.csv", "derived_profile_1d.csv"]
    for filename in filenames:
        with (reference / filename).open() as first, (current / filename).open() as second:
            a, b = csv.DictReader(first), csv.DictReader(second)
            if a.fieldnames != b.fieldnames:
                raise RuntimeError("field columns differ")
            columns = a.fieldnames
            a, b = list(a), list(b)
        if start_time is not None:
            # Restart metadata uses accumulated runtime time; regular output
            # uses step*dt. Match the existing CLI ValidateRestart clock gate.
            a = [row for row in a if float(row["time"]) >= start_time
                 or abs(float(row["time"])-start_time) <= 1e-12*max(1.0, abs(start_time))]
        if not a and not b and filename in ("species_profile_1d.csv", "derived_profile_1d.csv"):
            continue
        if not a or len(a) != len(b):
            raise RuntimeError("field row counts differ")
        for key in columns:
            if key in ("species", "field"):
                if [row[key] for row in a] != [row[key] for row in b]:
                    raise RuntimeError("species/derived field identities differ")
                continue
            av, bv = [float(row[key]) for row in a], [float(row[key]) for row in b]
            if not all(math.isfinite(value) for value in av + bv):
                raise RuntimeError("nonfinite output")
            if key == "time" and start_time is not None:
                if any(abs(x-y) > 1e-12*max(1.0, abs(x), abs(y)) for x, y in zip(av, bv)):
                    raise RuntimeError("restart clocks differ")
            elif key in ("time", "parent_id", "child_id", "cell", "x") and av != bv:
                raise RuntimeError("field identities differ")
            if key == "area" and min(av + bv) <= 0:
                raise RuntimeError("nonpositive area")
            norm = math.hypot(*av)
            error = math.hypot(*(x-y for x, y in zip(av, bv)))
            exceeds_tolerance = error / norm > 1e-6 if norm else error > 1e-12
            if exceeds_tolerance:
                raise RuntimeError(f"{filename}:{key} exceeds existing CPU tolerance")
            errors[f"{filename}:{key}"] = error / norm if norm else error
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--launcher", default="mpiexec --oversubscribe")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    directory = args.output_dir.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    binary = repo / "solvers/one_d/iga_1d"
    wrapper = repo / "scripts/hpc_failure_regression.py"
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               PETSC_OPTIONS="-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps")
    cases = [(f"{kind}-{ranks}", kind, ranks, "")
             for kind in ("rigid", "explicit", "implicit", "species") for ranks in (1, 3)]
    cases += [(name, "implicit", 3, stage) for name, stage in (
        ("arguments", "1d arguments: rank 1:"),
        ("check-control", "1d execution controls: rank 1:"),
        ("stop-control", "1d execution controls: rank 1:"),
        ("config-missing", "1d configuration input: rank 1:"),
        ("config-different", "1d configuration agreement: rank 1:"),
        ("geometry-missing", "asset content read: rank 1:"),
        ("petsc-flag", "PETSc option agreement: rank 1:"),
        ("output-setup", "1d output setup: rank 0:"),
        ("output-initial", "1d initial output: rank 0:"),
        ("output-step", "1d step output: rank 0:"),
        ("output-final", "1d final output: rank 0:"),
        ("output-buffered", "1d initial output: rank 0:"))]
    cases.append(("step-inlet-assets", "species", 3, "1d asset temporal oxygen_table: rank 1:"))
    cases.append(("step-inlet", "species", 3, "1d step input: rank 0:"))
    summary = {"status": "running", "binary": str(binary), "binary_sha256": digest(binary),
               "cases": [], "notes": ["Native CLI executions; no fault knob in production.",
               "MPI input copies use different absolute directories; bytes agree except deliberate mutations.",
               "Small correctness fixtures, not scalability or process-loss evidence."]}
    try:
        for name, kind, ranks, stage in cases:
            case = directory / name
            records, output = case / "ranks", case / "results"
            records.mkdir(parents=True)
            source = repo / "examples/one_d" / ("rigid_straight" if kind == "rigid" else
                "multispecies_physiology" if kind == "species" else "compliant_bifurcation")
            configuration = json.loads((source / "simulation_config.json").read_text())
            if name in ("step-inlet", "step-inlet-assets"):
                condition = next(item for item in configuration["boundaries"][0]["conditions"]
                                 if item["field"] == "oxygen")
                condition.pop("value")
                condition["waveform"] = "oxygen_table"
                configuration["temporal_functions"].append({"name": "oxygen_table", "kind": "periodic_table",
                    "units": "concentration", "period": 1.0, "file": "oxygen.csv", "interpolation": "linear"})
            if kind == "implicit":
                configuration["equation_systems"][0]["scheme"] = "implicit_petsc"
                configuration["time"] = {"dt": 0.001, "steps": 2, "output_every": 1}
            if name == "output-setup":
                output.write_text("deliberate non-directory")
            elif name.startswith("output-"):
                output.mkdir()
                blocked = {"output-initial": "profile_1d_000000.vtp",
                           "output-step": "profile_1d_000001.vtp",
                           "output-final": "profile_1d.pvd"}
                if name == "output-buffered":
                    (output / "flow_timeseries.csv").symlink_to("/dev/full")
                else:
                    (output / blocked[name]).mkdir()
            command = ["timeout", "--kill-after=5s", "90s", *shlex.split(args.launcher)]
            inputs = {}
            for rank in range(ranks):
                fixture = case / f"input-{rank}"
                fixture.mkdir()
                config = json.loads(json.dumps(configuration))
                if name in ("step-inlet", "step-inlet-assets"):
                    # Differing files now reject before initialization. A separate
                    # common input still verifies the actual step-time inlet gate.
                    value = -0.14 if name == "step-inlet" or rank == 1 else 0.14
                    (fixture / "oxygen.csv").write_text(f"time,value\n0,0.14\n0.002,{value}\n0.9,0.14\n")
                if rank == 1 and name == "config-different":
                    config["time"]["dt"] *= 2
                if not (rank == 1 and name == "config-missing"):
                    (fixture / "simulation_config.json").write_text(json.dumps(config) + "\n")
                if not (rank == 1 and name == "geometry-missing"):
                    shutil.copy2(source / config["geometry"]["file"], fixture)
                for path in fixture.iterdir():
                    inputs[str(path)] = digest(path)
                child = [str(binary), str(fixture), "--output-dir", str(output)]
                if rank == 1:
                    if name == "arguments": child += ["--stop-after-step", "invalid"]
                    if name == "check-control": child += ["--check"]
                    if name == "stop-control": child += ["--stop-after-step", "1"]
                    if name == "petsc-flag": child += ["-ksp_monitor"]
                if rank:
                    command += [":"]
                if ranks == 3:
                    worker = [sys.executable, str(wrapper), "--rank-worker", "--output-dir", str(records)]
                else:
                    worker = [sys.executable, str(repo / "scripts/hpc_rank_run.py"),
                              "--output-dir", str(records), "--expected-ranks", "1", "--timeout", "60"]
                command += ["-np", "1", *worker, "--", *child]
            expected = 1 if stage else 0
            evidence = {"case": name, "command_argv": command, "input_sha256": inputs,
                        "expected_returncode": expected, "status": "running", "ranks": []}
            summary["cases"].append(evidence)
            with (case / "launcher.log").open("x") as log:
                result = subprocess.run(command, cwd=repo, env=env, stdout=log, stderr=subprocess.STDOUT)
            evidence["launcher_returncode"] = result.returncode
            if result.returncode != expected:
                raise RuntimeError(f"{name}: unexpected launcher exit {result.returncode}")
            for rank in range(ranks):
                record_dir = records / f"rank-{rank}"
                report = json.loads((record_dir / "run.json").read_text())
                if (report["returncode"] != expected or report["timed_out"] or not report["resource"]
                        or report["rank"] != rank or report["ranks"] != ranks):
                    raise RuntimeError(f"{name}: rank {rank} did not exit as expected")
                for log, checksum in report["logs"].items():
                    if digest(record_dir / log) != checksum:
                        raise RuntimeError("log hash mismatch")
                evidence["ranks"].append(report)
            if stage:
                if stage not in (records / "rank-0/stderr.log").read_text():
                    raise RuntimeError(f"{name}: wrong failure boundary")
                if (output / "summary.json").exists():
                    raise RuntimeError(f"{name}: unexpected success summary")
                if name == "step-inlet" and (not (output / "profile_1d_000000.vtp").exists()
                        or (output / "profile_1d_000001.vtp").exists()):
                    raise RuntimeError("inlet fault did not occur after initialization and before step 1 output")
            else:
                completed = json.loads((output / "summary.json").read_text())
                if not completed["converged"] or completed["completed_step"] != configuration["time"]["steps"]:
                    raise RuntimeError("incomplete healthy run")
                evidence["fields"] = compare(directory / f"{kind}-1/results", output, include_species=kind == "species")
            evidence["status"] = "passed"
            print(f"{name}: expected exit {expected}, {ranks} rank reports verified", flush=True)
        summary["status"] = "passed"
    except Exception as error:
        summary.update(status="failed", error=str(error))
        raise
    finally:
        (directory / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
