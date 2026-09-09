#!/usr/bin/env python3
"""Run the three-rank controlled-failure gates, retaining nonzero job evidence."""

import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time

from hpc_inventory import digest
from hpc_rank_run import rank_identity, run as run_rank


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--launcher", default="mpiexec --oversubscribe",
        help="Open MPI default; override for another MPI installation")
    parser.add_argument("--rank-worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("command", nargs=argparse.REMAINDER, help=argparse.SUPPRESS)
    options = parser.parse_args()
    if options.rank_worker:
        command = options.command[1:] if options.command[:1] == ["--"] else options.command
        if not command:
            parser.error("rank worker requires a command")
        status = run_rank(command, options.output_dir, 3, 60)
        rank, _ = rank_identity(os.environ)
        # Open MPI can terminate peer wrappers as soon as one exits nonzero.
        # Publish readiness only after the measured child's logs/report close,
        # then retain the ORIGINAL nonzero status once every report is closed.
        # This test-only rendezvous requires the shared, fresh output directory.
        (options.output_dir / f"rank-{rank}" / "report-ready").touch(exist_ok=False)
        deadline = time.monotonic() + 15
        while not all((options.output_dir / f"rank-{peer}" / "report-ready").exists()
                      for peer in range(3)):
            if time.monotonic() >= deadline:
                print("timed out waiting for peer measurement reports", file=sys.stderr)
                return 2
            time.sleep(0.01)
        return status
    if options.command:
        parser.error("unexpected positional arguments")
    repository = Path(__file__).resolve().parent.parent
    directory = options.output_dir.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    wrapper = Path(__file__).resolve()
    unit = repository / "solvers/cpu/collective_failure_test"
    options_unit = repository / "solvers/cpu/collective_petsc_options_test"
    graph = repository / "solvers/coupling/multidomain_failure_test"
    environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                       PETSC_OPTIONS="-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps")
    for name in ("TUBULARFLOWIGA_INJECT_EXPLICIT_COUPLING_FAILURE_STEP",
                 "TUBULARFLOWIGA_INJECT_BIFURCATION_FAILURE_STEP"):
        environment.pop(name, None)
    cases = [("unit", unit, 0, ""),
             ("options-unit", options_unit, 0, ""),
             ("arguments", graph, 1, "graph arguments: rank 1:"),
             ("input", graph, 1, "graph input: rank 1:"),
             ("database", graph, 1, "graph asset junction / database: rank 1:"),
             ("output", graph, 1, "graph output: rank 0:"),
             ("stop", graph, 1, "graph execution controls: rank 1:"),
             ("newton", graph, 1, "graph execution controls: rank 1:"),
             ("injection", graph, 1, "graph execution controls: rank 1:"),
             ("manifest", graph, 1, "graph manifest: rank 1:"),
             ("one-d-config", graph, 1, "domain configuration source: rank 1:"),
             ("three-d-config", graph, 1, "domain configuration junction: rank 1:"),
             ("zero-d-config", graph, 1, "0D model source_0d: rank 1:"),
             ("healthy", graph, 0, ""),
             ("replica", graph, 0, ""),
             ("zero-d-healthy", graph, 0, ""),
             ("petsc-ksp", graph, 1, "PETSc option agreement: rank 1:"),
             ("petsc-prefix", graph, 1, "PETSc option agreement: rank 1:"),
             ("petsc-flag", graph, 1, "PETSc option agreement: rank 1:"),
             ("petsc-used", graph, 0, ""),
             ("petsc-file", graph, 0, ""),
             ("petsc-file-different", graph, 1, "PETSc option agreement: rank 1:"),
             ("asset-database", graph, 1, "graph asset junction / database: rank 1:"),
             ("asset-network", graph, 1, "graph asset source / network: rank 1:"),
             ("asset-mesh", graph, 1, "graph asset junction / control mesh: rank 1:"),
             ("asset-velocity", graph, 1, "graph asset junction / initial velocity: rank 1:"),
             ("asset-waveform-1d", graph, 1, "graph asset source / temporal inlet_flow: rank 1:"),
             ("asset-waveform-3d", graph, 1, "graph asset junction / temporal inlet_scale: rank 1:"),
             ("asset-missing", graph, 1, "graph asset catalog: rank 1:"),
             ("asset-fifo", graph, 1, "graph asset catalog: rank 1:"),
             ("asset-replica", graph, 0, "")]
    summary = {"schema_version": 1, "kind": "hpc_controlled_failure_regression",
               "status": "running", "cases": [],
               "binaries": {str(path): digest(path) for path in (unit, options_unit, graph)},
               "notes": ["Expected nonzero exits remain failed process records; this report verifies their failure contract.",
                         "These tiny fixtures establish correctness, not scaling or process-loss recovery."]}
    try:
        for name, binary, expected, diagnostic in cases:
            case = directory / name
            records = case / "ranks"
            records.mkdir(parents=True)
            arguments = [str(binary), str(case / "fixture")]
            if name == "options-unit":
                arguments = [str(binary)]
            elif name != "unit":
                arguments.append(name)
            command = ["timeout", "--kill-after=5s", "90s", *shlex.split(options.launcher),
                       "-np", "3", sys.executable, str(wrapper), "--rank-worker",
                       "--output-dir", str(records), "--", *arguments]
            with (case / "launcher.stdout").open("x") as stdout, (case / "launcher.stderr").open("x") as stderr:
                result = subprocess.run(command, cwd=repository, env=environment,
                                        stdout=stdout, stderr=stderr, check=False)
            evidence = {"case": name, "command_argv": command, "expected_returncode": expected,
                        "launcher_returncode": result.returncode, "status": "checking", "ranks": []}
            summary["cases"].append(evidence)
            if result.returncode != expected:
                raise RuntimeError(f"{name}: launcher exit {result.returncode}, expected {expected}")
            for rank in range(3):
                rank_directory = records / f"rank-{rank}"
                report = json.loads((rank_directory / "run.json").read_text())
                stdout = (rank_directory / "stdout.log").read_text()
                token = (f"collective failure rank={rank} passed" if name == "unit" else
                         f"PETSc options rank={rank} passed" if name == "options-unit" else
                         f"failure test mode={name} rank={rank} runner_status={expected} verified")
                if (report["rank"] != rank or report["ranks"] != 3
                        or report["returncode"] != expected or report["timed_out"]
                        or not report["resource"] or token not in stdout):
                    raise RuntimeError(f"{name}: rank {rank} did not verify its expected result")
                for log in ("stdout.log", "stderr.log"):
                    if digest(rank_directory / log) != report["logs"][log]:
                        raise RuntimeError(f"{name}: rank {rank} log hash mismatch")
                evidence["ranks"].append({"rank": rank, "returncode": report["returncode"],
                                          "wall_s": report["wall_s"], "resource": report["resource"]})
            if diagnostic and diagnostic not in (records / "rank-0/stderr.log").read_text():
                raise RuntimeError(f"{name}: failure occurred outside the intended production stage")
            evidence["status"] = "passed"
            print(f"{name}: expected job exit {expected}, all three ranks verified", flush=True)
        summary["status"] = "passed"
    except Exception as error:
        summary["status"] = "failed"
        summary["error"] = str(error)
        raise
    finally:
        (directory / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
