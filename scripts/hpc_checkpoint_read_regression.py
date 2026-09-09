#!/usr/bin/env python3
"""Reject damaged native 3D VCA restarts using retained vca_3d_smoke_test fixtures."""

import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import time

from hpc_inventory import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--launcher", default="mpiexec --oversubscribe")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    fixture = args.case_dir.resolve()
    output = args.output_dir.resolve()
    binary = repo / "solvers/cpu/iga_navier_stokes"
    sources = [(1, "fixture.ntiga", "split"), (2, "fixture-2.ntiga", "two-rank-split")]
    suffixes = (".json", ".state", ".vca.json", ".vca_transport.state")
    required = [fixture / "simulation_config.json", fixture / "controlmesh.vtk"]
    for _, database, directory in sources:
        required.append(fixture / database)
        required.extend(fixture / directory / ("checkpoint" + suffix) for suffix in suffixes)
    for path in required:
        if not path.is_file():
            parser.error(f"missing retained smoke fixture: {path}")
    output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               PETSC_OPTIONS="-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -ksp_rtol 1e-12")
    summary = {"status": "running", "binary_sha256": digest(binary),
               "inputs": {str(p): digest(p) for p in required},
               "environment": {k: env[k] for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "PETSC_OPTIONS")},
               "cases": []}
    modes = (
        ("missing-flow-metadata", ".json", "missing", "flow checkpoint metadata"),
        ("malformed-flow-metadata", ".json", "malformed", "flow checkpoint metadata"),
        ("truncated-flow", ".state", "truncate", "checkpoint local read"),
        ("missing-vca-metadata", ".vca.json", "missing", "VCA checkpoint metadata"),
        ("malformed-vca-metadata", ".vca.json", "malformed", "VCA checkpoint metadata"),
        ("truncated-transport", ".vca_transport.state", "truncate", "checkpoint local read"),
        ("trailing-transport", ".vca_transport.state", "trailing", "checkpoint local read"),
    )
    try:
        for ranks, database, directory in sources:
            for name, suffix, damage, stage in modes:
                case = output / f"ranks-{ranks}-{name}"
                checkpoint = case / "input"
                results = case / "results"
                checkpoint.mkdir(parents=True)
                results.mkdir()
                for item in suffixes:
                    shutil.copyfile(fixture / directory / ("checkpoint" + item),
                                    checkpoint / ("checkpoint" + item))
                target = checkpoint / ("checkpoint" + suffix)
                if damage == "missing":
                    target.unlink()
                elif damage == "malformed":
                    target.write_text("{invalid checkpoint metadata\n")
                elif damage == "truncate":
                    target.write_bytes(target.read_bytes()[:-1])
                else:
                    with target.open("ab") as stream:
                        stream.write(b"unexpected trailing data")
                command = ["timeout", "--kill-after=5s", "90s", *shlex.split(args.launcher),
                           "-np", str(ranks), str(binary), str(fixture / database), str(fixture),
                           "--max-newton", "12", "--restart", str(checkpoint / "checkpoint"),
                           "--checkpoint", str(results / "checkpoint"), "--output", str(results / "flow.txt")]
                entry = {"case": case.name, "ranks": ranks, "command_argv": command,
                         "required_stage": stage, "timeout_s": 90, "status": "running"}
                summary["cases"].append(entry)
                started = time.monotonic()
                log = case / "launcher.log"
                with log.open("x") as stream:
                    run = subprocess.run(command, cwd=repo, env=env, stdout=stream, stderr=subprocess.STDOUT)
                text = log.read_text()
                entry.update(returncode=run.returncode, elapsed_s=time.monotonic()-started,
                             log_sha256=digest(log))
                if run.returncode != 1 or stage not in text:
                    raise RuntimeError(f"{case.name}: expected exit 1 and diagnostic {stage!r}")
                if "restart=" in text or (results / "flow.txt").exists() or list(results.glob("checkpoint*")):
                    raise RuntimeError(f"{case.name}: failed restart published success/output")
                entry["status"] = "passed"
                print(case.name, "passed", flush=True)
        if any(digest(Path(p)) != checksum for p, checksum in summary["inputs"].items()):
            raise RuntimeError("source checkpoint fixture changed during regression")
        summary["status"] = "passed"
    except BaseException:
        summary["status"] = "failed"
        raise
    finally:
        (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
