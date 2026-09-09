#!/usr/bin/env python3
"""Exercise native 1D initial/step filename failures and healthy retries."""

import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--launcher", default="mpiexec --oversubscribe")
    parser.add_argument("--expect-swallowed", action="store_true",
                        help="negative control: require the pre-fix CLI to publish false success")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    binary = (args.binary or repo / "solvers/one_d/one_d_filename_failure_test").resolve()
    fixture = root / "fixture"
    shutil.copytree(repo / "examples/one_d/rigid_straight", fixture)
    config_path = fixture / "simulation_config.json"
    config = json.loads(config_path.read_text())
    config["time"].update(steps=2, output_every=1)
    config_path.write_text(json.dumps(config, indent=2) + "\n")
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               PETSC_OPTIONS="-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps")
    summary = {"status": "running", "negative_control": args.expect_swallowed, "cases": []}

    def run(name, ranks, mode, step):
        case = root / name
        reports, output = case / "ranks", case / "results"
        reports.mkdir(parents=True)
        expected = 0 if mode < 0 or args.expect_swallowed else 1
        wrapper = ([str(repo / "scripts/hpc_failure_regression.py"), "--rank-worker"]
                   if ranks == 3 else [str(repo / "scripts/hpc_rank_run.py"), "--timeout", "30"])
        command = ["timeout", "--kill-after=5s", "90s", *shlex.split(args.launcher),
                   "-np", str(ranks), sys.executable, *wrapper,
                   "--output-dir", str(reports), "--", str(binary),
                   str(mode), str(step), str(fixture), "--output-dir", str(output)]
        with (case / "launcher.log").open("x") as log:
            result = subprocess.run(command, cwd=repo, env=env, stdout=log, stderr=subprocess.STDOUT)
        entry = {"name": name, "command": command, "expected_exit": expected,
                 "launcher_exit": result.returncode, "ranks": []}
        summary["cases"].append(entry)
        if result.returncode != expected:
            raise RuntimeError(f"{name}: unexpected launcher exit {result.returncode}")
        for rank in range(ranks):
            report = json.loads((reports / f"rank-{rank}/run.json").read_text())
            entry["ranks"].append(report)
            if report["returncode"] != expected or report["timed_out"]:
                raise RuntimeError(f"{name}: wrong rank exit")
        stderr = (reports / "rank-0/stderr.log").read_text()
        if mode >= 0:
            if f"filename formatting injected mode={mode} step={step}" not in stderr:
                raise RuntimeError(f"{name}: missed filename formatter")
            if expected:
                stage = "1d initial output" if step == 0 else "1d step output"
                if stage + ": rank 0:" not in stderr:
                    raise RuntimeError(f"{name}: wrong failure stage")
                if (output / "summary.json").exists() or (output / "profile_1d_").exists():
                    raise RuntimeError(f"{name}: published false success or incomplete filename")
            elif not (output / "profile_1d_").is_file():
                raise RuntimeError(f"{name}: negative control did not demonstrate truncated filename")
        if expected == 0:
            final = json.loads((output / "summary.json").read_text())
            if not final["converged"] or final["completed_step"] != 2:
                raise RuntimeError(f"{name}: incomplete healthy run")
        entry["status"] = "passed"
        print(f"{name}: {ranks} rank exits {expected}, formatter and publication verified", flush=True)
        return output

    try:
        for ranks in (1, 3):
            baseline = run(f"healthy-{ranks}", ranks, -1, 0)
            for step in (0, 1):
                for mode in range(3):
                    run(f"fault-{ranks}-{step}-{mode}", ranks, mode, step)
                    if not args.expect_swallowed:
                        retry = run(f"retry-{ranks}-{step}-{mode}", ranks, -1, step)
                        files = {p.name for p in baseline.iterdir() if p.name != "summary.json"}
                        if files != {p.name for p in retry.iterdir() if p.name != "summary.json"}:
                            raise RuntimeError("retry file set differs")
                        for name in files:
                            if (baseline / name).read_bytes() != (retry / name).read_bytes():
                                raise RuntimeError(f"retry changed {name}")
        summary["status"] = "passed"
    except Exception as error:
        summary.update(status="failed", error=str(error))
        raise
    finally:
        (root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
