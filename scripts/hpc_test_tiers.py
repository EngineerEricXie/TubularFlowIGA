#!/usr/bin/env python3
"""Run one documented test tier with deadlines and a machine-readable result."""

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import time


ROOT = Path(__file__).resolve().parents[1]


def save(path, report):
    path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


def run(command, directory, environment, timeout, name):
    stdout = directory / f"{name}.stdout"
    stderr = directory / f"{name}.stderr"
    started = time.perf_counter()
    timed_out = False
    with stdout.open("x") as out, stderr.open("x") as err:
        process = subprocess.Popen(command, cwd=ROOT, env=environment, stdout=out,
                                   stderr=err, start_new_session=True)
        try:
            returncode = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            returncode = 124
    return {"name": name, "argv": command, "returncode": returncode,
            "timed_out": timed_out, "wall_s": time.perf_counter() - started,
            "stdout": str(stdout), "stderr": str(stderr),
            "status": "passed" if returncode == 0 and not timed_out else "failed"}


def unit_commands(_args):
    return [
        ["python3", "-m", "unittest", "discover", "-s", "scripts/tests"],
        ["make", "mesh-test"],
        ["make", "-C", "solvers/coupling", "test"],
        ["make", "-C", "solvers/cuda", "execution-test"],
    ], None


def mpi_commands(args):
    if not shutil.which(args.mpiexec):
        return [], f"MPI launcher is unavailable: {args.mpiexec}"
    petsc = os.environ.get("PETSC_DIR")
    if not petsc:
        return [], "PETSC_DIR is unset"
    build = ["make", "-C", "solvers/cpu", "distributed_immersed_extension_test",
             f"PETSC_DIR={petsc}"]
    if os.environ.get("PETSC_ARCH"):
        build.append(f"PETSC_ARCH={os.environ['PETSC_ARCH']}")
    commands = [build]
    launcher = shlex.split(args.launcher or args.mpiexec)
    binary = str(ROOT / "solvers/cpu/distributed_immersed_extension_test")
    for ranks in (1, 2, 4):
        commands.append(launcher + ["--bind-to", "core", "--oversubscribe",
                                    "-np", str(ranks), binary])
    return commands, None


def gpu_commands(_args):
    binary = ROOT / "solvers/cuda/iga_cuda"
    contract = ROOT / "solvers/cuda/cuda_execution_test"
    missing = [str(path) for path in (binary, contract) if not path.is_file()]
    if missing:
        return [], "GPU tier was not built: " + ", ".join(missing)
    if not Path("/dev/nvidiactl").exists() and not shutil.which("nvidia-smi"):
        return [], "no NVIDIA device interface is visible"
    return [[str(contract)], [str(binary), "device-info"]], None


def scheduled_commands(args):
    try:
        nodes = int(os.environ.get("SLURM_JOB_NUM_NODES", "0"))
        tasks = int(os.environ.get("SLURM_NTASKS", "0"))
    except ValueError:
        return [], "invalid Slurm allocation metadata"
    if nodes < 2 or tasks < 2:
        return [], "scheduled tier requires a Slurm allocation with at least two nodes and two tasks"
    binary = ROOT / "solvers/cpu/distributed_immersed_extension_test"
    if not binary.is_file():
        return [], "scheduled MPI test binary is missing"
    launcher = shlex.split(args.launcher or "mpiexec")
    return [launcher + ["--map-by", "ppr:1:node", "--bind-to", "core", "-np", "2",
                        str(binary)]], None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tier", choices=("unit", "mpi", "gpu", "scheduled"), required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--timeout", type=float, default=300)
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("--launcher", help="launcher prefix, parsed with shell-like quoting")
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("timeout must be positive")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    report_path = output / "result.json"
    report = {"schema_version": 1, "kind": "hpc_test_tier", "tier": args.tier,
              "created_utc": datetime.now(timezone.utc).isoformat(), "status": "running",
              "skip_reason": None, "commands": [],
              "source_commit": subprocess.check_output(
                  ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
              "allocation": {name: os.environ.get(name) for name in (
                  "SLURM_JOB_ID", "SLURM_JOB_NUM_NODES", "SLURM_NTASKS",
                  "SLURM_CPUS_PER_TASK")}}
    save(report_path, report)
    builders = {"unit": unit_commands, "mpi": mpi_commands,
                "gpu": gpu_commands, "scheduled": scheduled_commands}
    commands, reason = builders[args.tier](args)
    if reason:
        report["status"] = "skipped"
        report["skip_reason"] = reason
        save(report_path, report)
        print(f"{args.tier} tier skipped: {reason}")
        return
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS="1", OMP_THREAD_LIMIT="1", OMP_DYNAMIC="FALSE",
                       OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", BLIS_NUM_THREADS="1")
    for index, command in enumerate(commands):
        result = run(command, output, environment, args.timeout, f"command-{index:02d}")
        report["commands"].append(result)
        save(report_path, report)
        if result["status"] != "passed":
            report["status"] = "failed"
            save(report_path, report)
            raise SystemExit(result["returncode"] or 1)
    report["status"] = "passed"
    save(report_path, report)
    print(f"{args.tier} tier passed ({len(commands)} commands)")


if __name__ == "__main__":
    main()
