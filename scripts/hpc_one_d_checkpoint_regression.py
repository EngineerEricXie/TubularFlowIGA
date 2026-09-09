#!/usr/bin/env python3
"""Native checkpoint faults, restarted trajectories and PETSc binary compatibility."""

import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import struct
import subprocess
import sys

from hpc_inventory import digest
from hpc_one_d_cli_regression import compare


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--launcher", default="mpiexec --oversubscribe")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    binary = repo / "solvers/one_d/iga_1d"
    format_test = repo / "solvers/one_d/one_d_checkpoint_format_test"
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               PETSC_OPTIONS="-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps")
    summary = {"status": "running", "binaries": {str(p): digest(p) for p in (binary, format_test)}, "cases": []}

    def run(name, children, expected=0, stage=""):
        case = root / name
        records = case / "ranks"
        records.mkdir(parents=True)
        command = ["timeout", "--kill-after=5s", "90s", *shlex.split(args.launcher)]
        for rank, child in enumerate(children):
            if rank: command += [":"]
            command += ["-np", "1", sys.executable, str(repo / "scripts/hpc_failure_regression.py"),
                        "--rank-worker", "--output-dir", str(records), "--", *map(str, child)]
        entry = {"case": name, "command_argv": command, "expected_returncode": expected,
                 "status": "running", "ranks": []}
        summary["cases"].append(entry)
        with (case / "launcher.log").open("x") as log:
            result = subprocess.run(command, cwd=repo, env=env, stdout=log, stderr=subprocess.STDOUT)
        entry["launcher_returncode"] = result.returncode
        if result.returncode != expected:
            raise RuntimeError(f"{name}: unexpected job exit {result.returncode}")
        for rank in range(3):
            record_dir = records / f"rank-{rank}"
            report = json.loads((record_dir / "run.json").read_text())
            if (report["returncode"] != expected or report["timed_out"] or not report["resource"]
                    or report["rank"] != rank or report["ranks"] != 3):
                raise RuntimeError(f"{name}: rank {rank} exit/report mismatch")
            for log, checksum in report["logs"].items():
                if digest(record_dir / log) != checksum: raise RuntimeError("log hash mismatch")
            entry["ranks"].append(report)
        if stage and stage not in (records / "rank-0/stderr.log").read_text():
            raise RuntimeError(f"{name}: unexpected failure stage")
        if expected and (case / "results/summary.json").exists():
            raise RuntimeError(f"{name}: false success summary")
        entry["status"] = "passed" if expected else "process_verified"
        print(f"{name}: job and all rank exits {expected}, reports verified", flush=True)
        return entry

    def cli(name, fixture, extra):
        return [binary, fixture, "--output-dir", root / name / "results", *extra]

    try:
        for kind, source in (("rigid", "rigid_straight"), ("explicit", "compliant_bifurcation"),
                             ("species", "multispecies_physiology")):
            fixture = repo / "examples/one_d" / source
            config = json.loads((fixture / "simulation_config.json").read_text())
            stop = config["time"]["steps"] // 2
            prefix = root / f"{kind}-save/checkpoint"
            run(f"{kind}-full", [cli(f"{kind}-full", fixture, [])] * 3)
            run(f"{kind}-save", [cli(f"{kind}-save", fixture,
                ["--stop-after-step", stop, "--checkpoint", prefix, "--checkpoint-every", stop])] * 3)
            entry = run(f"{kind}-resume", [cli(f"{kind}-resume", fixture, ["--restart", prefix])] * 3)
            metadata = json.loads(prefix.with_suffix(".json").read_text())
            result = json.loads((root / f"{kind}-resume/results/summary.json").read_text())
            if (metadata["completed_step"] != stop or result["completed_step"] != config["time"]["steps"]
                    or not result["converged"]):
                raise RuntimeError("restart did not reach requested accepted step")
            entry["restart_field_errors"] = compare(root / f"{kind}-full/results",
                root / f"{kind}-resume/results", metadata["physical_time"], include_species=True)
            entry["checkpoint_sha256"] = {str(prefix.with_suffix(suffix)): digest(prefix.with_suffix(suffix))
                                          for suffix in (".json", ".state")}
        fixture = repo / "examples/one_d/multispecies_physiology"
        source = root / "species-save/checkpoint"
        legacy = root / "legacy"
        legacy.mkdir()
        run("legacy-format", [[format_test, source.with_suffix(".state"), legacy / "checkpoint.state"]] * 3)
        if digest(source.with_suffix(".state")) != digest(legacy / "checkpoint.state"):
            raise RuntimeError("legacy distributed VecView changed binary bytes")
        shutil.copy2(source.with_suffix(".json"), legacy / "checkpoint.json")
        entry = run("legacy-resume", [cli("legacy-resume", fixture, ["--restart", legacy / "checkpoint"])] * 3)
        metadata = json.loads(source.with_suffix(".json").read_text())
        entry["restart_field_errors"] = compare(root / "species-full/results", root / "legacy-resume/results",
                                               metadata["physical_time"], include_species=True)
        faults = [
            ("metadata-missing", "1d checkpoint metadata read: rank 1:"),
            ("metadata-malformed", "1d checkpoint metadata read: rank 1:"),
            ("metadata-different", "1d checkpoint metadata agreement: rank 1:"),
            ("metadata-overflow", "1d checkpoint metadata read: rank 1:"),
            ("state-missing", "1d checkpoint state read: rank 1:"),
            ("state-truncated", "1d checkpoint state read: rank 1:"),
            ("state-header", "1d checkpoint state read: rank 1:"),
            ("state-nonfinite", "1d checkpoint state read: rank 1:"),
            ("state-area", "1d checkpoint state read: rank 1:"),
            ("state-different", "1d checkpoint state agreement: rank 1:"),
            ("write-parent", "1d checkpoint preparation: rank 1:"),
            ("write-state", "1d checkpoint state write: rank 0:"),
            ("write-metadata", "1d checkpoint metadata write: rank 0:"),
            ("write-buffered", "1d checkpoint metadata write: rank 0:")]
        for name, stage in faults:
            case = root / name
            children = []
            for rank in range(3):
                local = case / f"replica-{rank}"
                local.mkdir(parents=True)
                prefix = local / "checkpoint"
                if name.startswith("write-"):
                    extra = ["--checkpoint", prefix, "--checkpoint-every", "1", "--stop-after-step", "1"]
                    if name == "write-parent" and rank == 1:
                        (local / "blocked").write_text("not a directory")
                        extra[1] = local / "blocked/checkpoint"
                    if rank == 0:
                        if name == "write-state": prefix.with_suffix(".state").mkdir()
                        if name == "write-metadata": prefix.with_suffix(".json").mkdir()
                        if name == "write-buffered": prefix.with_suffix(".json").symlink_to("/dev/full")
                else:
                    for suffix in (".json", ".state"):
                        shutil.copy2(source.with_suffix(suffix), prefix.with_suffix(suffix))
                    extra = ["--restart", prefix]
                    if rank == 1:
                        path = prefix.with_suffix(".json" if name.startswith("metadata-") else ".state")
                        if name.endswith("missing"): path.unlink()
                        elif name == "metadata-malformed": path.write_text("{")
                        elif name.startswith("metadata-"):
                            changed = json.loads(path.read_text())
                            if name == "metadata-overflow": changed["cells"] = 2147483647
                            else: changed["inlet_flow"] *= 2
                            path.write_text(json.dumps(changed))
                        else:
                            data = bytearray(path.read_bytes())
                            # PETSc real/double, 32-bit-index fixture: class ID and
                            # length precede big-endian scalar values. Reject other builds.
                            class_id, count = struct.unpack(">ii", data[:8])
                            if class_id != 1211214 or len(data) != 8+8*count:
                                raise RuntimeError("fault fixture requires PETSc real/double with 32-bit indices")
                            if name == "state-truncated": data = data[:-5]
                            elif name == "state-header": data[:4] = b"BAD!"
                            else:
                                old = struct.unpack(">d", data[8:16])[0]
                                value = float("nan") if name == "state-nonfinite" else 0.0 if name == "state-area" else old*1.01
                                data[8:16] = struct.pack(">d", value)
                            path.write_bytes(data)
                children.append(cli(name, fixture, extra))
            run(name, children, 1, stage)
        for entry in summary["cases"]:
            entry["status"] = "passed"
        summary["status"] = "passed"
    except Exception as error:
        summary.update(status="failed", error=str(error))
        raise
    finally:
        (root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
