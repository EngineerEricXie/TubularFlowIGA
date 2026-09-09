#!/usr/bin/env python3
"""Verify native 1D asset identity, local replicas, and unchanged numerical fields."""

import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

from hpc_inventory import digest
from hpc_one_d_cli_regression import compare


ROOT = Path(__file__).resolve().parents[1]


def fixture(kind):
    source = ROOT / "examples/one_d" / ("multispecies_physiology" if kind == "species" else "rigid_straight")
    config = json.loads((source / "simulation_config.json").read_text())
    config["time"]["steps"] = 2
    network = config["geometry"]["file"]
    files = {network: (source / network).read_text()}
    if kind == "obj":
        rows = [line.split() for line in files.pop(network).splitlines() if line and not line.startswith("#")]
        files["network.obj"] = "".join("v " + " ".join(row[2:6]) + " 0 0\n" for row in rows) + "l 1 2 3\n"
        config["geometry"].update(kind="obj_network", file="network.obj", root_node_id=1)
    if kind == "species":
        condition = next(c for b in config["boundaries"] for c in b["conditions"]
                         if b["role"] == "inlet" and c["field"] == "oxygen")
        condition.pop("value")
        condition["waveform"] = "oxygen_table"
        config["temporal_functions"].append(dict(name="oxygen_table", kind="periodic_table",
            units="concentration", period=1.0, file="oxygen.csv", interpolation="linear"))
        files["oxygen.csv"] = "time,value\n0,0.14\n0.5,0.14\n"
    else:
        config["temporal_functions"][0] = dict(name="inlet_flow", kind="periodic_table",
            units="m3/s", period=1.0, file="flow.csv", interpolation="linear")
        files["flow.csv"] = "time,value\n0,1e-8\n0.5,1e-8\n"
    if kind.startswith("replay-"):
        extension = kind.split("-", 1)[1]
        config["simulation_scope"] = dict(mode="vca_replay")
        config["coupling"] = dict(scheme="explicit_staggered", replay_file="replay." + extension)
        files["replay." + extension] = (
            "time_s,flow_m3_s\n0,1e-8\n0.03,1.2e-8\n" if extension == "csv" else
            json.dumps(dict(states=[dict(time_s=0, flow_m3_s=1e-8), dict(time_s=0.03, flow_m3_s=1.2e-8)])) + "\n")
    return config, files


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--reference-binary", type=Path, required=True)
    parser.add_argument("--launcher", default="mpiexec --map-by core --bind-to core")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    reference = args.reference_binary.resolve()
    binary = ROOT / "solvers/one_d/iga_1d"
    binaries = {str(p): digest(p) for p in (reference, binary)}
    environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1", BLIS_NUM_THREADS="1",
        PETSC_OPTIONS="-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps")
    cases = []
    for kind in ("flow", "species", "replay-csv", "replay-json", "obj"):
        cases.extend([(kind + "-before", kind, 1, "", False, True),
                      (kind + "-one", kind, 1, "", False, False),
                      (kind + "-three", kind, 3, "", False, False)])
    for name, kind, stage in (
        ("network-swc-different", "flow", "1d asset network: rank 1:"),
        ("network-obj-different", "obj", "1d asset network: rank 1:"),
        ("flow-different", "flow", "1d asset temporal inlet_flow: rank 1:"),
        ("species-different", "species", "1d asset temporal oxygen_table: rank 1:"),
        ("replay-csv-different", "replay-csv", "1d asset replay: rank 1:"),
        ("replay-json-different", "replay-json", "1d asset replay: rank 1:"),
        ("network-missing", "flow", "asset content read: rank 1:"),
        ("network-fifo", "flow", "asset content read: rank 1:"),
        ("network-directory", "flow", "asset content read: rank 1:"),
        ("configuration-fifo", "flow", "1d configuration input: rank 1:"),
        ("configuration-directory", "flow", "1d configuration input: rank 1:"),
        ("flow-missing", "flow", "asset content read: rank 1:"),
        ("flow-fifo", "flow", "asset content read: rank 1:"),
        ("replay-missing", "replay-csv", "asset content read: rank 1:"),
        ("replay-fifo", "replay-csv", "asset content read: rank 1:"),
        ("replay-malformed", "replay-csv", "1d replay input: rank 0:")):
        cases.append((name, kind, 3, stage, False, False))
    cases.extend([("unused-table", "flow", 3, "", False, False),
                  ("check-good", "flow", 3, "", True, False),
                  ("check-different", "flow", 3, "1d asset network: rank 1:", True, False),
                  ("check-fifo", "flow", 3, "asset content read: rank 1:", True, False)])
    summary = dict(status="running", binaries=binaries, cases=[])
    try:
        for name, kind, ranks, stage, check, before in cases:
            case = output / name
            records = case / "ranks"
            records.mkdir(parents=True)
            config, files = fixture(kind)
            if name == "unused-table":
                config["temporal_functions"].append(dict(name="unused", kind="periodic_table",
                    units="m3/s", period=1, file="deliberately-absent.csv", interpolation="linear"))
            command = ["timeout", "--kill-after=5s", "90s", *shlex.split(args.launcher)]
            input_hashes = {}
            for rank in range(ranks):
                local = case / f"input-{rank}"
                local.mkdir()
                payload = dict(files)
                network = config["geometry"]["file"]
                if rank == 1:
                    if name in ("network-swc-different", "check-different"):
                        payload[network] = payload[network].replace("0.001 -1", "0.002 -1", 1)
                    if name == "network-obj-different":
                        payload[network] = payload[network].replace("0.001 0 0", "0.002 0 0", 1)
                    if name == "flow-different": payload["flow.csv"] = payload["flow.csv"].replace("0.5,1e-8", "0.5,2e-8")
                    if name == "species-different": payload["oxygen.csv"] = payload["oxygen.csv"].replace("0.5,0.14", "0.5,0.15")
                    if name == "replay-csv-different": payload["replay.csv"] = payload["replay.csv"].replace("1.2e-8", "1.3e-8")
                    if name == "replay-json-different":
                        replay = json.loads(payload["replay.json"])
                        replay["states"][-1]["flow_m3_s"] = 1.3e-8
                        payload["replay.json"] = json.dumps(replay) + "\n"
                if name == "replay-malformed": payload["replay.csv"] = "time_s,flow_m3_s\n0,invalid\n"
                for filename, text in payload.items(): (local / filename).write_text(text)
                (local / "simulation_config.json").write_text(json.dumps(config) + "\n")
                if rank == 1 and name in ("network-missing", "network-fifo", "network-directory",
                        "flow-missing", "flow-fifo", "replay-missing", "replay-fifo", "check-fifo",
                        "configuration-fifo", "configuration-directory"):
                    target = local / ("simulation_config.json" if name.startswith("configuration-") else
                        network if name.startswith("network-") else
                        "replay.csv" if name.startswith("replay-") else "flow.csv")
                    target.unlink()
                    if name.endswith("fifo"): os.mkfifo(target, 0o600)
                    if name.endswith("directory"): target.mkdir()
                for path in local.iterdir():
                    if path.is_file(): input_hashes[str(path)] = digest(path)
                child = [str(reference if before else binary), str(local), "--output-dir", str(case / "result")]
                if check: child.append("--check")
                wrapper = ([sys.executable, str(ROOT / "scripts/hpc_failure_regression.py"), "--rank-worker",
                            "--output-dir", str(records)] if ranks == 3 else
                           [sys.executable, str(ROOT / "scripts/hpc_rank_run.py"), "--output-dir", str(records),
                            "--expected-ranks", "1", "--timeout", "60"])
                if rank: command.append(":")
                command += ["-np", "1", *wrapper, "--", *child]
            entry = dict(name=name, command_argv=command, input_sha256=input_hashes,
                         expected_returncode=1 if stage else 0, expected_stage=stage, status="running", ranks=[])
            summary["cases"].append(entry)
            with (case / "launcher.log").open("x") as log:
                result = subprocess.run(command, cwd=ROOT, env=environment, stdout=log, stderr=subprocess.STDOUT)
            entry["launcher_returncode"] = result.returncode
            if result.returncode != entry["expected_returncode"]: raise RuntimeError(f"{name}: unexpected job exit")
            for rank in range(ranks):
                d = records / f"rank-{rank}"
                report = json.loads((d / "run.json").read_text())
                if report["returncode"] != entry["expected_returncode"] or report["timed_out"]:
                    raise RuntimeError(f"{name}: unexpected rank exit or timeout")
                if report["rank"] != rank or report["ranks"] != ranks: raise RuntimeError("rank identity differs")
                for filename, checksum in report["logs"].items():
                    if digest(d / filename) != checksum: raise RuntimeError("rank log hash differs")
                entry["ranks"].append(report)
            if stage:
                if stage not in (records / "rank-0/stderr.log").read_text(): raise RuntimeError(f"{name}: wrong failure stage")
                if (case / "result").exists(): raise RuntimeError(f"{name}: rejected assets published output")
            elif check:
                if "schema_version=3 dimension=1d" not in (records / "rank-0/stdout.log").read_text():
                    raise RuntimeError("check did not publish its success summary")
                if (case / "result").exists(): raise RuntimeError("check unexpectedly ran a simulation")
            else:
                completed = json.loads((case / "result/summary.json").read_text())
                if not completed["converged"] or completed["completed_step"] != 2: raise RuntimeError("incomplete solution")
                entry["field_comparisons"] = compare(output / (kind + "-before/result"), case / "result", include_species=kind == "species")
            if any(digest(Path(p)) != h for p, h in input_hashes.items()): raise RuntimeError("inputs changed during execution")
            entry["status"] = "passed"
            print(name, "passed", flush=True)
        if any(digest(Path(p)) != h for p, h in binaries.items()): raise RuntimeError("binaries changed during regression")
        summary["status"] = "passed"
    except BaseException:
        summary["status"] = "failed"
        raise
    finally:
        (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
