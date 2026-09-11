#!/usr/bin/env python3
"""Measure body-fitted flow strong and weak scaling on prepared MPI cases."""

import argparse
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import shlex
import signal
import statistics
import subprocess
import sys
import time

from hpc_compare_fields import compare
from hpc_inventory import digest
from hpc_profile_summary import summarize


ROOT = Path(__file__).resolve().parents[1]


def write(path, value):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def case_path(pattern, ranks):
    return Path(pattern.format(ranks=ranks)).resolve(strict=True)


def samples(values):
    if len(values) < 3 or any(not math.isfinite(value) or value < 0 for value in values):
        raise ValueError("scaling statistics require at least three finite observations")
    return {"values": values, "median": statistics.median(values), "min": min(values),
            "max": max(values), "population_stdev": statistics.pstdev(values)}


def iterations(directory):
    total = 0
    configured = 0
    for line in (directory / "rank-0/stdout.log").read_text().splitlines():
        match = re.search(r"(?:^|\s)linear_iterations=(\d+)(?:\s|$)", line)
        if match:
            total += int(match.group(1))
        if line.startswith("solver_configuration "):
            values = dict(token.split("=", 1) for token in line.split()[1:] if "=" in token)
            if int(values.get("reason", "0")) <= 0:
                raise ValueError("solver configuration reports a nonconverged reason")
            configured += 1
    if not configured or total < 1:
        raise ValueError("missing positive solver iteration evidence")
    return total


def launch(command, directory, environment, timeout):
    started = time.perf_counter()
    with (directory / "launcher.stdout").open("x") as stdout, \
            (directory / "launcher.stderr").open("x") as stderr:
        process = subprocess.Popen(command, cwd=ROOT, env=environment, stdout=stdout,
                                   stderr=stderr, start_new_session=True)
        timed_out = False
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
    return {"argv": command, "returncode": returncode, "timed_out": timed_out,
            "wall_s": time.perf_counter() - started}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strong-case", required=True,
                        help="fixed case directory; may contain {ranks}")
    parser.add_argument("--strong-database", required=True,
                        help="database path relative to its case; may contain {ranks}")
    parser.add_argument("--weak-case", required=True,
                        help="scaled case directory pattern containing {ranks}")
    parser.add_argument("--weak-database", required=True,
                        help="database path relative to each weak case; may contain {ranks}")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--ranks", type=int, nargs="+", required=True)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=3600)
    parser.add_argument("--launcher", default="mpiexec --bind-to core --map-by core")
    parser.add_argument("--binary", type=Path,
                        default=ROOT / "solvers/cpu/iga_navier_stokes")
    parser.add_argument("--checker", type=Path, default=ROOT / "solvers/cpu/iga_mesh_check")
    parser.add_argument("--validator", type=Path,
                        default=ROOT / "solvers/cpu/iga_flow_validate")
    parser.add_argument("--rtol", type=float, default=1e-6)
    parser.add_argument("--zero-atol", type=float, default=1e-12)
    parser.add_argument("--weak-workload-relative-tolerance", type=float, default=0.15)
    args = parser.parse_args()
    ranks = sorted(args.ranks)
    if (not ranks or ranks[0] != 1 or len(ranks) != len(set(ranks)) or
            args.repetitions < 3 or args.timeout <= 0):
        parser.error("require unique ranks including 1, at least three repetitions, and positive timeout")
    if "{ranks}" not in args.weak_case:
        parser.error("--weak-case must contain {ranks}")
    for path in (args.binary, args.checker, args.validator):
        if not path.resolve().is_file():
            parser.error(f"missing executable: {path}")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    report_path = output / "summary.json"
    launcher = shlex.split(args.launcher)
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS="1", OMP_THREAD_LIMIT="1", OMP_DYNAMIC="FALSE",
                       OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", BLIS_NUM_THREADS="1",
                       IGA_PROFILE="1", PETSC_OPTIONS=(
                           "-ksp_type fgmres -ksp_rtol 1e-8 -pc_type bjacobi "
                           "-sub_ksp_type preonly -sub_pc_type ilu"))
    report = {"schema_version": 1, "kind": "hpc_cross_node_scaling",
              "created_utc": datetime.now(timezone.utc).isoformat(), "status": "running",
              "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"],
                                                        cwd=ROOT, text=True).strip(),
              "source_status": subprocess.check_output(
                  ["git", "status", "--porcelain=v1", "--untracked-files=all"],
                  cwd=ROOT, text=True).splitlines(),
              "ranks": ranks, "repetitions": args.repetitions, "launcher": launcher,
              "allocation": {name: os.environ.get(name) for name in (
                  "SLURM_JOB_ID", "SLURM_JOB_NUM_NODES", "SLURM_NODELIST",
                  "SLURM_NTASKS", "SLURM_CPUS_PER_TASK")},
              "inputs": {}, "runs": [], "summary": {}}
    provenance_paths = [args.binary.resolve(), args.checker.resolve(), args.validator.resolve(),
                        ROOT / "scripts/hpc_cross_node_scaling.py",
                        ROOT / "scripts/hpc_rank_run.py", ROOT / "scripts/hpc_profile_summary.py",
                        ROOT / "scripts/hpc_compare_fields.py", ROOT / "scripts/hpc_inventory.py"]
    report["provenance_sha256"] = {str(path): digest(path) for path in provenance_paths}
    write(report_path, report)
    failure = None
    try:
        cases = {}
        for mode in ("strong", "weak"):
            for count in ranks:
                case = case_path(args.strong_case if mode == "strong" else args.weak_case, count)
                if output.is_relative_to(case) or case.is_relative_to(output):
                    raise ValueError("scaling cases and output directory must be separate")
                database_pattern = args.strong_database if mode == "strong" else args.weak_database
                database = case / database_pattern.format(ranks=count)
                if not database.is_file():
                    raise ValueError(f"missing {mode} database for {count} ranks: {database}")
                cases[(mode, count)] = (case, database)
                report["inputs"][f"{mode}-{count}"] = {"case": str(case),
                    "database": str(database), "database_sha256": digest(database),
                    "case_sha256": {str(path): digest(path) for path in sorted(case.rglob("*"))
                                    if path.is_file()}}
        write(report_path, report)
        for mode in ("strong", "weak"):
            for count in ranks:
                case, database = cases[(mode, count)]
                check_dir = output / f"{mode}-np{count}-geometry"
                check_dir.mkdir()
                check = launch(launcher + ["-np", str(count), str(args.checker), str(database)],
                               check_dir, environment, args.timeout)
                text = (check_dir / "launcher.stdout").read_text()
                match = re.search(r"^elements=(\d+)\s+minimum_detJ=([^ ]+)\s+bad_elements=(\d+)",
                                  text, re.MULTILINE)
                if check["returncode"] or not match or int(match.group(3)):
                    raise RuntimeError(f"{mode} geometry preflight failed for {count} ranks")
                report["inputs"][f"{mode}-{count}"].update(
                    elements=int(match.group(1)), elements_per_rank=int(match.group(1)) / count,
                    minimum_det_j=float(match.group(2)))
        weak_work = [report["inputs"][f"weak-{count}"]["elements_per_rank"] for count in ranks]
        if max(weak_work) / min(weak_work) - 1 > args.weak_workload_relative_tolerance:
            raise ValueError("weak cases do not keep elements per rank within the configured tolerance")
        reference = None
        for repetition in range(args.repetitions):
            for mode in ("strong", "weak"):
                order = ranks[repetition % len(ranks):] + ranks[:repetition % len(ranks)]
                for count in order:
                    case, database = cases[(mode, count)]
                    directory = output / f"{mode}-np{count}-repeat{repetition}"
                    directory.mkdir()
                    field = directory / "field.txt"
                    command = launcher + ["-np", str(count), sys.executable,
                        str(ROOT / "scripts/hpc_rank_run.py"), "--output-dir", str(directory),
                        "--expected-ranks", str(count), "--timeout", str(args.timeout), "--",
                        str(args.binary.resolve()), str(database), str(case), "--output", str(field),
                        "--nonlinear-rtol", "1e-8", "--nonlinear-atol", "1e-12", "--mass-rtol", "1e-6"]
                    result = launch(command, directory, environment, args.timeout + 30)
                    row = {"mode": mode, "ranks": count, "repetition": repetition,
                           "directory": str(directory), "launch": result, "status": "failed"}
                    report["runs"].append(row)
                    write(report_path, report)
                    if result["returncode"]:
                        raise RuntimeError(f"{mode} solver failed for {count} ranks")
                    validation = subprocess.run([str(args.validator.resolve()), str(database), str(field)],
                                                cwd=ROOT, text=True, stdout=subprocess.PIPE,
                                                stderr=subprocess.STDOUT, timeout=args.timeout)
                    row["physical_validation"] = {"returncode": validation.returncode,
                                                   "output": validation.stdout.rstrip()}
                    if validation.returncode:
                        raise RuntimeError(f"{mode} physical validation failed for {count} ranks")
                    row["profile"] = summarize(directory)
                    row["solver_iterations"] = iterations(directory)
                    if mode == "strong":
                        if count == 1 and reference is None:
                            reference = (field, field.with_suffix(field.suffix + ".pressure"))
                        comparison = {
                            "velocity": compare(reference[0], field, args.rtol, args.zero_atol),
                            "pressure": compare(reference[1], field.with_suffix(field.suffix + ".pressure"),
                                                args.rtol, args.zero_atol),
                        }
                        row["field_comparison"] = comparison
                        if not all(value["passed"] for value in comparison.values()):
                            raise RuntimeError("strong-scaling field comparison failed")
                    row["status"] = "passed"
                    write(report_path, report)
        for mode in ("strong", "weak"):
            baseline = None
            mode_summary = {}
            for count in ranks:
                chosen = [row for row in report["runs"] if row["mode"] == mode and row["ranks"] == count]
                wall = samples([row["profile"]["process_wall_s"]["max"] for row in chosen])
                if count == 1:
                    baseline = wall["median"]
                phases = sorted(chosen[0]["profile"]["phases"])
                mode_summary[str(count)] = {
                    "elements": report["inputs"][f"{mode}-{count}"]["elements"],
                    "elements_per_rank": report["inputs"][f"{mode}-{count}"]["elements_per_rank"],
                    "max_rank_process_wall_s": wall,
                    "max_rank_peak_rss_bytes": samples([row["profile"]["peak_rss_bytes"]["max"] for row in chosen]),
                    "sum_individual_peak_rss_bytes": samples([row["profile"]["sum_individual_peak_rss_bytes"] for row in chosen]),
                    "solver_iterations": [row["solver_iterations"] for row in chosen],
                    "phases_max_rank_exclusive_s": {phase: samples([
                        row["profile"]["phases"][phase]["exclusive_s"]["max"] for row in chosen])
                        for phase in phases},
                }
            for count in ranks:
                item = mode_summary[str(count)]
                item["speedup_vs_one_rank"] = baseline / item["max_rank_process_wall_s"]["median"]
                item["parallel_efficiency"] = item["speedup_vs_one_rank"] / count
            report["summary"][mode] = mode_summary
        for item in report["inputs"].values():
            if item["database_sha256"] != digest(Path(item["database"])):
                raise RuntimeError("scaling database changed during measurement")
            if item["case_sha256"] != {str(path): digest(path)
                    for path in sorted(Path(item["case"]).rglob("*")) if path.is_file()}:
                raise RuntimeError("scaling case changed during measurement")
        if report["provenance_sha256"] != {str(path): digest(path) for path in provenance_paths}:
            raise RuntimeError("scaling executable or collector changed during measurement")
        report["status"] = "passed"
    except BaseException as error:
        failure = str(error)
        report["status"] = "failed"
        report["failure"] = failure
        raise
    finally:
        write(report_path, report)


if __name__ == "__main__":
    main()
