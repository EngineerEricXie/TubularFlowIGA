#!/usr/bin/env python3
"""Compare shared and per-domain MPI communicators on a two-duct graph."""

import argparse
import csv
import json
import math
import os
from pathlib import Path
import shutil
import statistics
import subprocess
import sys

from hpc_inventory import digest
from hpc_profile_summary import summarize


def graph(resources=None):
    port = lambda domain, name, locator, required: {
        "id": name, "locator_kind": "boundary_label" if domain.startswith("duct") else
        ("zero_d_port" if domain.startswith("rcr") or domain == "source" else "runtime_port"),
        "locator": locator, "provides": (["mean_pressure", "flow_rate"] if
        domain == "source" or domain.startswith("rcr") else ["area", "flow_rate", "mean_pressure"]),
        "requires": required}
    domains = [
        {"id": "source", "dimension": "0d", "kind": "zero_d_flow", "case": "source",
         "zero_d_model": "zero_d_model.json", "ports": [port("source", "port", "port", ["mean_pressure"])]},
        {"id": "duct_a", "dimension": "3d", "kind": "body_fitted_iga_flow", "case": "duct_a",
         "database": "duct_a.ntiga", "ports": [port("duct_a", "inlet", "1", ["flow_rate"]),
         port("duct_a", "outlet", "2", ["mean_pressure"]),
         {**port("duct_a", "wall", "0", []), "provides": ["flow_rate"]}]},
        {"id": "tree_a", "dimension": "1d", "kind": "network_flow", "case": "tree_a",
         "inlet_policy": "coupled_root", "ports": [port("tree_a", "root", "root", ["flow_rate"]),
         port("tree_a", "terminal_a", "outlet:2", ["mean_pressure"]),
         port("tree_a", "terminal_b", "outlet:3", ["mean_pressure"])]},
        {"id": "duct_b", "dimension": "3d", "kind": "body_fitted_iga_flow", "case": "duct_b",
         "database": "duct_b.ntiga", "ports": [port("duct_b", "inlet", "1", ["flow_rate"]),
         port("duct_b", "outlet", "2", ["mean_pressure"]),
         {**port("duct_b", "wall", "0", []), "provides": ["flow_rate"]}]},
        {"id": "tree_b", "dimension": "1d", "kind": "network_flow", "case": "tree_b",
         "inlet_policy": "coupled_root", "ports": [port("tree_b", "root", "root", ["flow_rate"]),
         port("tree_b", "terminal_a", "outlet:2", ["mean_pressure"]),
         port("tree_b", "terminal_b", "outlet:3", ["mean_pressure"]) ]},
        {"id": "tree_c", "dimension": "1d", "kind": "network_flow", "case": "tree_c",
         "inlet_policy": "coupled_root", "ports": [port("tree_c", "root", "root", ["flow_rate"]),
         port("tree_c", "terminal_a", "outlet:2", ["mean_pressure"]),
         port("tree_c", "terminal_b", "outlet:3", ["mean_pressure"])]}]
    for name in ("rcr_a", "rcr_b", "rcr_c", "rcr_d"):
        domains.append({"id": name, "dimension": "0d", "kind": "zero_d_flow", "case": name,
                        "zero_d_model": "zero_d_model.json",
                        "ports": [port(name, "port", "port", ["flow_rate"])]})
    links = [("source_tree", "source", "port", "tree_a", "root"),
             ("tree_duct_a", "tree_a", "terminal_a", "duct_a", "inlet"),
             ("tree_duct_b", "tree_a", "terminal_b", "duct_b", "inlet"),
             ("duct_a_tree", "duct_a", "outlet", "tree_b", "root"),
             ("duct_b_tree", "duct_b", "outlet", "tree_c", "root"),
             ("tree_b_rcr_a", "tree_b", "terminal_a", "rcr_a", "port"),
             ("tree_b_rcr_b", "tree_b", "terminal_b", "rcr_b", "port"),
             ("tree_c_rcr_c", "tree_c", "terminal_a", "rcr_c", "port"),
             ("tree_c_rcr_d", "tree_c", "terminal_b", "rcr_d", "port")]
    result = {"schema_version": 5, "time": {"dt": .01, "steps": 2}, "start_domain": "source",
              "execution": {"kind": "explicit", "maximum_iterations": 1,
              "pressure_relative_tolerance": 1e-6, "pressure_reference_pa": 1,
              "flow_relative_tolerance": 1e-10, "relaxation_factor": .5,
              "minimum_relaxation": .5, "maximum_relaxation": .5}, "domains": domains,
              "couplings": [{"id": edge, "a": {"domain": a, "port": ap},
              "b": {"domain": b, "port": bp}, "mode": "pressure_flow", "initial_pressure_pa": 0}
              for edge, a, ap, b, bp in links]}
    if resources:
        result["resources"] = resources
    return result


def prepare(repo, fixture, root, ranks, grouped, transverse, axial):
    subprocess.run([str(fixture), str(root), str(transverse), str(axial), str(ranks)],
                   cwd=repo, check=True)
    shutil.copytree(root / "root3d", root / "duct_a")
    shutil.copytree(root / "root3d", root / "duct_b")
    shutil.copy2(root / "duct.ntiga", root / "duct_a.ntiga")
    shutil.copy2(root / "duct.ntiga", root / "duct_b.ntiga")
    shutil.copytree(root / "tree", root / "tree_a")
    shutil.copytree(root / "tree", root / "tree_b")
    shutil.copytree(root / "tree", root / "tree_c")
    shutil.copytree(root / "rcr_a", root / "rcr_c")
    shutil.copytree(root / "rcr_a", root / "rcr_d")
    resources = None
    if grouped:
        resources = {"mode": "domain_groups", "groups": [
            {"id": "small", "ranks": 1,
             "domains": ["source", "tree_a", "tree_b", "tree_c", "rcr_a", "rcr_b", "rcr_c", "rcr_d"]},
            {"id": "duct_a", "ranks": 2, "domains": ["duct_a"]},
            {"id": "duct_b", "ranks": 2, "domains": ["duct_b"]}]}
    (root / "simulation_config.json").write_text(json.dumps(graph(resources), indent=2) + "\n")


def compare(reference, current):
    maximum = 0.0
    for name in ("pressure_flow_edges.csv", "pressure_flow_ports.csv",
                 "pressure_flow_domain_balances.csv", "zero_d_flow_history.csv"):
        with (reference / name).open() as stream:
            expected = list(csv.DictReader(stream))
        with (current / name).open() as stream:
            actual = list(csv.DictReader(stream))
        if len(expected) != len(actual):
            raise ValueError(f"row count differs for {name}")
        for first, second in zip(expected, actual):
            if first.keys() != second.keys():
                raise ValueError(f"columns differ for {name}")
            for key in first:
                try:
                    a, b = float(first[key]), float(second[key])
                except ValueError:
                    if first[key] != second[key]:
                        raise ValueError(f"identity differs for {name}.{key}")
                    continue
                error = abs(a-b) / max(1.0, abs(a), abs(b))
                if math.isfinite(error):
                    maximum = max(maximum, error)
                else:
                    raise ValueError(f"nonfinite comparison for {name}.{key}")
    if maximum > 1e-6:
        raise ValueError(f"field/history difference {maximum} exceeds 1e-6")
    return maximum


def main():
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=repo / "solvers/coupling/iga_multidomain_flow")
    parser.add_argument("--fixture", type=Path, default=repo / "solvers/coupling/hpc_duct_solver_fixture")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--transverse", type=int, default=4)
    parser.add_argument("--axial", type=int, default=8)
    parser.add_argument("--timeout", type=int, default=1200)
    parser.add_argument("--launcher", default="mpiexec --map-by core --bind-to core")
    args = parser.parse_args()
    if min(args.repeats, args.transverse, args.axial, args.timeout) < 1:
        parser.error("repeat, mesh and timeout values must be positive")
    root, binary, fixture = args.output_dir.resolve(), args.binary.resolve(), args.fixture.resolve()
    root.mkdir(parents=True, exist_ok=False)
    cases = root / "cases"
    cases.mkdir()
    prepare(repo, fixture, cases / "shared", 5, False, args.transverse, args.axial)
    prepare(repo, fixture, cases / "grouped", 2, True, args.transverse, args.axial)
    tracked = [binary, fixture, Path(__file__).resolve(), repo / "scripts/hpc_rank_run.py",
               repo / "scripts/hpc_profile_summary.py"]
    report = {"schema_version": 1, "kind": "hpc_domain_group_scaling", "status": "running",
              "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo,
              text=True).strip(), "source_dirty": bool(subprocess.check_output(
              ["git", "status", "--porcelain"], cwd=repo, text=True).strip()),
              "files": {str(path): digest(path) for path in tracked},
              "parameters": {"ranks": 5, "three_d_domains": 2,
              "elements_per_three_d_domain": args.transverse ** 2 * args.axial,
              "repeats": args.repeats, "launcher": args.launcher}, "runs": []}
    summary_path = root / "summary.json"
    def save():
        summary_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    try:
        save()
        environment = os.environ.copy()
        environment.update({"IGA_PROFILE": "1", "OMP_NUM_THREADS": "1",
                            "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"})
        launcher = args.launcher.split()
        for repeat in range(args.repeats):
            for mode in ("shared", "grouped"):
                run_dir = root / f"{mode}-{repeat}"
                run_dir.mkdir()
                rank_dir = run_dir / "ranks"
                rank_dir.mkdir()
                output = run_dir / "result"
                command = launcher + ["-np", "5", sys.executable,
                    str(repo / "scripts/hpc_rank_run.py"), "--output-dir", str(rank_dir),
                    "--expected-ranks", "5", "--timeout", str(args.timeout), "--",
                    str(binary), "--graph-case", str(cases / mode), "--output-dir", str(output)]
                completed = subprocess.run(command, cwd=repo, env=environment)
                profile = summarize(rank_dir)
                comparison = None
                if mode == "grouped":
                    comparison = compare(root / f"shared-{repeat}" / "result", output)
                row = {"mode": mode, "repeat": repeat, "returncode": completed.returncode,
                       "profile": profile, "numerical_maximum": comparison,
                       "output_manifest_sha256": digest(output / "graph_binding_manifest.json")}
                report["runs"].append(row)
                save()
                if completed.returncode:
                    raise RuntimeError(f"{mode} repeat {repeat} failed")
        report["summary"] = {}
        for mode in ("shared", "grouped"):
            selected = [row["profile"] for row in report["runs"] if row["mode"] == mode]
            report["summary"][mode] = {
                "median_max_process_wall_s": statistics.median(x["process_wall_s"]["max"] for x in selected),
                "median_max_application_s": statistics.median(x["application_elapsed_s"]["max"] for x in selected),
                "median_sum_individual_peak_rss_bytes": statistics.median(
                    x["sum_individual_peak_rss_bytes"] for x in selected),
                "median_max_rank_peak_rss_bytes": statistics.median(x["peak_rss_bytes"]["max"] for x in selected),
                "median_max_communication_s": statistics.median(
                    x["phases"]["communication"]["exclusive_s"]["max"] for x in selected),
                "median_mean_unscoped_s": statistics.median(x["unscoped_s"]["mean"] for x in selected)}
        shared, grouped = report["summary"]["shared"], report["summary"]["grouped"]
        report["comparison"] = {
            "wall_ratio_grouped_over_shared": grouped["median_max_process_wall_s"] /
            shared["median_max_process_wall_s"],
            "rss_sum_ratio_grouped_over_shared": grouped["median_sum_individual_peak_rss_bytes"] /
            shared["median_sum_individual_peak_rss_bytes"],
            "decision": "shared remains the default; domain_groups is opt-in for graphs whose per-domain memory or concurrency justifies coordination overhead"}
        if any(row["numerical_maximum"] is not None and row["numerical_maximum"] > 1e-6
               for row in report["runs"]):
            raise RuntimeError("numerical gate failed")
        report["status"] = "passed"
    except BaseException:
        report["status"] = "failed"
        raise
    finally:
        save()


if __name__ == "__main__":
    main()
