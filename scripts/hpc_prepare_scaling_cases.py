#!/usr/bin/env python3
"""Generate nontrivial fixed and constant-elements-per-rank duct cases."""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from hpc_inventory import digest


ROOT = Path(__file__).resolve().parents[1]


def weak_dimensions(ranks, elements_per_rank, initial_transverse):
    total = ranks * elements_per_rank
    transverse = initial_transverse
    while transverse <= 128:
        square = transverse * transverse
        if total % square == 0 and 1 <= total // square <= 128:
            return transverse, total // square
        transverse *= 2
    raise ValueError("cannot factor weak workload into supported duct dimensions")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--ranks", type=int, nargs="+", required=True)
    parser.add_argument("--strong-transverse", type=int, default=16)
    parser.add_argument("--strong-axial", type=int, default=64)
    parser.add_argument("--weak-elements-per-rank", type=int, default=256)
    parser.add_argument("--weak-initial-transverse", type=int, default=8)
    parser.add_argument("--fixture-builder", type=Path,
                        default=ROOT / "solvers/coupling/hpc_duct_solver_fixture")
    args = parser.parse_args()
    ranks = sorted(args.ranks)
    values = [*ranks, args.strong_transverse, args.strong_axial,
              args.weak_elements_per_rank, args.weak_initial_transverse]
    if (not ranks or ranks[0] != 1 or len(ranks) != len(set(ranks)) or
            any(value < 1 for value in values) or
            max(args.strong_transverse, args.strong_axial,
                args.weak_initial_transverse) > 128):
        parser.error("require unique positive ranks including 1 and dimensions in [1,128]")
    builder = args.fixture_builder.resolve()
    if not builder.is_file():
        parser.error(f"fixture builder is missing: {builder}")
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    report = {"schema_version": 1, "kind": "hpc_scaling_cases",
              "created_utc": datetime.now(timezone.utc).isoformat(),
              "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"],
                                                        cwd=ROOT, text=True).strip(),
              "fixture_builder": str(builder), "fixture_builder_sha256": digest(builder),
              "ranks": ranks, "strong": {}, "weak": {}, "status": "running"}
    manifest = root / "manifest.json"
    manifest.write_text(json.dumps(report, indent=2) + "\n")
    try:
        strong_elements = args.strong_transverse ** 2 * args.strong_axial
        if strong_elements < max(ranks) * 32:
            raise ValueError("strong case must provide at least 32 elements per maximum rank")
        for count in ranks:
            configurations = {
                "strong": (args.strong_transverse, args.strong_axial),
                "weak": weak_dimensions(count, args.weak_elements_per_rank,
                                        args.weak_initial_transverse),
            }
            for mode, (transverse, axial) in configurations.items():
                case = root / f"{mode}-{count}"
                result = subprocess.run([str(builder), str(case), str(transverse),
                                         str(axial), str(count)], cwd=ROOT, text=True,
                                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                if result.returncode:
                    raise RuntimeError(f"{mode}-{count} fixture failed: {result.stdout}")
                files = {str(path.relative_to(root)): digest(path)
                         for path in sorted(case.rglob("*")) if path.is_file()}
                report[mode][str(count)] = {
                    "case": str(case), "native_case": str(case / "root3d"),
                    "database": str(case / "duct.ntiga"), "transverse": transverse,
                    "axial": axial, "elements": transverse * transverse * axial,
                    "elements_per_rank": transverse * transverse * axial / count,
                    "stdout": result.stdout.rstrip(), "sha256": files,
                }
                manifest.write_text(json.dumps(report, indent=2) + "\n")
        if len({row["elements"] for row in report["strong"].values()}) != 1:
            raise RuntimeError("strong cases do not have one fixed element count")
        if len({row["elements_per_rank"] for row in report["weak"].values()}) != 1:
            raise RuntimeError("weak cases do not have constant elements per rank")
        if digest(builder) != report["fixture_builder_sha256"]:
            raise RuntimeError("fixture builder changed during case generation")
        for mode in ("strong", "weak"):
            for row in report[mode].values():
                current = {str(path.relative_to(root)): digest(path)
                           for path in sorted(Path(row["case"]).rglob("*")) if path.is_file()}
                if current != row["sha256"]:
                    raise RuntimeError("generated scaling case changed before publication")
        report["status"] = "passed"
    except BaseException as error:
        report["status"] = "failed"
        report["failure"] = str(error)
        raise
    finally:
        manifest.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
