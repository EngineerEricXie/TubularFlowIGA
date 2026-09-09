#!/usr/bin/env python3
"""Run a sequential workstation MPI matrix on prepacked body-fitted inputs.

Compilation and packing are deliberately outside the measured interval. Every
configuration gets a first process followed by at least three fresh processes.
The first one-rank field is the common numerical reference, never a timing sample.
"""

import argparse
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import signal
import statistics
import subprocess
import sys
import time

from hpc_compare_fields import compare
from hpc_inventory import digest, probe
from hpc_profile_summary import summarize


ROOT = Path(__file__).resolve().parents[1]


def write_json(path, value):
    with path.open("x") as stream:
        stream.write(json.dumps(value, indent=2, allow_nan=False) + "\n")


def samples(values):
    if len(values) < 3 or any(not math.isfinite(v) or v < 0 for v in values):
        raise ValueError("need at least three finite, nonnegative observations")
    return {"values": values, "count": len(values), "median": statistics.median(values),
            "min": min(values), "max": max(values), "population_stdev": statistics.pstdev(values)}


def aggregate(records, ranks, repetitions):
    """Fail closed: missing, duplicate, or failed observations cannot be dropped."""
    expected = {(n, i) for n in ranks for i in range(repetitions + 1)}
    if (len(records) != len(expected) or
            {(r["ranks"], r["repetition"]) for r in records} != expected or
            any(r["status"] != "field_comparison_passed" for r in records)):
        raise ValueError("matrix is incomplete or contains failed observations")
    result = {}
    for n in ranks:
        group = sorted((r for r in records if r["ranks"] == n), key=lambda r: r["repetition"])
        repeated = group[1:]
        names = set(group[0]["profile"]["phases"])
        if any(set(r["profile"]["phases"]) != names for r in group):
            raise ValueError("phase schemas differ between repetitions")
        result[str(n)] = {
            "first_run": group[0]["directory"],
            "launcher_wall_s": samples([r["launch"]["wall_s"] for r in repeated]),
            "max_rank_process_wall_s": samples([r["profile"]["process_wall_s"]["max"] for r in repeated]),
            "max_rank_application_s": samples([r["profile"]["application_elapsed_s"]["max"] for r in repeated]),
            "max_rank_peak_rss_bytes": samples([r["profile"]["peak_rss_bytes"]["max"] for r in repeated]),
            "sum_individual_peak_rss_bytes": samples([r["profile"]["sum_individual_peak_rss_bytes"] for r in repeated]),
            "phases_max_rank_exclusive_s": {
                name: samples([r["profile"]["phases"][name]["exclusive_s"]["max"] for r in repeated])
                for name in sorted(names)},
        }
    baseline = result["1"]["launcher_wall_s"]["median"]
    memory = result["1"]["sum_individual_peak_rss_bytes"]["median"]
    for n in ranks:
        item = result[str(n)]
        elapsed = item["launcher_wall_s"]["median"]
        item["observed_speedup"] = baseline / elapsed if elapsed else None
        item["observed_parallel_efficiency"] = baseline / elapsed / n if elapsed else None
        item["sum_peak_rss_ratio_to_one_rank"] = item["sum_individual_peak_rss_bytes"]["median"] / memory if memory else None
    return result


def launch(argv, directory, environment, timeout):
    """Preserve failures as observations, with an outer deadline for MPI startup."""
    start = time.perf_counter()
    timed_out = False
    with (directory / "launcher.stdout").open("x") as stdout, (directory / "launcher.stderr").open("x") as stderr:
        process = subprocess.Popen(argv, cwd=ROOT, env=environment, stdout=stdout,
                                   stderr=stderr, start_new_session=True)
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
    result = {"argv": argv, "returncode": returncode, "timed_out": timed_out,
              "wall_s": time.perf_counter() - start,
              "logs": {name: digest(directory / name) for name in ("launcher.stdout", "launcher.stderr")}}
    write_json(directory / "launch.json", result)
    return result


def snapshot(case, binaries):
    paths = sorted(p for p in case.rglob("*") if p.is_file()) + binaries
    return {str(p): digest(p) for p in paths}


def thread_plan(ranks, total_cores):
    if not ranks or any(n < 1 for n in ranks) or len(set(ranks)) != len(ranks):
        raise ValueError("require unique positive rank counts")
    if total_cores is None:
        return {n: 1 for n in ranks}
    if total_cores < 1 or total_cores not in ranks or any(total_cores % n for n in ranks):
        raise ValueError("fixed core count must be divisible by every rank count and include a pure MPI configuration")
    return {n: total_cores // n for n in ranks}


def hybrid_execution(directory, ranks, threads, topology_root=Path("/sys/devices/system/cpu")):
    """Check actual native teams and disjoint physical cores, not launcher intent."""
    result = []
    occupied = set()
    for rank in range(ranks):
        local = directory / f"rank-{rank}"
        run = json.loads((local / "run.json").read_text())
        if (run["rank"] != rank or run["ranks"] != ranks or run["returncode"] != 0
                or run["timed_out"] or run["status"] != "process_passed"):
            raise ValueError("invalid hybrid rank result")
        cores = set()
        for cpu in run["affinity_cpus"]:
            topology = topology_root / f"cpu{cpu}/topology"
            cores.add((int((topology / "physical_package_id").read_text()),
                       int((topology / "core_id").read_text())))
        if len(cores) != threads or occupied & cores:
            raise ValueError("hybrid rank bindings overlap or differ from requested physical core count")
        occupied.update(cores)
        reports = []
        for line in (local / "stdout.log").read_text().splitlines():
            if not line.startswith("body_fitted_element_assembly "):
                continue
            report = {k: int(v) for k, v in (token.split("=", 1) for token in line.split()[1:])}
            if (report["rank"] != rank or report["threads_requested"] != threads
                    or report["team_size"] != threads or report["elements"] < threads
                    or report["batches"] < 1 or not 1 <= report["maximum_resident_items"] <= min(threads, 8)):
                raise ValueError("actual body-fitted assembly team or batch differs from configuration")
            reports.append(report)
        if not reports:
            raise ValueError("missing native body-fitted assembly reports")
        result.append(dict(rank=rank, affinity_cpus=run["affinity_cpus"],
                           physical_cores=sorted(cores), assembly_reports=reports))
    return result


def run(options):
    case = options.case_dir.resolve(strict=True)
    output = options.output_dir.resolve()
    if output.is_relative_to(case) or case.is_relative_to(output):
        raise ValueError("case and output directories must be separate")
    ranks = sorted(set(options.ranks))
    if ranks[0] != 1 or len(ranks) != len(options.ranks) or options.repetitions < 3:
        raise ValueError("require unique positive ranks including 1, and at least three repetitions")
    if not math.isfinite(options.timeout) or options.timeout <= 0:
        raise ValueError("timeout must be positive and finite")
    total_cores = getattr(options, "total_cores", None)
    threads_by_rank = thread_plan(ranks, total_cores)
    if total_cores is not None and options.system != "flow":
        raise ValueError("fixed-core hybrid matrix currently supports body-fitted flow")
    if Path(options.database_stem).name != options.database_stem:
        raise ValueError("database stem must be a filename without directories")
    for n in ranks:
        if not (case / f"{options.database_stem}-{n}.ntiga").is_file():
            raise ValueError(f"missing prepacked database for {n} ranks")
    binary_name = "iga_navier_stokes_openmp" if total_cores is not None else "iga_navier_stokes"
    binary = ROOT / "solvers/cpu" / (binary_name if options.system == "flow" else "iga_solve")
    checker = ROOT / "solvers/cpu/iga_mesh_check"
    provenance_files = [binary, checker, ROOT / "benchmarks/hpc_baselines.json"] + [
        ROOT / "scripts" / name for name in ("hpc_cpu_matrix.py", "hpc_rank_run.py",
                                             "hpc_profile_summary.py", "hpc_compare_fields.py", "hpc_inventory.py")]
    provenance_files += list((ROOT / "include").glob("*.hpp")) + list((ROOT / "solvers/cpu/include").glob("*.hpp"))
    provenance_files += [ROOT / "solvers/cpu/Makefile", ROOT / "solvers/cpu/src/iga_navier_stokes.cpp"]
    before = snapshot(case, provenance_files)
    policy = json.loads((ROOT / "benchmarks/hpc_baselines.json").read_text())["comparison_policy"]
    environment = os.environ.copy()
    overrides = {"OMP_NUM_THREADS": "1", "OMP_THREAD_LIMIT": "1", "OMP_DYNAMIC": "FALSE",
                 "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "BLIS_NUM_THREADS": "1",
                 "IGA_PROFILE": "1", "PETSC_OPTIONS": ""}
    if options.system == "flow":
        overrides["PETSC_OPTIONS"] = "-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps"
    environment.update(overrides)
    output.mkdir(parents=True, exist_ok=False)
    metadata = {"schema_version": 1, "kind": "hpc_cpu_matrix", "created_utc": datetime.now(timezone.utc).isoformat(),
                "system": options.system, "case_dir": str(case), "ranks": ranks,
                "repetitions": options.repetitions, "threads_by_rank_count": threads_by_rank,
                "threads_per_rank": 1 if total_cores is None else None,
                "fixed_total_physical_cores": total_cores,
                "environment_overrides": overrides, "input_and_binary_sha256": before,
                "comparison_policy": policy, "controller_affinity": sorted(os.sched_getaffinity(0)),
                "cpu_topology": probe(["lscpu", "-e=CPU,CORE,SOCKET,ONLINE"]),
                "mpi_version": probe([options.mpiexec, "--version"]),
                "modes": {"mpi": "measured by this run",
                          "openmp_assembly": "measured by this run" if total_cores is not None else "not requested",
                          "hybrid_assembly": "measured by this run" if total_cores is not None else "not requested",
                          "single_gpu": "separate matrix"},
                "notes": ["First runs are separate fresh processes, not controlled cold filesystem caches.",
                          "Repetitions run sequentially, rotating configuration order each round.",
                          "Launcher wall includes MPI/Python startup, solver, output, and shutdown; preparation and field checks are excluded.",
                          "Phase maxima may occur on different ranks; their sum is not a critical-path wall time.",
                          "Sum of per-rank RSS peaks is not simultaneous aggregate RSS.",
                          "Field agreement and native process success do not replace independent physical validation.",
                          "Observed speedup is descriptive; isolation and physical gates must be audited before making a performance claim."]}
    write_json(output / "matrix.json", metadata)
    records = []
    preflights = []
    reference = output / "np1-first"
    failure = None
    try:
        for n in ranks:
            directory = output / f"np{n}-geometry"
            directory.mkdir()
            command = [options.mpiexec, "--bind-to", "core", "--map-by", "core", "--nooversubscribe", "-np", str(n),
                       str(checker), str(case / f"{options.database_stem}-{n}.ntiga")]
            check = launch(command, directory, environment, options.timeout)
            preflights.append({"ranks": n, "directory": str(directory), "launch": check})
            if check["returncode"] != 0:
                raise ValueError(f"geometry preflight failed for {n} ranks")
        for repetition in range(options.repetitions + 1):
            offset = (repetition - 1) % len(ranks) if repetition else 0
            order = ranks[offset:] + ranks[:offset]
            for n in order:
                threads = threads_by_rank[n]
                run_environment = environment.copy()
                run_environment.update(OMP_NUM_THREADS=str(threads), OMP_THREAD_LIMIT=str(threads),
                                       OMP_PROC_BIND="close", OMP_PLACES="cores",
                                       IGA_ASSEMBLY_THREADS=str(threads), IGA_ASSEMBLY_BATCH_SIZE=str(min(threads, 8)))
                name = f"np{n}-" + ("first" if repetition == 0 else f"repeat{repetition}")
                directory = output / name
                directory.mkdir()
                record = {"ranks": n, "threads": threads, "repetition": repetition,
                          "directory": str(directory), "status": "failed"}
                records.append(record)
                command = [str(binary), str(case / f"{options.database_stem}-{n}.ntiga"), str(case),
                           "--output", str(directory / "field.txt")]
                if options.system == "flow":
                    command += ["--nonlinear-rtol", "1e-8", "--nonlinear-atol", "1e-12", "--mass-rtol", "1e-6"]
                else:
                    command += ["--system", "neuron_transport"]
                mapping = f"slot:PE={threads}" if total_cores is not None else "core"
                argv = [options.mpiexec, "--bind-to", "core", "--map-by", mapping, "--nooversubscribe", "-np", str(n),
                        sys.executable, str(ROOT / "scripts/hpc_rank_run.py"), "--output-dir", str(directory),
                        "--expected-ranks", str(n), "--timeout", str(options.timeout), "--"] + command
                print(f"running {name}", flush=True)
                record["launch"] = launch(argv, directory, run_environment, options.timeout + 30)
                if record["launch"]["returncode"] != 0:
                    raise ValueError(f"solver launch failed: {name}")
                record["profile"] = summarize(directory)
                if total_cores is not None:
                    record["hybrid_execution"] = hybrid_execution(directory, n, threads)
                comparisons = {}
                record["fields"] = comparisons
                for field in (["field.txt", "field.txt.pressure"] if options.system == "flow" else ["field.txt"]):
                    comparisons[field] = compare(reference / field, directory / field,
                                                 policy["cpu_rank_field_relative_l2_max"],
                                                 policy["zero_reference_absolute_l2_max"], options.system != "flow")
                if not all(c["passed"] for c in comparisons.values()):
                    raise ValueError(f"field comparison failed: {name}")
                record["status"] = "field_comparison_passed"
                write_json(directory / "observation.json", record)
                print(f"passed {name}: launcher {record['launch']['wall_s']:.3f} s", flush=True)
        if snapshot(case, provenance_files) != before:
            raise ValueError("inputs or binaries changed during measurement")
        totals = aggregate(records, ranks, options.repetitions)
        if total_cores is not None:
            pure_mpi = totals[str(total_cores)]["launcher_wall_s"]["median"]
            pure_mpi_memory = totals[str(total_cores)]["sum_individual_peak_rss_bytes"]["median"]
            if pure_mpi <= 0 or pure_mpi_memory <= 0:
                raise ValueError("pure MPI reference time and memory must be positive")
            for n in ranks:
                item = totals[str(n)]
                if item["launcher_wall_s"]["median"] <= 0:
                    raise ValueError("hybrid observation time must be positive")
                item["threads_per_rank"] = threads_by_rank[n]
                item["total_physical_cores"] = total_cores
                item["observed_speedup_vs_pure_mpi"] = pure_mpi / item["launcher_wall_s"]["median"]
                item["sum_peak_rss_ratio_vs_pure_mpi"] = item["sum_individual_peak_rss_bytes"]["median"] / pure_mpi_memory
                # The serial aggregate's division by rank count is not a
                # parallel efficiency when every mode uses the same cores.
                item.pop("observed_parallel_efficiency")
                item["speedup_reference"] = "one MPI rank using all configured physical cores"
    except (OSError, ValueError, KeyError, TypeError) as error:
        failure = str(error)
        totals = None
    result = {"schema_version": 1, "kind": "hpc_cpu_matrix_results",
              "status": "failed" if failure else "field_comparisons_passed", "error": failure,
              "preflights": preflights, "observations": records, "repeated_statistics": totals}
    write_json(output / "summary.json", result)
    if failure:
        print(failure, file=sys.stderr)
    return 1 if failure else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-dir", required=True, type=Path)
    parser.add_argument("--database-stem", required=True)
    parser.add_argument("--system", choices=("flow", "neuron_transport"), required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--ranks", nargs="+", type=int, default=[1, 2, 4, 8])
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--total-cores", type=int,
                        help="fixed physical core budget for flow; ranks must include this count and divide it")
    parser.add_argument("--timeout", type=float, default=600)
    parser.add_argument("--mpiexec", default="mpiexec", help="Open MPI launcher executable")
    return run(parser.parse_args())


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError) as error:
        print(f"hpc_cpu_matrix: {error}", file=sys.stderr)
        sys.exit(1)
