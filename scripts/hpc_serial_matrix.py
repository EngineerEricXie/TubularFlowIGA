#!/usr/bin/env python3
"""Repeat the catalog's serial immersed/FSI and single-GPU benchmarks."""

import argparse
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import sys

from hpc_compare_fields import compare
from hpc_cpu_matrix import ROOT, launch, samples, write_json
from hpc_inventory import digest
from hpc_profile_summary import summarize


def finite_values(line):
    values = {key: float(value) for key, value in re.findall(r"([\w-]+)=([^\s]+)", line)
              if key != "converged"}
    if not values or any(not math.isfinite(value) for value in values.values()):
        raise ValueError("missing or nonfinite native diagnostics")
    return values


def native_diagnostics(kind, log):
    """Check complete diagnostic coverage; native code enforces unrounded gates."""
    lines = log.splitlines()
    if kind == "fsi":
        final = [line for line in lines if line.startswith("compliant_channel_fsi converged=true ")]
        history = [finite_values(line) for line in lines if line.startswith("fsi_iteration=")]
        if len(final) != 1 or not history:
            raise ValueError("missing FSI convergence diagnostics")
        values = finite_values(final[0])
        if (values["iterations"] != len(history) or
                [v["fsi_iteration"] for v in history] != list(range(len(history)))):
            raise ValueError("incomplete FSI iteration history")
        return {"kind": "native_fixture_gates", "final": values, "history": history,
                "scope": "Successful native fixture exit enforces full-precision force/moment, conservation, and coupling assertions; printed diagnostics are rounded."}
    fd = [line for line in lines if line.startswith("aneurysm zero-state centered-FD ")]
    expected = {(direction, block) for direction in ("velocity", "pressure", "controller") for block in range(3)}
    actual = []
    for line in fd:
        match = re.fullmatch(r"aneurysm zero-state centered-FD (\w+) output-block (\d+) relative-defect=(\S+)", line)
        if not match:
            raise ValueError("invalid immersed FD diagnostic")
        actual.append((match[1], int(match[2])))
        finite_values(line)
    solve = [line for line in lines if line.startswith("aneurysm static solve ")]
    conservation = [line for line in lines if line.startswith("aneurysm conservation ")]
    if len(fd) != 9 or set(actual) != expected or len(solve) != 1 or len(conservation) != 1:
        raise ValueError("incomplete immersed FD/solve/conservation diagnostics; solve-only is not this baseline")
    return {"kind": "native_fixture_gates", "fd": fd, "solve": finite_values(solve[0]),
            "conservation": finite_values(conservation[0]),
            "scope": "Successful native fixture exit enforces full-precision FD, Newton, geometry and conservation gates; timings include FD assemblies."}


def device_allocation(log):
    lines = [line for line in log.splitlines() if line.startswith("cuda_allocations ")]
    if len(lines) != 1:
        raise ValueError("expected exactly one CUDA allocation record")
    match = re.fullmatch(r"cuda_allocations scope=project_device_buffers requested_peak_bytes=(\d+) requested_live_bytes=(\d+)", lines[0])
    if not match or int(match[1]) <= 0 or int(match[2]) != 0:
        raise ValueError("invalid CUDA peak or project buffers still live at shutdown")
    return {"scope": "project_device_buffers", "requested_peak_bytes": int(match[1]),
            "requested_live_bytes": int(match[2]),
            "note": "Excludes driver/library allocations and allocator rounding; not whole-device memory."}


def aggregate_serial(records, repetitions, gpu):
    if (repetitions < 3 or len(records) != repetitions + 1 or
            sorted(r["repetition"] for r in records) != list(range(repetitions + 1)) or
            any(r["status"] != "accepted" for r in records)):
        raise ValueError("incomplete or failed serial matrix")
    records = sorted(records, key=lambda r: r["repetition"])
    repeated = records[1:]
    phases = set(records[0]["profile"]["phases"])
    if any(r["profile"]["ranks"] != 1 or set(r["profile"]["phases"]) != phases for r in records):
        raise ValueError("serial profile ranks or phase schemas differ")
    result = {"first_run": records[0]["directory"],
              "launcher_wall_s": samples([r["launch"]["wall_s"] for r in repeated]),
              "process_wall_s": samples([r["profile"]["process_wall_s"]["max"] for r in repeated]),
              "application_elapsed_s": samples([r["profile"]["application_elapsed_s"]["max"] for r in repeated]),
              "peak_rss_bytes": samples([r["profile"]["peak_rss_bytes"]["max"] for r in repeated]),
              "phases_exclusive_s": {name: samples([r["profile"]["phases"][name]["exclusive_s"]["max"] for r in repeated])
                                     for name in sorted(phases)}}
    if gpu:
        result["project_device_buffer_peak_bytes"] = samples([r["device_allocation"]["requested_peak_bytes"] for r in repeated])
    return result


def file_hashes(paths):
    return {str(path): digest(path) for path in sorted(set(paths))}


def cpu_reference(directory, system):
    """Require an intact accepted CPU matrix and its recorded input/binary snapshot."""
    metadata = json.loads((directory / "matrix.json").read_text())
    results = json.loads((directory / "summary.json").read_text())
    if (metadata["kind"] != "hpc_cpu_matrix" or metadata["system"] != system or
            results["status"] != "field_comparisons_passed"):
        raise ValueError("CPU reference must be an accepted matrix for the same system")
    for path, expected in metadata["input_and_binary_sha256"].items():
        if digest(Path(path)) != expected:
            raise ValueError("CPU reference inputs, binaries, or collection scripts changed: " + path)
    first = [r for r in results["observations"] if r["ranks"] == 1 and r["repetition"] == 0]
    if len(first) != 1 or first[0]["status"] != "field_comparison_passed":
        raise ValueError("missing accepted one-rank first reference")
    first = first[0]
    reference = Path(first["directory"])
    if summarize(reference) != first["profile"]:
        raise ValueError("CPU reference profile differs from recorded observation")
    fields = ["field.txt", "field.txt.pressure"] if system == "flow" else ["field.txt"]
    for name in fields:
        if digest(reference / name) != first["fields"][name]["candidate_sha256"]:
            raise ValueError("CPU reference field changed")
    return metadata, reference, fields


def run(options):
    gpu = options.case.startswith("cuda-")
    if options.repetitions < 3 or not math.isfinite(options.timeout) or options.timeout <= 0:
        raise ValueError("require at least three repetitions and a positive finite timeout")
    if options.cpu not in os.sched_getaffinity(0):
        raise ValueError("selected logical CPU is outside controller affinity")
    output = options.output_dir.resolve()
    policy_path = ROOT / "benchmarks/hpc_baselines.json"
    catalog = json.loads(policy_path.read_text())
    paths = [policy_path] + [ROOT / "scripts" / name for name in
                            ("hpc_serial_matrix.py", "hpc_cpu_matrix.py", "hpc_rank_run.py",
                             "hpc_profile_summary.py", "hpc_compare_fields.py", "hpc_inventory.py")]
    reference = None
    fields = []
    case_dir = None
    if gpu:
        if options.cpu_matrix is None or options.case_dir is not None:
            raise ValueError("GPU mode requires --cpu-matrix and no --case-dir")
        cpu_dir = options.cpu_matrix.resolve(strict=True)
        system = "flow" if options.case == "cuda-flow" else "neuron_transport"
        cpu_metadata, reference, fields = cpu_reference(cpu_dir, system)
        case_dir = Path(cpu_metadata["case_dir"])
        stem = "straight_tube" if system == "flow" else "straight_neurite"
        binary = ROOT / "solvers/cuda/iga_cuda"
        base_command = [str(binary), "navier-stokes" if system == "flow" else "solve",
                        str(case_dir / f"{stem}-1.ntiga"), str(case_dir)]
        if system == "flow":
            base_command += ["--max-newton", "30", "--nonlinear-rtol", "1e-8", "--nonlinear-atol", "1e-12", "--mass-rtol", "1e-6"]
        else:
            base_command += ["--system", "neuron_transport"]
        paths += [cpu_dir / "matrix.json", cpu_dir / "summary.json"] + [reference / name for name in fields]
        paths += [Path(path) for path in cpu_metadata["input_and_binary_sha256"]]
    elif options.case == "immersed":
        if options.case_dir is None or options.cpu_matrix is not None:
            raise ValueError("immersed mode requires --case-dir and no --cpu-matrix")
        case_dir = options.case_dir.resolve(strict=True)
        geometry = json.loads((case_dir / "immersed_geometry.json").read_text())
        if geometry["volume_quadrature"]["max_depth"] != 2:
            raise ValueError("this catalog benchmark requires an isolated depth-2 fixture")
        binary = ROOT / "solvers/cpu/immersed_aneurysm_jacobian_test"
        base_command = [str(binary), str(case_dir)]
        paths += [ROOT / "solvers/cpu/tests/test_immersed_aneurysm_jacobian.cpp"]
    else:
        if options.case_dir is not None or options.cpu_matrix is not None:
            raise ValueError("FSI uses the in-memory fixture; no input directory arguments")
        binary = ROOT / "solvers/cpu/compliant_channel_fsi_test"
        base_command = [str(binary)]
        paths += [ROOT / "solvers/cpu/tests" / name for name in
                  ("test_compliant_channel_fsi.cpp", "CompliantChannelFsiFixture.hpp")]
    if case_dir is not None:
        if output.is_relative_to(case_dir) or case_dir.is_relative_to(output):
            raise ValueError("case and output directories must be separate")
        paths += [p for p in case_dir.rglob("*") if p.is_file()]
    paths += [binary]
    before = file_hashes(paths)
    environment = os.environ.copy()
    overrides = {"OMP_NUM_THREADS": "1", "OMP_THREAD_LIMIT": "1", "OMP_DYNAMIC": "FALSE",
                 "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "BLIS_NUM_THREADS": "1",
                 "IGA_PROFILE": "1", "PETSC_OPTIONS": ""}
    if gpu:
        overrides["CUDA_VISIBLE_DEVICES"] = options.gpu
    environment.update(overrides)
    output.mkdir(parents=True, exist_ok=False)
    metadata = {"schema_version": 1, "kind": "hpc_serial_matrix", "case": options.case,
                "created_utc": datetime.now(timezone.utc).isoformat(), "repetitions": options.repetitions,
                "logical_cpu": options.cpu, "environment_overrides": overrides,
                "cuda_library_path": environment.get("LD_LIBRARY_PATH") if gpu else None,
                "input_binary_and_tool_sha256": before, "comparison_policy": catalog["comparison_policy"],
                "cpu_reference": str(reference) if reference else None,
                "notes": ["First run is a separate process, not a controlled cold-cache run; excluded from repeated statistics.",
                          "All runs in this matrix are sequential, with one host logical CPU and one requested BLAS/OpenMP thread.",
                          "Native fixture success enforces compiled assertions; field agreement alone is not independent physical validation.",
                          "CUDA peak covers project device buffers; host RSS and library/driver device allocations are different scopes.",
                          "Build and input preparation are outside the measured interval."]}
    write_json(output / "matrix.json", metadata)
    records = []
    failure = None
    totals = None
    try:
        if gpu:
            probe_dir = output / "device-probe"
            probe_dir.mkdir()
            probe_result = launch([str(binary), "device-info"], probe_dir, environment, 30)
            if probe_result["returncode"] != 0:
                raise ValueError("CUDA device probe failed; inspect device-probe logs")
        for repetition in range(options.repetitions + 1):
            directory = output / ("first" if repetition == 0 else f"repeat{repetition}")
            directory.mkdir()
            record = {"repetition": repetition, "directory": str(directory), "status": "failed"}
            records.append(record)
            command = base_command + (["--output", str(directory / "field.txt")] if gpu else [])
            argv = ["taskset", "--cpu-list", str(options.cpu), sys.executable, str(ROOT / "scripts/hpc_rank_run.py"),
                    "--output-dir", str(directory), "--expected-ranks", "1", "--timeout", str(options.timeout), "--"] + command
            print(f"running {options.case} {directory.name}", flush=True)
            record["launch"] = launch(argv, directory, environment, options.timeout + 30)
            if record["launch"]["returncode"] != 0:
                raise ValueError("fixture failed: " + directory.name)
            record["profile"] = summarize(directory)
            log = (directory / "rank-0/stdout.log").read_text()
            if gpu:
                record["device_allocation"] = device_allocation(log)
                record["fields"] = {name: compare(reference / name, directory / name,
                                                catalog["comparison_policy"]["cpu_cuda_field_relative_l2_max"],
                                                catalog["comparison_policy"]["zero_reference_absolute_l2_max"],
                                                options.case == "cuda-transport") for name in fields}
                if not all(c["passed"] for c in record["fields"].values()):
                    raise ValueError("CPU/CUDA field comparison failed: " + directory.name)
            else:
                record["native_diagnostics"] = native_diagnostics(options.case, log)
            record["status"] = "accepted"
            write_json(directory / "observation.json", record)
            print(f"accepted {directory.name}: process {record['profile']['process_wall_s']['max']:.3f} s", flush=True)
        if file_hashes(paths) != before:
            raise ValueError("inputs, reference, binaries or collection scripts changed during measurement")
        totals = aggregate_serial(records, options.repetitions, gpu)
    except (OSError, ValueError, KeyError, TypeError) as error:
        failure = str(error)
    write_json(output / "summary.json", {"schema_version": 1, "kind": "hpc_serial_matrix_results",
                                        "case": options.case, "status": "failed" if failure else "accepted",
                                        "error": failure, "observations": records, "repeated_statistics": totals})
    if failure:
        print(failure, file=sys.stderr)
    return 1 if failure else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=("immersed", "fsi", "cuda-flow", "cuda-transport"), required=True)
    parser.add_argument("--case-dir", type=Path)
    parser.add_argument("--cpu-matrix", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=1500)
    parser.add_argument("--cpu", type=int, default=min(os.sched_getaffinity(0)))
    parser.add_argument("--gpu", default="0", help="CUDA_VISIBLE_DEVICES selector for GPU mode")
    return run(parser.parse_args())


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError, TypeError) as error:
        print(f"hpc_serial_matrix: {error}", file=sys.stderr)
        sys.exit(1)
