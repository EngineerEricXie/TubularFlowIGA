#!/usr/bin/env python3
"""Capture reproducible HPC inputs and environment without running a simulation."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys


RESOURCE_ENV = (
    "OMP_NUM_THREADS", "OMP_PROC_BIND", "OMP_PLACES", "OMP_DYNAMIC",
    "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
    "SLURM_JOB_ID", "SLURM_NTASKS", "SLURM_CPUS_PER_TASK", "SLURM_JOB_NUM_NODES",
    "PETSC_DIR", "PETSC_ARCH", "PETSC_OPTIONS", "EIGEN_DIR", "CUDA_ARCHS", "IGA_PROFILE",
)


def probe(argv, cwd=None):
    """Failures are observations, never evidence of an absent installation."""
    try:
        result = subprocess.run(argv, cwd=cwd, capture_output=True, text=True,
                                errors="replace", timeout=15, check=False)
        return {"argv": argv, "status": "ok" if result.returncode == 0 else "failed",
                "returncode": result.returncode, "stdout": result.stdout.strip(),
                "stderr": result.stderr.strip()}
    except (OSError, subprocess.TimeoutExpired) as error:
        return {"argv": argv, "status": "unavailable", "error": str(error)}


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def file_records(root, names):
    result = []
    for name in sorted(set(names)):
        path = root / name
        if path.is_file():
            result.append({"path": name, "size_bytes": path.stat().st_size,
                           "sha256": digest(path)})
        else:
            result.append({"path": name, "missing": True})
    return result


def records_digest(records):
    payload = json.dumps(records, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def input_records(root, paths):
    names = []
    for relative in paths:
        path = (root / relative).resolve()
        if not path.is_relative_to(root.resolve()) or not path.exists():
            raise ValueError("fixture input must exist inside repository: " + relative)
        if path.is_dir():
            for item in path.rglob("*"):
                if item.is_file():
                    if not item.resolve().is_relative_to(root.resolve()):
                        raise ValueError("fixture symlink escapes repository: " + str(item))
                    names.append(str(item.relative_to(root)))
        else:
            names.append(str(path.relative_to(root)))
    return file_records(root, names)


def petsc_configuration(prefix):
    if not prefix:
        return {"status": "not_selected"}
    prefix = Path(prefix).resolve()
    # An explicit prefix must name the configured build, not guess among builds.
    conf = prefix / "include/petscconf.h"
    version = prefix / "include/petscversion.h"
    if not conf.is_file():
        return {"status": "unresolved", "prefix": str(prefix),
                "reason": "include/petscconf.h missing; supply configured PETSc prefix"}
    defines = {}
    for line in conf.read_text().splitlines():
        fields = line.split()
        if len(fields) >= 3 and fields[0] == "#define":
            defines[fields[1]] = " ".join(fields[2:])
    selected = ["PETSC_USE_64BIT_INDICES", "PETSC_USE_COMPLEX", "PETSC_USE_REAL_DOUBLE",
                "PETSC_HAVE_MUMPS", "PETSC_HAVE_HYPRE", "PETSC_HAVE_OPENMP",
                "PETSC_HAVE_HDF5", "PETSC_HAVE_CUDA", "PETSC_HAVE_MPIUNI"]
    return {"status": "ok", "prefix": str(prefix), "configuration_sha256": digest(conf),
            "version_header": version.read_text() if version.is_file() else None,
            "features": {name: defines.get(name) for name in selected}}


def collect(root, catalog_path, ranks, threads, petsc_prefix=None, nvcc=None, command=None):
    root = root.resolve()
    catalog = json.loads(catalog_path.read_text())
    tracked = subprocess.check_output(["git", "ls-files", "-z"], cwd=root).decode().split("\0")
    untracked = subprocess.check_output(
        ["git", "ls-files", "--others", "--exclude-standard", "-z"], cwd=root).decode().split("\0")
    suffixes = {".cpp", ".hpp", ".h", ".cu", ".cuh", ".py", ".sh", ".json", ".md"}
    added_sources = [name for name in untracked if name and
                     (Path(name).suffix in suffixes or Path(name).name in ("Makefile", "makefile"))
                     and not name.startswith(("results/", "outputs/"))]
    source = file_records(root, [name for name in tracked if name] + added_sources)
    fixtures = []
    for case in catalog["cases"]:
        files = input_records(root, case["source_inputs"])
        fixtures.append({**case, "input_files": files, "input_sha256": records_digest(files)})
    compiler = shlex.split(os.environ.get("CXX", "g++"))
    mpi_compiler = shlex.split(os.environ.get("MPICXX", "mpicxx"))
    nvcc_command = nvcc or os.environ.get("NVCC") or shutil.which("nvcc")
    if not nvcc_command:
        conda_nvcc = Path.home() / "anaconda3/envs/tubularflow-cuda/bin/nvcc"
        if conda_nvcc.is_file():
            nvcc_command = str(conda_nvcc)
    probes = {
        "compiler": probe(compiler + ["--version"]),
        "mpi": probe(["mpiexec", "--version"]),
        "mpi_compile": probe(mpi_compiler + ["--showme:compile"]),
        "mpi_link": probe(mpi_compiler + ["--showme:link"]),
        "cpu": probe(["lscpu", "--json"]),
        "hdf5": probe(["pkg-config", "--modversion", "hdf5"]),
        "hdf5_configuration": probe(["h5cc", "-showconfig"]),
        "gpu": probe(["nvidia-smi", "--query-gpu=name,driver_version,memory.total,compute_cap",
                      "--format=csv,noheader,nounits"]),
        "cuda_compiler": probe(shlex.split(nvcc_command) + ["--version"])
        if nvcc_command else {"status": "not_located"},
    }
    return {
        "schema_version": 1, "kind": "hpc_inventory", "simulation_executed": False,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "repository": {"root": str(root), "head": probe(["git", "rev-parse", "HEAD"], root),
                       "status": probe(["git", "status", "--short"], root),
                       "files": source, "content_sha256": records_digest(source)},
        "catalog_sha256": digest(catalog_path), "cases": fixtures,
        "execution": {"requested_mpi_ranks": ranks, "requested_omp_threads": threads,
                      "command_argv": command or [], "affinity_cpus": sorted(os.sched_getaffinity(0))
                      if hasattr(os, "sched_getaffinity") else None,
                      "environment": {name: os.environ.get(name) for name in RESOURCE_ENV}},
        "petsc": petsc_configuration(petsc_prefix), "probes": probes,
        "notes": ["Requested resources are not proof of actual rank/thread execution.",
                  "Probe failure may indicate sandbox restrictions; retain stderr and retry with permission.",
                  "HDF5 probe describes installed tools, not the solver binary's effective linkage.",
                  "Untracked source files are included; generated results are not reference solutions."],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--catalog", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--petsc-prefix")
    parser.add_argument("--nvcc")
    parser.add_argument("--binary", type=Path, action="append", default=[],
                        help="Built project executable to hash and inspect with ldd")
    parser.add_argument("--artifact", type=Path, action="append", default=[],
                        help="Explicit generated input or reference output to hash")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    options = parser.parse_args()
    if options.ranks < 1 or options.threads < 1:
        parser.error("ranks and threads must be positive")
    catalog = options.catalog or options.repo / "benchmarks/hpc_baselines.json"
    command = options.command[1:] if options.command[:1] == ["--"] else options.command
    result = collect(options.repo, catalog, options.ranks, options.threads,
                     options.petsc_prefix, options.nvcc, command)
    result["binaries"] = [{"path": str(path.resolve()), "sha256": digest(path),
                           "linked_libraries": probe(["ldd", str(path.resolve())])}
                          for path in options.binary]
    result["artifacts"] = [{"path": str(path.resolve()), "sha256": digest(path),
                            "size_bytes": path.stat().st_size}
                           for path in options.artifact]
    options.output.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive creation preserves previous evidence instead of silently replacing it.
    with options.output.open("x") as output:
        json.dump(result, output, indent=2, allow_nan=False)
        output.write("\n")
    print(f"HPC inventory: {options.output}; {len(result['cases'])} cases; no simulation executed")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"hpc_inventory: {error}", file=sys.stderr)
        sys.exit(1)
