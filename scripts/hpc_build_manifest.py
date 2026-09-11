#!/usr/bin/env python3
"""Record compiler, MPI, PETSc, HDF5, and CUDA build compatibility data."""

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess


ROOT = Path(__file__).resolve().parents[1]


def probe(argv):
    try:
        result = subprocess.run(argv, cwd=ROOT, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, timeout=20)
        return {"argv": argv, "returncode": result.returncode,
                "output": result.stdout.rstrip()}
    except (OSError, subprocess.TimeoutExpired) as error:
        return {"argv": argv, "returncode": 127, "output": str(error)}


def executable(name):
    value = shutil.which(name)
    return str(Path(value).resolve()) if value else None


def macros(path):
    if not path.is_file():
        return {"path": str(path), "status": "missing", "values": {}}
    selected = {}
    names = {
        "PETSC_USE_64BIT_INDICES", "PETSC_USE_COMPLEX", "PETSC_USE_REAL_SINGLE",
        "PETSC_USE_REAL_DOUBLE", "PETSC_USE_DEBUG", "PETSC_HAVE_HYPRE",
        "PETSC_HAVE_MUMPS", "PETSC_HAVE_HDF5", "PETSC_HAVE_MPIUNI",
    }
    for line in path.read_text(errors="replace").splitlines():
        match = re.match(r"#define\s+(PETSC_[A-Z0-9_]+)(?:\s+(.*))?$", line)
        if match and match.group(1) in names:
            selected[match.group(1)] = (match.group(2) or "1").strip()
    return {"path": str(path), "status": "read", "values": selected}


def petsc_prefix():
    root = os.environ.get("PETSC_DIR", "")
    arch = os.environ.get("PETSC_ARCH", "")
    if not root:
        return None
    return Path(root) / arch if arch else Path(root)


def cuda_architectures():
    configured = os.environ.get("CUDA_ARCHS", "").split()
    if configured:
        return {"source": "CUDA_ARCHS", "values": configured}
    makefile = ROOT / "solvers/cuda/Makefile"
    match = re.search(r"^CUDA_ARCHS\s*\?=\s*(.*)$", makefile.read_text(), re.MULTILINE)
    return {"source": "solvers/cuda/Makefile default",
            "values": match.group(1).split() if match else []}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--component", choices=("preprocessing", "cpu", "one-d", "cuda", "all"),
                        default="all")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dependency-status", choices=("passed", "failed", "not_run"),
                        default="not_run")
    parser.add_argument("--binary", action="append", type=Path, default=[])
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output already exists")

    cxx = os.environ.get("CXX", "g++")
    mpicxx = os.environ.get("MPICXX", "mpicxx")
    nvcc = os.environ.get("NVCC", "nvcc")
    prefix = petsc_prefix()
    report = {
        "schema_version": 1,
        "kind": "hpc_build_manifest",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "component": args.component,
        "dependency_status": args.dependency_status,
        "source": {
            "commit": probe(["git", "rev-parse", "HEAD"]),
            "status_porcelain": probe(["git", "status", "--porcelain=v1", "--untracked-files=all"]),
        },
        "host": {"node": platform.node(), "platform": platform.platform(),
                 "machine": platform.machine()},
        "environment": {name: os.environ.get(name) for name in (
            "PETSC_DIR", "PETSC_ARCH", "HDF5_CFLAGS", "HDF5_LIBS", "CUDA_ARCHS",
            "SLURM_JOB_ID", "SLURM_JOB_NUM_NODES", "SLURM_NTASKS",
            "SLURM_CPUS_PER_TASK", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS", "BLIS_NUM_THREADS")},
        "tools": {
            "cxx": {"path": executable(cxx), "version": probe([cxx, "--version"])},
            "mpicxx": {"path": executable(mpicxx), "version": probe([mpicxx, "--version"]),
                       "compile": probe([mpicxx, "--showme:compile"]),
                       "link": probe([mpicxx, "--showme:link"])},
            "mpiexec": {"path": executable("mpiexec"), "version": probe(["mpiexec", "--version"])},
            "hdf5": {"pkg_version": probe(["pkg-config", "--modversion", "hdf5"]),
                     "pkg_cflags": probe(["pkg-config", "--cflags", "hdf5"]),
                     "pkg_libs": probe(["pkg-config", "--libs", "hdf5"]),
                     "configuration": probe(["h5cc", "-showconfig"])},
            "cuda": {"nvcc_path": executable(nvcc), "version": probe([nvcc, "--version"]),
                     "target_architectures": cuda_architectures()},
        },
        "petsc": {"prefix": str(prefix) if prefix else None,
                  "configuration": macros(prefix / "include" / "petscconf.h") if prefix else None,
                  "variables": str(prefix / "lib/petsc/conf/petscvariables") if prefix else None},
        "binaries": [],
    }
    for requested in args.binary:
        path = requested.resolve()
        item = {"path": str(path), "exists": path.is_file()}
        if path.is_file():
            item["file"] = probe(["file", str(path)])
            item["dynamic_libraries"] = probe(["ldd", str(path)])
        report["binaries"].append(item)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    main()
