# Quick start

On a shared cluster, run these commands inside a compute allocation. If Codex or
VS Code is already running under a `SLURM_JOB_ID`, reuse that allocation for
lightweight builds and smoke tests after loading the required modules. Use
`sbatch` for large cases, long runs, GPU work without a GPU allocation, or
multiple independent jobs. Bridges-2 commands are in
[`BRIDGES2.md`](BRIDGES2.md).

Choose the path matching the backend that will be executed:

| Mode | Quick verification | Additional requirement |
|---|---|---|
| Preprocessing-only | Build and prepare the public smoke database | Eigen and `mpmetis` |
| CPU-only | Run mesh check and MPI transport on the smoke database | MPI and optimized PETSc |
| CUDA-only | Run CUDA mesh check or solver on the smoke database | CUDA Toolkit and NVIDIA GPU |

## Preprocessing-only

The fastest way to verify a fresh clone is:

```bash
git clone YOUR_TUBULARFLOWIGA_URL
cd TubularFlowIGA
./scripts/check_dependencies.sh preprocessing
./scripts/prepare_smoke.sh
```

This creates a temporary straight-tube case, then runs mesh generation, spline
extraction, METIS, database packing, inspection, and boundary validation. It
does not require PETSc and prints the generated case and database paths.

## CPU-only

Install or build an optimized PETSc with the same MPI compiler used at runtime:

```bash
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-linux-c-opt   # omit for an installed prefix
./scripts/check_dependencies.sh cpu
make cpu-petsc PETSC_DIR="$PETSC_DIR" PETSC_ARCH="$PETSC_ARCH"
```

Prepare the smoke case with a partition count matching the MPI process count:

```bash
export RANKS=2
./scripts/prepare_smoke.sh /path/to/empty/smoke-work
mpiexec -np "$RANKS" ./solvers/cpu/iga_mesh_check \
  /path/to/empty/smoke-work/smoke-$RANKS.ntiga
mpiexec -np "$RANKS" ./solvers/cpu/iga_solve \
  /path/to/empty/smoke-work/smoke-$RANKS.ntiga \
  /path/to/empty/smoke-work tracer_transport /path/to/empty/smoke-work/tracer.txt
```

Run MPI commands on a compute node when using a cluster. Navier-Stokes can then
produce a velocity field for transport as shown in the main README.

## CUDA-only

CUDA does not use PETSc:

```bash
./scripts/check_dependencies.sh cuda
make cuda CUDA_ARCHS="70 80 89 90"
./solvers/cuda/iga_cuda mesh-check /path/to/smoke-2.ntiga
```

Run CUDA commands on a GPU node. See the dependency and boundary-condition
guides for installation variables and case configuration.
