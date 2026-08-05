# TubularFlowIGA

TubularFlowIGA is a C++ isogeometric-analysis pipeline for stabilized steady
Navier-Stokes flow and transient two-field transport in tubular and branching
geometries. It combines MATLAB control-mesh generation, C++ Bezier extraction,
an MPI/PETSc CPU backend, and a single-GPU CUDA backend. The IGA basis,
quadrature, weak forms, assembly, nonlinear iteration, time integration, and
CUDA sparse kernels are project code; PETSc and cuBLAS provide low-level linear
algebra services.

This is research software specialized for tube-like networks. It is not a
general-purpose CFD package.

## Pipeline

```text
SWC skeleton
  -> MATLAB smoothing and hexahedral control mesh
  -> controlmesh.vtk
  -> C++ spline/Bezier extraction
  -> bzmeshinfo.txt + cmat.txt + bzpt.txt
  -> METIS partition + iga_pack
  -> partition-aware .ntiga database
  -> CPU (MPI/PETSc) or CUDA (single GPU)
  -> velocity, pressure, and transported fields
```

Repository layout:

- `meshgeneration/`: corrected MATLAB skeleton smoothing and mesh generation.
- `preprocessing/spline/`: spline construction and Bezier extraction.
- `solvers/cpu/`: packer, validators, Navier-Stokes, and transport with MPI/PETSc.
- `solvers/cuda/`: FP64 single-GPU solver consuming the same `.ntiga` database.
- `docs/`: pipeline details, validation conditions, and benchmark evidence.

Large cases and generated results are intentionally not versioned.

## Requirements

- MATLAB plus the external TREES Toolbox for SWC processing. TREES is not
  vendored; install it separately and add it to `MATLABPATH`.
- C++11, Eigen 3, and OpenMP for spline preprocessing.
- C++17, MPI, optimized PETSc, and METIS/`mpmetis` for CPU simulation.
- CUDA 12.x and cuBLAS for the optional GPU backend.

## Build on PSC Bridges-2

Compilation is allowed on a login node; simulations are not.

```bash
export TUBULARFLOWIGA_ROOT=/ocean/projects/mch260002p/thsieh1/TubularFlowIGA
cd "$TUBULARFLOWIGA_ROOT"

make spline
make cpu

module load openmpi/4.0.5-gcc10.2.0
make cpu-petsc

module load cuda/12.4.0
make cuda
```

Override `PETSC_DIR`, `PETSC_ARCH`, `CXX`, or `CUDA_ARCHS` when the
site defaults do not apply.

## Prepare and run a case

Keep case data outside the repository and preserve a trailing slash when calling
the legacy spline preprocessor:

```bash
export CASE_DIR=/path/to/case
export DATABASE=/path/to/work/case-16.ntiga

./preprocessing/spline/spline "$CASE_DIR/"
mpmetis "$CASE_DIR/bzmeshinfo.txt" 16
./solvers/cpu/iga_pack "$CASE_DIR" 16 "$DATABASE"
./solvers/cpu/iga_inspect "$DATABASE"
```

Before timing a solver, obtain a compute node and reject invalid geometry:

```bash
interact -A mch260002p -p RM-shared -t 00:30:00
module load anaconda3
module load openmpi/4.0.5-gcc10.2.0

mpiexec -np 16 ./solvers/cpu/iga_mesh_check "$DATABASE"
mpiexec -np 16 ./solvers/cpu/iga_navier_stokes \
  "$DATABASE" "$CASE_DIR" 8 velocity.txt
mpiexec -np 16 ./solvers/cpu/iga_transport \
  "$DATABASE" "$CASE_DIR" 300 concentration.txt velocity.txt
```

Use a GPU allocation for CUDA, then run `iga_cuda mesh-check` followed by the
same physics stages. See [the pipeline guide](docs/PIPELINE.md) for MATLAB,
Slurm, and file-interface details.

## Measured performance

On the corrected 35,949-node `NMO_54499_new` case, one V100 completed coupled
Navier-Stokes plus 300-step transport numerical work in 279.40 s versus
563.68 s on 16 CPU ranks (2.02x). Combined host and device peak allocation was
about 46% lower for Navier-Stokes and 51% lower for transport. CPU and CUDA
velocity, pressure, and transport relative L2 differences were
`5.33e-6`, `6.07e-6`, and `7.56e-6`.

On the 4,221-node cylinder, the optimized CPU first nonlinear update was about
16.6x faster and used 77.5% less memory than the legacy solver. The legacy
measurement was not a converged solve, so it is historical evidence rather than
a complete solver-to-solver baseline.

See [BENCHMARKS.md](docs/BENCHMARKS.md) and the backend
[CPU validation](solvers/cpu/VALIDATION.md) and
[CUDA validation](solvers/cuda/VALIDATION.md) for timings, job IDs, accuracy
gates, and limitations.
