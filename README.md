# TubularFlowIGA

TubularFlowIGA is a C++ isogeometric-analysis pipeline for stabilized steady
Navier-Stokes flow and transient two-field transport in tubular and branching
geometries. It combines dependency-free C++ control-mesh generation, C++
Bezier extraction, an MPI/PETSc CPU backend, and a single-GPU CUDA backend. The IGA basis,
quadrature, weak forms, assembly, nonlinear iteration, time integration, and
CUDA sparse kernels are project code; PETSc and cuBLAS provide linear-algebra
services.

This is research software specialized for tube-like networks, not a
general-purpose CFD package.

## Pipeline

```text
SWC skeleton
  -> C++ smoothing and hexahedral control mesh
  -> controlmesh.vtk
  -> C++ spline and Bezier extraction
  -> bzmeshinfo.txt + cmat.txt + bzpt.txt
  -> METIS partition + iga_pack
  -> partition-aware .ntiga database
  -> CPU (MPI/PETSc) or CUDA (single GPU)
  -> velocity, pressure, and transported fields
```

Repository layout:

- `preprocessing/mesh/`: standalone C++ skeleton smoothing and control-mesh
  generation.
- `meshgeneration/`: legacy MATLAB reference and template geometry assets.
- `preprocessing/spline/`: spline construction and Bezier extraction.
- `solvers/cpu/`: packer, validators, Navier-Stokes, and transport with MPI/PETSc.
- `solvers/cuda/`: FP64 single-GPU solver using the same `.ntiga` database.
- `docs/`: pipeline, platform notes, and benchmark evidence.

Large cases and generated results are intentionally not versioned.

## Requirements

- A C++17 compiler for control-mesh generation (no geometry library required).
- A C++11 compiler, Eigen 3, and OpenMP for spline preprocessing.
- A C++17 compiler, MPI, PETSc with C++ support, and METIS/`mpmetis` for CPU
  simulation.
- CUDA and cuBLAS for the optional GPU backend.

MATLAB and the external TREES Toolbox are optional and needed only to reproduce
the legacy reference workflow. TREES is not vendored.

## Build

From the repository root:

```bash
make mesh
make mesh-test
make spline EIGEN_DIR=/path/to/eigen3
make cpu

make cpu-petsc \
  PETSC_DIR=/path/to/petsc \
  PETSC_ARCH=your-petsc-arch

make cuda CUDA_ARCHS="70 80 89 90"
```

`PETSC_ARCH` may be omitted for an installed PETSc tree whose configuration
is directly under `PETSC_DIR/lib/petsc/conf`. Override `CXX`, `CPPFLAGS`,
`CXXFLAGS`, `NVCC`, or `CUDA_ARCHS` for the local toolchain.

## Prepare a case

Keep case data outside the repository. The spline preprocessor treats its
argument as a directory prefix, so include the trailing slash:

```bash
export CASE_DIR=/path/to/case
export DATABASE=/path/to/work/case-8.ntiga
export RANKS=8

./preprocessing/mesh/tubular_mesh pipeline \
  "$CASE_DIR" meshgeneration/template
./preprocessing/spline/spline "$CASE_DIR/"
mpmetis "$CASE_DIR/bzmeshinfo.txt" "$RANKS"
./solvers/cpu/iga_pack "$CASE_DIR" "$RANKS" "$DATABASE"
./solvers/cpu/iga_inspect "$DATABASE"
```

The METIS partition count, packed database rank count, and CPU MPI process count
must agree.

## Validate and run

Run the mesh check before either solver:

```bash
mpiexec -np "$RANKS" ./solvers/cpu/iga_mesh_check "$DATABASE"
mpiexec -np "$RANKS" ./solvers/cpu/iga_navier_stokes \
  "$DATABASE" "$CASE_DIR" 8 velocity.txt
mpiexec -np "$RANKS" ./solvers/cpu/iga_transport \
  "$DATABASE" "$CASE_DIR" 300 concentration.txt velocity.txt
```

For a CUDA-capable host:

```bash
./solvers/cuda/iga_cuda device-info
./solvers/cuda/iga_cuda mesh-check "$DATABASE"
./solvers/cuda/iga_cuda navier-stokes \
  "$DATABASE" "$CASE_DIR" 8 velocity.txt
./solvers/cuda/iga_cuda transport \
  "$DATABASE" "$CASE_DIR" 300 concentration.txt velocity.txt
```

Use the job scheduler required by the target cluster. Platform-specific setup
is documented separately:

- [PSC Bridges-2](docs/BRIDGES2.md)

The case directory must initially contain `skeleton_initial.swc` and
`mesh_parameter.txt`. See the [mesh generator guide](preprocessing/mesh/README.md),
[validation report](docs/MESH_CPP_VALIDATION.md), and
[pipeline guide](docs/PIPELINE.md) for assumptions and file-interface details.

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

See [BENCHMARKS.md](docs/BENCHMARKS.md), the
[CPU validation](solvers/cpu/VALIDATION.md), and the
[CUDA validation](solvers/cuda/VALIDATION.md) for hardware, timings, accuracy
gates, and limitations.
