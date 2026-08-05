# Repository Guidelines

## Project Structure & Module Organization

- `meshgeneration/`: MATLAB SWC smoothing and validated hexahedral control-mesh generation.
- `preprocessing/spline/`: C++11 spline construction and Bezier extraction.
- `solvers/cpu/`: C++17 MPI/PETSc packer, mesh checks, Navier-Stokes, and transport.
- `solvers/cuda/`: FP64 single-GPU implementation using the CPU database format.
- `docs/`: pipeline and benchmark evidence. Keep large cases and outputs outside Git.

Preserve the interfaces `controlmesh.vtk`, `bzmeshinfo.txt`, `cmat.txt`,
`bzpt.txt`, METIS partitions, and `.ntiga` databases.

## Build, Test, and Development Commands

```bash
make spline       # build Bezier preprocessing
make cpu          # build database utilities
make cpu-petsc    # build MPI/PETSc solvers
make cuda         # build the CUDA backend
```

On PSC Bridges-2, login nodes are for editing, compilation, and submission only.
Run CPU smoke tests interactively:

```bash
interact -A mch260002p -p RM-shared -t 00:30:00
module load anaconda3
module load openmpi/4.0.5-gcc10.2.0
mpiexec -np 8 ./solvers/cpu/iga_mesh_check CASE.ntiga
```

Request a GPU node before running `iga_cuda`; use `sbatch` for production.

## Coding Style & Naming Conventions

Use tabs in C++, braces on separate lines, `PascalCase` for classes and
methods, and descriptive lowercase locals. Keep CUDA kernels in `.cuh` files
and host orchestration in `.cu`. Preserve established MATLAB function names.
Compile with warnings enabled and resolve new warnings.

## Testing Guidelines

Rebuild only affected stages, then run the smallest valid case on a compute
node. Verify exit status, rank/partition agreement, positive Jacobians,
convergence reasons, output sizes, and numerical norms. For performance changes,
report assembly and solve times separately, plus host RSS and CUDA peak memory.
Compare CPU and CUDA fields using relative L2 error.

## Commit & Pull Request Guidelines

Use short imperative commits such as `Improve owned-row assembly`. Do not
commit executables, objects, `.ntiga`, partitions, VTK results, case data, or
Slurm logs. Pull requests should identify the pipeline stage, describe
numerical/file-format effects, list validation commands and hardware, and
include ParaView images for geometry changes.
