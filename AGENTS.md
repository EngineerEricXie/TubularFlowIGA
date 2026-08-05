# Repository Guidelines

## Project Structure & Module Organization

- `preprocessing/mesh/`: dependency-free C++ SWC smoothing and control meshes.
- `meshgeneration/`: legacy MATLAB reference and mesh template assets.
- `preprocessing/spline/`: C++11 spline construction and Bezier extraction.
- `solvers/cpu/`: C++17 MPI/PETSc packer, checks, Navier-Stokes, and transport.
- `solvers/cuda/`: FP64 single-GPU implementation using the CPU database format.
- `docs/`: pipeline, platform instructions, and benchmark evidence.

Preserve `controlmesh.vtk`, `bzmeshinfo.txt`, `cmat.txt`, `bzpt.txt`,
METIS partition, velocity-field, and `.ntiga` interfaces.

## Build, Test, and Development Commands

```bash
make mesh
make mesh-test
make spline EIGEN_DIR=/path/to/eigen3
make cpu
make cpu-petsc PETSC_DIR=/path/to/petsc PETSC_ARCH=your-arch
make cuda CUDA_ARCHS="70 80"
```

Run simulations on an allocated compute resource, not a shared cluster login
node. Follow the local scheduler policy. Bridges-2 users should use the
interactive CPU/GPU and `sbatch` examples in
[docs/BRIDGES2.md](docs/BRIDGES2.md).

## Coding Style & Naming Conventions

Use tabs in C++, function braces on separate lines, `PascalCase` for classes
and methods, and descriptive lowercase locals. Keep CUDA kernels in `.cuh`
files and host orchestration in `.cu`. Preserve established MATLAB function
names. Compile with warnings enabled and resolve new warnings.

## Testing Guidelines

Run `make mesh-test` after geometry changes. Rebuild affected stages and run the
smallest valid case. Verify exit status,
rank/partition agreement, positive Jacobians, convergence reasons, output
sizes, and numerical norms. For performance work, report assembly and solve
times separately, host peak RSS, and CUDA peak allocation. Compare CPU and CUDA
fields with relative L2 error. Run large validation cases through a scheduler.
For spline changes, compare all four generated files against a known case and
validate the packed `.ntiga` database.

## Commit & Pull Request Guidelines

Use short imperative commits such as `Improve owned-row assembly`. Do not
commit executables, objects, `.ntiga`, partitions, VTK results, case data, or
scheduler logs. Pull requests should identify the pipeline stage, describe
numerical and file-format effects, list validation commands and hardware, link
issues, and include geometry images when applicable.
