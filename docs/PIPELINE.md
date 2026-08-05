# Pipeline and Case Preparation

## 1. Generate the control mesh

A case directory starts with `skeleton_initial.swc` and
`mesh_parameter.txt`. The parameter file contains, in order:

1. noise-smoothing iterations;
2. bifurcation-node smoothing ratio;
3. noise-smoothing ratio;
4. target Bezier segment length;
5. bifurcation refinement ratio.

Install the MATLAB TREES Toolbox separately. Add both TREES and this repository
recursively to the MATLAB path, set `io_path` near the top of
`TreeSmooth.m` and `Hexmesh_main.m`, then run:

```matlab
cd('/path/to/TubularFlowIGA/meshgeneration')
addpath(genpath('/path/to/TREES'))
TreeSmooth
Hexmesh_main
```

The scripts write `skeleton_smooth.swc`, `controlmesh.vtk`, and an initial
velocity field. The generator rejects invalid lengths and diameters, limits
unsafe point motion, checks element corners and 4x4x4 quadrature samples, and
stops before writing a mesh below its scaled-Jacobian floor. Inspect the result
in ParaView as a second gate. A positive determinant is mandatory; a minimum
scaled Jacobian above 0.1 is a practical production target.

## 2. Extract the IGA representation

Build with `make spline`. The argument is a directory prefix, so include its
trailing slash:

```bash
./preprocessing/spline/spline "$CASE_DIR/"
```

It reads `controlmesh.vtk` and writes:

- `bzmeshinfo.txt`: Bezier element connectivity;
- `cmat.txt`: extraction coefficients in the legacy interchange format;
- `bzpt.txt`: 64 Bezier points per element;
- `bzmesh.vtk`: visualization output.

## 3. Partition and pack

CPU rank count, METIS partition suffix, and `mpiexec -np` must agree.

```bash
RANKS=8
mpmetis "$CASE_DIR/bzmeshinfo.txt" "$RANKS"
./solvers/cpu/iga_pack "$CASE_DIR" "$RANKS" "$DATABASE"
./solvers/cpu/iga_inspect "$DATABASE"
```

The packer validates legacy text once and creates a binary database with direct
element offsets, sparse extraction rows, adjacency, and per-rank
touching-element indices. Repack when changing the CPU rank count. CUDA ignores
ownership records and may reuse any valid packed database.

## 4. Validate and solve

Use a workstation only for cases that fit its resources. On a shared cluster,
request a compute allocation through its scheduler. Validate with
`iga_mesh_check` or `iga_cuda mesh-check` before timing either physics
stage. CPU runs use `mpiexec`; CUDA execution is one process on one GPU.

Platform-specific examples:

- [PSC Bridges-2](BRIDGES2.md)

Record the case snapshot, compiler, MPI/PETSc or CUDA versions, rank count or
GPU, scheduler job ID when applicable, solver tolerances, numerical norms,
stage timers, host RSS, and CUDA peak allocation. Compare output fields, not
runtime alone.
