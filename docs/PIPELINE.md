# Pipeline and Case Preparation

## 1. Generate the control mesh

A case directory starts with `skeleton_initial.swc` and
`mesh_parameter.txt`. The parameter file contains, in order:

1. noise-smoothing iterations;
2. bifurcation-node smoothing ratio;
3. noise-smoothing ratio;
4. target Bezier segment length;
5. bifurcation refinement ratio.

Install the MATLAB TREES Toolbox separately. From MATLAB, add both TREES and
this repository recursively to the path, set `io_path` near the top of
`TreeSmooth.m` and `Hexmesh_main.m`, then run:

```matlab
cd('/path/to/TubularFlowIGA/meshgeneration')
addpath(genpath('/path/to/TREES'))
TreeSmooth
Hexmesh_main
```

The scripts write `skeleton_smooth.swc`, `controlmesh.vtk`, and an initial
velocity field. The mesh generator rejects invalid lengths/diameters, limits
unsafe point motion, checks element corners and 4x4x4 quadrature samples, and
stops before writing a mesh below its scaled-Jacobian floor. Inspect the result
in ParaView as a second gate; a positive determinant is mandatory and a
minimum scaled Jacobian above 0.1 is a practical production target.

## 2. Extract the IGA representation

Build with `make spline`. The argument is a directory prefix, so include its
trailing slash:

```bash
./preprocessing/spline/spline "$CASE_DIR/"
```

It reads `controlmesh.vtk` and writes:

- `bzmeshinfo.txt`: Bezier element connectivity;
- `cmat.txt`: sparse extraction coefficients in the legacy interchange format;
- `bzpt.txt`: 64 Bezier points per element;
- `bzmesh.vtk`: visualization output.

## 3. Partition and pack

CPU rank count, METIS partition suffix, and `mpiexec -np` must agree.

```bash
RANKS=16
mpmetis "$CASE_DIR/bzmeshinfo.txt" "$RANKS"
./solvers/cpu/iga_pack "$CASE_DIR" "$RANKS" "$DATABASE"
./solvers/cpu/iga_inspect "$DATABASE"
```

The packer validates the legacy text once and creates a binary database with
direct element offsets, sparse extraction rows, adjacency, and per-rank
touching-element indices. Repack when changing the CPU rank count. CUDA ignores
ownership records and may reuse any valid packed database.

## 4. Validate and solve

Run expensive work only under Slurm. For CPU, request `RM-shared` with
`interact` for smoke tests or use `sbatch` for production. Load
`anaconda3` and the matching OpenMPI module, then run `iga_mesh_check`
before both physics stages.

For a V100 interactive test:

```bash
interact -A mch260002p -p GPU-shared --gres=gpu:v100-32:1 -t 00:30:00
module load cuda/12.4.0
./solvers/cuda/iga_cuda device-info
./solvers/cuda/iga_cuda mesh-check "$DATABASE"
./solvers/cuda/iga_cuda navier-stokes "$DATABASE" "$CASE_DIR" 8 velocity.txt
./solvers/cuda/iga_cuda transport \
  "$DATABASE" "$CASE_DIR" 300 concentration.txt velocity.txt
```

Record the case snapshot, compiler, PETSc/CUDA versions, rank count or GPU,
Slurm job ID, solver tolerances, numerical norms, stage timers, host RSS, and
CUDA peak allocation. Compare output fields, not runtime alone.
