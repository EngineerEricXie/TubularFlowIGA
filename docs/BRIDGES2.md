# PSC Bridges-2 Guide

This page contains optional site-specific instructions for PSC Bridges-2. The
solver and file formats do not depend on PSC.

## Environment

Login nodes are for editing, lightweight compilation, and job submission. Run
MPI simulations and all GPU commands on allocated compute nodes. Replace the
allocation below if you use another PSC project.

```bash
export PROJECT_ACCOUNT=mch260002p
module load anaconda3
```

Use an optimized PETSc installation with C++ and the same MPI toolchain used at
runtime:

```bash
module load openmpi/4.0.5-gcc10.2.0
make cpu-petsc \
  PETSC_DIR=/path/to/your/petsc \
  PETSC_ARCH=arch-linux-c-opt
```

The supplied `solvers/cpu/slurm/build_petsc_opt.sbatch` builds an optimized
PETSc source tree. Set `PETSC_DIR` and override its account with `sbatch -A`
when necessary.

## Interactive CPU smoke test

```bash
interact -A "$PROJECT_ACCOUNT" -p RM-shared -t 00:30:00
module load anaconda3
module load openmpi/4.0.5-gcc10.2.0

mpiexec -np 8 ./solvers/cpu/iga_mesh_check CASE.ntiga
mpiexec -np 8 ./solvers/cpu/iga_navier_stokes \
  CASE.ntiga CASE_DIR 8 velocity.txt
```

Use `mpiexec` for these repository MPI binaries. Plain `srun` may fail in
`MPI_Init_thread` with the validated OpenMPI configuration.

## Interactive GPU smoke test

```bash
interact -A "$PROJECT_ACCOUNT" -p GPU-shared \
  --gres=gpu:v100-32:1 -t 00:30:00
module load anaconda3
module load cuda/12.4.0

./solvers/cuda/iga_cuda device-info
./solvers/cuda/iga_cuda mesh-check CASE.ntiga
```

V100, H100, and L40S builds are supported by the default fat binary. For
production FP64 runs, select hardware based on availability, memory, and
double-precision throughput rather than consumer-GPU peak FLOPS.

## Batch validation

The scripts under `solvers/cpu/slurm/` and `solvers/cuda/slurm/` are
Bridges-2 examples. They contain the project allocation used for the published
measurements; override it with `sbatch -A "$PROJECT_ACCOUNT"`. CUDA validation
scripts require an external case root:

```bash
export IGA_CASE_ROOT=/path/to/case-data
sbatch -A "$PROJECT_ACCOUNT" \
  --export=ALL,IGA_CASE_ROOT="$IGA_CASE_ROOT" \
  solvers/cuda/slurm/validate_v100.sbatch
```

Keep case data, `.ntiga` databases, VTK output, and Slurm logs outside Git.
Record module versions, node or GPU type, task count, job ID, stage timings,
peak RSS, and CUDA peak allocation with every benchmark.
