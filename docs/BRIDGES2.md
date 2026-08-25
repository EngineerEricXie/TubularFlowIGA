# PSC Bridges-2 Guide

This page contains optional site-specific instructions for PSC Bridges-2. The
solver and file formats do not depend on PSC.

Complete module, optimized PETSc build, dependency-check, and public example
instructions are in the [dependency installation guide](DEPENDENCIES.md#psc-bridges-2-installation).

## Environment

Treat login nodes as coordination hosts only: edit, inspect, and submit jobs
there, but do not compile large targets or run simulations. Run builds, tests,
MPI programs, and GPU programs inside an allocated compute resource.

Codex or VS Code may already have been launched inside an allocation. Check
before requesting another one:

```bash
hostname
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  squeue -j "$SLURM_JOB_ID" -o '%.18i %.9P %.24j %.2t %.10M %.6D %R'
fi
```

Inside an existing CPU allocation, load the toolchain and run lightweight
serial builds, unit tests, dependency checks, or a small public example locally.
Do not start a nested `interact` session. Submit large cases, long benchmarks,
or several independent validations with `sbatch`. GPU execution still requires
an allocation containing a GPU.

The shell that launched Codex determines its inherited modules. Modules loaded
later in a different terminal are not automatically visible to that process.
Replace the allocation below if you use another PSC project.

```bash
export PROJECT_ACCOUNT=YOUR_PSC_PROJECT
export PROJECT_ROOT=/ocean/projects/${PROJECT_ACCOUNT}/${USER}
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

If already inside a CPU allocation, begin at the `module load` lines. Adjust
the PETSc path if the installation is elsewhere:

```bash
module load anaconda3
module load openmpi/4.0.5-gcc10.2.0
export PETSC_DIR="$PROJECT_ROOT/petsc"
export PETSC_ARCH=arch-linux-c-opt
./scripts/check_dependencies.sh cpu
make cpu-petsc PETSC_DIR="$PETSC_DIR" PETSC_ARCH="$PETSC_ARCH"
make cpu-test
```

To request a new allocation from a login node:

```bash
interact -A "$PROJECT_ACCOUNT" -p RM-shared -N 1 -n 8 -t 00:30:00
module load anaconda3
module load openmpi/4.0.5-gcc10.2.0

mpiexec -np 8 ./solvers/cpu/iga_mesh_check CASE.ntiga
mpiexec -np 8 ./solvers/cpu/iga_navier_stokes \
  CASE.ntiga CASE_DIR --max-newton 8 --output velocity.txt
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

Use batch jobs for large cases, long simulations, published measurements, or
multiple independent runs. A CPU allocation does not authorize GPU execution;
submit a GPU job when no GPU was requested for the current job.

The scripts under `solvers/cpu/slurm/` and `solvers/cuda/slurm/` are
Bridges-2 examples. They intentionally omit a project allocation; provide one
with `sbatch -A "$PROJECT_ACCOUNT"`. CUDA validation scripts require an
external case root:

```bash
export IGA_CASE_ROOT=/path/to/case-data
sbatch -A "$PROJECT_ACCOUNT" \
  --export=ALL,IGA_CASE_ROOT="$IGA_CASE_ROOT" \
  solvers/cuda/slurm/validate_v100.sbatch
```

The transient validation job additionally requires one configured case
directory and its packed database. Build both `make cpu` (for
`iga_flow_validate`) and `make cuda` before submission:

```bash
export IGA_TRANSIENT_CASE_DIR=/path/to/transient-case
export IGA_TRANSIENT_DATABASE=/path/to/transient-case/case-1.ntiga
sbatch -A "$PROJECT_ACCOUNT" \
  --export=ALL,IGA_TRANSIENT_CASE_DIR="$IGA_TRANSIENT_CASE_DIR",IGA_TRANSIENT_DATABASE="$IGA_TRANSIENT_DATABASE" \
  solvers/cuda/slurm/validate_transient_v100.sbatch
```

The case must contain backward-Euler Navier–Stokes and configured transport
systems with at least two time steps. The job compares uninterrupted and
checkpoint/restarted fields, reports boundary mass balance, and can compare
against CPU references through `IGA_CPU_FLOW_REFERENCE` and
`IGA_CPU_TRANSPORT_REFERENCE`.
Final-flow gates default to `1e-2` relative mass imbalance and `1e-6` relative
divergence-theorem error; override them with
`IGA_MASS_BALANCE_TOLERANCE` and `IGA_DIVERGENCE_THEOREM_TOLERANCE`. To gate a
multi-cycle run, also export `IGA_CARDIAC_PERIOD`, `IGA_FLOW_DT`, and
`IGA_FLOW_STEPS`. The script then checks every time-indexed field and defaults
to `1e-3` maximum cycle-to-cycle relative L2, configurable through
`IGA_CYCLE_REL_L2_TOLERANCE`.
Set `IGA_WOMERSLEY_CONFIG` to an analytical configuration whose selected sample
times occur in the numerical manifest. The job evaluates both fields at volume
quadrature points and defaults to a `5e-2` maximum physical relative L2 gate.
Alternatively, `IGA_WOMERSLEY_CASE_DIR` and `IGA_WOMERSLEY_MANIFEST` may identify
an independently projected coefficient reference; do not use raw point samples
for that comparison. Override the tolerance with
`IGA_WOMERSLEY_REL_L2_TOLERANCE`.

Keep case data, `.ntiga` databases, VTK output, and Slurm logs outside Git.
Record module versions, node or GPU type, task count, job ID, stage timings,
peak RSS, and CUDA peak allocation with every benchmark.
