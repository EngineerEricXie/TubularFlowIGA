# Dependencies and installation

TubularFlowIGA does not vendor Eigen, PETSc, MPI, METIS, or CUDA. Keeping these
platform-specific dependencies external makes the repository small and lets
each system use its native compiler, MPI, GPU driver, and optimized libraries.

A clone without PETSc can still build the mesh generator, spline preprocessor,
METIS database packer, inspector, and boundary-condition validator. PETSc is
required only by the MPI CPU solvers. CUDA does not use PETSc.

## Component matrix

| Mode | Required software | PETSc | GPU |
|---|---|---:|---:|
| Preprocessing-only | GNU Make, C++ compiler, Eigen 3, OpenMP, METIS/`mpmetis` | No | No |
| CPU-only solver | Preprocessing requirements, MPI, optimized PETSc with C++ | Yes | No |
| CUDA-only solver | Preprocessing requirements, CUDA Toolkit and cuBLAS | No | Yes at runtime |

The CPU and CUDA solvers consume the same packed `.ntiga` database. Preparing
that database requires Eigen and `mpmetis`, regardless of the selected solver.

## Automatic dependency check

Run the checker from the repository root:

```bash
./scripts/check_dependencies.sh preprocessing
./scripts/check_dependencies.sh cpu
./scripts/check_dependencies.sh cuda
./scripts/check_dependencies.sh all
```

The checker reports the C++ compiler, Eigen headers, `mpmetis`, MPI compiler,
PETSc configuration, and CUDA compiler needed by the selected mode. A missing
optional backend does not prevent checking another mode.

Recognized environment variables are:

```bash
export CXX=g++
export MPICXX=mpicxx
export EIGEN_DIR=/path/to/eigen3       # directory containing Eigen/
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-linux-c-opt     # omit for an installed PETSc prefix
export NVCC=nvcc
```

## General Linux installation

The following commands install the non-CUDA prerequisites using common package
names. Package names can differ on older distributions or installations with
restricted repositories. Run the dependency checker afterward instead of
assuming that a package installation supplied every command.

### Ubuntu or Debian

```bash
sudo apt update
sudo apt install \
  build-essential git \
  libeigen3-dev metis \
  openmpi-bin libopenmpi-dev \
  libblas-dev liblapack-dev
```

On these systems Eigen is normally under `/usr/include/eigen3`, so the default
spline Makefile works without setting `EIGEN_DIR`. The `metis` package must
provide `mpmetis`; verify it with `command -v mpmetis`.

### RHEL, Rocky Linux, or AlmaLinux

Enable the repositories used by your site for development packages, then run:

```bash
sudo dnf install \
  gcc-c++ make git \
  eigen3-devel metis \
  openmpi openmpi-devel \
  blas-devel lapack-devel
```

Some RHEL-family installations expose OpenMPI through Environment Modules. If
`mpicxx` is not initially in `PATH`, inspect the available MPI module and load
it, for example:

```bash
module avail mpi openmpi
module load mpi/openmpi-x86_64
```

Do not mix MPI implementations between PETSc compilation and solver runtime.

### If Eigen or `mpmetis` is unavailable

Prefer a distribution package. Otherwise install Eigen 3 headers and METIS 5
from their upstream projects, and point `EIGEN_DIR` at the directory containing
`Eigen/`. The required METIS executable is `mpmetis`, not only the METIS
library.

- Eigen: <https://gitlab.com/libeigen/eigen>
- METIS: <https://github.com/KarypisLab/METIS>
- OpenMPI: <https://docs.open-mpi.org/>

Verify the preprocessing environment:

```bash
./scripts/check_dependencies.sh preprocessing
make mesh mesh-test
make spline EIGEN_DIR="${EIGEN_DIR:-/usr/include/eigen3}"
make cpu cpu-test
./scripts/prepare_smoke.sh
```

`make cpu` in the last block builds only the PETSc-free packer, inspector, and
case validator. It does not build the MPI solvers.

## PETSc for the CPU solver

An optimized PETSc build is recommended. Use the same MPI compiler wrappers to
build PETSc and TubularFlowIGA, and use that MPI implementation at runtime.
The commands below follow the official PETSc release installation workflow.

```bash
git clone -b release https://gitlab.com/petsc/petsc.git petsc
cd petsc

export PETSC_DIR="$PWD"
export PETSC_ARCH=arch-linux-c-opt

./configure \
  PETSC_ARCH="$PETSC_ARCH" \
  --with-debugging=0 \
  --with-cc=mpicc \
  --with-cxx=mpicxx \
  --with-fc=0 \
  --with-x=0 \
  --with-shared-libraries=1 \
  COPTFLAGS=-O3 \
  CXXOPTFLAGS=-O3

make PETSC_ARCH="$PETSC_ARCH" -j"$(nproc)" all
make PETSC_ARCH="$PETSC_ARCH" check
```

If the system has no usable BLAS/LAPACK, add PETSc's
`--download-f2cblaslapack` configuration option. Consult the
[official PETSc installation tutorial](https://petsc.org/release/install/install_tutorial/)
for supported alternatives.

Return to TubularFlowIGA and build the CPU solver:

```bash
cd /path/to/TubularFlowIGA
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-linux-c-opt

./scripts/check_dependencies.sh cpu
make cpu-petsc PETSC_DIR="$PETSC_DIR" PETSC_ARCH="$PETSC_ARCH"
```

For an installed PETSc prefix whose configuration is directly under
`$PETSC_DIR/lib/petsc/conf`, leave `PETSC_ARCH` unset:

```bash
unset PETSC_ARCH
make cpu-petsc PETSC_DIR="$PETSC_DIR"
```

## CUDA for the CUDA solver

Install a CUDA Toolkit supported by the host compiler and NVIDIA driver. Enable
NVIDIA's repository for your Linux distribution first; then the toolkit package
is normally installed with one of:

```bash
# Ubuntu or Debian, after enabling NVIDIA's CUDA repository
sudo apt update
sudo apt install cuda-toolkit

# RHEL, Rocky Linux, or AlmaLinux, after enabling NVIDIA's CUDA repository
sudo dnf install cuda-toolkit
```

Driver setup is system-specific and separate from this repository. Follow the
[official NVIDIA CUDA Linux installation guide](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/)
rather than mixing runfile and package-manager installations.

Verify and build:

```bash
nvcc --version
nvidia-smi                         # requires a machine with an NVIDIA GPU
./scripts/check_dependencies.sh cuda
make cuda CUDA_ARCHS="70 80 89 90"
```

Compiling requires the toolkit; running `iga_cuda` additionally requires a
compatible NVIDIA GPU and driver. Set `CUDA_ARCHS` to the compute capabilities
that must be supported by the resulting binary.

## PSC Bridges-2 installation

Login nodes are coordination hosts for editing, inspection, and Slurm job
submission. Run builds, tests, MPI simulations, and GPU commands on allocated
compute resources. If Codex or VS Code is already inside an allocation, load
the modules there for lightweight validation; use `sbatch` for large, long,
GPU, benchmark, or parallel job sets. Do not request a nested allocation.

### Base environment

```bash
export PROJECT_ACCOUNT=mch260002p
export PROJECT_ROOT=/ocean/projects/${PROJECT_ACCOUNT}/${USER}

module load anaconda3
module load openmpi/4.0.5-gcc10.2.0

cd "$PROJECT_ROOT/TubularFlowIGA"
./scripts/check_dependencies.sh preprocessing
```

The validated Bridges-2 environment provides Eigen under
`/usr/include/eigen3` and `mpmetis` under `/usr/bin`.

### Build optimized PETSc on Bridges-2

Clone PETSc into project storage rather than adding it to the Git repository:

```bash
mkdir -p "$PROJECT_ROOT/software"
git clone -b release https://gitlab.com/petsc/petsc.git \
  "$PROJECT_ROOT/software/petsc"

export PETSC_DIR="$PROJECT_ROOT/software/petsc"
export PETSC_ARCH=arch-linux-c-opt

sbatch -A "$PROJECT_ACCOUNT" \
  --export=ALL,PETSC_DIR="$PETSC_DIR",PETSC_ARCH="$PETSC_ARCH" \
  solvers/cpu/slurm/build_petsc_opt.sbatch
```

Wait for the Slurm job to complete successfully, then enter or reuse a CPU
compute allocation to build and check the backend:

```bash
module load openmpi/4.0.5-gcc10.2.0
export PETSC_DIR="$PROJECT_ROOT/software/petsc"
export PETSC_ARCH=arch-linux-c-opt

./scripts/check_dependencies.sh cpu
make cpu-petsc PETSC_DIR="$PETSC_DIR" PETSC_ARCH="$PETSC_ARCH"
```

Only `arch-linux-c-opt` is needed for normal simulations. Do not commit PETSc
source, build directories, or libraries to TubularFlowIGA.

### CPU-node smoke test

Request enough tasks for the selected partition count:

```bash
interact -A "$PROJECT_ACCOUNT" -p RM-shared -N 1 -n 2 -t 00:30:00
module load anaconda3
module load openmpi/4.0.5-gcc10.2.0

cd "$PROJECT_ROOT/TubularFlowIGA"
export PETSC_DIR="$PROJECT_ROOT/software/petsc"
export PETSC_ARCH=arch-linux-c-opt
export RANKS=2

./scripts/prepare_smoke.sh "$PROJECT_ROOT/cases/tubular-smoke"
mpiexec -np "$RANKS" ./solvers/cpu/iga_mesh_check \
  "$PROJECT_ROOT/cases/tubular-smoke/smoke-$RANKS.ntiga"
mpiexec -np "$RANKS" ./solvers/cpu/iga_solve \
  "$PROJECT_ROOT/cases/tubular-smoke/smoke-$RANKS.ntiga" \
  "$PROJECT_ROOT/cases/tubular-smoke" tracer_transport \
  "$PROJECT_ROOT/cases/tubular-smoke/tracer.txt"
```

`prepare_smoke.sh` requires an empty target directory. Choose another path or
remove/archive an earlier generated smoke directory before rerunning it.

### GPU build and smoke test

CUDA compilation is lightweight but should still be performed in an allocated
CPU or GPU compute resource:

```bash
module load cuda/12.4.0
cd "$PROJECT_ROOT/TubularFlowIGA"
./scripts/check_dependencies.sh cuda
make cuda CUDA_ARCHS="70 80 89 90"
```

Run the executable only after requesting a GPU node:

```bash
interact -A "$PROJECT_ACCOUNT" -p GPU-shared \
  --gres=gpu:v100-32:1 -t 00:30:00
module load anaconda3
module load cuda/12.4.0

cd "$PROJECT_ROOT/TubularFlowIGA"
./solvers/cuda/iga_cuda device-info
./solvers/cuda/iga_cuda mesh-check \
  "$PROJECT_ROOT/cases/tubular-smoke/smoke-2.ntiga"
```

The CUDA reader ignores CPU ownership records, so a database packed for two CPU
ranks can also be used by the single-GPU backend.

Additional Bridges-2 scheduler and benchmarking notes are in
[`docs/BRIDGES2.md`](BRIDGES2.md).

## Final verification

After installation, the expected checks are:

```bash
make mesh-test
make cpu-test
./scripts/prepare_smoke.sh

# With MPI/PETSc configured:
./scripts/check_dependencies.sh cpu

# With the CUDA module/toolkit configured:
./scripts/check_dependencies.sh cuda
```

Successful dependency checks confirm that tools and headers are visible. The
public smoke case additionally confirms mesh generation, spline extraction,
METIS partitioning, database packing, and boundary-condition resolution.
See the [Bridges-2 guide](BRIDGES2.md) for allocation-aware validation and the
[development roadmap](ROADMAP.md) for larger numerical gates.
