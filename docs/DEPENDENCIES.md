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
| CPU-only solver | Preprocessing requirements, MPI, optimized PETSc with C++, HDF5 | Yes | No |
| Native 1D solver | C++17, OpenMP, MPI, PETSc, HDF5; MUMPS recommended for multi-rank nonlinear solves | Yes | No |
| CUDA-only solver | Preprocessing requirements, CUDA Toolkit, cuBLAS, HDF5; cuDSS for GPU flow/FSI | No | Yes at runtime |

The CPU and CUDA solvers consume the same packed `.ntiga` database. Preparing
that database requires Eigen and `mpmetis`, regardless of the selected solver.

For the optional liver DICOM SEG geometry workflow, the project-owned
`dicom_seg_to_surface.py`, `audit_dicom_seg_region_overlap.py`, and
`dicom_seg_to_multiregion_tet.py` additionally
need Python `pydicom`, NumPy, and SciPy; surface extraction also needs VTK.
`surface_to_fem_volume.py`, `dicom_seg_to_multiregion_tet.py`, and the
template-free `preprocessing/tet/skeleton_to_tet.py` path need the Gmsh Python
module. Skeleton-to-tet and smooth dual-tree geometry additionally use the VTK
Python module for non-shrinking vessel-surface smoothing, while the fTetWild
route invokes a separately installed [fTetWild](https://github.com/wildmeshing/fTetWild)
binary. These tools generate or audit geometry only—none supplies FEM elements,
weak forms, assembly, or physical coupling. The WSL patient-label checks used
`pydicom 3.0.1`, NumPy 1.21.5, SciPy 1.8.0, and VTK 9.3.20240617; pydicom was
installed only in a case-specific `/tmp` directory, not as a repository
dependency. See `docs/LIVER_GEOMETRY_CANDIDATES.md` for source and license
checks before processing medical images.

## Automatic dependency check

Run the checker from the repository root:

```bash
./scripts/check_dependencies.sh preprocessing
./scripts/check_dependencies.sh cpu
./scripts/check_dependencies.sh one-d
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
export HDF5_CFLAGS="-I/path/to/hdf5/include"  # optional override
export HDF5_LIBS="-L/path/to/hdf5/lib -lhdf5"
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
  build-essential cmake git \
  libeigen3-dev metis \
  openmpi-bin libopenmpi-dev \
  libblas-dev liblapack-dev libhdf5-dev pkg-config
```

On these systems Eigen is normally under `/usr/include/eigen3`, so the default
spline build finds it without setting `EIGEN_DIR`. The `metis` package must
provide `mpmetis`; verify it with `command -v mpmetis`.

### RHEL, Rocky Linux, or AlmaLinux

Enable the repositories used by your site for development packages, then run:

```bash
sudo dnf install \
  gcc-c++ make cmake git \
  eigen3-devel metis \
  openmpi openmpi-devel \
  blas-devel lapack-devel hdf5-devel pkgconf-pkg-config
```

Some RHEL-family installations expose OpenMPI through Environment Modules. If
`mpicxx` is not initially in `PATH`, inspect the available MPI module and load
it, for example:

```bash
module avail mpi openmpi
module load mpi/openmpi-x86_64
```

Do not mix MPI implementations between PETSc compilation and solver runtime.
The temporal VTKHDF writer runs only on rank zero and uses serial HDF5 calls.
Prefer a serial HDF5 development package. If only parallel HDF5 is available,
it must use the same MPI implementation as PETSc. `HDF5_CFLAGS` and
`HDF5_LIBS` can override automatic detection; this is useful on systems exposing
more than one HDF5 installation. The CPU build first uses a complete HDF5
installation under `CONDA_PREFIX` or `CONDA_ROOT`, with a runtime library path,
then falls back to system libraries (including `/usr/lib64` on Bridges-2) and
`pkg-config`. `HDF5_CFLAGS` applies to both tools and PETSc solver builds.

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
./scripts/generate_case.sh examples/vascular_flow/straight_tube --ranks 2
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

`make cpu-petsc` also builds the native `solvers/one_d/iga_1d` executable. To
build or test only the 1D subsystem, use:

```bash
./scripts/check_dependencies.sh one-d
make one-d-petsc
make one-d-test
```

The build reads PETSc's compile and link flags from its pkg-config file,
`$PETSC_DIR/$PETSC_ARCH/lib/pkgconfig/PETSc.pc`. For a system PETSc package, set
`PETSC_DIR` to the installed prefix (for example
`/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real`) and leave `PETSC_ARCH`
empty. Multi-rank nonlinear tests use distributed MUMPS LU by default;
include MUMPS in the PETSc build or provide alternate PETSc KSP/PC options.
The iterative Schur fieldsplit test also needs hypre; tests whose PETSc
packages are missing are skipped with the reason. A PETSc build that runs every
test:

```bash
./configure PETSC_ARCH=arch-linux-c-opt --with-debugging=0 \
  --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
  --download-mumps --download-scalapack --download-metis --download-parmetis \
  --download-hypre
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

Without root access, install the toolkit (not a driver) in a Conda environment:

```bash
conda create -n tubularflow-cuda -c nvidia cuda-toolkit=12.6
conda run -n tubularflow-cuda nvcc --version
make cuda CUDA_ARCHS=89 NVCC="$(conda run -n tubularflow-cuda which nvcc)"
```

The GPU flow and FSI solvers factorize their Newton systems with cuDSS
(NVIDIA's sparse direct solver, 0.8 or newer, CUDA 12 build). Install it into
its own prefix and pass that prefix as `CUDSS_DIR`:

```bash
conda create -n tubularflow-cudss -c conda-forge libcudss-dev
make cuda CUDA_ARCHS=89 NVCC="$(conda run -n tubularflow-cuda which nvcc)" \
  CUDSS_DIR="$(conda run -n tubularflow-cudss printenv CONDA_PREFIX)"
./solvers/cuda/iga_cuda device-info
```

Without cuDSS, `native_tet_flow_cuda` and `native_tet_fsi_cuda` are skipped with
that reason; the other CUDA solvers do not need it. The CUDA executables carry
an RPATH to the toolkit and cuDSS libraries, so they start without setting
`LD_LIBRARY_PATH`. The build pins the CUDA host compiler to the system C++
compiler, so the Conda environment does not need to be active while building;
see [Building and testing](BUILD.md#cuda).

This is also appropriate for WSL when `nvidia-smi` already works but `nvcc`
does not: the Windows NVIDIA driver is exposed to WSL, while `nvcc` is supplied
by a separate Linux toolkit. Do not install a second NVIDIA driver inside WSL.

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
compute resources. If your development session or IDE is already inside an
allocation, load the modules there for lightweight validation; use `sbatch` for large, long,
GPU, benchmark, or parallel job sets. Do not request a nested allocation.

### Base environment

```bash
export PROJECT_ACCOUNT=YOUR_PSC_PROJECT
export PROJECT_ROOT=/ocean/projects/${PROJECT_ACCOUNT}/${USER}

module load anaconda3
module load openmpi/4.0.5-gcc10.2.0
module load cmake        # any CMake 3.20 or newer; check with cmake --version

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

./scripts/generate_case.sh examples/vascular_flow/straight_tube \
  --output "$PROJECT_ROOT/cases/tubular-straight" --ranks "$RANKS"
mpiexec -np "$RANKS" ./solvers/cpu/iga_mesh_check \
  "$PROJECT_ROOT/cases/tubular-straight/database/straight_tube-$RANKS.ntiga"
mpiexec -np "$RANKS" ./solvers/cpu/iga_navier_stokes \
  "$PROJECT_ROOT/cases/tubular-straight/database/straight_tube-$RANKS.ntiga" \
  "$PROJECT_ROOT/cases/tubular-straight/preprocessing" \
  --output "$PROJECT_ROOT/artifacts/runs/tubular-straight/velocity.txt"
```

The generator preserves a nonempty target directory. Choose another path or
rerun with `--clean` when replacing an earlier generated case is intentional.

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
  "$PROJECT_ROOT/cases/tubular-straight/database/straight_tube-2.ntiga"
```

The CUDA reader ignores CPU ownership records, so a database packed for two CPU
ranks can also be used by the single-GPU backend.

Additional Bridges-2 scheduler and benchmarking notes are in
[`docs/hpc/BRIDGES2.md`](hpc/BRIDGES2.md).

## Final verification

After installation, the expected checks are:

```bash
make mesh-test
make cpu-test
./scripts/generate_case.sh examples/vascular_flow/straight_tube --ranks 2

# With MPI/PETSc configured:
./scripts/check_dependencies.sh cpu

# With the CUDA module/toolkit configured:
./scripts/check_dependencies.sh cuda
```

Successful dependency checks confirm that tools and headers are visible. The
public straight-vessel example additionally confirms mesh generation, spline
extraction, METIS partitioning, database packing, and boundary-condition
resolution.
See the [Bridges-2 guide](hpc/BRIDGES2.md) for allocation-aware validation and the
[CPU](../solvers/cpu/VALIDATION.md) and
[CUDA](../solvers/cuda/VALIDATION.md) reports for larger numerical gates.
