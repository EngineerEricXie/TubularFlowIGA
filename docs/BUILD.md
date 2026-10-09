# Building and testing

CMake (3.20 or newer) builds every executable and CTest runs every test. Each
directory's `CMakeLists.txt` lists its targets and tests; shared helpers live in
[`cmake/`](../cmake). Executables are still written beside their sources (for
example `solvers/cpu/iga_solve`), so every script and command that names those
paths keeps working. Objects and generated files stay in the build directory.

## Using make

The historical `make` commands are unchanged. The root `Makefile` and each
directory `Makefile` forward to CMake and CTest through
[`cmake/forward.py`](../cmake/forward.py), which configures `build/` on first use:

```bash
make mesh mesh-test
make cpu cpu-test
make cpu-petsc PETSC_DIR=/path/to/petsc PETSC_ARCH=arch-linux-c-opt
make -C solvers/cpu native_tet_hydraulic_graph PETSC_DIR=... PETSC_ARCH=...
make -C solvers/cpu native-tet-ale-petsc-runtime-test PETSC_DIR=... PETSC_ARCH=...
make cuda CUDA_ARCHS=89 NVCC=/path/to/nvcc
```

- An executable name (`native_tet_flow`) builds that executable.
- A test name (`native-tet-fem-test`, `test`, `petsc-test`) builds what the test
  needs, then runs the CTest tests labelled `make=<name>` in that directory.
- `all`, `petsc`, and `zero-d` build the same executable sets as before.
- `clean` removes that directory's executables. Delete `build/` for a full reset.

`PETSC_DIR`/`PETSC_ARCH`, `CUDA_ARCHS`/`NVCC`, `EIGEN_DIR`, and
`HDF5_CFLAGS`/`HDF5_LIBS` are passed to the CMake configuration. Once set, the
build directory remembers them; pass them again only to change them. Builds use
`-j` from make when given (sharing make's job slots), otherwise up to eight
jobs. `make -n` prints the CMake and CTest commands without running them. Set `TFI_BUILD_DIR` to
keep separate build directories, for example one per PETSc architecture.

## Using CMake directly

```bash
cmake -S . -B build -DPETSC_DIR=/path/to/petsc -DPETSC_ARCH=arch-linux-c-opt
cmake --build build -j 8                    # every configured target
ctest --test-dir build -L unit -j 8         # PETSc-free unit tests
ctest --test-dir build -L petsc             # PETSc tests (including MPI ones)
ctest --test-dir build/solvers/cpu -L '^make=test$'
```

Targets whose dependencies are absent (PETSc, HDF5, CUDA, Eigen) are skipped
during configuration rather than failing it; `build/tfi-forward.json` lists them
with the missing dependency.

Every target is compiled with `-O3` and **without** `NDEBUG`, as the former
Makefiles did; the tests rely on `assert`.

### Test labels

| Label | Meaning |
|---|---|
| `unit` | No PETSc, MPI launcher, GPU, or ParaView needed |
| `petsc` | Needs the configured PETSc |
| `mpi` | Launches with `mpiexec` (the launcher recorded in PETSc's configuration) |
| `gpu` | Needs an NVIDIA GPU |
| `paraview` | Needs `pvpython` (override with `PVPYTHON`) |
| `manual` | Needs a user-supplied input, for example `T4_ALE_MESH` |
| `diagnostic` | Convergence studies and documented negative results, not acceptance gates |
| `slow` | Runs for more than about a minute; left out of CI |
| `make=<name>` | The historical make target that ran this test |
| `cpu`, `coupling`, `one_d`, `mesh`, `cuda` | Source directory |

Tests run from their source directory with the original shell commands, so
variables such as `MPI_TEST_RANKS` or `T2_BIFURCATION_MESH` can be exported
before `ctest` or passed on the `make` command line. Each test gets a private
`TMPDIR` ([`cmake/run_test.sh`](../cmake/run_test.sh)), so tests can run in
parallel (`ctest -j`, or `make -j` for a test goal). The directory is deleted
afterwards unless `TUBULARFLOWIGA_KEEP_TEST_OUTPUT` is set.

## CUDA

```bash
cmake -S . -B build -DTFI_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=89 \
  -DCMAKE_CUDA_COMPILER=/path/to/nvcc
cmake --build build --target iga_cuda native_tet_flow_cuda
```

The CUDA host compiler is pinned to the C++ compiler, so an `nvcc` from a Conda
environment can be used without activating that environment's compilers or
linker. `native_tet_darcy_cuda` and `native_tet_fsi_cuda` include PETSc headers
(without linking PETSc) and therefore need `PETSC_DIR`. To change the CUDA
compiler, remove the build directory or use another `TFI_BUILD_DIR`.

## Native graph checkpoint identity

`iga_multidomain_flow`, `iga_1d_3d_bifurcation`, and their checkpoint tests
embed a hash of the graph sources, `solvers/coupling/CMakeLists.txt`, and
`cmake/TubularFlowHelpers.cmake`. The build regenerates it whenever one of those
files changes; a checkpoint written by a different source state is rejected.

## Continuous integration

[`.github/workflows/ci.yml`](../.github/workflows/ci.yml) runs on every pull
request and on pushes to `main`:

| Job | What it checks |
|---|---|
| Build and unit tests | Configures without PETSc, builds every remaining target, and runs `ctest -L unit -LE slow` |
| Python tests and evidence cards | `scripts/evidence.py check` and the `unittest` modules in `scripts/tests` |

PETSc, MPI-launched, GPU, ParaView, and `slow` tests are not run by CI; run them
before merging changes that affect them:

```bash
make -C solvers/cpu petsc-test PETSC_DIR=... PETSC_ARCH=...
make -C solvers/coupling petsc-test hpc-regression-test
python3 scripts/evidence.py run <card>    # for the evidence a change touches
```
