# Quick start

This walkthrough prepares one vascular case and one neuron-transport case from
a fresh clone, then shows the matching CPU and CUDA commands. Run commands from
the repository root.

On a shared cluster, build and simulate inside an allocated compute resource.
Use the local scheduler for large cases and GPU work; PSC-specific examples are
in [BRIDGES2.md](BRIDGES2.md).

## 1. Clone and check preprocessing dependencies

```bash
git clone https://github.com/EngineerEricXie/TubularFlowIGA.git
cd TubularFlowIGA
./scripts/check_dependencies.sh preprocessing
```

Preprocessing needs GNU Make, a C++ compiler, Eigen 3, OpenMP, and METIS with
`mpmetis`. Installation commands for Ubuntu/Debian, RHEL-family systems, WSL,
PETSc, and CUDA are in [DEPENDENCIES.md](DEPENDENCIES.md).

## 2. Prepare a vascular example

```bash
VASCULAR_SOURCE=examples/vascular_flow/straight_tube
./scripts/generate_case.sh "$VASCULAR_SOURCE" --ranks 2
VASCULAR_WORK="$VASCULAR_SOURCE/generated"
VASCULAR_CASE="$VASCULAR_WORK/preprocessing"
VASCULAR_DB="$VASCULAR_WORK/database/straight_tube-2.ntiga"
VASCULAR_RESULTS="$VASCULAR_WORK/results"
```

This step does not require PETSc or a GPU. It builds the dependency-free tools,
generates the control mesh, extracts the spline, creates a two-way METIS
partition, packs the `.ntiga` database, and validates Jacobians, configuration,
and boundary labels.

Use this small straight tube for installation checks. The source-only
[NMO_06840 neuron transport](../examples/neuron_transport/nmo_06840_bifurcation/README.md)
is a separate 41,097-node regression whose preparation and solve can take many
minutes and create more than 700 MiB of cache/database files.

Successful output includes `bad_elements=0` with a positive control-mesh
`min_detJ`, a successfully validated packed element count, `schema_version=4`,
and resolved wall/inlet/outlet conditions.

## 3. Prepare a neuron example

```bash
NEURON_SOURCE=examples/neuron_transport/straight_neurite
./scripts/generate_case.sh "$NEURON_SOURCE" --ranks 2
NEURON_WORK="$NEURON_SOURCE/generated"
NEURON_CASE="$NEURON_WORK/preprocessing"
NEURON_DB="$NEURON_WORK/database/straight_neurite-2.ntiga"
NEURON_RESULTS="$NEURON_WORK/results"
```

The vascular configuration contains only `blood_flow`. The neuron
configuration contains only the two-field `neuron_transport` system. Keeping
the applications separate prevents accidentally running an unrelated equation
system.

## 4A. Run with MPI/PETSc on CPU

Set the PETSc installation used with the active MPI implementation, then build:

```bash
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-linux-c-opt

./scripts/check_dependencies.sh cpu
make cpu-petsc PETSC_DIR="$PETSC_DIR" PETSC_ARCH="$PETSC_ARCH"
```

For an installed PETSc prefix with configuration directly under
`$PETSC_DIR/lib/petsc/conf`, leave `PETSC_ARCH` unset and omit it from the make
command.

Run vascular flow:

```bash
mpiexec -np 2 ./solvers/cpu/iga_mesh_check "$VASCULAR_DB"
mpiexec -np 2 ./solvers/cpu/iga_navier_stokes \
  "$VASCULAR_DB" "$VASCULAR_CASE" \
  --output "$VASCULAR_RESULTS/velocity-cpu.txt"

./solvers/cpu/iga_flow_validate \
  "$VASCULAR_DB" "$VASCULAR_RESULTS/velocity-cpu.txt"
```

Run neuron transport:

```bash
mpiexec -np 2 ./solvers/cpu/iga_solve \
  "$NEURON_DB" "$NEURON_CASE" \
  --system neuron_transport \
  --output "$NEURON_RESULTS/neuron-cpu.txt"
```

The MPI process count must equal the partition count used to pack the database.
The vascular example should converge and report bounded mass imbalance. The
neuron output should contain one row per node; its `.fields` file should list
`N0` followed by `Nplus`.

## 4B. Run on one CUDA GPU

CUDA compilation needs `nvcc`; runtime additionally needs a compatible NVIDIA
driver and GPU:

```bash
./scripts/check_dependencies.sh cuda
make cuda CUDA_ARCHS=YOUR_GPU_ARCH
./solvers/cuda/iga_cuda device-info
```

Run vascular flow and neuron transport using the same packed databases:

```bash
./solvers/cuda/iga_cuda mesh-check "$VASCULAR_DB"
./solvers/cuda/iga_cuda navier-stokes \
  "$VASCULAR_DB" "$VASCULAR_CASE" \
  --output "$VASCULAR_RESULTS/velocity-cuda.txt"

./solvers/cuda/iga_cuda solve \
  "$NEURON_DB" "$NEURON_CASE" \
  --system neuron_transport \
  --output "$NEURON_RESULTS/neuron-cuda.txt"
```

CUDA is single-GPU and ignores CPU ownership records in the packed database.
Its flow solver writes a VTK file next to the requested text output.

## 5. Understand the outputs

| File | Meaning |
|---|---|
| `preprocessing/skeleton_normalized.swc` | Validated, explicitly rooted skeleton |
| `preprocessing/skeleton.vtp` | ParaView centerline with radius and topology arrays |
| `preprocessing/mesh_diagnostics.json` | Machine-readable geometry risks, limits, and failures |
| `preprocessing/skeleton_diagnostics.vtp` | ParaView centerline colored by geometry risk |
| `preprocessing/controlmesh.vtk` | Labeled hexahedral control mesh |
| `preprocessing/mesh_quality.json` | Final Jacobian and surface-intersection results |
| `visualization/bzmesh.vtkhdf` | Cubic Bezier geometry preview for ParaView |
| `visualization/bzmesh.bezier_geometry.json` | Validation report for the preview geometry |
| `preprocessing/geometry_transform.json` | Source-to-normalized coordinate transform |
| `database/*.ntiga` | Packed binary database consumed by CPU and CUDA |
| `results/velocity-cpu.txt` | Three velocity coefficients per node |
| `results/velocity-cpu.txt.pressure` | One pressure coefficient per node |
| `results/neuron-cpu.txt` | `node_id N0 Nplus` |
| `results/neuron-cpu.txt.fields` | Ordered transport field names |
| `results/velocity-cuda.txt.vtk` | CUDA flow result for visualization |

Generated work directories are reproducible and intentionally outside Git.
Keep them while inspecting results and delete only the exact work directories
you created when they are no longer needed.

## Common first-run failures

| Symptom | Likely cause | Check |
|---|---|---|
| `mpmetis: command not found` | METIS executable is missing | `command -v mpmetis` |
| Eigen headers not found | `EIGEN_DIR` does not contain `Eigen/` | `./scripts/check_dependencies.sh preprocessing` |
| PETSc compile or runtime link failure | `PETSC_DIR`, `PETSC_ARCH`, or MPI implementation mismatch | `./scripts/check_dependencies.sh cpu` |
| Packed rank mismatch | `mpiexec -np` differs from the database partition count | Reprepare with matching `RANKS` |
| `nvidia-smi` works but `nvcc` is missing | Driver is visible but CUDA Toolkit is absent | Follow the WSL/Conda section in `DEPENDENCIES.md` |
| CUDA runtime cannot see a device | Driver/GPU is unavailable in the current node or container | Run `iga_cuda device-info` on a GPU resource |
| Work directory is rejected | It already contains files | Create a fresh empty directory |

## Next steps

- Browse all cases in the [examples catalog](../examples/README.md).
- Learn the file stages in [PIPELINE.md](PIPELINE.md).
- Modify equations using [PDE_CONFIGURATION.md](PDE_CONFIGURATION.md).
- Modify inlet, wall, and outlet conditions using
  [BOUNDARY_CONDITIONS.md](BOUNDARY_CONDITIONS.md).
- Review current scope limits in the root [README](../README.md) and numerical
  evidence in the CPU and CUDA validation documents.
