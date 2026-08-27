<p align="center">
  <img src="docs/images/iga-wordmark-flow.png"
       alt="ParaView rendering of velocity magnitude through a connected pipe spelling IGA"
       width="100%">
</p>

# TubularFlowIGA

TubularFlowIGA is a native C++ simulation toolkit for tubular and branching
networks. It provides a direct SWC-to-1D vascular flow/transport path and a
three-dimensional isogeometric-analysis (IGA) pipeline that generates a
hexahedral control mesh, constructs the spline and Bezier representation,
packs a partition-aware database, and solves on MPI/PETSc CPUs or one CUDA GPU.

This is research software specialized for tube-like networks. It is not a
general-purpose CFD package.

## What can it simulate?

| Application | Available now | Important boundary |
|---|---|---|
| Vascular flow | Native 1D rigid Poiseuille and compliant A/Q networks; CPU 3D rigid-wall steady/transient Navier--Stokes; native 1D and CPU 3D `vca_closed_loop` vascular coupling | 3D VCA requires backward-Euler CPU flow; species-coupled runs support one in-memory transport system. CUDA VCA, 3D replay/open-loop, and 3D FSI are not included |
| Neuron transport | Configurable two-field `N0`/`Nplus` axonal transport on straight and branching neurites | This is material transport, not membrane voltage, action potentials, synapses, or network electrophysiology |
| Generic biological transport | Config-selected 1D and 3D multispecies transport with reaction, source, wall exchange, metabolism, oxygen capacity, and blood-gas derived fields | The physiology layer is a configurable reduced model; 3D vasodilation is disabled because rigid-wall flow has no FSI |

CPU and CUDA use the same `simulation_config.json` schema and packed `.ntiga`
database. CUDA configured transport supports one through eight scalar fields.

## Example results

The image above shows velocity magnitude from the connected
`vascular_flow/iga_wordmark` example.

| 3D pulse multispecies physiology | 1D pulse multispecies physiology |
|---|---|
| ![3D pulse flow with oxygen, glucose, and lactate transport](docs/images/multispecies-3d-pulse.gif) | ![1D pulse flow with oxygen, glucose, and lactate transport](docs/images/multispecies-1d-pulse.gif) |
| Transient Navier--Stokes velocity plus six transported species on a curved vessel. | Compliant pressure-network flow plus six conservative network species. |

| Steady-state 3D vascular flow | Neuron material transport |
|---|---|
| ![ParaView center slice of velocity magnitude in the steady-state Y-bifurcation vascular example](docs/images/vascular-y-bifurcation-velocity.png) | ![Animated Nplus transport through the branched-neurite example](docs/images/neuron-branched-transport.gif) |
| Steady-state rigid-wall flow through a Y-bifurcation. | Time-dependent two-field material transport through a branched neurite. |

Reproduction commands are in the [examples catalog](examples/README.md), and
numerical checks are recorded in the
[public-example validation report](examples/VALIDATION.md).

## Pipeline

```text
SWC or radius-annotated line-OBJ centerline
  -> C++ smoothing and hexahedral control mesh
  -> controlmesh.vtk
  -> C++ spline and Bezier extraction
  -> bzmeshinfo.txt + spline_cache.igacache
  -> METIS partition + iga_pack
  -> partition-aware .ntiga database
  -> CPU (MPI/PETSc) or CUDA (single GPU)
  -> velocity, pressure, and transported fields
```

## Choose a first example

Every 3D runnable example contains `skeleton_initial.swc` or
`skeleton_initial.obj` plus a schema-v4 `simulation_config.json`; its
`geometry` and `mesh` blocks configure preprocessing in the same validated
document as the physics. A native 1D example needs only the skeleton and
schema-v3 configuration. Generated meshes, databases, and results are written
to a separate work directory.

| Goal | Recommended case | Command |
|---|---|---|
| First vascular run | Straight rigid vessel | `./scripts/prepare_example.sh vascular_flow/straight_tube` |
| Curved vascular geometry | Planar bend | `./scripts/prepare_example.sh vascular_flow/bent_tube` |
| Branching vascular flow | Symmetric bifurcation | `./scripts/prepare_example.sh vascular_flow/y_bifurcation` |
| Large neuron regression | NMO_06840 transport; long-running | `RANKS=12 ./scripts/prepare_example.sh neuron_transport/nmo_06840_bifurcation` |
| README showcase | Connected IGA wordmark | `./scripts/prepare_example.sh vascular_flow/iga_wordmark` |
| First neuron run | Straight neurite | `./scripts/prepare_example.sh neuron_transport/straight_neurite` |
| Branching neuron transport | Branched neurite | `./scripts/prepare_example.sh neuron_transport/branched_neurite` |
| First native 1D run | Straight Poiseuille vessel | `./solvers/one_d/iga_1d examples/one_d/rigid_straight --check` |
| Compliant 1D flow | Pulsatile Y-bifurcation | `./solvers/one_d/iga_1d examples/one_d/compliant_bifurcation` |
| 1D multispecies physiology | Six-species pulse network | `./solvers/one_d/iga_1d examples/one_d/multispecies_physiology` |
| 3D multispecies pulse | Navier--Stokes plus six species | `./scripts/prepare_example.sh vascular_flow/multispecies_pulse` |
| 3D VCA closed loop | Two-outlet vascular coupling smoke case | `RANKS=2 ./scripts/prepare_example.sh vascular_flow/vca_bifurcation` |

See the [examples catalog](examples/README.md) for the input contract and case
descriptions.

## Install dependencies

On Ubuntu or WSL Ubuntu, install the preprocessing and CPU build prerequisites
with:

```bash
sudo apt update
sudo apt install \
  build-essential git \
  libeigen3-dev metis \
  openmpi-bin libopenmpi-dev \
  libblas-dev liblapack-dev
```

After cloning the repository and entering its directory, check the selected
backend:

```bash
./scripts/check_dependencies.sh preprocessing
./scripts/check_dependencies.sh cpu       # requires PETSC_DIR
./scripts/check_dependencies.sh one-d    # accepts PETSc pkg-config or PETSC_DIR
./scripts/check_dependencies.sh cuda      # requires the CUDA Toolkit and nvcc
```

The CPU simulation executables additionally require PETSc built with the same
MPI implementation used at runtime. CUDA compilation requires the CUDA Toolkit,
while an NVIDIA GPU and compatible driver are needed only at runtime. See the
[dependency and installation guide](docs/DEPENDENCIES.md) for PETSc setup,
CUDA/Conda on WSL, RHEL-family systems, and Bridges-2.
ParaView/`pvbatch` and Pillow (`python3-pil`) are optional and are needed only
to inspect PVD/VTU results or regenerate the README animations.

## Native 1D quick start

The 1D path consumes only a rooted SWC or radius-annotated line-OBJ centerline
and schema-v3 `simulation_config.json`; it does not run mesh generation,
Python, FEniCS, or HexSim. Coordinates and radii are converted to SI using
`geometry.length_scale_to_m`.

```bash
./scripts/check_dependencies.sh one-d
make one-d-petsc

./solvers/one_d/iga_1d examples/one_d/rigid_straight --check
./solvers/one_d/iga_1d examples/one_d/rigid_straight \
  --output-dir /tmp/tubularflowiga-1d-rigid
```

Open `/tmp/tubularflowiga-1d-rigid/profile_1d.pvd` in ParaView. For MPI and
PETSc options, compliant formulations, transport, checkpoint/restart, SI units,
and the Hex-to-schema-v3 field map, see the [native 1D guide](docs/ONE_D.md).

## Five-minute preprocessing check

The shortest first run needs GNU Make, a C++ compiler, Eigen 3, OpenMP, and
METIS with `mpmetis`. PETSc and CUDA are not required merely to generate and
validate a database.

```bash
git clone https://github.com/EngineerEricXie/TubularFlowIGA.git
cd TubularFlowIGA

./scripts/check_dependencies.sh preprocessing

VASCULAR_WORK="$(mktemp -d /tmp/tubularflowiga-vascular.XXXXXX)"
RANKS=2 ./scripts/prepare_example.sh \
  vascular_flow/straight_tube "$VASCULAR_WORK"
```

The preparation script builds the dependency-free tools and runs mesh
generation, spline extraction, two-way METIS partitioning, database packing,
inspection, configuration validation, and boundary-label validation. A
successful run prints the work directory, database path, and matching CPU/CUDA
solver commands.

The prepared directory contains, among other generated files:

- `skeleton_normalized.swc`: validated, rooted canonical skeleton;
- `skeleton.vtp`: centerline, radius, topology, and branch data for ParaView;
- `mesh_diagnostics.json` and `skeleton_diagnostics.vtp`: geometry feasibility,
  dimensionless quality metrics, and collision candidates;
- `controlmesh.vtk`: labeled hexahedral control mesh;
- `mesh_quality.json`: final Jacobian and surface-intersection results;
- `bzmesh.vtk` and `bzmeshinfo.txt`: Bezier visualization and extraction data;
- `geometry_transform.json`: the source-to-normalized coordinate transform;
- `initial_velocityfield.txt`: generated spatial velocity profile;
- `straight_tube-2.ntiga`: packed solver database.

`prepare_example.sh` requires an empty work directory. Its default `RANKS=2`
is intentional because `mpmetis` does not create a one-partition output.

## Run on CPU with MPI/PETSc

PETSc must be built with the same MPI implementation used at runtime. After
setting `PETSC_DIR` and, when applicable, `PETSC_ARCH`:

```bash
./scripts/check_dependencies.sh cpu
make cpu-petsc PETSC_DIR="$PETSC_DIR" PETSC_ARCH="${PETSC_ARCH:-}"
```

Run the prepared vascular example with exactly the rank count used during
packing:

```bash
VASCULAR_DB="$VASCULAR_WORK/straight_tube-2.ntiga"

mpiexec -np 2 ./solvers/cpu/iga_mesh_check "$VASCULAR_DB"
mpiexec -np 2 ./solvers/cpu/iga_navier_stokes \
  "$VASCULAR_DB" "$VASCULAR_WORK" \
  --output "$VASCULAR_WORK/velocity-cpu.txt"

./solvers/cpu/iga_flow_validate \
  "$VASCULAR_DB" "$VASCULAR_WORK/velocity-cpu.txt"
```

For neuron transport, prepare the neuron case and run its named equation
system:

```bash
NEURON_WORK="$(mktemp -d /tmp/tubularflowiga-neuron.XXXXXX)"
RANKS=2 ./scripts/prepare_example.sh \
  neuron_transport/straight_neurite "$NEURON_WORK"
NEURON_DB="$NEURON_WORK/straight_neurite-2.ntiga"

mpiexec -np 2 ./solvers/cpu/iga_solve \
  "$NEURON_DB" "$NEURON_WORK" \
  --system neuron_transport \
  --output "$NEURON_WORK/neuron-cpu.txt"
```

The flow solver writes three velocity columns to the requested path and one
pressure column to the neighboring `.pressure` file. Configured transport
writes `node_id` followed by fields in configured order; the neighboring
`.fields` file records their names. Both solvers also write a ParaView-ready
final `.vtu`; `--output-every N` additionally creates step-indexed `.vtu`
snapshots and a `.pvd` time collection.

## Run on one CUDA GPU

The CUDA backend does not require PETSc, but compilation requires the CUDA
Toolkit and runtime requires a compatible NVIDIA driver and GPU.

```bash
./scripts/check_dependencies.sh cuda
make cuda CUDA_ARCHS=89

./solvers/cuda/iga_cuda device-info
./solvers/cuda/iga_cuda mesh-check "$VASCULAR_DB"
./solvers/cuda/iga_cuda navier-stokes \
  "$VASCULAR_DB" "$VASCULAR_WORK" \
  --output "$VASCULAR_WORK/velocity-cuda.txt"

./solvers/cuda/iga_cuda solve \
  "$NEURON_DB" "$NEURON_WORK" \
  --system neuron_transport \
  --output "$NEURON_WORK/neuron-cuda.txt"
```

Set `CUDA_ARCHS` to the compute capability of the target GPU. For example,
V100 is `70` and RTX 4080 SUPER is `89`. In WSL, `nvidia-smi` confirms driver
access but does not install `nvcc`; see the [dependency guide](docs/DEPENDENCIES.md)
for native and Conda CUDA Toolkit options. The CUDA flow path also writes a VTK
file next to the requested text output.

## Build and test targets

```bash
make mesh
make mesh-test
make spline EIGEN_DIR=/path/to/eigen3
make cpu
make cpu-test
make one-d-petsc
make one-d-test

make cpu-petsc \
  PETSC_DIR=/path/to/petsc \
  PETSC_ARCH=your-petsc-arch

make cuda CUDA_ARCHS="70 80 89 90"
```

`make cpu` builds the PETSc-free packer, inspectors, validators, and reference
utilities. `make cpu-petsc` builds the MPI simulation executables. MATLAB and
the external TREES Toolbox are optional and are needed only to reproduce the
legacy reference workflow.

## Create or modify a case

For 3D, copy one complete directory from `examples/neuron_transport/` or
`examples/vascular_flow/` to a work directory, then:

1. Edit the SWC or OBJ named by `geometry.file` for the centerline and radii.
2. Edit the schema-v4 `geometry` and `mesh` blocks in `simulation_config.json`
   for smoothing, adaptive spacing, junction clearance, and quality gates.
3. Edit the remaining blocks for fields, equations, time integration, and named
   boundary conditions.
4. Run the preparation pipeline, inspect positive Jacobians and boundary labels,
   then execute the matching solver.

For 1D, copy a directory from `examples/one_d/`, edit its SWC and schema-v3
configuration, then run `iga_1d CASE_DIR --check` before simulation. There is no
control-mesh generation or `.ntiga` packing step.

The mesh generator assigns wall label 0, inlet label 1, and terminal outlet
labels starting at 2. Do not assume a branch label without checking
`iga_case_check` output.

The spline stage translates coordinates by the domain minima and divides them
by the smallest domain-axis extent before writing the IGA representation.
`geometry_transform.json` records that affine map, and version-5 `.ntiga`
databases also store it together with `geometry.length_scale_to_m`.
Consequently, the packed database uses normalized coordinates. Example
viscosity, density, time, velocity, and transport coefficients are internally
consistent numerical values, not automatic SI or patient-specific parameters.
Document and apply a complete nondimensionalization when interpreting a case
physically.

## Repository layout

- `examples/`: source-only neuron, vascular, and validation cases.
- `preprocessing/mesh/`: dependency-free C++ SWC smoothing and control meshes.
- `meshgeneration/`: legacy MATLAB reference and template assets.
- `preprocessing/spline/`: C++11 spline construction and Bezier extraction.
- `solvers/cpu/`: C++17 packer, checks, MPI/PETSc flow, and transport.
- `solvers/one_d/`: native C++17 SWC-network flow, transport, physiology, and PETSc solvers.
- `solvers/cuda/`: FP64 single-GPU backend using the CPU database format.
- `docs/`: installation, pipeline, configuration, and validation.

Large generated meshes, databases, caches, partitions, and results are
intentionally not versioned. The source-only NMO_06840 regression is committed
under `examples/`, but preparing it creates hundreds of MiB of work files.

## Documentation map

| Need | Document |
|---|---|
| Fresh-clone walkthrough | [Quick start](docs/QUICKSTART.md) |
| Linux, WSL, PETSc, MPI, CUDA, and Bridges-2 setup | [Dependencies](docs/DEPENDENCIES.md) |
| Files produced at every pipeline stage | [Pipeline](docs/PIPELINE.md) |
| Fields, operators, time stepping, and solver CLI | [PDE configuration](docs/PDE_CONFIGURATION.md) |
| Native 1D schema, solvers, units, outputs, and Hex field map | [Native 1D guide](docs/ONE_D.md) |
| Run and validate native CPU 3D VCA | [VCA bifurcation case](examples/vascular_flow/vca_bifurcation/README.md) |
| Run the large morphology-derived neuron regression | [NMO_06840 transport](examples/neuron_transport/nmo_06840_bifurcation/README.md) |
| SWC and radius-annotated line-OBJ inputs | [Skeleton formats](docs/SKELETON_FORMATS.md) |
| Boundary labels and supported conditions | [Boundary conditions](docs/BOUNDARY_CONDITIONS.md) |
| CPU solver details | [CPU solver README](solvers/cpu/README.md) |
| CUDA solver details | [CUDA solver README](solvers/cuda/README.md) |
| PSC Bridges-2 scheduler workflow | [Bridges-2](docs/BRIDGES2.md) |
| Numerical evidence and limitations | [Benchmarks](docs/BENCHMARKS.md) |

## Validation status and limitations

The public straight, bent, Y-shaped, and NMO_06840 meshes have positive sampled
Jacobians. Representative CPU/CUDA transport and steady-flow comparisons,
mass-balance results, restart checks, and a rigid straight-tube Womersley gate
are recorded in the [vascular example validation](examples/vascular_flow/VALIDATION.md),
[CPU validation](solvers/cpu/VALIDATION.md), and
[CUDA validation](solvers/cuda/VALIDATION.md).

Current scope limits are important when interpreting results:

- 3D vessel and neurite walls are geometrically rigid; 1D compliant wall laws
  do not constitute 3D FSI;
- native 1D and CPU 3D `vca_closed_loop` coupling are available, but CUDA VCA,
  3D replay/open-loop, arbitrary multi-system coupling, and 3D FSI are not;
- 1D pressure/R/RCR and 3D pressure/R/RC/RCR outlets are reduced terminal-bed
  models, not tissue-resolved circulation;
- the 1D physiology layer is configurable and reduced, not automatically a
  patient-validated model;
- the packed IGA geometry is normalized, so physical interpretation requires a
  documented dimensional scaling;
- reduced metabolism and blood-gas derived fields are available in 1D and 3D,
  but full RBC/PFC chemistry, tissue calibration, and parameter fitting are not;
- neuron transport is not neuron electrophysiology.

For shared clusters, run simulations on allocated compute resources rather
than login nodes and follow the local scheduler policy.

## Performance evidence

On the corrected 35,949-node `NMO_54499_new` case, one V100 completed coupled
Navier--Stokes plus 300-step transport numerical work in 279.40 s versus
563.68 s on 16 CPU ranks (2.02x). CPU and CUDA velocity, pressure, and transport
relative L2 differences were `5.33e-6`, `6.07e-6`, and `7.56e-6`.

Here **legacy** refers to the original
[NeuronTransportIGA](https://github.com/EngineerEricXie/NeuronTransportIGA).
See [BENCHMARKS.md](docs/BENCHMARKS.md) for hardware, timings, comparison scope,
and interpretation limits.

## License

TubularFlowIGA is distributed under the [BSD 3-Clause License](LICENSE).
