<p align="center">
  <img src="docs/images/iga-wordmark-flow.png"
       alt="ParaView rendering of velocity magnitude through a connected pipe spelling IGA"
       width="100%">
</p>

# CoupledFlow

**CoupledFlow** was originally developed as **TubularFlowIGA**. The repository
URL remains unchanged so existing links in CVs and application materials keep
their original context.

The [shared solver entry](docs/SOLVER_ENTRY.md) runs existing 0D, 1D, IGA,
and native tetra FEM cases with `python3 scripts/solver.py CONFIG.json`.
Editable outer-config examples are in [examples/solver](examples/solver).

CoupledFlow is a native C++ framework for coupled flow and multiphysics
simulation across 0D, 1D, and 3D domains. It combines lumped circuits and
line-network models with tetrahedral finite elements and three-dimensional
isogeometric analysis (IGA) on MPI/PETSc CPUs or one CUDA GPU. The body-fitted
3D pipeline generates hexahedral control meshes, constructs splines and Bezier
extraction, and packs a partition-aware database. Additional solvers provide
Darcy flow, species transport, immersed and moving domains, solid mechanics,
and fluid--structure interaction (FSI).

This is research software for mixed-dimensional coupled simulation. It is not
a general-purpose CFD package.

The [artificial cube dual-tree FSI case](docs/CUBE_DUAL_TREE_FSI.md) now builds
two centerline/radius trees, five conforming tetra regions, conservative
artery-to-Darcy-to-vein flow, and a small-strain vessel-wall feedback example.
An optional [refined artificial input](cases/idealized_cube_dual_tree_refined.json)
uses a 0.10 mm test wall and denser pipe circumference; independent 2/8-rank
functional reruns pass, but mesh convergence and physiology are not claimed.
The [multi-step passive-oxygen functional case](docs/CUBE_DUAL_TREE_OXYGEN.md)
then uses the frozen post-FSI flow to transport a project-owned FEM tracer
through artery, fixed Darcy tissue, and vein, with full tissue-mesh ParaView
time series and a conservative venous-outlet breakthrough check. Its artificial
porosity-one, zero-consumption concentration is not physiological oxygenation.
Its hydraulic coupling currently transfers terminal flux one way; vessel–tissue
pressure continuity, transient FSI accounting, and physiological validation
are not claimed.
The [geometry audit](docs/LIVER_GEOMETRY_CANDIDATES.md) and
[versioned functional case](cases/liver_roi_functional.json) record provenance,
labels, hashes, limits, and the raw-SEG-to-result local rerun command. DOLFIN is not used as the FEM core;
PETSc supplies algebra/MPI and fTetWild may generate tetrahedra.

## What can it simulate?

| Application | Available now | Important boundary |
|---|---|---|
| Vascular flow | Native 0D R/RC/RLC circuits; native 1D rigid Poiseuille/inertance and compliant A/Q networks; CPU/CUDA body-fitted 3D rigid-wall Navier--Stokes; native C++/PETSc tetrahedral P2/P1 fixed-wall FEM vertical slice with manufactured temporal and spatial convergence; native network and CPU 3D `vca_closed_loop` vascular coupling | The tetrahedral FEM route has local idealized tube/Y evidence but not a general production case schema, scalable large-mesh Schur solver, completed fine Poiseuille gate, or physiological validation. 3D VCA requires backward-Euler CPU flow; CUDA VCA and 3D VCA replay/open-loop are unavailable |
| Multiscale circulation | CPU 0D/1D/3D pressure/flow graphs with explicit or strong coupling; conservative 1D/body-fitted-3D species transfer; native tetra ALE FEM functional graph tests with 1D and generic 0D source/terminal RCR species; a configured surface/volume→native tetra hydraulic workflow with explicit or fixed-point coupling, separate outlet RCRs, and file restart; optional disjoint MPI groups for independent domains | Executable graphs require supported acyclic topology. The native hydraulic workflow is one 0D source→tetra→one or more RCRs, not the general IGA graph runner or physiological validation. A local two-outlet Y example passes on Gmsh/fTetWild with bounded loads; high-impedance Newton robustness remains open. Native source→FEM→RCR species has only a fixed-case five-shard file-restart coordinator, not an arbitrary production graph restart; a full closed-loop 0D heart is deferred |
| Immersed and moving flow | CPU closed-surface immersed IGA with cut-cell integration, Nitsche wall conditions, ghost stabilization, distributed PETSc fields, and prescribed moving geometry; native tetrahedral ALE operator, mesh-motion, rollback, and prescribed-motion runtime vertical slice | Immersed motion uses a fixed Eulerian background. The idealized LV has shared-geometry/motion parity and 16-step functional cycles for native ALE (96, 384, and 672 tetrahedra) and single-process immersed with the same optional backflow port model. Full-cycle common-point velocity L2 peaks at 71.37%/58.15%/54.22% across the three ALE meshes, with every point in the shared domain; the errors remain large and the mesh series is not asymptotically converged. Native ALE pressure–volume and kinetic-energy QoI are recorded; an additional second geometry cycle leaves a 27.40%/37.22% velocity-field return error on the 96/384 meshes, so the fluid state is not yet periodic. Spatial/time convergence, full energy balance, wall-enforcement attribution, distributed immersed backflow, remeshing, physiological validation, and a production moving-domain graph CLI remain deferred |
| Thin structures | Native single-patch NURBS Kirchhoff--Love shell element/static vertical slice with membrane, bending, tangent, convergence, traction, constraint, and patch-interface rejection tests; legacy pre-tensioned P1 membrane | The current shell runtime is a small dense reference implementation with isotropic St. Venant--Kirchhoff material; multipatch coupling, anisotropic artery/myocardium material, distributed algebra, contact, and fluid coupling remain deferred |
| Solid mechanics | Native total-Lagrangian tetrahedral compressible neo-Hookean and stabilized mixed P1/P1 displacement-pressure formulations; explicit inverse-elastostatics prestress initialization | Element tangents, dense global Newton, near-incompressible bending/Poisson-sweep locking gates, and loaded-to-unloaded reference reconstruction are verified locally. These are bounded dense reference runtimes; distributed mixed algebra, fibres, active contraction, parameter identification, contact, and patient-specific validation remain deferred |
| Foundational FSI | Distributed moving immersed flow coupled to a bounded single-owner pre-tensioned membrane with owned surface transfer, global strong Dirichlet--Neumann convergence, and dynamic Aitken relaxation; native matching tetra ALE–solid adapters | Immersed route is validated for a small-displacement compliant-channel benchmark at 1/2/4 ranks, including paired restart. The autonomous native ALE-flow/solid path passes transfer conservation, zero-state invariant, nonzero traction-driven response, an idealized one-step compliant-channel functional smoke, independent acceptance gates, rollback, paired commit, single-partition paired restart, and coordinator/Aitken tests. It still lacks a native compliant-channel/elastic-tube physics benchmark, production PETSc fluid adapter, and distributed interface ownership. Nonmatching transfer, monolithic FSI, and valve/contact models are deferred |
| Neuron transport | Configurable two-field `N0`/`Nplus` axonal transport on straight and branching neurites | This is material transport, not membrane voltage, action potentials, synapses, or network electrophysiology |
| Generic biological transport | Config-selected 1D and 3D multispecies transport with reaction, source, wall exchange, metabolism, oxygen capacity, and blood-gas derived fields; a separate native moving-tetra P1 species weak form with `u−w`, previous/current inventories, and a prescribed-velocity external-mesh CLI; native well-mixed 0D species balance with source/RCR staged graph adapters | The physiology layer is a configurable reduced model. Native moving-tetra scalar transport has dense/PETSc-MPI tests, a configurable prescribed-velocity CLI with P1 VTU output, and controlled source→FEM→RCR moving-graph and five-shard restart tests. The CLI does not solve fluid flow or couple graph ports; the fixed-case graph restart is not a generic production coordinator. High-Péclet positivity for arbitrary data and physiological validation remain open. 3D physiology-driven vasodilation is disabled in the rigid-wall transport path |

Immersed transient, moving-flow and FSI runtimes distribute PETSc rows and
cell work across MPI ranks, and also have optional OpenMP volume assembly
within each rank. Body-fitted flow has an optional MPI/OpenMP volume assembly
CLI. Thread configuration and worker validation are documented in the
[CPU guide](solvers/cpu/README.md#optional-openmp-volume-assembly).

Standalone body-fitted CPU and CUDA solvers share the configuration format
and packed `.ntiga` database within their supported feature sets. CUDA
configured transport supports one through eight scalar fields. Multidomain,
immersed, moving-domain, and FSI execution use CPU runtimes.

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

The body-fitted 3D workflow preserves the original file interfaces:

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

Native 1D runs read a centerline and configuration directly. Immersed 3D runs
build a Cartesian cubic B-spline background and cut-cell quadrature from a
closed triangulated surface. They do not require the SWC-to-control-mesh
pipeline. Multidomain runs connect supported native domains through named
pressure/flow and species ports, using SI units and outward-positive flow.

| Configuration | Purpose | Entry point |
|---|---|---|
| Schema v3 | Standalone native 0D circuit or 1D flow and transport | `iga_0d`, `iga_1d` |
| Schema v4 | Standalone body-fitted 3D geometry, mesh, flow, and transport | `prepare_example.sh`, then CPU or CUDA solver |
| Schema v5 graph | Heterogeneous flow-only coupling, including supported 0D and steady/fixed-transient immersed domains | `iga_multidomain_flow --graph-case ROOT --output-dir DIR` |
| Schema v6 graph | Conservative species transport across native 1D and body-fitted 3D domains | `iga_multidomain_flow --graph-case ROOT --output-dir DIR` |

Graph manifests reference each domain's native inputs. Schema v6 rejects 0D
and immersed domains. Moving-anatomy and FSI benchmarks have dedicated runtime
and validation targets described in their architecture guides below.

## Choose a first example

Each standalone body-fitted 3D example contains `skeleton_initial.swc` or
`skeleton_initial.obj` plus a schema-v4 `simulation_config.json`; its
`geometry` and `mesh` blocks configure preprocessing in the same validated
document as the physics. A native 1D example needs only the skeleton and
schema-v3 configuration. Generated meshes, databases, and results are written
to a separate work directory.

| Goal | Recommended case | Command |
|---|---|---|
| First native 0D run | Straight Poiseuille R network | `./solvers/one_d/iga_0d examples/zero_d/steady_resistive_straight --check` |
| Native 0D RC/RCR | Lumped transient bifurcation | `./solvers/one_d/iga_0d examples/zero_d/rc_rcr_bifurcation --check` |
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
descriptions. The [immersed aneurysm chain](examples/vascular_flow/immersed_aneurysm_chain/README.md)
has separate surface/native-domain inputs and dedicated validation targets.

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
MPI implementation used at runtime. CPU and CUDA visualization builds require
HDF5. CUDA compilation requires the CUDA Toolkit, while an NVIDIA GPU and
compatible driver are needed only at runtime. See the
[dependency and installation guide](docs/DEPENDENCIES.md) for PETSc setup,
CUDA/Conda on WSL, RHEL-family systems, and Bridges-2.
ParaView/`pvbatch` and Pillow (`python3-pil`) are optional and are needed only
to inspect PVD/VTU/VTKHDF results or regenerate the README animations.

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

The preparation script builds the preprocessing and PETSc-free CPU tools and
runs mesh generation, spline extraction, two-way METIS partitioning, database packing,
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
result: transient solvers default to temporal Bézier VTKHDF, while steady flow
uses VTU. See the [visualization output guide](docs/VISUALIZATION.md).

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

## Multidomain coupling, immersed flow, and FSI

Build the CPU coupling executables with the same PETSc/MPI installation used
for the native solvers:

```bash
make coupling PETSC_DIR="$PETSC_DIR" PETSC_ARCH="${PETSC_ARCH:-}"
make coupling-test
```

This builds `iga_1d_3d_explicit`, `iga_1d_3d_bifurcation`, and
`iga_multidomain_flow`. The generic runner accepts a schema-v5 or schema-v6
graph root containing `simulation_config.json` and its referenced native case
assets. Its output directory must not already exist. For body-fitted domains,
each packed database must match the MPI rank count.

For native 0D/1D/body-fitted 3D graphs, `--checkpoint-dir ROOT` saves accepted
steps and `--restart-dir ROOT` restores a complete checkpoint in a new MPI job.
Use the same inputs, build, rank membership, and numerical options. See
[graph restart](docs/COUPLED_RESTART.md) for save intervals, complete history
output, and `SIGUSR1` handling at an accepted-step boundary.

An optional schema-v5/v6 `resources` object assigns domains to disjoint MPI
rank groups. Independent 3D domains can solve concurrently while port exchange,
rollback, commit, and species donor reversal remain globally coordinated.
Shared communicator execution remains the default. See
[multidomain resources](docs/MULTIDOMAIN_RESOURCES.md) for the manifest format,
mapping reconstruction, measured tradeoffs, and current limitations.

`SimulationGraph` validates topology and port capabilities;
`DomainRuntimeRegistry` owns the native runtimes. Component executors exchange
boundary data, converge trial states, and prepare every domain before
finalizing the accepted step. Failed trials can be rolled back without
advancing committed physical state. See the
[coupling architecture](docs/architecture/COUPLING_ARCHITECTURE.md) for graph
configuration and runtime contracts.

The numerical milestones include 1D--3D pulsatile coupling, conservative
species transfer, a quasi-static immersed aneurysm chain, a prescribed
idealized left-ventricle cycle, a compliant-channel FSI benchmark, and a
five-domain `0D source -> 3D -> 1D -> two 0D RCR outlets` circulation case.
Focused reproduction targets include:

```bash
make -C solvers/cpu phase6-aneurysm-depth2-regression PETSC_DIR="$PETSC_DIR"
make phase7-lv-closure-test PETSC_DIR="$PETSC_DIR"
make compliant-channel-fsi-test PETSC_DIR="$PETSC_DIR"
make -C solvers/coupling phase9-multiscale-closure-test PETSC_DIR="$PETSC_DIR"
```

Add `PETSC_ARCH` when required by your installation. These are numerical
validation workloads; use an appropriate compute allocation. The depth-2
aneurysm target checks the local immersed Jacobian and conservation; the full
chain has a separate closure target. See the
	[moving-domain architecture](docs/architecture/MOVING_DOMAIN_ARCHITECTURE.md),
	[FSI architecture](docs/architecture/FSI_ARCHITECTURE.md), and
	[benchmark summary](docs/BENCHMARKS.md) for scope, prerequisites, and recorded
	results.

## Build and test targets

```bash
make mesh
make mesh-test
make spline EIGEN_DIR=/path/to/eigen3
make cpu
make cpu-test
make one-d-petsc
make one-d-test
make coupling-test

# Machine-readable build compatibility and bounded test tiers
make hpc-build-manifest HPC_BUILD_MANIFEST=/path/to/cpu-build.json
make hpc-test-unit HPC_TEST_OUTPUT=/path/to/test-results
make hpc-test-mpi PETSC_DIR="$PETSC_DIR" HPC_TEST_OUTPUT=/path/to/test-results
make hpc-test-gpu HPC_TEST_OUTPUT=/path/to/test-results

make cpu-petsc \
  PETSC_DIR=/path/to/petsc \
  PETSC_ARCH=your-petsc-arch

make cuda CUDA_ARCHS="70 80 89 90"
```

`make cpu` builds the PETSc-free packer, inspectors, validators, and reference
utilities. `make cpu-petsc` builds the CPU 3D and native 1D simulation
executables; `make coupling` builds the separate PETSc coupling runners.
`make coupling-test` exercises the PETSc-free graph, runtime, 0D, species,
surface, and FSI contracts. MATLAB and the external TREES Toolbox are optional
and are needed only to reproduce the legacy reference workflow.

## Create or modify a case

For standalone body-fitted 3D, copy one complete case directory from
`examples/neuron_transport/` or `examples/vascular_flow/` to a work directory,
then:

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

For a multidomain case, define a schema-v5 flow graph or schema-v6 species
graph and supply the native assets for each domain. Use the coupling
architecture guide for supported topology, port capabilities, and unit
conversion. Immersed cases supply a closed surface and background-grid
configuration through their dedicated case loader.

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

- `examples/`: source-only applications and editable shared-solver configurations.
- `cases/`: small versioned cross-domain case definitions; local geometry and
  patient-derived data remain ignored.
- `benchmarks/`: compact validation contracts and numerical evidence; raw runs
  and solver output remain ignored.
- `artifacts/`: ignored local workspaces grouped into `cases/`, `examples/`,
  `benchmarks/`, `visualization/`, `evidence/`, and `profiles/`.
- `include/`: shared configuration, I/O, physiology, domain adapters, and coupling contracts.
- `preprocessing/mesh/`: dependency-free C++ SWC smoothing and control meshes.
- `meshgeneration/`: legacy MATLAB reference and template assets.
- `preprocessing/spline/`: C++11 spline construction and Bezier extraction.
- `solvers/cpu/`: C++17 packer, checks, MPI/PETSc flow and transport, immersed/moving flow, and FSI runtimes.
- `solvers/one_d/`: native C++17 SWC-network flow, transport, physiology, and PETSc solvers.
- `solvers/coupling/`: PETSc multidomain runners and coupling validation harnesses.
- `solvers/cuda/`: FP64 single-GPU backend using the CPU database format.
- `scripts/`: the shared solver entry plus geometry, workflow, validation, and
  rendering tools; Python regression tests are under `scripts/tests/`.
- `docs/`: installation, pipeline, configuration, architecture, and validation reports.

Large generated meshes, databases, caches, partitions, and results belong under
`artifacts/` and are intentionally not versioned. The source-only NMO_06840
regression is committed under `examples/`, but preparing it creates hundreds of
MiB of work files.

## Documentation map

- [T1 surface-to-FEM volume mesh](docs/T1_FEM_VOLUME_MESH.md): bounded Gmsh
  baseline plus the fTetWild robust tetrahedralization adapter,
  label/quality/provenance contract, and local regression.
- [T6 native vessel–0D wall reservoir exchange](docs/T6_NATIVE_WALL_RESERVOIR_EXCHANGE.md):
  conservative, fixed-volume tissue-storage computation coupled to native
  tetra species wall labels; one or multiple regions are exposed through the
  standalone CLI with paired restart. This is not a 3D perfusion model.
- [T6 native tetra Darcy foundation](docs/T6_NATIVE_TET_DARCY.md):
  autonomous 3D P1 pressure/flux FEM, standalone labelled-mesh CLI and
  surface/volume fTetWild/Gmsh workflow with functional gates; not yet a
  coupled or calibrated organ-perfusion model.
- [T9 native tetra hydraulic workflow](docs/T9_NATIVE_TET_WORKFLOW.md):
  explicit surface/volume and mesher selection, source→native FEM→RCR graph execution
  with explicit or fixed-point coupling, provenance, safe output policy, and current
  functional-only scope.
- [T9 native tetra prescribed-velocity species route](docs/T9_NATIVE_TET_SPECIES_CLI.md):
  versioned surface/volume workflow and standalone external-mesh P1 ALE transport with explicit diffusion, source, optional implicit first-order decay, labelled wall exchange and finite 0D wall storage,
  inflow donors, conservative budgets, rank-owned concentration VTU, and
  standalone per-step checkpoint/restart;
  this separate route does not solve flow or graph coupling.
- [Generic 0D species reservoir contract](docs/T7_GENERIC_ZERO_D_SPECIES.md):
  autonomous well-mixed multi-species balance and reverse-flow unit tests;
  source/terminal RCR staged graph coupling and a bounded native tetra species
  chain with fixed-test five-shard restart; arbitrary-case production restart
  remains open.
- [Geometry data contract](docs/GEOMETRY_DATA_CONTRACT.md): shared
  reference/current, unit, stable-ID, region-role, and boundary-label semantics.
- [T2 fixed-wall validation contract](docs/T2_FIXED_FLOW_VALIDATION.md): frozen
  Poiseuille physics, signs, analytic QoI, convergence series, and result schema.
- [T3 native IGA shell](docs/T3_NATIVE_IGA_SHELL.md): frozen autonomous
  Kirchhoff--Love contract and the native NURBS/kinematics/patch-system status.

| Need | Document |
|---|---|
| Fresh-clone walkthrough | [Quick start](docs/QUICKSTART.md) |
| Linux, WSL, PETSc, MPI, CUDA, and Bridges-2 setup | [Dependencies](docs/DEPENDENCIES.md) |
| Files produced at every pipeline stage | [Pipeline](docs/PIPELINE.md) |
| Fields, operators, time stepping, and solver CLI | [PDE configuration](docs/PDE_CONFIGURATION.md) |
| Native 1D schema, solvers, units, outputs, and Hex field map | [Native 1D guide](docs/ONE_D.md) |
| Native 0D R/RC/RLC circuits, parameters, and outputs | [Native 0D guide](docs/ZERO_D.md) |
| 0D/1D/3D graphs, ports, species routing, and runtime lifecycle | [Coupling architecture](docs/architecture/COUPLING_ARCHITECTURE.md) |
| Prescribed moving immersed anatomy | [Moving-domain architecture](docs/architecture/MOVING_DOMAIN_ARCHITECTURE.md) |
| Membrane coupling, traction transfer, and foundational FSI limits | [FSI architecture](docs/architecture/FSI_ARCHITECTURE.md) |
| Immersed aneurysm inputs and focused validation | [Immersed aneurysm chain](examples/vascular_flow/immersed_aneurysm_chain/README.md) |
| Repeat CPU MPI, serial fixture, and single-GPU timing, memory, and field comparisons | [HPC benchmark guide](docs/HPC_BENCHMARKS.md) |
| Reproducible builds, test tiers, scheduler staging/requeue, and cross-node scaling | [HPC deployment](docs/HPC_DEPLOYMENT.md) |
| Run and validate native CPU 3D VCA | [VCA bifurcation case](examples/vascular_flow/vca_bifurcation/README.md) |
| Run the large morphology-derived neuron regression | [NMO_06840 transport](examples/neuron_transport/nmo_06840_bifurcation/README.md) |
| SWC and radius-annotated line-OBJ inputs | [Skeleton formats](docs/SKELETON_FORMATS.md) |
| Boundary labels and supported conditions | [Boundary conditions](docs/BOUNDARY_CONDITIONS.md) |
| Per-domain PETSc options for 1D and body-fitted graphs | [Solver options](docs/SOLVER_OPTIONS.md) |
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

- standalone body-fitted 3D vessel and neurite walls remain rigid; two-way FSI
  is currently limited to the small-displacement immersed membrane benchmark;
- immersed and moving flow distribute active cell work and PETSc field rows;
  FSI surface publications use unique owned nodes and global reductions. The
  bounded membrane solve is still centralized on one selected rank, and the
  current moving/FSI checkpoint path is a library/test integration rather than
  a native graph CLI;
- moving anatomy uses a fixed Eulerian background. ALE/remeshing, nonmatching
  FSI transfer, monolithic FSI, and advanced valves/leaflet contact are deferred;
- multidomain execution supports validated acyclic pressure/flow graphs and
  staged 1D/body-fitted-3D species transfer. Optional domain groups support
  0D, 1D, and body-fitted 3D, while grouped checkpoint/restart and grouped
  immersed execution remain deferred. 0D species, a full closed-loop 0D heart,
  multirate coupled clocks, and moving/FSI native graph restart remain deferred;
- native 1D and CPU 3D `vca_closed_loop` coupling remain separate from generic
  domain coupling; CUDA VCA and 3D VCA replay/open-loop are unavailable;
- 0D/1D pressure/R/RC/RCR and 3D pressure/R/RC/RCR outlets are reduced terminal-bed
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
The Bridges-2 workflow supports cross-node graph checkpoint/requeue/restart,
strong FSI paired restart, and 1/64/128/256-rank scaling runs. Reproduce these
checks with the commands in the [HPC benchmark guide](docs/HPC_BENCHMARKS.md)
and [deployment guide](docs/HPC_DEPLOYMENT.md).

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

CoupledFlow is distributed under the [BSD 3-Clause License](LICENSE).
