# Documentation index

Start with the root [README](../README.md) for capabilities and the shortest
workflows. Documents are grouped by purpose:

```text
docs/
├── *.md            user guides and input/output reference (this level)
├── cases/          walkthroughs of specific showcase cases
├── validation/     numerical contracts, evidence, and benchmark results
├── hpc/            cluster deployment, scaling, and provenance
├── architecture/   runtime design contracts for developers
└── images/         figures used by the README and guides
```

## Getting started

| Task | Document |
|---|---|
| Install Linux, WSL, PETSc, MPI, METIS, or CUDA dependencies | [DEPENDENCIES.md](DEPENDENCIES.md) |
| Follow a fresh-clone vascular and neuron walkthrough | [QUICKSTART.md](QUICKSTART.md) |
| See which physics and discretizations are available | [TECHNOLOGY_MATRIX.md](TECHNOLOGY_MATRIX.md) |
| Run any case through the shared `scripts/solver.py` entry | [SOLVER_ENTRY.md](SOLVER_ENTRY.md) |

## User guides

| Task | Document |
|---|---|
| Understand generated files and stage interfaces | [PIPELINE.md](PIPELINE.md) |
| Configure fields, operators, time stepping, checkpointing, and solver CLI | [PDE_CONFIGURATION.md](PDE_CONFIGURATION.md) |
| Configure wall, inlet, outlet, waveform, and Windkessel conditions | [BOUNDARY_CONDITIONS.md](BOUNDARY_CONDITIONS.md) |
| Select PETSc solvers and option prefixes | [SOLVER_OPTIONS.md](SOLVER_OPTIONS.md) |
| Configure and run native lumped 0D R/RC/RLC networks | [ZERO_D.md](ZERO_D.md) |
| Configure spatially distributed compliant 1D A/Q networks | [ONE_D.md](ONE_D.md) |
| Run the native tetrahedral FEM workflow from a labelled surface or volume | [NATIVE_TET_WORKFLOW.md](NATIVE_TET_WORKFLOW.md) |
| Run the native tetrahedral hydraulic graph CLI | [NATIVE_TET_HYDRAULIC_CLI.md](NATIVE_TET_HYDRAULIC_CLI.md) |
| Run the native tetrahedral species CLI | [NATIVE_TET_SPECIES_CLI.md](NATIVE_TET_SPECIES_CLI.md) |
| Checkpoint and restart coupled graphs | [COUPLED_RESTART.md](COUPLED_RESTART.md) |
| Choose and read visualization output | [VISUALIZATION.md](VISUALIZATION.md) and [NATIVE_TET_VTKHDF.md](NATIVE_TET_VTKHDF.md) |
| Repair legacy midpoint-based 0D/1D VTKHDF files | [NETWORK_VTKHDF_REPAIR.md](NETWORK_VTKHDF_REPAIR.md) |

## Input reference

| Topic | Document |
|---|---|
| SWC and radius-annotated OBJ centerlines | [SKELETON_FORMATS.md](SKELETON_FORMATS.md) |
| Shared geometry vocabulary across IGA, immersed, and FEM routes | [GEOMETRY_DATA_CONTRACT.md](GEOMETRY_DATA_CONTRACT.md) |
| Supported inputs, units, and boundary labels | [SUPPORTED_INPUT_CONTRACT.md](SUPPORTED_INPUT_CONTRACT.md) |
| Public liver CT/SEG source data and licensing | [LIVER_GEOMETRY_CANDIDATES.md](LIVER_GEOMETRY_CANDIDATES.md) |

## Showcase cases

| Case | Document |
|---|---|
| Idealized cube with arterial/venous trees, Darcy tissue, and FSI | [cases/CUBE_DUAL_TREE_FSI.md](cases/CUBE_DUAL_TREE_FSI.md) |
| Passive tracer through the same dual-tree cube | [cases/CUBE_DUAL_TREE_OXYGEN.md](cases/CUBE_DUAL_TREE_OXYGEN.md) |
| Prescribed moving immersed flow | [cases/MOVING_IMMERSED_CASE.md](cases/MOVING_IMMERSED_CASE.md) |
| Native CPU 3D VCA coupling | [VCA bifurcation case](../examples/vascular_flow/vca_bifurcation/README.md) |
| Large morphology-derived neuron regression | [NMO_06840 transport](../examples/neuron_transport/nmo_06840_bifurcation/README.md) |

Runnable inputs and their three-file contract are indexed in
[examples/README.md](../examples/README.md).

## Validation and evidence

The `T1`–`T7` prefixes are validation milestones. Each matches a
`make t<N>-…` audit target in the root `Makefile` and, where applicable, an
evidence file under [`benchmarks/`](../benchmarks).

| Milestone | Topic | Document |
|---|---|---|
| — | Benchmark results and performance | [validation/BENCHMARKS.md](validation/BENCHMARKS.md) |
| — | C++ hexahedral mesher correctness | [validation/MESH_CPP_VALIDATION.md](validation/MESH_CPP_VALIDATION.md) |
| T1 | Surface-to-FEM volume meshing | [validation/T1_FEM_VOLUME_MESH.md](validation/T1_FEM_VOLUME_MESH.md) |
| T2 | Fixed-wall flow (Poiseuille, bifurcation) | [validation/T2_FIXED_FLOW_VALIDATION.md](validation/T2_FIXED_FLOW_VALIDATION.md) |
| T3 | Native IGA Kirchhoff–Love shell | [validation/T3_NATIVE_IGA_SHELL.md](validation/T3_NATIVE_IGA_SHELL.md) |
| T3 | Native tetrahedral hyperelastic solid | [validation/T3_NATIVE_TET_SOLID.md](validation/T3_NATIVE_TET_SOLID.md) |
| T4 | Native tetrahedral ALE | [validation/T4_NATIVE_ALE.md](validation/T4_NATIVE_ALE.md) |
| T5 | Matching ALE–solid interface | [validation/T5_NATIVE_MATCHING_INTERFACE.md](validation/T5_NATIVE_MATCHING_INTERFACE.md) |
| T6 | Prescribed left-ventricle QoI | [validation/T6_NATIVE_LV_QOI.md](validation/T6_NATIVE_LV_QOI.md) |
| T6 | Tetrahedral Darcy tissue flow | [validation/T6_NATIVE_TET_DARCY.md](validation/T6_NATIVE_TET_DARCY.md) |
| T6 | Vessel-to-tissue source map | [validation/T6_NATIVE_VESSEL_TISSUE_SOURCE_MAP.md](validation/T6_NATIVE_VESSEL_TISSUE_SOURCE_MAP.md) |
| T6 | Vessel–0D wall reservoir exchange | [validation/T6_NATIVE_WALL_RESERVOIR_EXCHANGE.md](validation/T6_NATIVE_WALL_RESERVOIR_EXCHANGE.md) |
| T7 | Generic 0D species reservoirs | [validation/T7_GENERIC_ZERO_D_SPECIES.md](validation/T7_GENERIC_ZERO_D_SPECIES.md) |
| T7 | Moving-tetra species transport | [validation/T7_NATIVE_MOVING_SPECIES.md](validation/T7_NATIVE_MOVING_SPECIES.md) |
| T7 | Tetra ALE flow with 0D source and RCR terminal | [validation/T7_NATIVE_ZERO_D_FLOW.md](validation/T7_NATIVE_ZERO_D_FLOW.md) |
| — | Imported VascularFM mesher/flow/transport patches (2026-09-29) | [validation/VASCULARFM_PATCHES_20260929.md](validation/VASCULARFM_PATCHES_20260929.md) |

Backend-specific validation reports live next to their code:
[CPU](../solvers/cpu/VALIDATION.md) and [CUDA](../solvers/cuda/VALIDATION.md).

## HPC and clusters

| Task | Document |
|---|---|
| Build, tier tests, stage scheduler inputs, requeue checkpoints, and collect scaling | [hpc/HPC_DEPLOYMENT.md](hpc/HPC_DEPLOYMENT.md) |
| Record HPC provenance and repeat CPU MPI, serial, and GPU baselines | [hpc/HPC_BENCHMARKS.md](hpc/HPC_BENCHMARKS.md) |
| Run on PSC Bridges-2 | [hpc/BRIDGES2.md](hpc/BRIDGES2.md) |
| Assign MPI resource groups to multidomain graphs | [hpc/MULTIDOMAIN_RESOURCES.md](hpc/MULTIDOMAIN_RESOURCES.md) |

The tools these guides call are in [`scripts/hpc/`](../scripts/hpc).

## Architecture

| Topic | Document |
|---|---|
| Multidomain graph, ports, and coupling iteration | [architecture/COUPLING_ARCHITECTURE.md](architecture/COUPLING_ARCHITECTURE.md) |
| Moving immersed-domain geometry | [architecture/MOVING_DOMAIN_ARCHITECTURE.md](architecture/MOVING_DOMAIN_ARCHITECTURE.md) |
| Fluid–structure interaction | [architecture/FSI_ARCHITECTURE.md](architecture/FSI_ARCHITECTURE.md) |
| Parallel node, element, halo, and surface ownership | [architecture/PARALLEL_OWNERSHIP.md](architecture/PARALLEL_OWNERSHIP.md) |
| MPI failure boundaries and collective error handling | [architecture/MPI_FAILURE_BOUNDARIES.md](architecture/MPI_FAILURE_BOUNDARIES.md) |
| Body-fitted runtime cleanup | [architecture/RUNTIME_CLEANUP.md](architecture/RUNTIME_CLEANUP.md) |
| Coupled checkpoint state contract | [architecture/COUPLED_CHECKPOINT_CONTRACT.md](architecture/COUPLED_CHECKPOINT_CONTRACT.md) |
| Coupled checkpoint bundle format | [architecture/COUPLED_CHECKPOINT_BUNDLE.md](architecture/COUPLED_CHECKPOINT_BUNDLE.md) |

Implementation-specific guides live next to their code:

- [Control-mesh generator](../preprocessing/mesh/README.md)
- [Template-free skeleton-to-tet mesher](../preprocessing/tet/README.md)
- [Spline and Bezier extraction](../preprocessing/spline/README.md)
- [MPI/PETSc CPU backend](../solvers/cpu/README.md) and its [architecture](../solvers/cpu/ARCHITECTURE.md)
- [CUDA backend](../solvers/cuda/README.md) and its [architecture](../solvers/cuda/ARCHITECTURE.md)
- [Workflow and maintenance scripts](../scripts/README.md)
