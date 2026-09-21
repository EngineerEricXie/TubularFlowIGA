# MPI failure boundaries

This document identifies the coordinated error boundaries used by production
MPI entry points. A listed boundary means the code and its regression test can
be located; it does not imply in-process recovery from an MPI or PETSc runtime
failure.

## Coordination contract

- `CollectiveLocalStage` wraps local work only. It catches local exceptions
  and then coordinates the failing rank and diagnostic message.
- PETSc collectives run outside the local callback. Their return codes are
  coordinated afterward with `RequireCollectivePetscSuccess` or the
  corresponding helper.
- Runtimes, adapters, and executors borrow a communicator. Every member must
  enter the same stages in the same order. Control values and complete input
  signatures are built before group agreement.
- A stage or returned-outcome wrapper does not replace coordination inside the
  called runtime. An outer catch cannot recover a peer already blocked in a
  different collective.
- The protocol assumes that all processes remain alive and MPI communication
  still works. MPI initialization failures, internal MPI faults, `SIGKILL`, OOM
  termination, and node loss require job-level recovery.

## Supported entry points

| Entry point | Coordinated stages | Regression coverage |
|---|---|---|
| [CPU flow CLI](../../solvers/cpu/src/iga_navier_stokes.cpp) | Arguments, assets, configuration, boundaries, runtime construction, initialization, steps, VCA, checkpoint, field output, VTKHDF, completion | Input, stepping, output, and injected failures |
| [Configured transport CLI](../../solvers/cpu/src/iga_solve.cpp) | Controls, assets, velocity series, runtime, steps, checkpoint, fields, VTKHDF open/close | CLI, initialization, output, and injected failures |
| [Legacy transport](../../solvers/cpu/src/iga_transport.cpp) | Input snapshots, PETSc options, local assembly, returned errors, cleanup, completion | Separate legacy CLI regression |
| [Mesh check](../../solvers/cpu/src/iga_mesh_check.cpp) | Database and owner checks, owned elements, geometry reduction, result logging | Asset regression, including empty ranks |
| [Assembly smoke](../../solvers/cpu/src/iga_assembly_smoke.cpp) | Input, PETSc options, owned-row construction/insertion, assembly, matrix information, destruction, logging | Assembly regression and returned PETSc errors |
| [Native 1D CLI](../../solvers/one_d/src/iga_1d.cpp) | Input/assets, runtime, implicit solve, initial/step/final output, checkpoint | CLI, checkpoint, and stream failures |
| [Multidomain runner](../../solvers/coupling/src/iga_1d_3d_bifurcation.cpp) | Catalog, registry, flow/species executors, histories, and manifest through a borrowed-communicator runner | Registry, pressure, and species executor regressions |
| [Sequential coupling CLI](../../solvers/coupling/src/iga_1d_3d_explicit.cpp) | Input/runtime, fixed-point or Aitken loop, global convergence, abort outcomes, output, completion | Initialization, strong coupling, and output regressions |
| [Immersed graph loader](../../include/ImmersedFlowCase.hpp) | Local preflight followed by collective steady/transient initialization, mode/clock agreement, rollback | Case and graph regressions; MPI callers use `InitializeDistributed` |
| [FSI ParaView exporter](../../solvers/cpu/src/phase8_compliant_channel_fsi_paraview.cpp) | Rejects multi-rank execution before changing output; single-rank FSI and export | Exporter regression |
| [CUDA CLI](../../solvers/cuda/src/iga_cuda.cu) | Single process and single GPU with top-level exception handling | CUDA CLI and stream regressions; multi-process launch is rejected |
| [Mesh](../../preprocessing/mesh/src/main.cpp) and [spline](../../preprocessing/spline/main.cpp) | Serial entry points with optional OpenMP | Local and OpenMP worker failures |

Packing, inspection, case/configuration checking, flow/transport validation,
and Womersley tools are serial utilities. Their format compatibility remains a
pipeline requirement, but they do not implement an MPI runtime contract.

## Shared runtimes and adapters

| Component | Coordinated behavior | Notes |
|---|---|---|
| `TransientFlowRuntime` | Initialization, `BeginStep`, and `MeasurePorts` construct checked local signatures before agreement | Newton and PETSc calls use the established runtime protocol |
| `TransientTransportRuntime` | `GatherState` and `IntegrateFields` complete their signatures before gather/reduction | Source preparation uses the same checked stream handling |
| `ThreeDBodyFittedFlowTransportDomainAdapter` | Plan/input signatures are checked; candidates are published only after collective success | Hydraulic and transport trial/rollback tests cover the lifecycle |
| `OneDRuntime` and adapter | Runtime-internal coordination is separate from outer stage outcomes; CLI output uses a collective local stage | Retry with a new process and output directory after a writer failure |
| Registry and pressure/species executors | Runtime ownership, stage outcomes, global convergence, and abort ordering are tested | Safety depends on internal runtime coordination, not only an outer wrapper |
