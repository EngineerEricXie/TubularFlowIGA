# TubularFlowIGA CPU Backend

A project-owned C++17 isogeometric-analysis pipeline for stabilized steady and
backward-Euler Navier-Stokes flow plus configurable transient multi-field transport. PETSc provides distributed
sparse matrices and Krylov solvers; Bezier extraction, basis derivatives,
quadrature, VMS/SUPG weak forms, boundary conditions, nonlinear iteration, and
time integration are implemented in this repository.

This backend replaces the legacy solver. The matching single-GPU implementation is in [../cuda](../cuda).

## Capabilities

The immersed MPI assembly layer now supports uniquely owned cell/ghost
integration, distributed sparse rows, scalar constraints, and required-state
exchange. `ImmersedStaticDistributedOperator` now uses the shared serial/MPI
geometry setup and existing physics kernels, with global port/wall/pressure
diagnostics. `ImmersedStaticDistributedRuntime` adds owned Newton/KSP updates,
global convergence and conservation checks, and transactional commit/rollback.
Its C++ interface and the quasi-static graph entry support distributed steady
solves. Fixed-geometry backward-Euler and moving immersed entry points remain
serial pending further HPC-03C/D work. See
the [assembly report](../../docs/progress/HPC_03A_ASSEMBLY_PROGRESS.md),
[static operator report](../../docs/progress/HPC_03A_STATIC_OPERATOR_PROGRESS.md),
and [static runtime report](../../docs/progress/HPC_03A_STATIC_RUNTIME_PROGRESS.md).

Static distributed constructors also accept
`ImmersedWorkPartition::WeightedContiguous` to distribute cell work using
quadrature and stabilization cost estimates. The default remains `CellCount`;
the measured small case showed no clear speedup. See the
[partition comparison](../../docs/progress/HPC_03B_WORK_PARTITION_PROGRESS.md)
for per-rank work, halo, assembly timings, and reproducible benchmark commands.

The [steady graph MPI report](../../docs/progress/HPC_03C_STATIC_GRAPH_PROGRESS.md)
covers collective port updates, case preflight/initialization, transactional
rollback, and native explicit/fixed/Aitken graph validation.

`ImmersedDistributedVelocityHistory` freezes owned committed velocities and
exchanges only required nodes. `ImmersedDistributedTransientVolume` couples that
snapshot to uniquely owned, frozen body-force quadrature and backward-Euler
volume integration. Both compact and expanded catalogs are covered by the
[transient volume regression](../../docs/progress/HPC_03C_TRANSIENT_VOLUME_PROGRESS.md).
`ImmersedTransientDistributedOperator` combines these inputs with stationary
material-aware walls (including optional inertial impedance), mixed traces,
ghost stabilization, ports, and the gauge. Its
[operator regression](../../docs/progress/HPC_03C_TRANSIENT_OPERATOR_PROGRESS.md)
compares five boundary modes against the serial transient runtime on 1/2/4 ranks.
`ImmersedTransientDistributedRuntime` adds shared MPI Newton solves and owned
committed/prepared fields with a transactional accepted clock. Frozen history
survives Newton rollback; `AbortTrial` releases it before changing controls.
The [runtime regression](../../docs/progress/HPC_03C_TRANSIENT_RUNTIME_PROGRESS.md)
covers two accepted steps, field/conservation parity, failures, split groups,
and weighted work partitions.
`ThreeDImmersedTransientDistributedFlowDomain` supplies graph transactions that
restore accepted controls and discard frozen trial inputs before coupling retry;
see the [adapter regression](../../docs/progress/HPC_03C_TRANSIENT_DOMAIN_PROGRESS.md).
Transient case/native-graph integration remains pending; entry restrictions still apply.

| Executable | Purpose |
| --- | --- |
| `iga_pack` | Validate sparse cache or legacy text and create indexed `.ntiga` |
| `iga_inspect` | Print database metadata and partition statistics |
| `iga_case_check` | Validate case inputs and resolved boundary roles without PETSc |
| `iga_mesh_check` | Evaluate every element at 4x4x4 quadrature points and reject non-positive Jacobians |
| `iga_config_check` | Strictly validate `simulation_config.json` without PETSc |
| `iga_flow_validate` | Integrate boundary flow, check divergence-theorem closure, mass balance, and cycle repeatability |
| `iga_womersley_reference` | Generate a time-aligned analytical straight-tube Womersley velocity series |
| `iga_solve` | Solve configured one-to-many-field transport/operator systems |
| `iga_assembly_smoke` | Exercise distributed owned-row sparse assembly |
| `iga_navier_stokes` | Solve the four-field stabilized steady or transient velocity-pressure system |
| `iga_transport` | Transitional old-input CLI, lowered through the generic assembler |

`iga_solve` writes fields in configured order plus a `.fields` name file.
`iga_navier_stokes` reads flow physics and boundary values from that same config.
`iga_solve` can also read a configured time-resolved velocity manifest and
linearly interpolate flow snapshots onto the transport time grid.

Before constructing its runtime, `iga_navier_stokes` coordinates CLI controls
and local input errors across MPI ranks. Database, control mesh, initial
velocity, selected configuration/legacy parameters, and referenced transient
boundary tables must have identical contents across rank-local copies. Inputs
must remain immutable during execution; missing or non-regular required files
are rejected before parsing. This adds a full startup scan of each asset per
rank. See the [flow input validation report](../../docs/progress/HPC_01C_FLOW_INPUT_PROGRESS.md)
for coverage, legacy numerical limitations, and remaining runtime boundaries.

Flow output now coordinates path preparation, PETSc gather errors, field writes
and output indexes across ranks. Existing non-regular output targets are rejected
before opening, and the final solve summary is printed after those writers return
successfully. Failed writes may leave partial files; atomic publication remains
separate work. See the
[flow output validation report](../../docs/progress/HPC_01C_FLOW_OUTPUT_PROGRESS.md).

CPU flow and configured transport now explicitly finalize VTKHDF before their
success summaries and coordinate close errors across ranks. The shared writer
also serves CUDA flow and transport. Geometry reports and coupling histories
check the stream after closing. See the
[I/O finalization report](../../docs/progress/HPC_01C_IO_FINALIZATION_PROGRESS.md)
for native close-failure injection, HDF5 roundtrips, and remaining limitations.

The native flow CLI also coordinates timestep configuration, waveform and VCA
inlet preparation, port-result processing, transport budgets, and circuit/history
updates. A failed VCA step ends the job; this does not provide whole-step rollback
for an in-process retry. See the
[flow step validation report](../../docs/progress/HPC_01C_FLOW_STEP_PROGRESS.md).

`iga_solve` also coordinates input and timestep errors. It validates the
configuration, labels, selected velocity source and referenced boundary tables
across ranks; snapshot files are checked when selected for interpolation.
Unused snapshots and temporal definitions remain unopened. Keep all case inputs
immutable during execution. Checkpoint reads accept verified identical local
replicas, while checkpoint writes still require shared storage. State and
metadata are published separately, so a complete crash-safe restart bundle
remains pending. See the [transport CLI validation report](../../docs/progress/HPC_01C_TRANSPORT_CLI_PROGRESS.md)
for failure, numerical, memory-report, VTKHDF and restart evidence.

The legacy `iga_transport` CLI also checks packed-database, parameter, mesh,
selected velocity and optional case-configuration contents across rank-local
replicas. It compares effective PETSc options, coordinates local assembly and
returned PETSc errors, and releases its matrix/vector/solver owners before the
success summary. Its positional arguments and numerical defaults are preserved.
See the [legacy transport report](../../docs/progress/HPC_01C_LEGACY_TRANSPORT_PROGRESS.md)
for native faults, cleanup checks, replica compatibility and remaining limits.

`iga_mesh_check` and `iga_assembly_smoke` require identical database contents
across rank-local paths. The mesh checker also rejects out-of-range owners and
disagreement between the ownership index and loaded element records; an empty
rank remains valid. Assembly smoke coordinates returned PETSc errors and checks
matrix destruction before printing its result. See the
[tool asset report](../../docs/progress/HPC_01C_TOOL_ASSET_PROGRESS.md).

Schema-v5 pressure-flow graph and sequential explicit executors coordinate
operation outcomes and require a common convergence decision before commit.
Schema-v6 species execution also checks donor ownership and transport order
before starting transport, and captures port observations before its local
precommit callback. Registry construction now coordinates storage and index
validation before adopting runtime owners; see the
[registry report](../../docs/progress/HPC_01C_REGISTRY_PROGRESS.md),
[pressure-flow report](../../docs/progress/HPC_01C_PRESSURE_EXECUTOR_PROGRESS.md)
and [species report](../../docs/progress/HPC_01C_SPECIES_EXECUTOR_PROGRESS.md).
Sequential CSV and manifest writers check explicit close outcomes before
publishing the graph completion marker; see the
[sequential output report](../../docs/progress/HPC_01C_SEQUENTIAL_OUTPUT_PROGRESS.md).
Sequential initialization now coordinates local preparation, checks effective
controls/options and selected input contents, and supports identical rank-local
replicas on a borrowed communicator; see the
[sequential initialization report](../../docs/progress/HPC_01C_SEQUENTIAL_INITIALIZATION_PROGRESS.md).
The sequential fixed/Aitken loop now coordinates operation outcomes, global
convergence and all three abort attempts. Native fault/retry tests cover
1/2/3-rank communicators and verify committed snapshots and PETSc ownership
release; see the [strong loop report](../../docs/progress/HPC_01C_SEQUENTIAL_STRONG_PROGRESS.md).
Configured transport Bezier VTKHDF initialization also coordinates standard and
nonstandard exceptions; see the [initialization report](../../docs/progress/HPC_01C_TRANSPORT_VISUALIZATION_INIT_PROGRESS.md).
Configuration and metadata readers now reject incomplete stream reads, and the
selected control, JSON, CSV and checkpoint serializers propagate stream errors;
see the [checked text report](../../docs/progress/HPC_01C_CHECKED_TEXT_PROGRESS.md).
The 1D checkpoint metadata/fingerprint, VTU filenames, VTKHDF array schema and
resource report also propagate intermediate stream failures. The allocation
sweep and native compatibility results are documented in the
[text helper report](../../docs/progress/HPC_01C_TEXT_HELPER_PROGRESS.md).
Run the focused MPI regression with `make text-helper-failure-test PETSC_DIR=...`.
Flow initialization, step controls, port measurements, transport integrals and
staged adapter signatures also reject stringstream formatting errors before
agreement. The [MPI boundary index](../../docs/architecture/MPI_FAILURE_BOUNDARIES.md)
maps supported entries to tests and identifies remaining gaps.
The remaining helper and supported-entry audit is still open.

## Optional OpenMP volume assembly

`ImmersedTransientFlowRuntime` can compute volume element systems with OpenMP.
Its moving-flow and FSI wrappers use the same implementation. PETSc reads,
matrix/vector insertion, wall/port terms and ghost penalties execute on the
MPI initialization thread; volume results are inserted in the original order.
`TransientFlowRuntime` also supports body-fitted volume assembly with MPI and
OpenMP. Each rank gathers its required nodal values on the initialization
thread, computes private element systems, and inserts them in the original
order. Pressure traction integration remains on the initialization thread.
Static immersed flow and body-fitted transport assembly remain serial within
each rank.

```bash
make compliant_channel_fsi_openmp_test PETSC_DIR=/path/to/petsc
OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4 OMP_DYNAMIC=FALSE \
OMP_PROC_BIND=close OMP_PLACES=threads \
OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 \
taskset -c 0,2,4,6 ./compliant_channel_fsi_openmp_test
```

Adapt CPU IDs to the allocated machine. The ordinary
`compliant_channel_fsi_test` target remains a build without OpenMP unless the
caller adds OpenMP flags. An OpenMP build defaults to `omp_get_max_threads()`;
`IGA_ASSEMBLY_THREADS=1` explicitly selects serial assembly. An explicit
request above one thread is rejected by a build without OpenMP. Multiple
threads require MPI to provide at least `MPI_THREAD_FUNNELED`.

The default batch capacity is `min(requested_threads, 8)` elements.
`IGA_ASSEMBLY_BATCH_SIZE` overrides it with a positive count; this is an element
count, not a byte budget. Both settings are captured at runtime construction.
`IGA_PROFILE=1` reports actual team size and maximum resident batch items.
Failures join workers before propagation; retry completes pending PETSc
insertions and clears them before reassembly. Runtime calls remain confined
to the initialization thread. Numerical evidence and remaining acceptance
work are recorded in [HPC-02 progress](../../docs/progress/HPC_02_VOLUME_PROGRESS.md).

For body-fitted flow, build `iga_navier_stokes_openmp` with the same PETSc
installation as the ordinary CLI. For example, on a single allocated node:

```bash
make iga_navier_stokes_openmp PETSC_DIR=/path/to/petsc
OMP_NUM_THREADS=2 OMP_THREAD_LIMIT=2 OMP_DYNAMIC=FALSE \
OMP_PROC_BIND=close OMP_PLACES=cores IGA_ASSEMBLY_THREADS=2 \
OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 \
mpiexec --map-by slot:PE=2 --bind-to core -np 2 \
  ./iga_navier_stokes_openmp DATABASE_FOR_2_RANKS.ntiga CASE_DIR --output flow.txt
```

This Open MPI example allocates two physical cores per rank; scheduler jobs
must request matching resources. A failed worker batch is joined before the
cross-rank error agreement. Retry flushes pending matrix insertions before
clearing them, retaining preallocation even if the first assembly failed.
See [body-fitted hybrid progress](../../docs/progress/HPC_02_HYBRID_PROGRESS.md)
for tested configurations and outstanding performance comparisons.

## Embedding MPI runtimes

Call `RequireExecutionResources(comm, &stream)` from the MPI initialization
thread before numerical work to validate numeric thread settings and the
real/double PETSc ABI. It reports communicator size, PetscInt width, MPI thread
support, compiled OpenMP, effective OpenMP thread limits, and BLAS thread requests.
Multiple OpenMP threads require at least `MPI_THREAD_FUNNELED`. Legal differences
between ranks are retained and reported as ranges; the report does not measure
active workers or certify CPU binding. The native flow, configured transport,
1D, and graph entrypoints call this preflight automatically. See the
[resource validation report](../../docs/progress/HPC_01D_PROGRESS.md).
`iga_mesh_check`, `iga_assembly_smoke`, and legacy `iga_transport` also perform
this preflight and coordinate input failures. Assembly `FIELDS` must be a full
positive integer within the PETSc row capacity; legacy `STEPS` must be a full
nonnegative integer, and ranks must agree on step count and whether to output.
See the [auxiliary tool report](../../docs/progress/HPC_01CD_TOOLS_PROGRESS.md)
for before/after field comparisons, bad-Jacobian exit codes, and failure coverage.

`OwnedRowAssembler`, `TransientFlowRuntime`, and `TransientTransportRuntime`
borrow the communicator supplied to their constructors. The owner must keep it
valid until these objects and their PETSc resources are destroyed, before
`PetscFinalize`. Their local ghost vectors use `COMM_SELF`; distributed solves,
reductions, and transport checkpoint operations use the supplied group.

Call flow `InitializeState` collectively while the runtime is committed. It
compares the resolved initial boundaries across the group, prepares candidate
fields in existing scratch vectors, and preserves the public state Vec handle.
Transport `GatherState` is also collective and validates the ordered field
layout before gathering. It retains raw values, including NaN and infinity for
diagnostics, and allocates the complete field on every rank; use the required-node
query for element work. See the
[initialization and gather report](../../docs/progress/HPC_01C_INITIALIZATION_PROGRESS.md)
for failure coverage, compatibility tests, and remaining constructor/I/O work.

Transport `ReadState` is collective and requires the committed phase. Native
flow restart and this transport reader validate the path, exact binary length,
Vec header, and finite coefficients before publishing a candidate vector. Each
rank reads its owned rows from the same immutable regular file on shared storage;
the reader accepts an ordinary uncompressed PETSc binary Vec matching the current
PetscInt/PetscScalar build and ignores adjacent `.info` options. A successful
transport load retains the existing warm-start convention (`Steps() == 1`).
The [checkpoint read report](../../docs/progress/HPC_01C_CHECKPOINT_READ_PROGRESS.md)
records failure and restart coverage. Coordinated checkpoint publication and
whole-coupled-state recovery remain HPC-05 work.

Transport `WriteState` is collective and accepts only committed state. The native
flow and transport writers copy owned coefficients, reject nonfinite input, and
write disjoint sections of a sibling temporary file on shared storage. After all
ranks finish and sync their writes, rank 0 renames the file onto the destination.
Caught errors before that rename preserve the previous Vec file and trigger
temporary-file cleanup. Use one writer group per destination; separate groups
must use distinct paths. New files have owner-only permissions; replacement
preserves an existing regular file's permission bits. Nonregular destinations,
including symlinks, are rejected. See the
[checkpoint write report](../../docs/progress/HPC_01C_CHECKPOINT_WRITE_PROGRESS.md)
for format compatibility and real I/O failure tests. Metadata and coupled state
are still separate publications, so this is not a complete restart transaction.

The synchronous [RunMultidomainFlow API](../../include/MultidomainRunner.hpp)
runs a complete graph on a borrowed communicator. Initialize PETSc and configure
its options before calling it; the function parses graph arguments and does not
initialize/finalize PETSc or insert a new PETSc options database. Compile the
coupling runner with `IGA_MULTIDOMAIN_NO_MAIN` when embedding it. Independent
graphs require distinct output paths. Existing CLIs continue to select world.
This is whole-graph isolation; assigning separate groups to individual domains
remains HPC-08. See the [communicator report](../../docs/progress/HPC_01A_PROGRESS.md)
for tests, numerical limits, and the currently serial immersed/FSI paths.

`OwnedRowAssembler::RequiredRows` validates the global IDs and PETSc row
arithmetic used by the runtime scatters. Its optional collective
`ValidateOwnership()` audit checks partition integration, minimum-node
integration and owned-row element contributions independently. Use
`make -C solvers/cpu parallel-ownership-test PETSC_DIR=...` from the repository
root for the three-rank fault-injection and empty-rank tests. Surface publication
ownership has a separate MPI adapter; see the
[ownership contract](../../docs/architecture/PARALLEL_OWNERSHIP.md).

The `OwnedRowAssembler` constructor and both `CreateMatrix` overloads now
require participation by every member of the supplied communicator. They
coordinate local database/preallocation exceptions before PETSc creation.
`CollectiveLocalStage` is for local callbacks only: all group members must call
the same stages in the same order, and callbacks must not enter MPI/PETSc
collectives. It propagates the lowest failing group rank's bounded diagnostic.
Body-fitted element assembly and graph preflight/output use this protocol;
coverage of other failure boundaries remains in progress. See the
[failure protocol report](../../docs/progress/HPC_01C_PROGRESS.md) for the exact
tested scope and the controlled-failure regression command.

Graph execution now compares captured manifest, 0D model and native 1D/3D
configuration bytes within each communicator before creating collective
runtimes. Identical copies may live at different paths; different formatting
also counts as different input bytes. Stop/Newton/injection controls must agree.
The assembler additionally checks global node, element and field counts.
These checks do not establish identity of all geometry, waveform or packed
database contents. See the
[configuration progress report](../../docs/progress/HPC_01C_CONFIGURATION_PROGRESS.md).

Graph startup also compares the visible PETSc option entries within its supplied
communicator, preserving unused-option tracking. Application arguments are checked
separately; different `-options_file` locations are allowed when the loaded entries
agree. Separate communicators may use different options. This checks the startup
database, not later changes or preconfigured PETSc objects. The native 1D CLI now
uses the same startup check; standalone CPU CLI coverage remains pending. See the
[PETSc options report](../../docs/progress/HPC_01C_PETSC_OPTIONS_PROGRESS.md) and
[1D CLI report](../../docs/progress/HPC_01C_ONE_D_CLI_PROGRESS.md).

## Why this version is faster and smaller

The legacy MPI programs made every rank scan large ASCII extraction files,
stored dense extraction data, preallocated 1,000 matrix entries per row, and
relied on PETSc off-process insertion. This version changes the data and
assembly path:

- Pack the sparse preprocessing cache, mesh, and partition data into a binary
  database with direct element offsets and per-rank touching-element indices.
- Store only nonzero Bezier extraction rows and seek directly to records needed
  by each rank.
- Assemble only owned matrix rows and derive exact diagonal/off-diagonal
  preallocation from symbolic node adjacency.
- Turn unexpected sparse allocation into an error with
  `MAT_NEW_NONZERO_ALLOCATION_ERR` instead of silently growing matrices.
- Keep element matrices in temporary contiguous buffers and release them after
  insertion.
- Assemble time-independent transport operators once, reuse KSP and
  preconditioners, and reuse the current state as the next initial guess.
- Validate input lengths, partition consistency, Jacobian signs, allocation,
  and `KSPConvergedReason` collectively before accepting results.
- Run the full Navier-Stokes nonlinear loop instead of the legacy hard-coded
  single update.

See [ARCHITECTURE.md](ARCHITECTURE.md) for implementation details.

## Measured improvement over the legacy solver

All measurements below were collected on PSC Bridges-2 with optimized PETSc and
OpenMPI 4.0.5. The cylinder contains 4,221 nodes and 3,600 elements.

| Cylinder measurement | Legacy release | TubularFlowIGA CPU |
| --- | ---: | ---: |
| Comparable first Newton update | 199 s | about 12 s |
| First-update speedup | 1.0x | about 16.6x |
| Peak memory | 6.72 GiB | 1.51 GiB |
| Memory reduction | baseline | about 77.5% |
| Complete converged solve | not implemented | 94 s, 7 completed updates |

The legacy timing is one hard-coded update and is not a converged nonlinear
reference. The CPU solver completes all seven updates in less than half the
legacy single-update time while using about 4.45x less peak memory.

Two cylinder transport steps reproduce the legacy fields to relative L2 errors
of `4.18e-7` for `N0` and `3.13e-7` for `Nplus`. The converged default
block-Jacobi/ILU Navier-Stokes run at the time of this benchmark used 94 s and
1.51 GiB. A different measured Schur field-split alternative used 180 s and
1.67 GiB. These historical measurements do not characterize the current
size-aware GAMG Schur default described below; rebenchmark it before making a
new Bridges-2 performance claim.

### Larger validated CPU case

`NMO_66748_subtree` contains 57,456 nodes and 50,940 elements. Its 16-way
database is 507 MiB and all quadrature samples are positive
(`min(detJ)=8.26e-9`).

| Stage, 16 MPI ranks | Result |
| --- | ---: |
| Navier-Stokes convergence | 2 updates, 1,929 Krylov iterations |
| Navier-Stokes solver time | 401 s |
| Navier-Stokes peak RSS | 10.0 GiB |
| Two-step transport assembly | 32.5 s |
| Two-step transport solves | 5.74 s, 372 iterations |
| Transport total wall time | 45 s |
| Transport peak RSS | 6.34 GiB |

`NMO_54499_new` was not accepted for this CPU validation: 6 of 31,500 elements
contained 43 non-positive quadrature samples with minimum determinant
`-5.67e-6`. The legacy solver silently assembled signed volumes; this solver
fails early. Detailed numerical checks are recorded in
[VALIDATION.md](VALIDATION.md).

## Repository layout

- `include/`: database, case input, element kernels, and owned-row assembly.
- `src/`: packer, inspector, mesh checker, smoke test, and production solvers.
- `slurm/`: Bridges-2 PETSc build workflow.
- `ARCHITECTURE.md`: data layout and numerical design.
- `VALIDATION.md`: measured regression and large-case evidence.

Generated executables, `.ntiga` databases, simulation results, and Slurm output
are intentionally ignored.

## Requirements

- A C++17 compiler.
- MPI; Bridges-2 validation used OpenMPI 4.0.5.
- PETSc with C++ and MPI support; the validated build uses optimized
  `arch-linux-c-opt`.
- HDF5 development headers and library for temporal VTKHDF output. The writer
  uses serial HDF5 calls from rank zero; do not mix a parallel HDF5 built for a
  different MPI implementation with PETSc.
- METIS/`mpmetis` to create a partition matching the MPI rank count.
- Spline outputs: `bzmeshinfo.txt` and `spline_cache.igacache`.
- Optional legacy outputs: `cmat.txt` and `bzpt.txt`.

## Build

From the repository root:

~~~bash
make cpu
make cpu-petsc PETSC_DIR=/path/to/petsc PETSC_ARCH=your-petsc-arch
~~~

The first command builds `iga_pack`, `iga_inspect`, `iga_case_check`,
`iga_config_check`, `iga_flow_validate`, and `iga_womersley_reference` with a
standard C++17 compiler.
PETSc targets require an MPI-aware PETSc configuration.
`PETSC_ARCH` may be omitted for an installed PETSc layout.

## Prepare a case

Keep large cases outside this repository:

~~~bash
export IGA_CPU_ROOT="$(pwd)/solvers/cpu"
export CASE_DIR=/path/to/case
export DATABASE=/path/to/work/case-8.ntiga
export RANKS=8

mpmetis "$CASE_DIR/bzmeshinfo.txt" "$RANKS"
"$IGA_CPU_ROOT/iga_pack" "$CASE_DIR" "$RANKS" "$DATABASE"
"$IGA_CPU_ROOT/iga_inspect" "$DATABASE"
~~~

The partition suffix and `mpiexec -np` count must match the rank count passed
to `iga_pack`. The packer prefers the sparse cache. Add `--legacy-text` after
`$DATABASE` to force the legacy reader for comparison. New databases include
per-element boundary-face labels; version 3 databases remain readable but do
not provide surface metadata.

## Run

Use an allocated compute resource on shared clusters:

~~~bash
mpiexec -np "$RANKS" "$IGA_CPU_ROOT/iga_mesh_check" "$DATABASE"
mpiexec -np "$RANKS" "$IGA_CPU_ROOT/iga_navier_stokes" \
  "$DATABASE" "$CASE_DIR" --max-newton 12 --output velocity.txt
mpiexec -np "$RANKS" "$IGA_CPU_ROOT/iga_solve" \
  "$DATABASE" "$CASE_DIR" --system neuron_transport --output neuron.txt
~~~

Use the Navier--Stokes command with a prepared `vascular_flow/*` case and the
transport command with a prepared `neuron_transport/*` case. The
application-specific source cases and preparation command are documented in
the [examples catalog](../../examples/README.md). Most small configurations
contain one matching equation system; `vascular_flow/multispecies_pulse`
intentionally runs flow first and then six-species transport from its generated
velocity snapshot manifest.

The canonical [`simulation_config.json`](../../docs/PDE_CONFIGURATION.md) names
fields, systems, operators, viscosity, time integration, and per-field boundary
conditions. A velocity Dirichlet condition can use a three-component value or
`initial_velocityfield.txt` profile plus scale. The old two-file input remains
accepted only as a transition adapter. Navier-Stokes uses nonlinear relative
tolerance `1e-5`, absolute residual RMS tolerance `1e-10` per equation, and
boundary mass-imbalance tolerance `1e-3`; configured transport uses relative
tolerance `1e-8`. Override the flow criteria with `--nonlinear-rtol R`,
`--nonlinear-atol A`, and `--mass-rtol R`. The solver accepts a Newton state
only when the nonlinear and mass criteria both pass and prints L2/RMS,
continuity, and boundary-flow diagnostics for every iteration.

The default Navier--Stokes preconditioner is size-aware. Cases below 1,000
control points retain block Jacobi, which is efficient and deterministic for
small validation problems. Larger mixed velocity-pressure systems use a full
Schur field split with the stabilized pressure block (`A11`) and PETSc GAMG for
both subproblems. Standard PETSc options still take precedence; for example,
`-pc_type bjacobi` explicitly selects the legacy preconditioner.

Backward-Euler flow supports time-indexed text output and a PETSc checkpoint
state with validated JSON metadata:

~~~bash
mpiexec -np "$RANKS" "$IGA_CPU_ROOT/iga_navier_stokes" \
  "$DATABASE" "$CASE_DIR" --output velocity.txt --output-every 10 \
  --checkpoint flow-checkpoint --checkpoint-every 10

mpiexec -np "$RANKS" "$IGA_CPU_ROOT/iga_navier_stokes" \
  "$DATABASE" "$CASE_DIR" --output resumed.txt \
  --restart flow-checkpoint
~~~

Transient output defaults to temporal Bézier VTKHDF; steady output defaults to
VTU. Use `--visualization-format auto|vtkhdf|vtu` to select explicitly. File
layout, geometry validation, restart behavior, and ParaView validation are in
the [visualization output guide](../../docs/VISUALIZATION.md).

The metadata validates node count, completed step, physical time, `dt`, density,
and viscosity before restart. A final checkpoint is written even when the final
step is not an interval boundary. Legacy positional `MAX_NEWTON` and `OUTPUT`
remain accepted; the default maximum is 12 nonlinear iterations.
For a reproducible interruption test, add `--stop-after-step N`; the checkpoint
retains the original configured final step, so the same case can resume without
editing its time configuration.

Configured transport uses equivalent flags and validates ordered field names,
system name, and velocity-source name in addition to time metadata. A
snapshot-series restart therefore resumes interpolation at the completed
physical time:

~~~bash
mpiexec -np "$RANKS" "$IGA_CPU_ROOT/iga_solve" \
  "$DATABASE" "$CASE_DIR" --system neuron_transport \
  --output neuron.txt --output-every 10 --stop-after-step 10 \
  --checkpoint neuron-checkpoint --checkpoint-every 10

mpiexec -np "$RANKS" "$IGA_CPU_ROOT/iga_solve" \
	"$DATABASE" "$CASE_DIR" --system neuron_transport \
	--output neuron-resumed.txt --restart neuron-checkpoint
~~~

Add `--memory-report transport-memory.jsonl` to collect one JSON object per
stage after database loading, required-element expansion, matrix
preallocation, first operator assembly, first KSP setup, and the completed
solve loop. The report
contains current and peak RSS, PETSc process memory, and PETSc allocator usage
for every rank plus maximum and aggregate values. Allocation tracking is
enabled before `PetscInitialize` only when this option is present, so normal
runs do not pay the PETSc malloc-debug overhead.

Configured transport derives separate left/previous field-coupling patterns
from the compiled terms and Robin/Dirichlet requirements. PETSc scalar AIJ
preallocation and reusable element `nen x nen` blocks include only those active
equation/trial pairs; `MAT_NEW_NONZERO_ALLOCATION_ERR` remains enabled.

Distributed transport diagnostics and Flux/Robin boundary preflight assign
each required element to the PETSc owner of its minimum global connectivity
node. This counts every packed element exactly once without loading a second
METIS-owned element copy; METIS `element.owner` is not the ownership rule for
the row-owner-loaded assembler set.

Validate one velocity snapshot, or a complete manifest with an optional cardiac
period for cycle-to-cycle comparison:

~~~bash
"$IGA_CPU_ROOT/iga_flow_validate" "$DATABASE" velocity.txt
"$IGA_CPU_ROOT/iga_flow_validate" "$DATABASE" \
  --compare velocity-reference.txt velocity-current.txt
"$IGA_CPU_ROOT/iga_flow_validate" "$DATABASE" \
  --manifest "$CASE_DIR" velocity_series.csv 0.8
"$IGA_CPU_ROOT/iga_flow_validate" "$DATABASE" \
  --compare-manifests cpu-output velocity_series.csv \
  "$CASE_DIR" velocity_series.csv
"$IGA_CPU_ROOT/iga_flow_validate" "$DATABASE" \
  --womersley examples/validation/womersley/womersley_reference.json \
  "$CASE_DIR" velocity_series.csv
~~~

The reported relative mass imbalance is
`2*abs(sum(Q))/sum(abs(Q))`, using outward flow on every packed boundary label.
The divergence-theorem error is also normalized by the larger of total
absolute boundary flow and absolute volume divergence, making its gate portable
across geometric and unit scales.
Manifest mode also reports the maximum velocity relative L2 between snapshots
separated by the supplied period. The utility accepts either CPU or CUDA
three-column velocity output. `--compare-manifests` requires every reference
time to appear in the numerical manifest (which may contain additional startup
snapshots) and reports the per-snapshot and maximum velocity relative L2,
suitable for CPU/CUDA coefficient parity or an independently projected
coefficient reference. Raw analytical point samples are not IGA control
coefficients. `--womersley` instead evaluates the numerical spline and the
analytical field at volume quadrature points and integrates their physical
relative L2 norm, which is the Womersley acceptance metric.
The generator configuration and pressure-gradient convention are documented in
[`examples/validation/womersley`](../../examples/validation/womersley/README.md).

For PSC module, interactive-node, and Slurm examples, see the
[Bridges-2 guide](../../docs/BRIDGES2.md).

## Reproducibility notes

- Record compiler, MPI, PETSc commit/configuration, rank count, node type, and
  PETSc options with every benchmark.
- Compare numerical fields as well as wall time; changing partition or
  preconditioner can alter reduction order.
- Run `iga_mesh_check` before timing assembly or solves.
- Benchmark from allocated compute nodes and report both stage timers and peak RSS.
- The large validation databases and outputs are not versioned; use matching
  case snapshots when reproducing the published numbers.

## Current limitations

- The CPU solver depends on PETSc for distributed sparse linear algebra.
- Scaling has been validated to 16 MPI ranks, not exhaustively across nodes.
- The binary database currently stores partition-specific touching-element
  indices and should be repacked when the rank count changes.
- CPU Navier-Stokes supports steady and backward-Euler rigid-wall flow,
  per-step velocity Dirichlet and pressure-traction waveforms,
  checkpoint/restart, and time-indexed output.
- Configured scalar flux/Robin surface assembly is available in `iga_solve`
  with a version 4-or-later `.ntiga` database.
- Natural pressure traction and R/RC/RCR outlets are available on CPU. The
  outlet models apply `-p n` without replacing continuity rows; compliant walls
  are not implemented.
- The stabilized formulations follow the legacy model; this repository does not
  claim a new physical model or discretization order.
