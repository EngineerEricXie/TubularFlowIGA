# Phase 0 Progress Report

Status: complete — implementation revision `dc05f8a` passes Gates A, B, and C.

## Objective

Prepare the existing 1D and body-fitted 3D solvers for reusable direct
coupling without changing validated standalone or VCA numerical behavior.

## Completed work

### PR 0.1 baseline and runtime-boundary audit

- Added the next-generation roadmap at its canonical path.
- Audited the 3D transient runtime, 3D transport/VCA orchestration, native 1D
  runtime, checkpoint boundaries, mutable physical state, and output side
  effects.
- Defined the proposed subsystem lifecycle, in-memory rollback semantics,
  explicit SI port quantities, global outward-positive sign convention,
  compatibility matrix, first 1D--3D benchmark, and planned Phase 0 file
  boundaries in `docs/architecture/COUPLING_ARCHITECTURE.md`.

No solver source, schema, database format, or runtime behavior changed in this
slice.

### PR 0.2 common coupling port values

- Added dependency-free C++17 port data in `include/CouplingPort.hpp`:
  `PortQuantity`, `PortOrientation`, `PortState`, `PortBoundaryData`,
  `CouplingPort`, and `CouplingResidual`.
- The fields state SI units in their names. Flow and species flux are explicitly
  positive outward from their subsystem; `PortOrientation` converts only a
  backend-native sign to that convention.
- Optional boundary values preserve the distinction between absent data and a
  supplied zero. Optional port-state measurements use the same representation,
  so a port that does not provide an area or pressure does not report a fake
  zero value.
- Added validation for finite values, positive supplied area, nonempty port
  identifiers and locators, orientation signs of exactly `+1` or `-1`, duplicate
  raw quantity declarations, contradictory provide/require declarations, and
  port IDs unique within each subsystem.
- Added the dependency-free `coupling_port_test` to `make cpu-test`. It covers
  validation, orientation conversion, zero-versus-missing boundary data, and
  conservative flow/species edge residual signs.

This slice does not alter solver runtimes, VCA adapters, CLIs, schemas, or
`.ntiga` handling.

### PR 0.3 generic 3D port measurement

- Added `TransientFlowRuntime::MeasurePorts` overload accepting explicit
  `CouplingPort` values and a physical time. It returns `PortState` values with
  area, outward flow, mean pressure, and optional outward species flux.
- The generic 3D adapter accepts only strict nonnegative `boundary_label`
  locators. Unsupported locator kinds, malformed labels, duplicate boundary
  locators, missing ports, and invalid port declarations fail before assembly.
- The implementation retains the existing owned-element traversal, ghost-state
  access, quadrature calls, and MPI reductions. It applies the declared local
  orientation only when converting the resulting flow and species flux to the
  common outward-positive contract.
- The existing VCA `MeasurePorts(ThreeDVascularPortDefinition, ...)` overload
  now creates `boundary_label` ports with the native outward orientation and
  adapts generic states back to the legacy maps. It does not duplicate surface
  integration or MPI reduction logic.
- Extended the focused PETSc VCA runtime test with a unit-cube flow state. It
  checks area, flow, mean pressure, species flux, physical time, strict locator
  validation, orientation conversion, and parity between generic and legacy VCA
  measurements.

No CLI, schema, database, 1D, CUDA, time-integration, history, or rollback
behavior changed in this slice.

### PR 0.4 rollback-safe 3D flow lifecycle

- Added explicit `BeginStep`, `SetTrialBoundaryConfiguration`, `SetPortInput`,
  `SolveTrial`, `RollbackTrial`, `CommitStep`, and `GetPortState` operations to
  `TransientFlowRuntime`, with clear rejection of invalid transitions.
- `BeginStep` takes an in-memory snapshot of the committed PETSc flow vector,
  resolved boundary values, pressure tractions, and full outlet-model state.
  Every trial restores the same flow vector and outlet snapshot and uses that
  vector as the backward-Euler history at `t_n`.
- Trial linear-iteration counts are provisional. Rollback discards them and
  commit publishes them to the physical cumulative total exactly once.
- The CLI now performs the explicit lifecycle before transport, VCA circuit
  advancement, histories, output, and checkpoints. The retained `Advance`
  method is a compatibility adapter using the same lifecycle and numerical
  implementation.
- The current 3D port input supports only mean static pressure or outward mean
  normal traction on a `boundary_label`, using the existing pressure-traction
  weak form. It rejects scalar prescribed flow because no scientifically valid
  velocity-profile reconstruction exists, and rejects conflicts with active
  resistance/RC/RCR outlet models.
- Extended the focused PETSc runtime test with exact trial--rollback--trial
  flow-vector and iteration-count equality, boundary-traction restoration,
  RC capacitor rollback/replay, single commit, invalid transitions, port-state
  timing, unsupported input rejection, and `Advance` adapter parity.
- The constructor-time velocity/pressure constraint masks are invariant:
  configured trial boundaries that add or remove a constrained row are rejected
  before assignment because PETSc `boundary_rows_` is not rebuilt mid-step.

This slice does not add 3D transport rollback or move the VCA circuit into the
flow lifecycle. Those objects remain CLI-owned and advance only after a
successful flow commit.

### PR 0.5 reusable rollback-safe native 1D runtime

- Added header-only C++17 `OneDFlowRuntime`, which owns configuration, network,
  `OneDFlowState` (including time, completed/internal substeps, and outlet/RCR
  state), and every transport/inlet mutable state.
- The runtime exposes `BeginStep`, `SetPortInput`, `SetCoupledInlet`,
  `SetOpenLoopInlet`, `SolveTrial`, `GetPortState`, `RollbackTrial`, and
  `CommitStep`, with legal-transition checks and exactly-once commit.
- Every trial snapshots and restores configuration/perfusate mutations,
  dynamic radius/area/resistance, flow/outlet state, and transport state. It
  calls the existing rigid, explicit, or injected PETSc implicit advance then
  the existing transport and vasodilation routines; no solver was copied.
- The CLI now drives this lifecycle. Writers, checkpoints, coupling history,
  prior mass, replay provider, and external circuit remain outside the runtime
  and execute only after commit. Checkpoint format and schema-v3 handling are
  unchanged.
- Generic 1D measurements use root and `outlet:<node-id>` ports. Flow and
  species flux are outward-positive SI: root uses explicit `-1` orientation,
  distal leaves use `+1`. Trial boundaries support root flow with optional
  concentration and static pressure at a terminal already configured with a
  pressure closure. Resistance/RCR replacement and unsupported traction inputs
  are rejected.
- Added `one_d_runtime_test` to the dependency-free `core-test` target. It
  checks invalid transitions, root/distal signs, exact rollback/re-solve,
  transport/configuration/dynamic-network/last-inlet rollback, time/step state,
  safe terminal pressure input, resistance-closure rejection, failed-solve
  rollback, and single commit.
- It also covers a true RCR capacitor trial mutation, rollback restoration,
  deterministic replay, exactly-once commit, and terminal-pressure rejection.

## Baseline evidence

On 2026-08-31, before Phase 0 implementation:

| Command | Result | Coverage |
|---|---|---|
| `make cpu-test` | pass | dependency-free CPU configuration, boundary, checkpoint, outlet, VCA, and visualization unit tests |
| `make -C solvers/one_d core-test` | pass | native 1D flow/transport and coupling unit tests |
| `make cpu-test` after PR 0.2 | pass | baseline CPU tests plus dependency-free common coupling-port tests |
| `make -C solvers/one_d core-test` after PR 0.2 | pass | unchanged native 1D fast test target |
| `make -C solvers/cpu vca_3d_runtime_test PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real` | pass | focused generic 3D-port and VCA-adapter PETSc test build |
| `./solvers/cpu/vca_3d_runtime_test` | pass | generic measurement values, locator rejection, and VCA parity |
| `make -C solvers/cpu vca_3d_smoke_test PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real` | pass | VCA flow/transport executable and smoke-test build |
| `./solvers/cpu/vca_3d_smoke_test` | pass | one-/two-rank VCA flow, transport, checkpoint, and restart regression |
| `make -C solvers/cpu vca_3d_runtime_test PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real` after PR 0.4 | pass | warning-clean focused lifecycle/port/outlet test build |
| `./solvers/cpu/vca_3d_runtime_test` after PR 0.4 | pass | deterministic flow replay, boundary/outlet rollback, transition guards, single commit, and adapter parity |
| `make -C solvers/cpu petsc-test PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real` after PR 0.4 | pass | PETSc kernels, two-rank boundary preflight, lifecycle runtime, and one-/two-rank VCA checkpoint/restart smoke |
| `make cpu-test && make -C solvers/one_d core-test` after PR 0.4 | pass | unchanged dependency-free CPU and native 1D fast baselines |
| `make -C solvers/one_d core-test` after PR 0.5 | pass | core, coupling, and dependency-free runtime lifecycle tests |
| `make one-d-test PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real` after PR 0.5 | pass | one- and two-rank PETSc native 1D baselines plus runtime lifecycle test |
| four documented `iga_1d --check` cases after PR 0.5 | pass | rigid, compliant, physiology, and VCA schema-v3 configuration compatibility |
| direct three-step rigid run versus one-step checkpoint plus restart to step three | pass | final `profile_1d_000003.vtp` byte-identical, SHA-256 `3cbb0781eaefce6df3c39e9437401871797bd4a0729d7af153280fb6426b7472` |
| two-step `vca_pfc_closed_loop` run after PR 0.5 | pass | external circuit path remains executable with advancement after runtime commit |
| base `9aa28d3` versus Phase 0 implementation `dc05f8a` configured transient 3D VCA fixture | pass | velocity, pressure, flow checkpoint, transport checkpoint, reservoir metadata, and coupling manifest are byte-identical |
| base `9aa28d3` versus Phase 0 implementation `dc05f8a` three-step compliant PETSc/multispecies 1D case | pass | time-series/profile outputs and final VTP are byte-identical |

Base-versus-implementation hashes for the final 3D velocity and pressure fields are
`9c6394c7b66cf2af9105986fdff21c607cf9ded795679272c9074d1c972f2f51`
and `407c77b51d739e08d9870c48fa1b82bcbd06c01e47150d7114e4efb59469e46e`.
The 3D flow checkpoint hash is
`ee488d0604880e36077bdbeb3c05ff1a16182eabcee18c24eac56a91e57cacbd`;
the coupling manifest hash is
`e8bf5f76dcff8500bce7a4ffd5843442924e2316ff025eb6b2afc6f48464f513`.
The final 1D VTP hash is
`55cf6a8efbab2d9a2d908f7265f4eba63824347543824c74a6bfebe7ea526d01`.

PETSc was discoverable locally through `pkg-config` as version 3.15.5. PETSc
runtime, multi-rank VCA, and numerical parity gates passed for Phase 0 and
remain required for later implementation slices that affect them.

### Reproduction commands

The final warning-clean and MPI/PETSc gate commands were:

```bash
make cpu-test
make -C solvers/one_d core-test
make -C solvers/cpu petsc-test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
make one-d-test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
```

The base/current executables were built from separate source trees at
`9aa28d3` and `dc05f8a`, then run against identical copied case directories.
The byte comparisons and recorded hashes used these commands, where
`BASE_3D`, `PHASE0_3D`, `BASE_1D`, and `PHASE0_1D` name the four output
directories:

```bash
cmp "$BASE_3D/flow.txt" "$PHASE0_3D/flow.txt"
cmp "$BASE_3D/flow.txt.pressure" "$PHASE0_3D/flow.txt.pressure"
cmp "$BASE_3D/checkpoint.state" "$PHASE0_3D/checkpoint.state"
cmp "$BASE_3D/checkpoint.vca_transport.state" \
  "$PHASE0_3D/checkpoint.vca_transport.state"
cmp "$BASE_3D/checkpoint.vca.json" "$PHASE0_3D/checkpoint.vca.json"
cmp "$BASE_3D/coupling_manifest.json" \
  "$PHASE0_3D/coupling_manifest.json"

for artifact in branch_flow_1d.csv derived_profile_1d.csv flow_1d.csv \
  profile_1d.csv species_profile_1d.csv profile_1d.pvd \
  physiology_manifest.json profile_1d_000003.vtp; do
  cmp "$BASE_1D/$artifact" "$PHASE0_1D/$artifact"
done

sha256sum "$PHASE0_3D/flow.txt" "$PHASE0_3D/flow.txt.pressure" \
  "$PHASE0_3D/checkpoint.state" "$PHASE0_3D/coupling_manifest.json" \
  "$PHASE0_1D/profile_1d_000003.vtp"
```

## Design decisions

- Coupled port flow and species flux are positive outward from each subsystem.
- Coupling edges, not label order or inlet/outlet names, perform sign changes.
- Coupling values use SI for time, area, velocity, flow, and pressure.
- A missing boundary quantity is distinct from a supplied zero value.
- Trial rollback is an in-memory physical-state snapshot; disk checkpoints
  remain persistence/restart artifacts.
- Writers, checkpoints, histories, and external circuits advance only after a
  successful commit.
- Schema v5 and `.ntiga` changes are deferred because Phase 0 needs neither.

## Backward compatibility

The Phase 0 slices preserve existing schema-v3 1D cases, schema-v4 3D cases,
older accepted configurations, `.ntiga` databases, standalone executables, VCA
workflows, and CUDA sources. PR 0.3 changes only the internal implementation
behind the existing VCA port-measurement overload; the focused adapter parity
test and VCA smoke regression pass.

## Known limitations and risks

- The trial API accepts only pressure/normal-traction port inputs; scalar flow
  input remains unsupported until a scientifically explicit profile model is
  designed.
- 3D transport and VCA circuit state have no rollback lifecycle yet.
- Native 1D accepts root flow-controlled input and static pressure at terminals
  already using a pressure closure; it does not silently replace R/RCR data.
- Runtime snapshots do not include the external VCA circuit, writer/history
  objects, or disk checkpoint persistence; the CLI advances those after commit.
- The current 1D checkpoint does not persist `last_inlet_`; immediately after
  restart, root species port measurement reflects the initialized inlet until
  the first post-restart solve. Flow state is unaffected. Coupled-transport
  restart must extend the persistence contract before that later phase.
- Full numerical regression cases can require PETSc, MPI, generated databases,
  or an allocated compute resource.

## Phase-gate result

| Gate | Result | Evidence |
|---|---|---|
| A: lifecycle correctness | pass | Both runtimes restore identical committed physical state before every trial, replay deterministically after rollback, publish provisional counters only on commit, and reject duplicate or out-of-order transitions. The 1D test includes real RCR capacitor mutation and restoration. |
| B: numerical and backward compatibility | pass | Dependency-free, PETSc, and MPI suites pass. Base `9aa28d3` and implementation `dc05f8a` produce byte-identical 3D VCA and compliant/multispecies 1D artifacts. Existing schemas and file formats are unchanged. |
| C: coupling contract readiness | pass | Ports use explicit SI quantities and outward-positive orientation, constraints and unsupported inputs fail early, 3D PETSc boundary-row topology cannot change during a trial, and subsystem side effects occur only after commit. |

The final Sol-level review found no unresolved lifecycle, SI-unit, sign,
topology, or compatibility defect. Extending rollback to 3D transport is
deliberately deferred until transport participates in rejected coupled trials;
it is not needed for the flow-only Phase 1 prototype.

## Phase 0 exit condition

The condition is satisfied at `dc05f8a`: both runtimes solve from identical
committed state after rollback, commit exactly once, reject invalid lifecycle
transitions, preserve standalone and VCA results, and expose port measurements
using the documented SI/outward-positive contract. Phase 1 PR 1.1 may begin
from this revision. Its first deliverable is the flow-only explicit-staggered
straight-vessel path `1D -> body-fitted 3D -> 1D`; strong coupling, schema v5,
and coupled transport remain out of scope until that path passes its own gate.
