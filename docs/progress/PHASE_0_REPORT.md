# Phase 0 Progress Report

Status: in progress

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

PETSc was discoverable locally through `pkg-config` as version 3.15.5. PETSc
runtime, multi-rank VCA, and numerical parity gates remain required when their
corresponding implementation slices begin.

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

- `TransientFlowRuntime::Advance` currently advances backward-Euler history
  based on its CLI step argument and is not trial-safe.
- 3D transport and VCA circuit state have no rollback lifecycle yet.
- The 1D CLI still owns runtime orchestration and mutable dynamic geometry.
- `ApplyOneDCoupledInlet` mutates configuration and transport inlet state.
- Native 1D boundary support does not yet expose generic terminal overrides.
- Full numerical regression cases can require PETSc, MPI, generated databases,
  or an allocated compute resource.

## Remaining Phase 0 work

1. Add 3D trial/rollback/commit and deterministic replay tests.
2. Extract a reusable 1D runtime with equivalent standalone behavior.
3. Include transport, dynamic-radius, circuit, time, and output state in the
   correct lifecycle boundary.
4. Run engineering, numerical, and architecture gates.

## Phase 0 exit condition

Phase 1 may begin only when both runtimes solve from identical committed state
after rollback, commit exactly once, reject invalid lifecycle transitions,
preserve standalone and VCA results, and expose port measurements using the
documented SI/outward-positive contract.
