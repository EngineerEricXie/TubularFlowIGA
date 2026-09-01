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

## Baseline evidence

On 2026-08-31, before Phase 0 implementation:

| Command | Result | Coverage |
|---|---|---|
| `make cpu-test` | pass | dependency-free CPU configuration, boundary, checkpoint, outlet, VCA, and visualization unit tests |
| `make -C solvers/one_d core-test` | pass | native 1D flow/transport and coupling unit tests |
| `make cpu-test` after PR 0.2 | pass | baseline CPU tests plus dependency-free common coupling-port tests |
| `make -C solvers/one_d core-test` after PR 0.2 | pass | unchanged native 1D fast test target |

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

The baseline slice is documentation-only. Existing schema-v3 1D cases,
schema-v4 3D cases, older accepted configurations, `.ntiga` databases,
standalone executables, VCA workflows, and CUDA sources are unchanged.

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

1. Add common coupling data types and unit tests.
2. Generalize 3D port measurements while preserving VCA values.
3. Add 3D trial/rollback/commit and deterministic replay tests.
4. Extract a reusable 1D runtime with equivalent standalone behavior.
5. Include transport, dynamic-radius, circuit, time, and output state in the
   correct lifecycle boundary.
6. Run engineering, numerical, and architecture gates.

## Phase 0 exit condition

Phase 1 may begin only when both runtimes solve from identical committed state
after rollback, commit exactly once, reject invalid lifecycle transitions,
preserve standalone and VCA results, and expose port measurements using the
documented SI/outward-positive contract.
