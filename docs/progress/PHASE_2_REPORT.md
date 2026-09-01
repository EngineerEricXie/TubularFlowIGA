# Phase 2 report

Status: **in progress**. PR 2.1 provides validated in-memory topology, and the
first PR 2.3 slices provide atomic domain transactions, backend adapters, exact
runtime binding, and a runnable sequential 1D--3D--1D graph executor. Schema
v5, branching execution, and multiple 3D islands remain open.

## Objective

Remove the assumption that a complete case has one global dimension or exactly
one hard-coded 1D/3D pair. Phase 2 will introduce heterogeneous domain nodes,
explicit coupling edges, multiple interfaces, a 3D bifurcation benchmark, and
multiple independent 3D islands while preserving standalone schema-v3 1D and
schema-v4 3D behavior.

## PR 2.1: immutable SimulationGraph topology

`include/CouplingEdge.hpp` defines typed domain kinds, port references, and the
initial `pressure_flow` law. `include/SimulationGraph.hpp` owns validated copies
of domain/port metadata and coupling edges. It contains no PETSc/MPI types,
runtime pointers, solver controls, trial state, parser, or output behavior.

Graph construction rejects empty or duplicate IDs, unknown endpoints, port
ownership and locator conflicts, endpoint reuse, self-edges, and ambiguous or
under-specified pressure/flow capabilities. A logical port may both accept and
report the same quantity. This is necessary for controlled boundaries such as
a 3D inlet that accepts a scalar flow target through the profile adapter and
also reports realized outward flow for conservation.

The graph core permits disconnected and branched topology. The first
`MakeSequentialPlan` is intentionally narrower: given an explicit start domain,
it produces a deterministic plan only for a connected, acyclic, nonbranching
heterogeneous chain. It rejects topology that needs a component or branch
executor instead of guessing an execution direction from edge insertion,
labels, or domain names.

The existing `iga_1d_3d_explicit` driver constructs this graph after all
backend-specific preflight checks and resolves its logical 1D terminal/root and
3D inlet/outlet ports through it. The established Phase 1 runtimes still own
initialization and the exact trial/rollback/commit transaction. Strong-fixed
and strong-Aitken still iterate both interface pressures as one vector; no
algorithm state was incorrectly moved onto individual edges.

## Design decisions and limitations

- The graph remains immutable topology metadata. Runtime ownership and
  transaction state live in the separate PR 2.3 registry and executor.
- Schema v5 waits for a runnable graph-authored case in PR 2.2/2.3. Existing
  v3/v4 dispatch and `.ntiga` files are unchanged.
- `DomainKind` currently lists native 1D flow and body-fitted 3D flow. Future
  0D and immersed 3D kinds will be added when their runtime contracts exist.
- The sequential plan rejects branches, cycles, disconnected components, and
  same-kind neighbors. The underlying graph representation accepts the
  topology required by later Phase 2 executors.

## PR 2.3 foundation: atomic domain lifecycle

The native 1D and body-fitted 3D runtimes now expose a graph-safe close path.
`AbortStep` restores the complete committed snapshot and closes an open,
solved, failed-solved, or commit-prepared step; it is idempotent in the idle
committed phase. `PrepareCommitStep` validates a successful trial without
publishing it, and `FinalizeCommitStep` performs only nonthrowing publication.
The existing `CommitStep` API remains as a compatibility wrapper.

This two-stage commit prevents a future multidomain executor from partially
committing a graph when a later domain rejects its trial. Native tests cover
open, solved, failed-solved, and prepared aborts, committed-state restoration,
trial-work reset, boundary-traction restoration, idempotence, and legacy
commit behavior. Runtime adapters, registry binding, and the sequential
multi-interface executor are supplied by the next PR 2.3 slice.

## PR 2.3: runtime adapters and multiple interfaces

`CoupledDomainRuntime` is a dependency-free lifecycle and port contract.
`DomainRuntimeRegistry` owns exactly one adapter for every graph domain and
rejects mismatched domain kinds or any difference in logical port metadata.
The 1D adapter makes the inlet policy explicit: an upstream domain uses its
configured open-loop waveform, while a downstream domain receives coupled
root flow. The body-fitted 3D adapter caches logical inputs and executes the
required order: materialize endpoint waveforms, scale reference flow profiles,
install the trial boundary configuration, apply pressure inputs, then solve.

`PressureFlowComponentExecutor` consumes the deterministic sequential plan and
supports any length alternating 1D/body-fitted-3D chain whose pressure
receivers precede its flow receivers. It transfers outward-positive flow with
the opposite sign at the peer, updates all interface pressures as one vector
with explicit, fixed, or dynamic Aitken iteration, checks every pressure and
flow residual, rolls retries back in reverse order, and prepares every domain
before any nonthrowing finalizer. Cleanup attempts every active domain and
reports both the primary and any abort failures. A precommit observation hook
lets output/history code inspect accepted trial states without taking ownership
of runtime state or weakening atomic commit. Every adapter state is validated
at the executor boundary before transfer or commit. A typed nonconvergence
error retains every iteration and the final pressure/flow residuals.

The production explicit 1D--3D--1D driver now binds its three native runtimes
to the graph and advances its two interfaces through this executor. Its CSV,
manifest, work counters, lagged-pressure behavior, and injected precommit
failure test remain unchanged. The existing strong-mode loop remains in place
for its established detailed iteration diagnostics; component-wide fixed and
Aitken behavior in the generic executor is covered independently before a
schema-v5 graph-authored case adopts it. Per-attempt backend work remains
adapter-specific, so the strong production loop will not migrate until the
generic diagnostic API can retain those counters as well.

## Verification

The dependency-free graph suite covers a valid chain in both directions,
insertion-order-independent planning, logical-port overlap, query behavior,
invalid identifiers/references/ownership, duplicate physical locators,
capability mismatches, endpoint reuse, self-edges, branches, disconnected
components, and same-kind sequential plans.

The executor suite covers explicit, fixed, and component-wide dynamic Aitken
execution over two interfaces, deterministic retry counts, complete port
observation, exact registry binding, direction rejection, solve and prepare
failures, and best-effort cleanup with combined error reporting. Native adapter
tests cover both 1D inlet policies and an actual PETSc body-fitted 3D reference
profile solve.

Validation commands for this snapshot:

```text
make coupling-test
make cpu-test
make coupling PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
make coupling-petsc-test PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
git diff --check
```

The long steady-spatial and pulsatile-temporal ladders need not be rerun because
the same native solvers, schemas, boundary values, and output definitions are
used. The final smoke run passed explicit, strong-fixed, and strong-Aitken
modes with
`N=4` subcycling on one and two MPI ranks, including invalid-input rejection,
injected precommit failure suppression, work accounting, and rank-parity
checks.

## Entry condition for PR 2.2

PR 2.2 may define schema-v5 domain and edge records only after PR 2.1 passes
the dependency-free and one-/two-rank MPI smoke gates. The parser must produce
the same validated `SimulationGraph` metadata and must not force existing v3/v4
standalone cases through graph dispatch. A graph-authored case is not considered
runnable until PR 2.3 supplies backend adapters and a multi-interface executor.
That runtime gate is now satisfied; schema-v5 parsing is the next active slice.
