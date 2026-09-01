# Phase 2 report

Status: **in progress**. PR 2.1 provides the validated in-memory domain and
coupling-edge topology for the existing sequential 1D--3D--1D case. Schema v5,
generic runtime ownership, branching execution, and multiple 3D islands remain
open.

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

- The graph is immutable topology metadata. Runtime-owning adapters wait for
  PR 2.3, where 1D/3D lifecycle mismatches and graph transaction atomicity can
  be handled deliberately.
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
multi-interface executor remain the next PR 2.3 slice.

## Verification

The dependency-free graph suite covers a valid chain in both directions,
insertion-order-independent planning, logical-port overlap, query behavior,
invalid identifiers/references/ownership, duplicate physical locators,
capability mismatches, endpoint reuse, self-edges, branches, disconnected
components, and same-kind sequential plans.

Validation commands for this snapshot:

```text
make coupling-test
make cpu-test
make coupling PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
make coupling-petsc-test PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
git diff --check
```

The long steady-spatial and pulsatile-temporal ladders need not be rerun for
this metadata-only refactor because the numerical loop, schemas, CSV fields,
manifest fields, boundary materialization, and runtime lifecycle are unchanged.
The final smoke run passed explicit, strong-fixed, and strong-Aitken modes with
`N=4` subcycling on one and two MPI ranks, including invalid-input rejection,
injected precommit failure suppression, work accounting, and rank-parity
checks.

## Entry condition for PR 2.2

PR 2.2 may define schema-v5 domain and edge records only after PR 2.1 passes
the dependency-free and one-/two-rank MPI smoke gates. The parser must produce
the same validated `SimulationGraph` metadata and must not force existing v3/v4
standalone cases through graph dispatch. A graph-authored case is not considered
runnable until PR 2.3 supplies backend adapters and a multi-interface executor.
