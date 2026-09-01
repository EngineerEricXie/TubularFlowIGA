# Phase 2 report

Status: **in progress**. PR 2.1 provides validated in-memory topology, the
first PR 2.3 slices provide atomic domain transactions, backend adapters, exact
runtime binding, and a runnable sequential 1D--3D--1D graph executor. PR 2.2
defines the schema-v5 graph manifest, and that manifest is now bound to the
production coupling driver. Branching execution and multiple 3D islands remain
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

- The graph remains immutable topology metadata. Runtime ownership and
  transaction state live in the separate PR 2.3 registry and executor.
- Schema v5 is runnable through the production coupling driver. Existing v3/v4
  standalone dispatch and `.ntiga` files are unchanged.
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

The production 1D--3D--1D driver now binds its three native runtimes to the
graph and advances its two interfaces through this executor. Its CSV,
manifest, work counters, lagged-pressure behavior, and injected precommit
failure test remain unchanged. The established strong-mode loop remains in
place for its detailed iteration diagnostics, while schema-v5 fixed and Aitken
controls route to that same verified path. Per-attempt backend work remains
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

## PR 2.2: schema-v5 multidomain manifest

`MultidomainConfig.hpp` parses schema v5 independently from the standalone
schema-v3 and schema-v4 parsers and constructs the same validated
`SimulationGraph` used by the runtime registry and executor. `iga_config_check`
dispatches v5 to this parser while retaining the existing direct v3/v4 paths.
The `.ntiga` database format is unchanged.

The root record contains `time`, an explicit `start_domain`, one graph-wide
`execution` policy, `domains`, and `couplings`. Each domain declares its stable
ID, dimension/kind, case directory, backend-specific assets, and complete
logical-port metadata: physical locator, optional native-to-outward sign, and
provided/required quantities. A 1D domain declares whether its root uses the
configured open-loop waveform or receives coupled flow. A body-fitted 3D
domain declares its database. All asset paths are relative to the graph-case
root, lexical `..` traversal is rejected, and the runner canonicalizes existing
case directories and database files against the canonical root so symlinks
cannot escape it. It also canonicalizes each native configuration, geometry,
mesh, initial-velocity, and periodic-table asset inside its declared case
directory. Graph preflight failures are reduced collectively before solver
construction so a rank-local filesystem failure cannot strand peer ranks.

Each coupling names two logical endpoints, the `pressure_flow` law, and a
finite initial interface pressure. Fixed and Aitken controls live at the
component execution level because all interface pressures are iterated as one
vector; putting an independent algorithm on each edge would conflict with the
validated executor semantics. Strict key validation rejects misspellings,
backend-inconsistent fields, invalid controls, duplicate IDs, reused ports,
unknown endpoints, unattached required inputs, multiple scalar inputs on one
3D flow port, and pressure/flow capability mismatches.

Schema v5 requires a single connected multidomain graph because its execution
policy applies to one component. Branches remain schema-valid for PR 2.4.
`MakeSequentialPressureFlowPlan` separately checks the narrower executor
contract: one heterogeneous, acyclic, nonbranching chain whose pressure
receivers precede its flow receivers. `iga_config_check` reports that runner
compatibility independently of schema validity.

## Schema-v5 production runner

`iga_1d_3d_explicit --graph-case ROOT --output-dir DIR` reads
`ROOT/simulation_config.json`, resolves the current sequential semantic roles
without relying on domain or port names, canonicalizes the declared assets,
and binds the native runtimes to the exact graph metadata. The schema supplies
the start domain, interface IDs and initial pressures, time grid, and explicit,
fixed, or Aitken execution controls. Backend preflight confirms that declared
1D locators, 3D boundary labels, time grid, inlet policy, and databases match
the loaded native cases.

The legacy positional invocation remains available and constructs the same
runtime graph. Successful graph runs additionally write
`graph_binding_manifest.json` with the canonical graph root, resolved case and
database paths, domain IDs, coupling IDs, initial pressures, and selected
execution kind. Coupling records retain both logical endpoints. Smoke coverage
runs graph-authored explicit, fixed, and Aitken
cases with deliberately non-legacy IDs and compares their complete physical
histories with the corresponding legacy invocations. The same gate retains
one- and two-rank parity, reversed native-port orientation and conservation,
subcycling, input rejection, and atomic failure rollback checks. The binding
manifest is published by same-filesystem rename after the other outputs and
serves as the graph-run completion marker. Graph output directories must not
preexist, preventing a failed rerun from leaving an older marker beside mixed
outputs.

## PR 2.4 foundation: acyclic branch execution

`MakeAcyclicPressureFlowPlan` directs every interface from its pressure
receiver/flow provider toward its flow receiver/pressure provider. It requires
one declared source, rejects cycles, disconnected or same-kind coupling, and
topologically orders ready domains by stable ID. Interfaces are ordered
lexically by edge ID, defining the component-wide fixed/Aitken residual vector
independently of manifest insertion or endpoint order.

`PressureFlowComponentExecutor` now consumes this general plan. Before each
domain solve it installs all pressure guesses owned by that domain; after the
single solve it transfers every provided outward flow with the conservative
peer sign. Thus a 3D junction with multiple pressure-controlled outlets is
solved once per component iteration, not once per branch. Reverse rollback,
prepare-all/finalize-all commit, and best-effort abort cover every domain in
the topological component.

Dependency-free branch tests cover deterministic planning under domain, edge,
and endpoint permutation; explicit, fixed, and vector Aitken execution over
three interfaces; distinct sibling flow routing; exactly one junction solve
per iteration; missing initial pressure; nonfinite state; callback, sibling
solve, and prepare failures; and all-domain abort without partial commit. The
production runtime builder, Y-junction fixture, and long-form output remain the
next PR 2.4 slices.

`ResolveOneDThreeDBifurcation` now supplies the topology-to-runtime boundary
for that production work. It accepts exactly one configured-open-loop 1D
source, one body-fitted 3D junction with one flow-receiving inlet and at least
two pressure-receiving outlets, and coupled-root 1D leaves. It derives all
roles from the directed plan, inlet policies, port capabilities, and physical
locators, sorts branch records by edge ID, and requires the root, wall, and
terminal observations needed by the benchmark. Multiple 3D domains, deeper
trees, missing observations, and incorrect inlet policies remain rejected.
