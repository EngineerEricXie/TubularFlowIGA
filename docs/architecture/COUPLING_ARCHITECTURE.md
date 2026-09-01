# Coupling Architecture

Status: Phase 1 complete; Phase 2 in progress. The verified 1D--3D--1D
lifecycle and coupling algorithms now consume the PR 2.1 in-memory multidomain
topology. PR 2.3 provides runtime-owned sequential graph execution, and PR 2.2
provides the strict schema-v5 graph manifest and its production runner binding.
Branch execution and multiple 3D islands remain subsequent Phase 2 work.

This document fixes the runtime and interface contracts that precede direct
1D--3D coupling. The implementation remains incremental: existing standalone
executables continue to own configuration and output while reusable subsystem
objects acquire explicit trial, rollback, and commit semantics.

## Current runtime boundaries

### Three-dimensional flow

`solvers/cpu/src/iga_navier_stokes.cpp` owns configuration materialization,
physical time and step counting, VCA circuit state, transport orchestration,
checkpoints, and output. `TransientFlowRuntime` owns the PETSc flow objects,
resolved boundary data, outlet-model state, assembly, nonlinear solves, and
boundary measurements.

The runtime now owns an explicit committed/trial lifecycle. `BeginStep`
snapshots `state_`, resolved boundary data, pressure tractions, and every outlet
model value, and sets backward-Euler `previous_` from the committed flow vector
exactly once for the macro-step. `SolveTrial` restores the committed vector and
outlet state before every attempt. Trial linear iterations remain provisional;
`CommitStep` adds them to the physical total exactly once, while
`RollbackTrial` discards them. The legacy `Advance` entry point is a
compatibility adapter over begin, configure, solve, and commit.

The PETSc boundary-row topology is fixed at construction. Every configured
boundary resolution must retain the constructor velocity and pressure
constraint masks; a trial that adds or removes either constraint is rejected
before its values are assigned, rather than rebuilding `boundary_rows_` during
a lifecycle step.

`MeasurePorts` already performs distributed outward-normal flow integration,
area-weighted pressure integration, and advective species-flux integration.
It is coupled to `ThreeDVascularPortDefinition`, exposes only VCA outlet
labels, and discards the integrated area after forming mean pressure. The
integration implementation is the basis of the generic 3D port adapter; it
must not be duplicated.

### One-dimensional flow and transport

`solvers/one_d/include/OneDRuntime.hpp` owns configuration, network, flow
state, transport state, dynamic-radius updates, physical time, completed step,
internal substeps, and outlet/RCR state. `iga_1d` is the thin CLI driver for
initialization, replay/circuit input selection, checkpoint persistence,
history, and output. Existing numerical routines remain the implementation;
the runtime does not duplicate a solver. Its implicit advance is injected by
the PETSc CLI, preserving a dependency-free common lifecycle header/test.

The native 1D inlet is flow-controlled. Outlet leaves use pressure,
resistance, or RCR closures. Direct 1D--3D coupling will eventually require
explicit per-port boundary overrides, but Phase 0 must not introduce a new
case schema or silently reinterpret existing boundaries.

### Existing VCA path

VCA is an explicit vascular-to-circuit bridge, not direct 1D--3D coupling.
The common implementation should reuse its SI validation and measurement
logic while keeping the generic port layer independent of reservoir-specific
types and the fixed inlet/outlet layout.

## Phase 1 explicit straight-chain prototype

`iga_1d_3d_explicit` is a deliberately narrow flow-only prototype for one
upstream native 1D terminal, one body-fitted 3D inlet/outlet pair, and one
downstream native 1D root. At each `t_n -> t_n+1` it solves upstream 1D with
the lagged 3D-inlet static pressure, scales the 3D reference inlet profile to
the negative measured upstream-terminal outward flow, solves 3D with the
lagged downstream-root static pressure, and supplies the negative measured 3D
outlet flow to the downstream root. It validates the complete provisional
history row before committing all three runtimes exactly once. CSV and
manifest output are written only after the complete run succeeds.

This is explicit staggered coupling: `iteration_count=1` and
`relaxation_factor=1`. The pressures exchanged across both interfaces remain
one-step lagged throughout; they are not an interface pressure solve and
should not be interpreted as continuity of instantaneous pressure. Constant
inlet flow can still exhibit the documented first-step lag. The prototype has
no restart/checkpoint option. Its explicit path now routes through the generic
sequential graph executor; the established strong modes retain their richer
driver-owned per-attempt diagnostics while the generic executor independently
supports fixed and Aitken iteration. A test-only environment variable,
`TUBULARFLOWIGA_INJECT_EXPLICIT_COUPLING_FAILURE_STEP=<positive-step>`, throws
after all trial solves and history validation but before a commit; it exists to
prove that rejected steps emit neither coupling history nor manifest output.

The same executable accepts a schema-v5 case through
`--graph-case ROOT --output-dir DIR`. This mode derives all domain, coupling,
and logical-port identities from the manifest, resolves declared assets inside
`ROOT`, and confines transitive native assets to their declared case
directories. MPI ranks collectively agree on graph preflight before solver
construction. The mode takes execution controls and initial interface
pressures from schema v5 and applies each port's native-to-outward orientation
to both measurements and reference-profile scaling. It publishes
`graph_binding_manifest.json` by same-filesystem rename as the final completion
marker, recording the canonical cases, database, graph IDs, and execution kind
used. A graph output directory must not preexist, so a failed rerun cannot
retain an older success marker. The positional Phase 1 CLI remains supported
for compatibility.

## Phase 2 bifurcation runner

`iga_1d_3d_bifurcation` is the schema-v5 production path for one native 1D
source, one body-fitted 3D junction, and two or more coupled-root 1D leaves.
Semantic roles come from pressure/flow capabilities and the acyclic component
plan, not domain names or manifest insertion order. Every pressure input is
installed before the junction's single trial solve; its outlet flows are then
routed independently to branch roots. All interface pressures form one
lexically edge-ordered fixed/Aitken vector.

The runtime registry owns one adapter per graph node, and the component uses
prepare-all/finalize-all commit with reverse abort. Each native 1D domain is
initialized from its own configured inlet waveform at time zero, preserving
independent branch state. Accepted data is serialized in long form so branch
count does not change the schema. A renamed `graph_binding_manifest.json` is
the final success marker; failure before commit or output completion cannot
publish it.

`iga_multidomain_flow` uses the same construction core without the
bifurcation topology restriction. Preflight parses every native domain,
validates rank-specific database records and all physical 3D boundary audit
ports, and collectively agrees on errors before constructing PETSc objects.
Instantiation follows deterministic component order and preserves ownership:
adapters are destroyed before native runtimes, and each 3D runtime is destroyed
before its database. Flow-receiving 3D ports each retain their own oriented
reference profile; graph-controlled pressure ports each override a native
pressure-traction boundary.

The generic runner supports one connected acyclic heterogeneous pressure-flow
component with one 1D source, including arbitrary alternating chain and branch
composition and multiple independently owned body-fitted 3D regions. Output
adds per-3D-domain boundary-flow balances to the edge, iteration, port, and
initialization records. Its final manifest records canonical assets and sorted
edge endpoints and is published only after every CSV closes successfully.
Cycles, disconnected components, same-kind edges, transport, moving/immersed
domains, FSI, 0D models, restart, per-edge relaxation, and active graph-driven
3D outlet models remain outside this phase.

`--coupling-mode strong-fixed` adds the equally narrow fixed-relaxation
alternative without changing the explicit default or its output names.  For a
single physical step all three runtimes enter `BeginStep` once.  Each trial
sweep uses `x=[upstream terminal pressure, 3D outlet pressure-traction parameter]` and
`G=[measured 3D inlet mean static pressure, measured downstream root mean
static pressure]`, checks conservative flow transfer, and then either rolls
back downstream, 3D, and upstream in that order or commits all three exactly
once.  The update is `x_next=x+omega*(G-x)` and convergence uses
`abs(G_i-x_i)/max(Pref,abs(G_i),abs(x_i))`.  In particular, the measured 3D
outlet pressure remains a diagnostic; it is not the right-side fixed-point
pressure because the applied 3D quantity is a pressure-traction parameter.  Rejected
attempts contribute to coordinator work diagnostics but never to runtime
committed counters.  A successful strong run writes separate
`strong_coupling_history.csv`, `strong_coupling_iterations.csv`, and
`strong_coupling_manifest.json`; restart remains unsupported.

PR 1.3 adds `strong-aitken` for the same two-pressure vector only. It resets
at every macro step and uses the signed raw pressure residual `r=G-x` in the
documented vector order. A proposed relaxation is clamped to its configured
interval and accepted only after reverse rollback, immediately before the
next trial uses it. It does not add restart or graph coupling.

PR 1.4 permits each 1D domain to use an integer subdivision of the 3D macro
step. The 3D grid remains the macro grid; a 1D horizon must match it and its
configured step count must be `macro_steps*N`. Interface data are zero-order
held through the N configured 1D substeps, while a configured upstream
waveform is sampled at each endpoint. A macro trial snapshots once and rolls
all N substeps back together. Configured-substep counts and explicit CFL work
are separate diagnostics; `internal_substeps` remains a native CFL counter.

## State ownership and rollback inventory

The committed snapshot is the complete physical state at `t_n`. Solver
matrices, Krylov work vectors, ghost scatters, topology, and immutable
configuration are rebuildable or reusable scratch and are not physical state.

| Owner | Physical or trial state | Required handling |
|---|---|---|
| 3D flow runtime | PETSc `state_` | Snapshot at `BeginStep`; restore on rollback; accept once on commit. |
| 3D flow runtime | backward-Euler `previous_` | Set from the committed `state_` once per macro-step; never from a rejected trial. |
| 3D flow runtime | `outlet_models_`, including flow, pressure, and capacitor pressure | Snapshot and restore with the flow vector. |
| 3D flow runtime | `boundaries_` and `pressure_tractions_` | Treat materialized trial inputs as trial state; restore the committed values. |
| 3D flow runtime | linear/nonlinear diagnostic counters | Trial diagnostics remain attributable to the attempt; physical cumulative counters update only on commit. |
| 3D transport runtime | PETSc `current_` and logical `steps_` | `BeginStep` snapshots both; rollback/abort restore both; prepare/finalize publish once. Do not advance on rejected flow trials. |
| 3D CLI/VCA | `VcaExternalCircuit` state, last arterial species, and previous species mass | Advance exactly once after the coupled step converges. |
| 1D flow runtime | area, flow, pressure, nodal pressure, segment flow, inlet flow, outlet/RCR states | Copy as one committed `OneDFlowState` snapshot. |
| 1D flow runtime | physical time, completed step, internal substep count | Advance in the trial copy; publish only on commit. |
| 1D network | dynamic `radius0`, `area0`, and resistance from vasodilation | Snapshot when dynamic geometry is enabled. |
| 1D transport | every species concentration and mutable inlet value/waveform | Snapshot with the trial state. |
| 1D configuration | hematocrit/hemoglobin fields changed by `ApplyOneDCoupledInlet` | `OneDFlowRuntime` snapshots its complete configuration with transport and network state. |
| CLI output | writers, coupling history, prior mass, checkpoints | Side effects occur after `CommitStep`, never during `SolveTrial`. |

Disk checkpoint/restart remains a persistence interface. It is not the
in-memory rollback mechanism.

The 1D internal-substep count is therefore a rollback-owned physical counter:
rejected strong sweeps restore it with the complete 1D trial state. The current
coordinator reports accepted/rejected 3D KSP work only; rejected 1D
computational work is intentionally not reported.

## Common port contract

The first common types will be dependency-free C++17 values under `include/`.
They will not depend on PETSc, MPI, JSON parsing, VCA, or either solver.

Conceptually:

```cpp
enum class PortQuantity {
	Area,
	FlowRate,
	MeanPressure,
	MeanNormalTraction,
	TotalPressure,
	SpeciesConcentration,
	SpeciesFlux
};

struct PortOrientation {
	// Converts the dimension-specific positive direction to positive outward.
	int native_to_outward_sign = 1; // must be +1 or -1
};

struct PortState {
	double time_s = 0.0;
	double area_m2 = 0.0;
	double outward_flow_m3_s = 0.0;
	double mean_pressure_pa = 0.0;
	double mean_normal_traction_pa = 0.0;
	double total_pressure_pa = 0.0;
	std::map<std::string, double> concentration;
	std::map<std::string, double> outward_species_flux;
};

struct PortBoundaryData {
	double time_s = 0.0;
	std::optional<double> outward_flow_m3_s;
	std::optional<double> mean_pressure_pa;
	std::optional<double> mean_normal_traction_pa;
	std::optional<double> total_pressure_pa;
	std::map<std::string, double> concentration;
	std::map<std::string, double> outward_species_flux;
};

struct CouplingPort {
	std::string id;
	std::string subsystem_id;
	std::string locator_kind;
	std::string locator;
	PortOrientation orientation;
	std::set<PortQuantity> provides;
	std::set<PortQuantity> requires;
};
```

Field presence in `PortBoundaryData` is explicit. A numeric zero is valid data
and must not mean "not provided." Species concentration units are defined by
each configured species; species flux uses concentration unit times m^3/s.
The first implementation validates finite values, positive area when area is
provided, a sign of exactly `+1` or `-1`, and unique non-empty identifiers.
Beginning with PR 2.1, `provides` means a quantity readable from a logical
physical port and `requires` means a boundary quantity accepted at that port.
The sets may overlap: a controlled 3D inlet, for example, both accepts a target
flow and reports its realized flow for conservation. Edge-law validation, not
an artificial split into input and measurement port IDs, determines whether
the bidirectional exchange is well posed.

## SI units and sign convention

All coupling interfaces use:

| Quantity | Unit |
|---|---|
| time | s |
| area | m^2 |
| flow rate | m^3/s |
| velocity | m/s |
| pressure and normal traction | Pa |
| species concentration | the species' explicitly configured concentration unit |
| species flux | concentration unit m^3/s |

`PortState::outward_flow_m3_s > 0` means material leaves the subsystem through
the port. In 3D this is `integral(u dot n dA)` with the geometric outward
normal. In a rooted 1D vessel, root-to-leaf flow is outward at a distal port
but inward at the root port. `PortOrientation::native_to_outward_sign`
performs only this local conversion.

A conservative edge joins ports A and B using

```text
Q_A,out + Q_B,out = 0
Phi_A,out(species) + Phi_B,out(species) = 0.
```

### Schema-v6 species metadata

Schema v6 adds the metadata foundation for conservative species coupling while
schema v5 remains the executable flow-only graph contract. A graph-level
registry assigns every logical species a stable ID and concentration unit;
its flux unit is derived explicitly as that concentration unit times `m^3/s`.
Each domain maps the logical ID to its native scalar field, so coupling code
never gives special meaning to names such as `oxygen`. Each participating port
and pressure-flow edge carries an explicit sorted logical-species set.

Species-capable edge endpoints must both provide and accept concentration and
total outward flux. This symmetric capability is required because physical
flow reversal swaps donor and receiver; hydraulic input roles remain separate
and retain the pressure/flow formulation. Routing uses the measured outward
flows, not endpoint ordering or the topological flow plan. Opposite-sign,
conservative flows select the positive-outward endpoint as donor. Near zero,
the last committed donor owns the tie; without committed ownership the state
is ambiguous and rejected. Same-sign or materially nonconservative pairs are
rejected.

PR 3.1 is intentionally metadata and dependency-free routing only. Production
runners reject schema v6 until transactional 1D/3D transport, total-flux
measurement, and time-integrated balance accounting are implemented.

The edge applies the sign reversal explicitly. Boundary labels, graph order,
and inlet/outlet names never imply a sign.

## Subsystem lifecycle

The reusable interface is deliberately small:

```cpp
class CoupledSubsystem {
public:
	virtual ~CoupledSubsystem() = default;
	virtual void Initialize() = 0;
	virtual void BeginStep(double time_s, double dt_s) = 0;
	virtual void SetPortInput(const std::string&, const PortBoundaryData&) = 0;
	virtual void SolveTrial() = 0;
	virtual PortState GetPortState(const std::string&) const = 0;
	virtual void RollbackTrial() = 0;
	virtual void CommitStep() = 0;
};
```

The implemented 3D legal sequence is:

```text
committed -> BeginStep -> trial-ready
trial-ready -> SetPortInput -> SolveTrial -> trial-solved
trial-solved -> RollbackTrial -> trial-ready
trial-solved -> CommitStep -> committed
```

Invalid transitions fail clearly. `SolveTrial` always starts from the same
committed `t_n` snapshot, including internal 1D subcycling and terminal-model
state. `CommitStep` succeeds exactly once. A failed solve leaves a rollback
path and cannot emit output or advance external circuit state.

`OneDFlowRuntime` implements the same state machine after either
`InitializeOpenLoop` or `InitializeCoupled`. `BeginStep(time_s, dt_s)` requires
the committed time; trial input is supplied through `SetPortInput` (common
ports), `SetCoupledInlet` (existing VCA/replay state), or `SetOpenLoopInlet`.
`SolveTrial` restores its complete in-memory image before every solve, then
uses the existing rigid/explicit/PETSc-implicit routine, transport advance, and
vasodilation update. `RollbackTrial` restores configuration, network, flow,
outlet/RCR, time/substep, and transport state. Output, checkpoints, histories,
and external VCA circuit advancement remain CLI work after `CommitStep`.
The focused runtime test includes an RCR capacitor mutation, rollback, exact
replay, and single commit.

The native 1D generic ports are `root` and `outlet:<node-id>`. They measure
area, mean pressure, concentration, outward flow, and outward species flux in
SI units. Root native root-to-leaf flow is converted with orientation `-1`;
terminal native flow uses `+1`. A root flow-controlled trial boundary accepts
optional species concentration. A distal port may override mean static pressure
only when its configured closure is already a pressure boundary; resistance
and RCR closures cannot be replaced. Traction, total-pressure, and imposed
species-flux inputs remain unsupported.

`TransientFlowRuntime::BeginStep` additionally receives the existing CLI step
index and nonlinear tolerances because those controls belong to the 3D solve,
not to a generic port value. Existing materialized 3D boundary conditions are
applied through `SetTrialBoundaryConfiguration`. The current overload
`SetPortInput(const CouplingPort&, const PortBoundaryData&)` accepts only a
`boundary_label` port requiring either mean static pressure or outward mean
normal traction. Static pressure is applied through the existing pressure-
traction weak form; outward normal traction uses `pressure = -traction`.
Prescribed flow, total pressure, species values, and flow profiles are rejected:
the current solver has no scientifically defined way to derive a velocity
profile from one scalar flow value. A trial pressure override also cannot
replace an active resistance/RC/RCR outlet model. This narrow interface is
intentional until Phase 1 defines those models.

For a configured runtime, `SetTrialBoundaryConfiguration` must precede
`SetPortInput` so pressure-Dirichlet and outlet-model conflicts are checked
against the current step rather than stale committed boundaries. During an
active step, backward-Euler `previous_` is deliberately the committed `t_n`
flow vector; rollback preserves that value for exact trial replay. The older
pre-step `previous_` vector is solver scratch once `BeginStep` establishes the
new macro-step and is not a separately committed physical state.

For fixed-point coupling, the driver performs:

```text
BeginStep(A, B)
  -> apply trial interface data
  -> SolveTrial(A), measure A
  -> SolveTrial(B), measure B
  -> compute dimensionless residuals and relax
  -> RollbackTrial(A, B) and repeat, or CommitStep(A, B)
```

## Smallest implementation boundaries

The 3D refactor is flow-only. Additions to `TransientFlowRuntime` make snapshot
ownership explicit and remove the `step > 0` history copy from trial execution.
Generic measurements accept
an explicit set of boundary-label locators and return area as well as flow and
pressure. The VCA adapter calls that implementation so its output remains
unchanged. Transport rollback follows only when flow lifecycle parity passes.

The 1D runtime collects numerical orchestration without duplicating routines.
The CLI drives its lifecycle and keeps persistence/circuit/history side effects
post-commit. The focused dependency-free runtime test covers legal transitions,
root/outlet signs, coupled configuration mutation, transport, dynamic network
mutation, and exact rollback/re-solve equality. PETSc checkpoint and broader
VCA numerical regression remain the compatibility gate for changes beyond this
runtime extraction.

## Straight-vessel benchmark definition

Phase 1 will compare two physically equivalent cases:

```text
reference: upstream 1D segment -> replacement 1D segment -> downstream 1D segment
coupled:   upstream 1D domain  -> body-fitted 3D straight tube -> downstream 1D domain
```

The executable smoke uses geometrically identical circular sections. The
spatial-verification case uses a unit-square body-fitted 3D duct because it has
an independent analytic resistance and velocity series on a regular
tensor-product spline grid. Its all-1D replacement has the same area and a
hydraulic-equivalent circular length, so the reference and coupled paths have
the same exact steady resistance even though their cross-section shapes differ.
This distinction is deliberate and is recorded by the test rather than hidden
as an apparent geometric-equivalence claim.

The spatial gate uses open-uniform cubic C2 bases and shape-regular refinements.
At each transverse resolution it solves 3D lengths 1, 1.5, and 2 m with the
same axial spacing, verifies that consecutive half-length pressure slopes
converge to one length-independent bulk coefficient, and compares that
coefficient with the analytic square-duct value. This cancels finite cap
effects without assuming they are constant: the third length directly gates
that assumption. The fine cases additionally compare the absolute external
drop with the all-1D reference and bound both interface pressure jumps.

All cases use identical SI density, viscosity, terminal model, inlet waveform,
macro-step interval, and outward-normal port convention. The temporal benchmark
advances four pulse periods, checks the final two for repeatability, and
measures the final period on nested macro grids.
It records P/Q waveforms at both interfaces, pressure drop, first-harmonic
amplitude and phase, per-step mass imbalance, coupling residual histories,
iteration counts, and tolerance sensitivity. Fieldwise pressure and flow
self-convergence supplements the completed spatial gate. An independently
advanced rigid resistance--inertance network uses the same outer segments and
a unit-area, unit-length middle element with exact square-duct drag. It thereby
matches both steady resistance and fluid inertance of the 3D replacement while
remaining independent of the coupling driver. Both paths use one 1D step per
3D macro step for this comparison, so their backward-Euler inertance increments
use the same consecutive endpoint flows; subcycling remains covered by the
separate N=4 MPI gate. The fixed-area incompressible
model has pressure phase but no finite-speed pulse propagation; the benchmark
therefore gates the physically expected zero flow-transit delay rather than
claiming compliant-wave behavior.

## Compatibility matrix

| Interface | Phase 0 requirement |
|---|---|
| 1D schema v3 | Parse and execute unchanged. |
| 3D schema v2/v3/v4 accepted by current readers | Parse and execute unchanged. |
| `.ntiga` databases | No format or reader change. |
| standalone `iga_1d` | Same CLI, numerical path, checkpoints, and output. |
| standalone `iga_navier_stokes` | Same CLI, numerical path, checkpoints, and output. |
| VCA replay/closed loop | Same signs, SI values, histories, and reservoir advancement. |
| CUDA | No Phase 0 behavior change; common interfaces remain backend-neutral. |

Schema v5 now describes connected multidomain graphs and produces the same
validated topology consumed by the runtime registry. It does not alter the
standalone schema-v3/schema-v4 dispatch paths or the `.ntiga` format.

## Phase 2 PR 2.1 topology boundary

`SimulationGraph` is initially an immutable, dependency-free topology object.
It owns copied `DomainNode` and `CouplingEdge` metadata, not solver runtimes.
Each domain has its own `DomainKind`; there is no case-global dimension. A
typed `pressure_flow` edge joins two logical ports and requires one flow
receiver opposite one mean-pressure receiver. Both endpoints must report flow
and pressure, so outward-flow conservation and pressure diagnostics remain
available independently of which boundary value each side accepts.

The existing production driver now constructs the in-memory
`upstream 1D -- body-fitted 3D -- downstream 1D` graph after its established
backend preflight. It obtains the runtime port locators and unified 3D inlet and
outlet descriptors from that graph, then retains the verified Phase 1 solve,
rollback, commit, diagnostics, and output sequence unchanged. A sequential
plan requires an explicit start domain and accepts only one connected,
heterogeneous, acyclic chain. The core graph itself permits branches and
disconnected components so PR 2.4 and PR 2.5 do not require a topology rewrite;
their executors will add the corresponding execution plans.

PR 2.3 now owns each native backend behind a typed domain adapter and binds one
runtime per graph domain. The adapters share dependency-free metadata
validators with the v5 parser, so a manifest accepted by `iga_config_check`
cannot defer locator, orientation, or 1D inlet-policy errors until runtime
binding. The component executor performs reverse rollback and prepare-all /
nonthrowing-finalize-all commit across the full sequential chain. Coupling
algorithm state belongs to that connected-component integrator rather than
individual edges because all current interface pressures form one fixed-point
or Aitken residual vector.

Schema v5 requires one connected multidomain graph but permits branch topology
needed by PR 2.4. `MakeSequentialPressureFlowPlan` is the separate current
runner-compatibility check: it additionally requires a heterogeneous,
nonbranching, acyclic chain ordered from pressure receivers to flow receivers.
The configuration checker reports this compatibility without conflating a
future-valid branched manifest with a malformed schema.

For the current sequential runner, `ResolveSequentialOneDThreeDOneD` derives
the upstream 1D, body-fitted 3D, and downstream 1D roles from the validated
plan and port capabilities rather than reserved IDs. It requires the upstream
configured-open-loop inlet policy, the downstream coupled-root policy, and the
observation ports needed by the established history contract. Native binding
then confirms every declared locator and boundary label before the first step.

PR 2.4 adds a separate acyclic pressure-flow plan for branch execution. Edge
direction comes only from port requirements: a pressure receiver provides the
flow transferred to its peer, while that flow receiver provides the measured
pressure residual. The plan requires a unique declared source and uses
topological domain order with lexical tie-breaking; its component residual
vector is always ordered by edge ID. The executor installs every pressure input
for a domain, solves it once, then transfers all of its outgoing flows. Fixed
and Aitken relaxation therefore remain one component-wide operation across all
interfaces. Cycles, disconnected components, and same-kind edges remain
unsupported.

The PR 2.4 bifurcation resolver narrows this general executor to one supported
production topology: a configured-open-loop 1D source, one body-fitted 3D
junction, and at least two coupled-root 1D leaves. Branch identity and ordering
come from directed interfaces and stable edge IDs; domain and port names carry
no semantics. Required source-root, 3D-wall, and branch-terminal observation
ports keep output construction explicit rather than inferring native state.

## Validation manifest

Fast dependency-free baseline:

```bash
make cpu-test
make -C solvers/one_d core-test
```

PETSc baseline when the local PETSc environment is configured:

```bash
make cpu-petsc PETSC_DIR=/path/to/petsc PETSC_ARCH=arch
make -C solvers/cpu petsc-test PETSC_DIR=/path/to/petsc PETSC_ARCH=arch
make one-d-test PETSC_DIR=/path/to/petsc PETSC_ARCH=arch
```

Configuration checks after building `iga_1d`:

```bash
./solvers/one_d/iga_1d examples/one_d/rigid_straight --check
./solvers/one_d/iga_1d examples/one_d/compliant_bifurcation --check
./solvers/one_d/iga_1d examples/one_d/multispecies_physiology --check
./solvers/one_d/iga_1d examples/one_d/vca_pfc_closed_loop --check
```

The two-rank 3D VCA regression follows
`examples/vascular_flow/vca_bifurcation/README.md` on an allocated compute
resource. Large numerical cases do not run on a shared cluster login node.

## Planned Phase 0 touch list

PR 0.2 is limited to common value types and dependency-free tests:

- `include/CouplingPort.hpp` (new)
- `solvers/cpu/tests/test_coupling_port.cpp` (new shared test location)
- `solvers/cpu/Makefile`

Later 3D lifecycle/measurement work is expected to touch:

- `solvers/cpu/include/TransientFlowRuntime.hpp`
- `solvers/cpu/src/iga_navier_stokes.cpp`
- `include/ThreeDVcaCoupling.hpp`
- focused CPU/PETSc tests and Makefile dependencies

Later 1D runtime extraction is expected to touch:

- `solvers/one_d/include/OneDRuntime.hpp` (new)
- `solvers/one_d/src/iga_1d.cpp`
- `solvers/one_d/include/OneDCoupling.hpp`
- focused 1D tests and Makefile dependencies

Changes outside these boundaries require a documented reason and renewed
compatibility review.
