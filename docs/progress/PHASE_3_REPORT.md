# Phase 3 report

Status: **in progress**. PR 3.1 provides the schema and dependency-free
species contract. Native conservative transfer remains open for PR 3.2, and
the two-species production gate remains open for PR 3.3.

## Objective

Conserve arbitrary transported quantities across 1D/3D interfaces using
outward-positive total species flux and report independently normalized
per-species global balances.

## PR 3.1: general species port data

Schema v6 extends the flow graph without changing standalone schema v3/v4 or
flow-only schema v5. It declares logical species and concentration units at
graph scope, maps each logical ID to a native scalar field in every domain,
and lists sorted species IDs on every participating port and pressure-flow
edge. Flux units are explicitly derived as concentration unit times `m^3/s`.
Arbitrary logical and native names are supported; no species is privileged.

Both endpoints must provide and accept concentration and total outward flux.
This is a deliberate reversal-safe capability contract, independent of the
hydraulic pressure/flow receiver roles. Validation rejects empty, duplicate,
or undeclared species, missing domain bindings, reused native fields,
endpoint/edge set differences, and incomplete species capabilities.

`ResolveSpeciesDonor` uses only measured outward flows. It requires a
conservative opposite-sign pair, supports forward and reversed flow, and uses
the last committed donor as the near-zero hysteresis owner. Same-sign flow,
excess flow residual, nonfinite input, and near-zero flow without committed
ownership fail explicitly. Existing native runners reject schema v6 until
PR 3.2 makes transport transactional.

Focused validation:

```text
make -C solvers/coupling test
make -C solvers/cpu coupling_port_test config-check-dispatch-test
./solvers/cpu/coupling_port_test
./solvers/cpu/vca_3d_runtime_test
```

The focused legacy VCA runtime regression passed with the local PETSc/MPI
toolchain. The PR-closing Sol review accepted schema ownership, unit and sign
semantics, reversal/hysteresis routing, graph binding invariants, and v3/v4/v5
compatibility.

## PR 3.2 numerical gate

### Transactional 3D transport foundation

`TransientTransportRuntime` now has the same trial lifecycle needed by a
composite flow/transport adapter: `BeginStep`, deterministic `SolveTrial`,
`RollbackTrial`, idempotent `AbortStep`, two-stage prepare/finalize, and the
existing `Advance` compatibility wrapper. The committed PETSc concentration
vector and logical step count are snapshotted and restored together. Failed or
replayed trials cannot increment physical time state, and finalization is
nonthrowing. The legacy VCA runtime regression covers exact rollback, exact
replay, one commit, open-step abort, and `Advance` parity.
The focused Sol review accepted physical-state ownership, failed-solve
recovery, prepared-step abort, nonthrowing finalization, checkpoint refresh,
and legacy wrapper semantics.

### Composite 3D flow/transport transaction

`ThreeDBodyFittedFlowTransportDomainAdapter` now holds one 3D flow trial open
while transport advances with that trial's ordered nodal velocity. Flow and
transport must expose identical required-node orderings. A joint trial may be
observed or prepared only after both solves succeed; both runtimes prepare
before either nonthrowing finalization. Rollback restores both committed
vectors and clears the input generation so a strong-coupling replay must
reapply, and may change, hydraulic and concentration data. A transport failure
after a successful flow solve leaves both runtimes rollback/abort capable.

Coupled scalar boundaries switch per trial: complete logical concentration
data maps to native Dirichlet values on measured inflow, while measured
outflow uses the natural `advective_outflow` condition. The transport runtime
therefore rebuilds its owned Dirichlet row list from each trial configuration
instead of retaining construction-time rows. Near-zero direction ownership is
represented by whether the graph executor supplies the complete concentration
set, consistent with edge hysteresis. Port output maps native scalar fields
back to logical species and reports the already-verified compiled total
advective/diffusive flux.

Peer total species flux is deliberately not a scalar boundary condition:
imposing it together with concentration would overconstrain the PDE. The
adapter rejects such input explicitly. The species-aware graph executor must
compare the two measured outward total fluxes and reject a nonconservative
trial before preparing any domain; that executor remains part of the open PR
3.2 gate.

The focused PETSc regression covers changed inputs after rollback, exact replay,
dynamic Dirichlet removal and recreation, logical/native mapping, reversal,
two-stage commit, and a real transport-solve failure after accepted flow. It
passes with the local PETSc/MPI toolchain. Focused Sol re-review accepted the
transaction ownership, dynamic row semantics, executor-owned flux residual,
and failure recovery.

### Staged 3D prerequisite

Species routing cannot safely use a guessed combined solve: under reversal the
transport dependency can oppose the hydraulic execution order, guesses have no
convergence contract, and provisional trials would corrupt exact integrated
amount semantics. The accepted architecture therefore adds a separate
`StagedFlowTransportDomainRuntime` beside the unchanged common runtime API.
Hydraulics converge first; accepted outward flows define a donor-to-receiver
transport graph; transport then solves once in its actual topological order.

The composite 3D adapter now implements this staged contract. Hydraulic retries
restore only flow while leaving transport open. Transport retries preserve the
accepted flow vector and velocity while restoring only concentration state.
Compatibility `SolveTrial` and full rollback remain wrappers over both stages.
Observation and prepare require both stages to have succeeded, and abort still
restores both native committed images.

Logical `SpeciesStepAccounting` reports initial/final mass, source amount,
per-port outward amounts, and
`M1-M0+sum(outward)-source`. For the current single-step backward-Euler 3D
transport, physical end-step compiled flux and source rates are multiplied by
the native timestep. `BeginStep` rejects a macro timestep that differs from the
compiled transport timestep before either native transaction opens. The
residual remains an explicit physical diagnostic; SUPG algebraic flux is not
silently folded into it or claimed conservative.

Focused PETSc tests cover changed hydraulic replay, fixed accepted velocity,
changed concentration replay, amount sign/time scaling, timestep mismatch,
compatibility parity, and failure recovery. The staged state and accounting
contract passed focused Sol review and re-review.

The remaining numerical gate is a species-aware graph transaction: route donor
concentrations after each hydraulic trial, compare both measured interface
fluxes, validate time-integrated global balances, and only then prepare every
domain before any finalization. The 3D port flux is

```text
integral((c*u - D*grad(c)) dot n dA).
```

The 1D runtime must report the exact finite-volume face flux used by its
conservative `A*C` update, including diffusion. Global balance must use
time-integrated boundary transfer and sources, not only end-step flux. A
species failure aborts the still-open flow transaction.

The first 1D flux slice now centralizes the numerical face formula
`Q*C_upwind - A*D*(C_right-C_left)/dx` for the conservative `A*C` update. The
runtime stores and reports the actual final-substep numerical boundary flux;
an immediate root bifurcation is aggregated across every root face before the
root outward sign is applied. Terminal reporting uses its stored numerical
face flux with the currently supported zero-gradient exterior trace. Focused
core and runtime regressions pass, and focused numerical re-review accepted
the aggregation, publication timing, and orientation. Time integration of
those face transfers and reversal-safe outlet boundary inputs remain part of
the PR 3.2 gate.

The next 1D conservation slice accumulates root, terminal, and volumetric/wall
source transfers over every internal finite-volume and configured runtime
substep. Per-species accounting reports initial and final mass, outward
boundary amounts, source amount, and the discrete residual
`M1-M0+sum(outward)-source`. Trial reset, failure invalidation, rollback, and
deterministic replay use the existing runtime snapshot. When flow or
vasodilation changes cross-sectional area, the transport update now preserves
the pre-flow `A*C` storage state instead of creating mass by rebuilding it with
the new area. The legacy VCA step result uses the exact accumulated residual
and time-averaged outlet/source rates when available.

Focused numerical review found no formulation or ownership defect and asked
for three additional acceptance checks. They now pass: a high-diffusivity,
nonzero-source case forces three internal finite-volume substeps and compares
independently accumulated root/outlet/source amounts; the macro test confirms
that vasodilation actually changes area while balance remains conserved; and
a partially failed trial rejects accounting before rollback. The affected
PETSc `iga_1d` application also compiles warning-free.

Reversal-safe 1D boundary transport is the next bounded slice. A distal port
may now receive an exterior concentration together with its pressure input;
that trace is used by the same advective-diffusive face operator. A materially
reversed, explicitly coupled terminal without a supplied concentration fails
before the parallel species update; local outlet closures retain their legacy
zero-gradient trace. At an outflowing root, omission of remote concentration
selects the interior trace, avoiding reuse of a stale inlet value. A focused
reversal regression covers missing-data failure, rollback, supplied terminal
upwinding, root and terminal outward signs, and the integrated balance.

The reversed-root donor state uses that same decision in both numerical flux
and `PortState`: without a current supplied trace, concentration is the
outward-flow-weighted interior value across root branches. Mixed materially
inward/outward root branches are rejected because one aggregate concentration
cannot represent both roles. The regression makes configured/stale inlet and
interior concentrations distinct and checks both pre-advance and solved-state
reporting. Focused numerical re-review accepted the direction epsilon,
outward weighting, mixed-flow rejection, fallback consistency, and reversal
signs with no remaining blocker.

The first 3D flux slice evaluates each transport equation's boundary flux from
its compiled advection and diffusion terms, including cross-field terms:
`sum_j(a_ij*c_j*u - D_ij*grad(c_j)) dot n`. The same face quadrature also
reports area-mean concentration. The legacy fields-only measurement remains
an advection-only compatibility path, while the VCA application passes its
compiled transport system and therefore consumes total flux. Dependency-free
unit-cube tests isolate advective and diffusive signs; the PETSc runtime test
checks a linear concentration field with nonzero diffusion.

Physical surface averages now use `2*det(J)*|row(J^-1)|`; flux keeps the
oriented cofactor `sign*2*det(J)*row(J^-1)`. A scaled 2-by-3-by-4 element
regression checks area, concentration integral, and both face-flux signs. VCA
carries measured nonnegative concentration separately from total flux instead
of dividing total flux by flow, and its species residual now sums measured
total outward flux at the inlet and every outlet. SUPG remains excluded from
the physical boundary flux; any discrete mismatch is reported as a residual.
Focused numerical re-review accepted the surface metric, cofactor orientation,
compiled cross-field mapping, independent concentration semantics, and the
`dM/dt + sum(Phi_outward) - source` diagnostic with no remaining blocker.
Measured concentration is retained even on a stagnant outlet: its total
diffusive flux still contributes to species transfer, while its zero flow
contributes no `Q*C` weight to the aggregate concentration. A focused
flowing-plus-stagnant regression prevents fallback contamination; narrow
numerical re-review accepted the fix.

### Staged 1D prerequisite

Native 1D now separates hydraulic advancement from conservative transport
replay. Every configured substep retains its pre-flow area, accepted flow/area
state, sampled inlet schedule, start time, and exact interval. Transport restores
the committed species image and replays those frames through the unchanged
finite-volume `A*C` update, so integrated boundary/source amounts and the final
state match the legacy combined path. Configured-open-loop replay uses each
frame's sampled species values rather than one macro-end value.

The staged adapter exposes logical/native species mapping and exact native
accounting through `StagedFlowTransportDomainRuntime`. A rejected scalar attempt
consumes its concentration-input generation but preserves accepted hydraulic
frames, allowing a changed concentration retry without another flow solve.
Hydraulic rollback and scalar rollback restore their respective committed
images, while abort restores both exactly. Concentration-to-flow feedback such
as vasodilation remains rejected before opening a staged native transaction.

One concentration ownership decision cannot represent a port whose material
flow direction changes within the stored hydraulic frames. The adapter therefore
checks every frame with the configured flow epsilon and rejects such a macro
step until per-frame graph routing exists. Focused tests cover two-substep
nonconstant schedule parity, changed-input scalar retry at fixed flow, direction
change rejection, exact replay and abort restoration, core transport, and the
PETSc 1D build. Focused Sol re-review accepted the replay formulation, ownership,
and transaction semantics.

### Species-aware graph transaction

`SpeciesPressureFlowComponentExecutor` now drives the staged runtimes without
changing the flow-only executor. It converges hydraulics through the existing
explicit/fixed/Aitken pressure-flow formulation, resolves each edge's donor
from accepted outward flows and committed near-zero hysteresis, rejects cyclic
transport dependencies, and advances transport once per domain in deterministic
donor-to-receiver topological order. Receivers obtain complete logical species
maps from already-solved donor ports; there is no predictor or fallback solve.

Conservation gates use native time-integrated physical amounts. For every
edge/species the executor checks the sum of the two outward amounts. It also
independently recomputes every domain residual from initial/final mass, all
outward port amounts, and sources, checks agreement with the native diagnostic,
and reports a component-global balance using a gross-activity normalization
scale. Per-species controls combine an absolute amount tolerance, reference
amount, and relative tolerance.

No domain prepares until routing, transport, accounting, all conservation
gates, and the optional callback succeed. All domains prepare before their
nonthrowing finalizations, and candidate donor ownership becomes committed only
after every finalization. Failure aborts domains in reverse deterministic order,
preserves the primary error, and leaves donor history unchanged. Dependency-free
tests cover forward, reversed, and near-zero routing; two-species complete-map
transfer; deterministic order; hydraulic replay; cycle rejection; independent
edge/domain/global gates; callback, transport, and partial-prepare failure; and
cleanup error aggregation. Focused Sol review and re-review accepted the
formulation, conservation scales, ownership, and transaction semantics.

The remaining PR 3.2 work is schema-v6 execution-control plumbing and native
runner integration. PR 3.3 then closes the phase with a two-species
1D-to-3D-to-1D production regression and verified per-species global balances.
