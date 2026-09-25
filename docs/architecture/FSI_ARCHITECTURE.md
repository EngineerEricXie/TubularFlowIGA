# FSI Architecture

The foundational compliant-channel path combines a bounded strong
Dirichlet--Neumann coordinator, a moving immersed-flow adapter, a
pre-tensioned membrane adapter, patch-to-closed-material composition, and
traction extraction/projection. The reference benchmark remains a sequential
`PETSC_COMM_SELF`, one-partition configuration.

## Scope and first benchmark

The first vertical slice is a small-displacement, pre-tensioned compliant
membrane patch over a rectangular immersed-flow channel. It is intentionally a
benchmark-sized coupling problem: full Cartesian interface fields are carried
even though the initial membrane relaxes only the normal generalized traction.

## Separate surface contracts

`MaterialSurfaceKinematics.hpp` is the immutable producer-neutral geometry
payload consumed by the moving cut geometry, moving-wall Nitsche path, and
moving immersed runtime.  One payload binds exact material/topology identities,
immutable reference material vertices kept separately from current canonical
geometry and current material vertices, wall-velocity provenance, and
the exact time/step interval.  It owns neither an FSI graph endpoint nor a
distributed field layout.  The public neutral factory recomputes and verifies
domain-separated material, topology, and complete current-content digests from
its owned fields.  The material digest binds immutable reference coordinates,
source connectivity, labels, and topology provenance; current displacement
cannot change it.  The factory accepts reference geometry separately from the
current state and recomputes every supplied digest, so a producer cannot claim
an unrelated material identity.  It rejects incomplete, non-bijective,
mislabeled, split-coincident, or orientation-reversed canonical/material
provenance (while accepting cyclic/even triangle permutations).  The content digest
binds current coordinates and wall velocity and is included in downstream cut
geometry/publication identities, so equal claimed producer strings cannot
alias numerically different Nitsche inputs.  Prescribed motion currently
produces the payload through the legacy `PrescribedSurfaceMotion::Evaluation`
spelling; its historical evaluation hash remains the separately verified
geometry-epoch identity byte-for-byte, while the new content digest provides
neutral integrity.  A later membrane producer can use the same public factory
by supplying its immutable reference material geometry independently of its
evaluated displacement state.
Consumers reject stale time/dt/content or changed material/topology before
trial geometry is built.

The implementation provides a producer-neutral `CreateFromSourceTopology` route and an
immutable `MaterialSurfacePatchMap`.  The map is initially deliberately narrow:
one fully-owned partition, one label, and an explicit conforming P1 subset of
the directed source triangles.  It derives a patch-scoped reference digest,
requires the layout and interface to carry it, accepts no coordinate-search
authority, and proves exact source/canonical correspondence, orientation,
connectivity, and the clamped one-sided seam.  Its reference digest is strictly
patch-scoped: sorted global IDs, their mapped immutable reference positions,
the directed/cyclic patch triangle layout, and exact selected triangle labels.
It intentionally excludes whole-material/topology identities and clamp/mapping
configuration, so an unchanged patch has a stable reference digest when an
unrelated closed-surface remainder changes.  The distinct map identity binds
that patch digest to complete material/topology identities, the full
distributed-interface identity, explicit source mapping, triangle map, and
clamps.  The endpoint is the exact structural role: it provides sorted
`Displacement, Velocity` and requires only `TractionOnStructure`.
`MaterialSurfacePatchKinematics` is stateless: it composes a context/stamp-exact
patch trial into a newly validated whole closed surface, fixes the remainder at
the immutable reference, requires zero seam fields and a roundoff-scaled
backward-Euler check using actual SI terms and the patch geometry scale (with
exact zero handled exactly, not a one-metre floor), and binds map, committed
content, patch publication, context, and target content in its composition
identity.

`DistributedSurfaceInterface.hpp` is the dependency-free contract for
field-valued interface exchange. It is separate from scalar P/Q
`CouplingPort`: a fluid publishes traction on the structure and a structure
publishes displacement and velocity. A surface layout identifies reference
material nodes, owned IDs, reference triangles, and positive reference lumped
areas. Every published field stamp also carries the immutable partition
identity derived from its exact owned IDs and weights, so equal-sized ownership
slices cannot cross-bind. Ghost values are local scratch only and never confer
publication ownership. Material-surface interfaces require nonempty,
canonical sorted boundary labels (label zero is valid). The first fluid and
structure interfaces are bidirectional: each must declare a nonempty
`provides` list and a nonempty `requires` list; both lists are canonical sorted,
unique, and disjoint. Interface identities delimit and count each list
explicitly.

The sign convention for the fluid load is fluid-on-structure Cauchy traction:

\[
t_{\mathrm{on\ structure}}=-\sigma_f n_f.
\]

The first slice requires identical material topology and global node IDs on
both sides. Nonmatching interpolation/projection is explicitly deferred.

Fluid-side extraction is patch-authoritative: it requires the
immutable `MaterialSurfacePatchMap`, accepts a fluid endpoint distinct from
the map's structural endpoint, and requires the exact fluid role
`TractionOnStructure <- Displacement, Velocity`.  The fluid interface and
layout bind the map's patch-scoped reference digest, labels, layout and
partition; full material/topology identities remain bound to the Cartesian
domain and quadrature catalog.  Selected cut cells and points are determined
only by mapped canonical-triangle membership.  Every retained point checks
its canonical/source/map labels, and P1 barycentric coordinates are routed
through the explicit mapped triangle/node permutation.  State and projection
identities include map identity and membership.

The immersed Cartesian path provides a read-only fluid-side kernel. It
uses the catalog normal exactly as the closed-surface outward fluid normal (the
catalog and moving-cut provenance audit compare it to
`ClosedTriangulatedSurface::outward_unit_normal`), evaluates
\(\sigma_f=-pI+\mu(\nabla u+\nabla u^T)\), then publishes
\(-\sigma_f n_f\).  Retained points are filtered by the exact material
boundary labels.  A bounded, single-partition consistent P1 surface-mass
projection produces nodal traction while its unmodified right-hand side is the
consistent nodal force; the kernel checks resultant and current-configuration
moment conservation with compensated sums, explicit SI absolute tolerances,
and a relative tolerance scaled by accumulated absolute force/moment
contributions.  The selected interface must explicitly provide
`TractionOnStructure`.  It rejects distributed layouts, incomplete/duplicate
or reordered state-cell ownership, inconsistent coefficients at shared global
IGA nodes, stale stamps, unrelated material/layout/current-content identity,
and caller-supplied producer-state hashes that do not recompute from the
supplied IGA state plus the bound Cartesian domain, element
connectivity/extraction, and quadrature catalog.  Projection identity is
derived from the projection form, exact field/material/layout identities,
state, labels, and retained quadrature content.  Runtime publication and
coupling remain deferred.

## Typed FSI edges and runtime capabilities

`FsiCouplingEdge.hpp` defines `FsiCouplingEdge`, a deliberately separate type
from scalar `CouplingEdge`/`PortRef`.  Its `fluid` and `structure` endpoints
are explicitly directional `SurfaceInterfaceRef` values: each binds domain,
subsystem, and interface ID. Declarations must repeat that exact subsystem;
the endpoints belong to distinct domains and currently permit only the
`FluidStructureTractionKinematics` law.  Its deterministic SHA-256 identity
includes the edge ID, both directional endpoints (including subsystem), and
law; subsystem mutation or swapping endpoints changes the identity and fails
directional catalog validation.

The edge helpers require a fluid surface to provide fluid-on-structure
traction and require displacement plus velocity, while a structure surface
has the inverse declaration.  They require exact reference-mesh identity,
exact boundary-label compatibility (including valid label zero), equal global
reference-layout identities, and equal rank-local partition identities when
layouts are available. The partition check binds rank, owned global IDs, and
reference lumped weights, so equal global layouts with different ownership
slices cannot bind. Catalog and binding helpers reject duplicate edge IDs,
duplicate catalog endpoints, and any attempt to bind one surface endpoint to
more than one FSI edge. These are rank-local contract checks:
they do not assert distributed collective coverage from rank-local catalogs or
layouts.

`FsiDomainRuntime.hpp` supplies pure capability mixins instead of extending
`CoupledDomainRuntime`: `FsiFluidDomainRuntime` catalogs surfaces, accepts
kinematics, and returns traction; `FsiStructureDomainRuntime` catalogs
surfaces, accepts traction, and returns kinematics. `FsiTrialLifecycle` is a
dependency-free embedded availability owner required by the first runtime
integration and its mocks. It binds exact edge/endpoints/layouts/partitions
and moves through `idle -> step-active -> iteration-awaiting-input ->
input-ready -> solved -> prepared -> idle`. `BeginStep` accepts the existing
nonnegative `int` `DomainStepContext::step_index`, validates it, then widens
the accepted value for uint64 storage and context hashing; it does not accept
the full uint64 macro-step input range. Iteration, time/dt/end time, and the
exact input stamp are also validated. Before solve, an output envelope fixes
the exact time, step, coupling iteration, reference/layout, and partition but
intentionally has no producer-state identity. The producer supplies that
identity with its actual solved field; the lifecycle validates it against the
envelope and retains that exact actual stamp for trial, prepared, and committed
output gates. Input can be
set once only while awaiting that iteration; output can be read only after
solve; rollback/reject clears trial input/output; abort returns idle and
invalidates trial state. Prepare does not make output committed; only finalize
records the distinct committed stamp. A production fluid runtime that also
supports scalar ports inherits both
`CoupledDomainRuntime` and `FsiFluidDomainRuntime`, exposing separate APIs
without scalar overloads or ownership ambiguity.  The focused contract test
uses that arrangement and validates exact field stamps, layouts, partitions,
subsystems, and lifecycle availability before accepting an input field.

This interface layer does **not** by itself wire FSI edges into `SimulationGraph`, instantiate a
production FSI runtime, or add a coordinator/executor. It made no claim of a
running FSI solve or collective transaction.

## Structure-side runtime adapter

`PretensionedMembraneFsiRuntime.hpp` is the first production implementation of
`FsiStructureDomainRuntime`. It composes exactly one `FsiTrialLifecycle` with
exactly one `PretensionedMembrane`; construction independently binds both
owners to the same directional edge, local structure endpoint/capability,
fluid peer endpoint/capability, reference layout, and rank-local partition.
The inherited membrane restriction remains explicit: this is a
single-partition rank-local runtime and makes no MPI or collective claim.

The public sequence is `BeginMacroStep`, `BeginCouplingIteration` with an exact
traction stamp and producer-neutral kinematics envelope, `SetSurfaceTraction`, `SolveMembraneTrial`,
`GetSurfaceKinematics`, then reject/abort or prepare/finalize. Inputs and
outputs are value snapshots. Rejected iterations clear both owners, so every
retry starts from the same committed membrane state. The generated kinematics
stamp must match the declared envelope; its producer-state identity is created
only by the membrane solve and becomes the exact stored output identity. The
runtime exposes no membrane owner alias, only immutable value snapshots and
diagnostic identity copies.

Prepare copies the prospective committed public field and completes all
fallible identity work. The lifecycle stores a prepared output stamp and the
membrane stores a prepared committed-state identity. Finalize prevalidates both
owners and then uses only noexcept exchanges for numerical state, committed
kinematics, and lifecycle availability; neither owner can expose or commit a
partial trial. The committed-field accessor therefore remains on its old
snapshot until finalize.

## Moving immersed-flow runtime adapter

`MovingImmersedTransientFlowFsiRuntime.hpp` is the rank-local fluid-side
implementation of `FsiFluidDomainRuntime`. It composes one moving immersed
flow owner, one `FsiTrialLifecycle`, an immutable patch map, and value
snapshots for accepted kinematics, target material state, traction, transition
conservation, and force/moment projection audits. The accepted patch is
composed into the closed material surface before moving cut geometry,
moving-wall Nitsche assembly, state transfer, and patch-authoritative traction
projection are evaluated.

Input acceptance is one-shot: phase/duplicate guards precede snapshot
mutation, all fallible validation and copying happen in a local candidate, and
the lifecycle gate is opened immediately before noexcept publication. Stale,
foreign, wrong-interface, or duplicate input therefore leaves an accepted
input and its pending solve intact. Conservation is exposed only with solved
trial traction or a committed transition; reject/abort retain committed flow,
geometry, traction, and diagnostics. Prepare completes fallible preparation
before publication, and finalize uses only prevalidated noexcept handoffs. The
lifecycle's final handoff is private and friended only to the two production
adapters.

The solved transition follows the same rule: traction, target material state,
composition identity, and value diagnostics are fully staged before
`MarkSolved`. After that call succeeds, the adapter performs only noexcept
swaps, so a projection/snapshot allocation failure cannot advertise a solved
field without its matching numerical epoch. The focused channel fixture uses
triangle-derived reference lumped areas `{.03,.045,.015,.045,.09,.045,.015,.045,.03}`
(sum `.36 m^2`), a nonzero moving wall speed, and a normalized moving-mass
acceptance bound below `.03`; these are a bounded regression check, not a
general conservation claim.

## In-memory graph integration

`SimulationGraph` now stores typed `FsiCouplingEdge` values separately from
scalar `CouplingEdge` values. `DomainNode` declares material-surface catalogs
and exact per-surface layouts; graph construction validates the directional
fluid/structure capabilities, endpoint subsystem/interface identity, mesh,
layout, rank-local partition, boundary labels, globally unique edge IDs, and
one-to-one endpoint use. Declared surfaces may not be orphaned. The initial
topology is deliberately limited to one `ThreeDImmersedFlow` domain paired
with one `SurfaceMembraneStructure` domain. A membrane has material topology
dimension two and embedding dimension three; it is not a scalar-flow domain.
`MakeFluidStructurePairPlan` returns the two directional field exchanges
(structure displacement/velocity to fluid, fluid traction to structure)
without changing scalar P/Q graph ordering or cycle checks. Parser and schema
formats remain unchanged.

## Weighted Aitken

Strong coupling uses one scalar dynamic Aitken factor over the flattened
interface vector. Each rank retains its positive, unnormalized owned
reference-area weights \(a_i\) and receives the explicit positive global area
\(A=\sum_i a_i\). For the accepted residual \(r_{k-1}\) and current residual
\(r_k\), each rank contributes

\[
N_{\mathrm{local}}=\sum_i a_i r_{k-1,i}(r_{k,i}-r_{k-1,i}),\qquad
D_{\mathrm{local}}=\sum_i a_i(r_{k,i}-r_{k-1,i})^2.
\]

The runtime must sum those two terms and max-reduce a common residual scale
\(s\). Before that reduction/proposal, it must allgather-and-compare or
broadcast the `ControlStateIdentitySha256()` value. The proposal API requires
that expected identity and fails closed when a rank has desynchronized control
state. The identity covers only globally identical state that can affect the
scalar: reset phase/generation, accepted-iteration count, previous relaxation,
history-present flag, controls, and global area. It also encodes whether a
proposal is pending; in that phase it binds the globally identical proposal
identity, so a pending rank and a proposal-ready rank compare unequal while
synchronized pending ranks compare equal. It deliberately excludes the
partition-local previous residual vector and local weights, which remain bound
to each rank's partition/history. With identical validated \(N\), \(D\), \(s\),
and control identity, every rank produces the same \(\omega\) and applies
\(\omega_{k+1}=-\omega_k N/D\). There is intentionally no hidden MPI
dependency or rank-local normalization. The positive finite
`reference_scale` is a lower bound, not a fixed scale: the common scale is
`max(reference_scale, abs(all current and accepted previous residual
entries))`. It is used only by the tiny-difference fallback, which tests
\(D/(A s^2)\) against the configured threshold; it does not alter the Aitken
candidate. The single-rank convenience route computes the same local
contributions and calls the reduced proposal route exactly. The partition and
unnormalized weight identity are fixed by the Aitken owner; any change requires
a fresh owner. A proposal has two identities: its global identity covers only
the collectively identical control and reduction facts and is the identity
bound into a pending control state; its local identity additionally covers the
bound partition/weight identity and exact-bit local iterate, residual, update,
and next vectors. `AcceptApplied` recomputes both from the supplied fields and
requires an exact-bit match to its locally pending proposal (including status
and diagnostics), rather than trusting a carried identity string. The proposal
is one-use, its applied scalar must match exactly, and a reset invalidates all
pending proposals. This contract defines the distributed reduction and
ownership requirements only; it makes no runtime FSI claim.

## Ownership and future transaction

Before a surface field can be published, its producer must be producer-neutral
with respect to the surface contract: it must establish the exact material
surface, reference/layout identities, owned values, and a producer-state
identity without depending on a future coupling owner. The lifecycle enforces
the field portion of a transaction shaped as `idle -> step-active -> iteration
-> solved -> prepared -> finalized`: validate unpublished fields and
convergence data first, then perform a non-allocating ownership exchange.
Abort discards trial data; ghosts remain scratch.

## Strong Dirichlet--Neumann coordinator

`StrongFluidStructureCoupling.hpp` owns one exact `FsiCouplingEdge`, one
single-partition layout, its convergence controls, and a
`DynamicWeightedAitkenRelaxation` owner. It publishes a producer-owned,
backward-Euler-consistent predictor, performs exact kinematics-to-fluid then
traction-to-structure trials, and measures scalar normal-displacement
residuals against immutable reference normals. RMS uses reference lumped area
weights and convergence is explicit: `rms <= absolute_m + relative *
max(reference_m, raw/current displacement scale)`. Clamped nodes stay exactly
zero. Rejected trials are discarded by both runtimes before the one-use Aitken
proposal is accepted and a new Cartesian displacement/BE velocity field gets a
derived identity binding raw field, committed structure state/model, context,
and Aitken facts; a raw stamp is never reused for a modified iterate.

Both runtimes prepare before coordinator-only prevalidation permits private
noexcept finalization handoffs.  The coordinator first obtains the exact
prepared final field/state/composition identities and allocates/hashes its
complete result; only then may it prevalidate and execute the two private
noexcept handoffs.  Thus preparation, result staging, or prevalidation failure
aborts both with no cross-runtime partial commit, while the post-first-finalize
tail contains only noexcept handoffs, scalar state updates, and a statically
verified noexcept result move. Exceptions and maximum-iteration
nonconvergence also abort both and reset Aitken pending state. Diagnostics are
value snapshots binding all convergence options, context, histories,
norms/scales/thresholds, exact fluid/raw/current/accepted identities, and
Aitken control/proposal/factor facts. A converged trial explicitly records that
no Aitken proposal was applied; it does not invent an initial status or a
factor of one. Immutable reference normals must already be unit length to
roundoff, so scalar projection and Cartesian reconstruction cannot rescale a
field.
Focused coverage includes deterministic mock adapters and the controlled real
compliant-channel benchmark.  The latter is closure evidence for this bounded
single-partition slice, not a distributed performance or production-anatomy
claim.

Performance measurements separate assembly and solve time, host peak RSS,
CUDA peak allocation, rank/partition
agreement, and CPU/CUDA field differences. Representative distributed solves
belong on allocated resources rather than login nodes.

## Distributed moving-FSI runtime and restart

The current CPU path uses owned surface partitions and collective execution.
`ImmersedMovingDistributedFsiRuntime`
borrows a distributed moving-flow runtime and material patch map. It composes
owned kinematics into the bounded replicated cut-surface geometry, solves only
owned PETSc field rows, publishes traction and consistent nodal force slices,
and checks five conservation gates before paired commit.

`SingleOwnerMembraneRuntime` deliberately keeps the small dense membrane solve
on one configured rank. Its input traction is gathered from unique surface
owners; its committed displacement and velocity are redistributed to owned
surface slices. This is an implemented distributed fluid/transfer path with a
bounded centralized structure cost, not a claim of distributed structural
assembly.

`SolveDistributedStrongFsiStep` computes weighted Aitken reductions, RMS and
maximum residuals globally. Every iteration either leaves both runtimes at the
same accepted state or reaches the paired no-throw finalize tail. Rank-local
errors, convergence failure and exhaustion collectively abort both candidates.

`MovingFsiCheckpointBundle` publishes the accepted moving-flow shards, one
traction slice per source rank, and owner-only membrane metadata/state under a
single manifest. Restore constructs fresh flow, adapter and membrane owners;
failure closes the unpublished flow candidate. For a changed rank count, all
authenticated source surface records must cover every stable node exactly once.
The target keeps only its owned records and receives new producer, projection
and partition identities derived from the source payload set. Same-partition
restore retains exact saved publication bytes. This bounded surface restore
reads every source surface shard on each target rank; it does not gather the
full fluid field and is not intended as a large distributed shell solver.

Local regression targets cover 1/2/4-rank execution and restart. Cross-node
execution, scaling measurements, nonmatching transfer, contact and
native graph CLI integration remain outside the currently accepted scope.

## Explicit exclusions

The bounded single-partition path excludes MPI/collective strong execution, nonmatching transfer,
contact, ALE/remeshing, monolithic coupling, and any distributed performance
claim.  The closed compliant-channel benchmark does not close the subsequent
compliant-tube, aneurysm-wall, thin-shell valve, valve opening/closing,
leaflet-contact, or patient-specific valve work.
The structure kernel is a bounded, dense, single-rank P1 pre-tensioned membrane
model with normal scalar displacement, reference normals, explicit Dirichlet IDs,
and backward-Euler trial/prepare/finalize state. Its local structure runtime and
the bounded local fluid adapter do not by themselves add a coordinator or
collective execution. Each successful solve issues
one membrane-instance-owned generation capability from the unchanged committed
state; it must be prepared and finalized, or explicitly rejected/aborted,
before another solve. Rejected, aborted, stale, foreign, superseded, or
modified capabilities cannot commit. Static coercivity and the dynamic SPD
solve use scale-relative Cholesky checks; the dynamic solve also checks finite
values and its backward-error residual. It excludes MPI structural
assembly/collectives, nonmatching transfer, contact, ALE/remeshing, a
monolithic FSI solve, a coordinator/executor, and any claim of a distributed
operating FSI benchmark.
The patch mapping excludes distributed patch ownership, multiple labels, nonmatching
transfer, patch remeshing/contact, and any change to fluid traction/runtime
paths.
