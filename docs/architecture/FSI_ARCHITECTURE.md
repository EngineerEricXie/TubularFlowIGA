# FSI Architecture

Status: **PR8.2a bounded rank-local membrane numerical kernel plus PR8.1b
fluid-side traction extraction/projection boundary**. No fluid--structure
runtime solve or moving-domain FSI execution exists in this revision.

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

PR8.1b adds a read-only fluid-side kernel for the immersed Cartesian path.  It
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

## Typed FSI edges and runtime capabilities (PR8.0b1)

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
the full uint64 macro-step input range. Iteration, time/dt/end time, and exact
input/output stamps are also validated. Input can be
set once only while awaiting that iteration; output can be read only after
solve; rollback/reject clears trial input/output; abort returns idle and
invalidates trial state. Prepare does not make output committed; only finalize
records the distinct committed stamp. A production fluid runtime that also
supports scalar ports inherits both
`CoupledDomainRuntime` and `FsiFluidDomainRuntime`, exposing separate APIs
without scalar overloads or ownership ambiguity.  The focused contract test
uses that arrangement and validates exact field stamps, layouts, partitions,
subsystems, and lifecycle availability before accepting an input field.

PR8.0b1 did **not** wire FSI edges into `SimulationGraph`, instantiate a
production FSI runtime, or add a coordinator/executor. It made no claim of a
running FSI solve or collective transaction.

## In-memory graph integration (PR8.0b2)

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

Performance work carries forward the Phase 7 requirement to measure assembly
and solve time separately, host peak RSS, CUDA peak allocation, rank/partition
agreement, and CPU/CUDA field differences. Representative distributed solves
belong on allocated resources rather than login nodes.

## Explicit exclusions

PR8.2a adds only a bounded, dense, single-rank P1 pre-tensioned membrane
kernel: normal scalar displacement, reference normals, explicit Dirichlet IDs,
and backward-Euler trial/prepare/finalize state.  Each successful solve issues
one membrane-instance-owned generation capability from the unchanged committed
state; it must be prepared and finalized, or explicitly rejected/aborted,
before another solve. Rejected, aborted, stale, foreign, superseded, or
modified capabilities cannot commit. Static coercivity and the dynamic SPD
solve use scale-relative Cholesky checks; the dynamic solve also checks finite
values and its backward-error residual. It excludes MPI structural
assembly/collectives, nonmatching transfer, contact, ALE/remeshing, a
monolithic FSI solve, production FSI runtime or coordinator/executor, and any
claim of an operating FSI benchmark.
