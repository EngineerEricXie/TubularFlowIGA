# FSI Architecture

Status: **PR8.0b1 typed edge and runtime capability contracts**. No
fluid--structure solve, membrane model, or moving-domain FSI execution exists
in this revision.

## Scope and first benchmark

The first vertical slice is a small-displacement, pre-tensioned compliant
membrane patch over a rectangular immersed-flow channel. It is intentionally a
benchmark-sized coupling problem: full Cartesian interface fields are carried
even though the initial membrane relaxes only the normal generalized traction.

## Separate surface contracts

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

PR8.0b1 does **not** wire FSI edges into `SimulationGraph`, instantiate a
production FSI runtime, or add a coordinator/executor.  It makes no claim of
a running FSI solve or collective transaction.

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

PR8.0b1 excludes a structural solver, traction integration implementation,
fluid-side surface extraction, MPI collectives, nonmatching transfer,
contact, ALE/remeshing, a monolithic FSI solve, production FSI runtime or
coordinator/executor, and any claim of an operating FSI benchmark.
