# FSI architecture

Fluid--structure interaction couples a native tetrahedral ALE fluid to a
native tetrahedral hyperelastic solid through a matching interface. A strong
Dirichlet--Neumann coordinator with dynamic weighted Aitken relaxation drives
the reference runtime, which runs on a single partition (one MPI rank).

The earlier immersed-flow and pre-tensioned-membrane path was removed on
2026-10-08. Its code, tests, and documentation are preserved at the git tag
`archive/immersed-shell-2026-10`.

## Components

| Layer | Implementation | Role |
|---|---|---|
| Surface contract | [`DistributedSurfaceInterface.hpp`](../../include/DistributedSurfaceInterface.hpp) | Interface IDs, boundary labels, reference mesh/layout/partition identities, and provided/required surface fields |
| Edge | [`FsiCouplingEdge.hpp`](../../include/FsiCouplingEdge.hpp) | Typed directional fluid/structure edge |
| Runtime capabilities | [`FsiDomainRuntime.hpp`](../../include/FsiDomainRuntime.hpp) | Fluid and structure capability mixins and the `FsiTrialLifecycle` owner |
| Coordinator | [`StrongFluidStructureCoupling.hpp`](../../solvers/cpu/include/StrongFluidStructureCoupling.hpp) | Strong Dirichlet--Neumann iteration and paired commit |
| Relaxation | [`DynamicWeightedAitkenRelaxation.hpp`](../../include/DynamicWeightedAitkenRelaxation.hpp) | Area-weighted dynamic Aitken factor |
| Interface transfer | [`NativeTetMatchingFsiInterface.hpp`](../../solvers/cpu/include/NativeTetMatchingFsiInterface.hpp) | Matching-vertex kinematics and traction transfer |
| Fluid adapter | [`NativeTetAleFsiRuntime.hpp`](../../solvers/cpu/include/NativeTetAleFsiRuntime.hpp) | Harmonic ALE mesh motion and P2/P1 flow |
| Structure adapter | [`NativeTetSolidFsiRuntime.hpp`](../../solvers/cpu/include/NativeTetSolidFsiRuntime.hpp) | Hyperelastic tetrahedral solid |
| Checkpoint | [`NativeTetFsiCheckpoint.hpp`](../../solvers/cpu/include/NativeTetFsiCheckpoint.hpp) | Single-payload checkpoint of the single-partition pair |

The CPU executable is `solvers/cpu/native_tet_fsi_steps` (one MPI rank). The
CUDA executable `solvers/cuda/native_tet_fsi_cuda` uses explicit sequential
coupling instead of the strong coordinator: the GPU fluid step supplies
traction, the host solves the static solid, and the host updates the harmonic
ALE mesh before the next fluid step.

## Matching interface

"Matching" means that fluid ALE vertices and solid reference vertices share
stable node IDs and the selected interface triangles have identical
connectivity. Solid displacement defines the current boundary positions;
solid velocity becomes the ALE boundary velocity and, under no slip, the fluid
wall velocity. A native P2 fluid edge node takes the average of its two P1
solid endpoint velocities, and the harmonic mesh-motion stage extends boundary
motion into the interior. The fluid publishes `traction_on_structure_pa`
(the negated area-averaged `sigma * n_fluid`), integrated as three equal P1
nodal forces per current triangle. The contract and its regressions are in
[T5 native matching interface](../validation/T5_NATIVE_MATCHING_INTERFACE.md).

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

## Graph integration

`SimulationGraph` stores typed `FsiCouplingEdge` values separately from scalar
`CouplingEdge` values. `DomainNode` declares surface catalogs and exact
per-surface layouts; graph construction validates directional fluid/structure
capabilities, endpoint subsystem/interface identity, mesh, layout, rank-local
partition, boundary labels, globally unique edge IDs, and one-to-one endpoint
use. Declared surfaces may not be orphaned. `MakeFluidStructurePairPlan`
returns the two directional field exchanges without changing scalar
pressure/flow ordering or cycle checks.

No graph domain kind implements an FSI endpoint yet:
`IsSupportedFsiDomainPair` accepts no pair, so every FSI edge is rejected at
graph construction. Adding native tetrahedral ALE-fluid and solid domain kinds
is the integration step that will enable graph FSI.

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
Focused coverage includes deterministic mock adapters and the native
tetrahedral compliant-channel benchmark (`make -C solvers/cpu
native-tet-compliant-channel-fsi-test`). The latter is closure evidence for this
bounded single-partition slice, not a distributed performance or
production-anatomy claim.

Performance measurements separate assembly and solve time, host peak RSS,
CUDA peak allocation, rank/partition
agreement, and CPU/CUDA field differences. Representative distributed solves
belong on allocated resources rather than login nodes.

## Current limitations

- The reference runtime is single-partition. Both the ALE fluid
  (`NativeTetAleDenseRuntime`) and the solid static solver use dense direct
  solves, so it is suitable only for small reference cases. The distributed
  PETSc ALE flow runtime (`NativeTetAlePetscRuntime`) is not yet connected to
  the FSI adapter.
- The multidomain graph has no FSI-capable domain kinds.
- Nonmatching interface transfer, contact, remeshing, and monolithic coupling
  are not implemented.
