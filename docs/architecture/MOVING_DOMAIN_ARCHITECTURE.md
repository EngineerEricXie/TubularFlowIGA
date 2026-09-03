# Moving immersed-domain geometry

Phase 7 keeps the Cartesian background fixed in space: this is an Eulerian
immersed method, not ALE.  Every evaluated prescribed-surface time owns a new,
immutable full-rebuild chain (surface, spatial classification, cut volume,
immersed surface, ghost penalty).  Prepared/trial states are owned separately
and are only committed by replacing an owner after successful construction;
published states are never mutated.

Immersed surface quadrature records canonical material-triangle and
barycentric provenance, so a prescribed wall velocity is evaluated from the
same material point that created each quadrature point.  The Phase 6 static
path remains unchanged and `.ntiga` is unchanged.

Clipping propagates those canonical barycentrics from the original triangle
through exact intersections and Dunavant interpolation; it never recovers
material coordinates from rounded physical coordinates.  Each evaluation owns
its evaluated time.  `GeometryIdentitySha256()` names only current immutable
geometry, while `PublicationIdentitySha256()` additionally names compatible
predecessor and transition context used for publication/retry decisions.

S7-A is GO for the fixed-geometry backward-Euler runtime (PR7.3b): it consumes
one immutable `MovingCutGeometry`, requires its exact current identity/layout
at both ends of a step, uses current velocity for convection and an
identity-only velocity history keyed by global Cartesian node ID.  Geometry
motion, extended history, ALE terms, and material wall motion are rejected;
the latter is certified from each retained quadrature point's canonical
provenance before assembly.  Publication is an explicit
idle → trial → prepared → finalized owner transaction: preparation may fail
without visibility, finalize only swaps already prepared owners, rollback
restores the frozen seed/history, and abort leaves the committed global state
and port measurements unchanged.  PETSc vectors are trial-local; canonical
publication remains `ImmersedGlobalFlowState`.

PR7.4a is committed as a standalone, non-integrated immersed velocity/history
extension.  It admits only exact positive-cell layouts on a fixed Cartesian
grid, grows a deterministic six-face time-slab band from the old positive
cells, and uses that forward-only containment check as its sole acceptance
gate.  A deterministic sorted-level, six-neighbor multi-source reverse BFS
from the new positive cells is retained only as an exact capped diagnostic: it
counts seeds, stops once every old positive cell is reached, and fail-closes on
its explicit resource cap or an exhausted unresolved frontier.  Its storage is
grown only from actual bounded seeds and visits, never from the configured cap.
A reverse
distance greater than the requested forward layers is reported, not rejected.
The
extension minimizes the unregularized cubic third-normal-derivative jump
energy on interior band faces, preserves old-layout anchor values bitwise, and
uses one bounded dense scaled-Cholesky factor for all velocity components and
optional scalar warm starts.  It publishes only a fully covered target
`ImmersedVelocityHistory`, with deterministic operator/reduced/extension
hashes and residual/pivot/coverage gates.  PR7.4b slice 1 adds only a
target-geometry runtime entry accepting an exact PR7.4a history, target-layout
warm-start seed, and immutable map identity.  The outer moving runtime,
publication orchestration, output/washout, and FSI remain out of scope.
Independent PR7.4a continuation, long-double LDLT, deterministic-storage,
directionality/cap, and fail-closed evidence is passing.

## Endpoint moving-domain conservation diagnostics

PR7.4b keeps the accepted PDE exactly as a fixed-grid Eulerian endpoint
backward-Euler solve on \(\Omega^{n+1}\): the old velocity is extended to the
target layout, spatial terms use the conservative mixed form with
current-velocity convection, and the material-velocity Nitsche wall term is
evaluated at target canonical material points.  It introduces no ALE
convection, mesh Jacobian, local swept-volume term, space-time cut integration,
or FSI coupling.

The diagnostics use outward normals and target-endpoint quadrature.  They
report the lab-frame fluid flux \(Q_u=\int_{\partial\Omega^{n+1}}u\cdot n\),
by boundary label and in total; port \(Q_u\) semantics remain absolute.  The
material flux is also evaluated by boundary label over the *entire* target
surface, \(Q_w=\int_{\partial\Omega^{n+1}}w\cdot n\), using immutable
canonical triangle and barycentric provenance.  A separate wall-only material
flux \(Q_{w,wall}=\int_{\Gamma_w^{n+1}}w\cdot n\) is retained for
\(Q_{u,wall}-Q_{w,wall}\) wall-relative leakage.  The fixed divergence
diagnostic is retained as \(R_{div}=\int_{\Omega^{n+1}}\nabla\cdot u-Q_u\).

For each moving transition, immutable old/new closed-surface audited volumes
and exact source/target epoch identity give
\(G_{BE}=(V^{n+1}-V^n)/\Delta t\),
\(R_{Reynolds}=G_{BE}-Q_w\), and
\(R_{moving}=G_{BE}+Q_u-Q_w\).  Positive values mean expansion or outward
flux.  Flux maps, the raw defects, and normalized defects use the explicit
scale \(\max(Q_{ref},|G_{BE}|,|Q_u|,|Q_w|)\).  In particular,
\(R_{Reynolds}\) need not vanish for a general finite deformation: it is a
diagnostic, not a hidden GCL or momentum correction.  Conservation records are
built before inner prepare and are retained after the noexcept final owner
exchange without exposing a trial as committed.
