# Moving immersed-domain geometry

The moving immersed path keeps the Cartesian background fixed in space: this is an Eulerian
immersed method, not ALE.  Every evaluated prescribed-surface time owns a new,
immutable full-rebuild chain (surface, spatial classification, cut volume,
immersed surface, ghost penalty).  Prepared/trial states are owned separately
and are only committed by replacing an owner after successful construction;
published states are never mutated.

Immersed surface quadrature records canonical material-triangle and
barycentric provenance, so a prescribed wall velocity is evaluated from the
same material point that created each quadrature point. The existing static
path remains unchanged and `.ntiga` is unchanged.

Clipping propagates those canonical barycentrics from the original triangle
through exact intersections and Dunavant interpolation; it never recovers
material coordinates from rounded physical coordinates.  Each evaluation owns
its evaluated time.  `GeometryIdentitySha256()` names only current immutable
geometry, while `PublicationIdentitySha256()` additionally names compatible
predecessor and transition context used for publication/retry decisions.

The fixed-geometry backward-Euler runtime consumes
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

## Transient immersed-wall impedance

The moving-wall formulation keeps the existing symmetric Nitsche traction,
adjoint-consistency, pressure-gap, and velocity-pressure blocks unchanged. It
only replaces the coefficient in the existing penalty residual and tangent:

\[
\eta=\eta_\mu+\eta_t,\qquad
\eta_\mu=16\,\gamma_\mu c_\alpha\mu/h_n,\qquad
\eta_t=16\,\gamma_t\rho h_n/\Delta t.
\]

For the cubic basis, `C_p=16`; `c_alpha=1` under the ghost-covered policy and
`1/alpha` under the legacy non-ghost policy. The inertial term is exactly zero
when `dt=0`, before any reciprocal is formed. Thus the default
`wall_inertial_gamma0=0.0` retains the steady and old low-level paths
bit-for-bit. `wall_gamma0` remains the positive viscous gamma (default 2), and
the inertial gamma is finite and nonnegative. Runtime diagnostics stage and
publish min/max viscous, inertial, and total penalties, component
nondimensional ratios, and component fractions.

This is impedance stabilization only: it adds no ALE or space-time term, mesh
motion, convective penalty, changed traction, FSI coupling, or snapshot-schema
change.

The standalone immersed velocity/history extension admits only exact
positive-cell layouts on a fixed Cartesian
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
hashes and residual/pivot/coverage gates. The target-geometry runtime accepts
an exact compatible history, target-layout
warm-start seed, and immutable map identity.  The outer moving runtime,
publication orchestration, output/washout, and FSI remain out of scope.
Continuation, long-double LDLT, deterministic-storage, directionality/cap,
and fail-closed cases are covered by regression tests.

## Endpoint moving-domain conservation diagnostics

The accepted PDE remains a fixed-grid Eulerian endpoint
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
\(Q_{u,wall}-Q_{w,wall}\) wall-relative leakage.  The first-class endpoint
continuity diagnostic is \(R_{cont}=Q_{port}+Q_{w,wall}\): every configured
port participates, while only configured wall labels contribute material flux.
Its scale is \(\max(Q_{ref},|Q_{port}|,|Q_{w,wall}|)\), and the record checks
\(R_{cont}+(Q_{u,wall}-Q_{w,wall})=Q_u\) to summation roundoff.  This is a
discrete target-endpoint continuity measurement, not a new PDE constraint.
The retained fixed divergence diagnostic
\(R_{div}=\int_{\Omega^{n+1}}\nabla\cdot u-Q_u\), including its historical
`normalized_open_balance` meaning, is explicitly a volume/surface quadrature
consistency measurement; it is not used as the moving-wall continuity gate.

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

## Pure moving immersed field snapshots

The field-snapshot sampler is a pure, immutable post-processing
operation over one exact `MovingCutGeometry`, `ImmersedActiveLayout`, and
`ImmersedGlobalFlowState`, including matching gauge topology.  It neither
captures runtime state nor writes VTU/JSON/PVD output.  It samples usable volume
quadrature in x-fast cell and catalog order, supports both expanded and compact
volume storage, and records velocity, pressure, vorticity, Q, enstrophy density,
and physical integration weight.  Compact coalescing may change stored points,
weights, order, and count while preserving its supported polynomial tensor
moments.  Accordingly, manufactured polynomial fields require expanded/compact
parity to relative-or-absolute tolerance `1e-12` for volume, enstrophy integral,
mean enstrophy, and mean Q.  Discontinuous pointwise indicator measures--the
Q-positive volume/fraction and stagnant volume/fraction--are instead
storage- and quadrature-sampling-dependent estimates.  They must be finite and
their fractions lie in `[0, 1]`, but are not required to agree across storage
modes; refinement/convergence evidence is required before interpreting them
physically.  Snapshot requests explicitly carry sorted configured wall labels
and sorted configured port labels; these disjoint sets cover every retained
surface boundary label and every label has positive area.  Port labels are not
inferred from `layout.PortIds()`: that list contains only flow-controller
multipliers, while pressure-like ports are part of the endpoint boundary too.
Endpoint port flows are explicitly available or unavailable.  Unavailable is
never converted to measured zero: turnover availability is false and the
canonical rate/replacement are zero with no turnover time.  Known zero inflow
is available (also with zero rate/replacement and no turnover time); positive
inflow additionally has the usual turnover time.  Snapshot content identity
carries storage mode, surface partition, and flow availability, so compact and
expanded content hashes intentionally differ; the separate snapshot identity
binds geometry/publication/layout/state/time/index.

`endpoint_turnover_rate_per_s`, `endpoint_turnover_time_s`,
`well_mixed_replacement_fraction_over_step`, and
`stagnant_volume_fraction` are endpoint, instantaneous/global well-mixed, or
threshold proxies.  They are not particle or tracer residence time, retained
blood fraction, a blood-age distribution, or a washout curve.

## Recoverable moving snapshot publication

Publication is a separate, immutable consumer of an already committed field
snapshot. It emits one `<stem>.epoch%020llu` directory with `fields.vtu` and
`metrics.json`, then atomically rebuilds the `<stem>.pvd` collection from all
valid complete epoch directories. The VTU is a quadrature-support point cloud:
every ordered volume support point is an independent `VTK_VERTEX` cell. It is
not a boundary-conforming moving volume mesh and must not be interpreted as
one. The version-3 manifest and VTU field metadata carry identities, hashes,
units, thresholds and aggregate metrics. The manifest also canonically records
transition availability and `dt_s`, wall and port label partitions,
endpoint-flow availability, endpoint-turnover availability, and (only when
available) label-ordered outward port flows. Thus unavailable flows, measured
zero flows, and no-port snapshots remain distinct during recovery and
idempotent retry. A complete orphaned epoch is recoverable after a crash
between the directory and PVD renames; schema or request-semantic mismatches
are rejected while rebuilding the unchanged PVD collection. This is the
current schema-v3 contract; it does not provide a broader artifact-integrity or
recovery protocol.
Static geometry remains on the Bézier VTKHDF visualization path.

Wall mismatch reports area, maximum and RMS relative velocity in m/s, plus
`wall_relative_velocity_squared_area_integral_m4_per_s2`, the dimensionally
explicit integral of squared relative velocity over wall area.  Each requested
wall label must independently resolve to finite, positive surface area.

## Focused idealized-chamber closure driver

`make phase7-lv-closure-test PETSC_DIR=/path/to/petsc` runs one deterministic,
16-step prescribed-motion benchmark.  It uses a closed 12-sector, four-ring
idealized prolate chamber with a wall label and two connected basal half-disk
ports.  The ED/ES/ED frames are source-periodic; ES is the affine `.90, .90,
.95` contraction about the basal plane.  The test uses a fixed `6x6x7`
Cartesian background and depth-2 cut quadrature with empty-rule rescue through
depth 9, density `1050 kg/m3`, dynamic
viscosity `0.012 Pa s`, period `0.8 s`, and an outward inlet controller target
of `-1e-6 m3/s` (inflow magnitude `1e-6 m3/s`) with a zero-pressure outlet.
These are numerical benchmark
settings, not a patient, FSI, ALE, tracer, or clinical model.  The `--one-step`
path additionally re-integrates the solved same endpoint field at depths 1/2/3
(with the same rescue cap) to audit only volume/surface divergence quadrature
consistency; it does not re-solve or change the production state.

By default it removes its temporary output after asserting transactional
publication/retry.  `solvers/cpu/phase7_lv_closure_test /path/to/output-base`
retains the committed VTU/JSON/PVD epochs for inspection.
