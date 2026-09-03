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

The moving extension remains NO-GO and unimplemented: it is the face-jump
energy minimization over a time-slab band, pending its own tests and Sol
review.
