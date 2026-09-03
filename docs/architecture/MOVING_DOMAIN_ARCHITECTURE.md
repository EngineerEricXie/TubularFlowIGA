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

Transient extension and any time-slab formulation remain gated by Sol S7-A;
this geometry oracle does not implement transient flow or FSI coupling.
