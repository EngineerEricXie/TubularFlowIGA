# Phase 5 report

Status: **in progress**.

## PR 5.2b2: strict VTP PolyData ingestion complete

The dependency-free CPU layer now reads bounded VTP PolyData buffers, text,
and paths into the canonical closed-surface validator.  It supports strict
ASCII arrays plus uncompressed VTK length-prefixed inline `binary` base64 and
`appended` base64 blocks, using only little-endian UInt32/UInt64 headers.
Appended base64 blocks use VTK encoded-character offsets and must cover that
encoded stream exactly without overlaps, gaps, or trailing bytes.  Raw appended
data and compression remain unsupported.
Its small XML tokenizer accepts a leading declaration and comments while
rejecting DTDs, entities, later processing instructions, malformed
nesting/attributes, and oversized input/tokens.  It accepts exactly one
PolyData Piece with points and triangular polygons, and optionally maps one
integral CellData boundary array (default `boundary_id`); an empty `CellData`
is treated as absent and uses the configured default label.  Counts, offsets,
indices, labels, finite numeric tokens, base64 padding, typed payload sizes,
and IEEE binary values are checked before the existing builder applies physical
scale once.  VTP reader work is complete for this Phase 5 scope.
Focused coverage is `make -C solvers/cpu surface_vtp_test`; it is included in
the dependency-free aggregate test target.

## PR 5.2a4: strict STL ingestion

The dependency-free CPU layer now reads strict ASCII and binary STL buffers or
paths into the canonical closed-surface validator.  Binary detection is based
solely on exact record-size arithmetic, including `solid` binary headers;
facet normals are advisory but finite, attributes must be zero, and malformed
or trailing content is rejected.  ASCII parsing accepts only complete STL
facet structure and finite full numeric tokens.  Physical scaling is applied
once by the existing canonical surface builder.

## PR 5.2a3: deterministic welding and self-intersection preflight

The in-memory surface preflight now supports exact or bounded deterministic
vertex welding, retaining a lexicographic exact input representative and
rejecting transitive clusters wider than the configured tolerance.  It records
source/canonical/welded counts and applied physical scale/tolerance.  A
dependency-free median-split triangle AABB BVH performs closed-shell
self-intersection preflight after topology/orientation, with candidate-count
diagnostics.  Readers, point classification, and cut quadrature remain
outside this PR.

## PR 5.3a: immutable surface spatial queries

`SurfaceSpatialIndex` owns a validated canonical surface and builds a
deterministic median-split immutable triangle BVH for point and closed-box
queries. Point queries first use exact dyadic point-on-triangle tests, then a
half-open exact `+x` parity ray; arithmetic-cap failures are reported as
`Ambiguous`. Box queries use exact dyadic triangle-box SAT axes and
distinguish no contact, boundary-only contact, and contact through the open
box interior. Query diagnostics report checked candidate, narrow-phase, and
crossing counts. The self-intersection preflight BVH remains separate here to
preserve its established validation behavior.

Focused coverage is `make -C solvers/cpu surface_spatial_index_test`; it is
included in `make -C solvers/cpu test`.

## PR 5.2 prerequisite: canonical hashing

The dependency-free CPU layer includes a self-contained incremental SHA-256
utility with FIPS known-answer and chunk-boundary coverage.

## PR 5.2a2: canonical closed-surface topology and orientation

The dependency-free CPU layer now accepts an in-memory raw triangle soup and
constructs an immutable canonical closed triangulated surface in physical
metres.  It exactly deduplicates scaled coordinates, validates closed manifold
edge use, consistent local winding, a single face component, and one cyclic
incident-face fan at every vertex.  It rejects degeneracy, duplicate facets,
open/non-manifold/pinched topology, and numerically zero enclosed volume.
Inward shells are flipped consistently.  Canonical lexicographic vertices and
oriented triangles produce a versioned SHA-256 identity that is invariant to
input vertex/triangle ordering, cyclic triangle starts, and global orientation.
Equivalent source units have the same identity only when conversion produces
identical canonical physical double values (there is no coordinate
quantization); boundary labels are retained and hash-sensitive.  Diagnostics
expose counts, bounds, area, volume, edge/area
minima, topology counts, flip state, label histogram, and hash.

Focused coverage is `make -C solvers/cpu surface_geometry_test`; it is also
included in the dependency-free aggregate `make -C solvers/cpu test`.

Strict STL and bounded VTP readers are complete; classification and cut
quadrature remain pending work.

## PR 5.1: cubic Cartesian background

The CPU backend now has a dependency-light immutable cubic Cartesian background
grid.  It stores only the grid specification and exact per-axis knot-insertion
Bezier extraction operators, materializing elements lazily with x-fast ids and
connectivity.  Volume coordinates are metres, Bezier points are affine, and
all boundary labels remain unset for later classification/ownership work.
Cell neighbor slots are explicitly ordered `ZMinus, YMinus, XPlus, YPlus,
XMinus, ZPlus`; absent box-edge neighbors use `kNoNeighbor`.

Focused dependency-free coverage verifies one-cell identity ordering and affine
geometry, anisotropic element/node ids and 48-basis face sharing, plus invalid
bounds rejection.  Evidence: `make -C solvers/cpu cartesian_background_test`
and the focused binary.
