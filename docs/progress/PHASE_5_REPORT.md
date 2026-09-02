# Phase 5 report

Status: **complete**.

## PR 5.10: aggregate manufactured-closure evidence (2026-09-02)

The aggregate command `./solvers/cpu/phase5_closure_test --all` completed with
exit status 0 against system PETSc 3.15 (`PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real`;
linked `/lib/x86_64-linux-gnu/libpetsc_real.so.3.15`).  It exercises an
arbitrary closed immersed box on a compact cubic background through a
converged incompressible-flow solve.

| cells | h | dofs | uL2 | uH1 | pL2 | div L2 | pressure mean |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | 1 | 865 | 1.9851558135e-3 | 6.1261210319e-3 | 1.4589680951e-2 | 2.4672897644e-3 | 2.4678194537e-18 |
| 6 | 0.5 | 2917 | 4.1436317970e-5 | 2.2766384238e-4 | 3.4918893201e-4 | 1.3381586806e-4 | 9.2160914576e-19 |
| 12 | 0.25 | 8789 | 2.4734497479e-6 | 1.2401881371e-5 | 2.4039081405e-5 | 3.4648080904e-6 | -3.1080232377e-18 |

The `(uL2,uH1,pL2,div)` observed orders are `(5.5822126268,4.7499967844,
5.3847966928,4.2046060512)` from 3→6 and
`(4.0662993167,4.1982742807,3.8605541585,5.2713299463)` from 6→12.  At the
finest level, normalized errors are `(1.4657479987e-4,3.7218720498e-4,
1.4428518143e-3,1.0398077520e-4)`.

Independent high-order integration changed the prior accepted level-6 evidence
by `0.03840%`, `0.000404%`, `0.000713%`, and `0%` respectively.  With the
requested outer tolerance tightened 100×, the error norms changed by less than
`1.2e-5%`.  The achieved residual was not asserted to be 100× smaller; the
request was a tolerance change, not a residual-ratio claim.

The centered finite-difference diagnostic keeps `J(0)d` frozen, then compares
it with `-(b(+eps*d)-b(-eps*d))/(2 eps)` for the stored negative residual.  Its
full relative defect was `0.213162`, stable over `eps=1e-4..1e-7`; volume-only
was `2.35395`; wall+ghost-only converged to the expected roundoff/truncation
range (`8.9e-13` at `1e-4` through `9.45e-10` at `1e-7`).  This is consistent
with the frozen-VMS-tangent hypothesis, not a claim that the full Jacobian is
exact.  The nonzero-state contribution-decomposition assertion also passed.

The solver's block residual gates passed.  Initially nonzero momentum and
continuity blocks must fall by `1e9`; an initially zero gauge instead has the
scale-aware absolute bound `atol + rtol*||F0||`, so no `1e9` gauge-reduction
claim is made.

### Predeclared gates (all passed)

- Each manufactured solve requires convergence, expected volume/surface/ghost
  counts, `total_dofs=4*active_nodes+1`, pressure measure within
  `5e-8` relative of the analytic box measure, `|pmean|<1e-10`, global final
  residual no larger than `rtol*initial` (which must exceed `atol`), and each
  independent error-integration change below `2e-3`.  Each initially nonzero
  residual block must be no larger than its initial value divided by `1e9`; an
  initially zero block must be no larger than `atol + rtol*initial`.
- The level triplet requires strictly decreasing four error norms for both
  pairs, strictly increasing dofs, 6→12 order floors `(2.5,2.0,2.0,1.5)`, and
  level-12 normalized-error caps `(0.02,0.05,0.05,0.05)`.
- The nonzero-state decomposition requires each matrix and residual difference
  to be below `2e-11` for `full_ungauged=volume+wall+ghost` and
  `full_gauged=full_ungauged+gauge`.
- Every sliver case requires the expected minimum fraction within `5e-5`
  relative, a covered nonempty ghost catalog, analytic gauge measure within
  `5e-8` relative, exact repeated-assembly/action agreement, a finite nonzero
  mixed spectrum, sign-flipped symmetry defect at most `2e-10`, and normalized
  constant-pressure defect below `2e-10`; aggregate condition growth must be
  finite and below `100`.

The sliver sequence had `dofs=865`, `faces=54`, `records=27`, and
`logical=1728` in every case; the gauge measure equalled the analytic box
measure.  Its results were:

| m | alpha | normalized pressure defect | min abs eig | condition |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 0.3828125 | 1.4847756290e-16 | 0.0025999412605 | 3248.396953 |
| 2 | 0.19140625 | 1.1682996957e-16 | 0.00078432965147 | 10754.912455 |
| 3 | 0.095703125 | 1.3122712676e-16 | 0.00043263506790 | 19521.898638 |
| 4 | 0.0478515625 | 1.1660008727e-16 | 0.00031946694595 | 26456.657277 |
| 5 | 0.02392578125 | 1.1659537158e-16 | 0.00027401546389 | 30846.322401 |
| 6 | 0.011962890625 | 1.1660958035e-16 | 0.00025365318610 | 33318.483234 |
| 7 | 0.0059814453125 | 1.1662177314e-16 | 0.00024401712163 | 34630.585697 |
| 8 | 0.00299072265625 | 1.1896193476e-16 | 0.00023932999445 | 35306.537317 |

The maximum eigenvalue was approximately `8.435`–`8.452`, all symmetry defects
were at most `2.04e-16`, and the aggregate condition-growth ratio was
`10.868910982`.

Phase 5 is complete: arbitrary closed immersed box → compact cubic background
→ converged incompressible flow.  The phase-closing Sol review returned **GO**
with no blocking findings.

## PR 5.9: serial global immersed static flow

`ImmersedStaticFlowRuntime.hpp` provides a deliberately serial PETSc global
assembly/runtime over the Phase 5 catalog APIs.  It certifies every classified
Inside/Cut volume rule before eligibility filtering (only a usable,
certified-empty Cut cell may be skipped), then compacts nodes used by
positive-volume cells into sorted dense four-field rows, checks every
cell/face mapping, and appends one pressure-gauge scalar.  Cells,
surface cells, and ghost faces are scattered in their catalog order using only
the existing bounded element/face blocks.  Expanded rules are consumed directly;
Compact rules stream through `ForEachVolumePoint(UsableCompactRule(...))` and
never form a logical-point temporary.  The Nitsche adapter accepts a supplied
streamed volume base before adding its surface terms, and ghost residuals use a
bounded face-local state accessor rather than a background-sized expansion.
Cut walls require the matching,
covered ghost catalog; mismatched bindings and uncovered cut cells fail before
assembly.  One deterministic lazy volume pass after catalog preflight and
active compaction caches the immutable gauge weight per active pressure node;
the gauge is `R=[F+g lambda; g^T p]`, with positive
matrix blocks and negative-residual signs in both pressure and gauge rows.
The constant-pressure defect is measured before gauge insertion.

The static Newton trial starts from an immutable committed vector, requires a
finite converged PETSc KSP step, and accepts a damped update only if the
assembled residual strictly decreases (halving through `1/128`), including the
candidate after the final allowed Newton update.  Exceptions,
KSP failure, and iteration/backtracking caps restore the exact committed vector;
only `Commit()` publishes a successful trial.  Diagnostics distinguish trial
and committed state and report active rows, volume/wall/ghost counts, gauge
measure/defect, nonlinear/KSP status, damping, commits, and rollbacks.

`NavierStokesElement.hpp` now has a backward-compatible physical body-force
density (N/m^3) evaluator overload.  It contributes to both the weak momentum residual and
the strong residual driving VMS/PSPG; its tangent remains zero for a
state-independent callback.  Focused coverage is
`make -C solvers/cpu immersed_static_flow_test PETSC_DIR=…`.

## PR 5.6: static immersed Nitsche wall

`ImmersedNitscheWall.hpp` is a catalog-bound element adapter: it obtains both
the volume and immersed-surface rules only through their bound `UsableRule`
interfaces, then evaluates cubic basis values and physical gradients at each
immersed (not face) parametric point.  The static wall uses
`sigma=-p I + mu (grad u + grad u^T)`, outward surface normals, and the
symmetric Nitsche residual/Jacobian with a prescribed velocity `g` that has no
degrees of freedom.  Surface quadrature weights are already physical `dA` and
are used directly.

For cubic degree `p=3`, the penalty is
`eta = 16 gamma0 mu / (alpha h_n)`, where `h_n=1/||J^{-1}n||` and `alpha` is
the catalog's validated estimated reference-volume fraction (the matching
physical fraction for an affine Cartesian cell).  The adapter rejects invalid
binding, labels, gamma, fraction, length, penalty, and wall velocity, and
publishes per-label selected/skipped area and point counts plus fraction,
length, penalty, and gap diagnostics.  It deliberately neither floors
`alpha` nor caps `eta`: PR 5.7 owns cut-independent conditioning and sliver
stabilization.  Focused evidence is `make -C solvers/cpu nitsche_wall_test`
with PETSc configured; it is included in `petsc-test`.

## PR 5.5: immersed surface quadrature

`ImmersedSurfaceQuadratureCatalog` stores immutable, x-fast physical boundary
rules bound to the exact Cartesian grid and canonical-surface hash.  It clips
canonical oriented triangles to closed cells with constraint-provenance
homogeneous dyadics (never recursive rational intersections),
uses deterministic fan triangulation and a positive degree-six 12-point rule,
and audits canonical triangle, boundary-label, and whole-surface emitted areas
against compensated canonical source sums.  Diagnostics publish source,
accumulated, signed residual, and absolute residual measures; tolerances scale
with measure and operation count without a one-square-metre floor.  Internal
coplanar faces have one normal-directed owner.  Cell construction is
transactional through rule validation: any non-cap clipping, representation,
mapping, normal, or weight failure discards the cell's retained rule and area
state and globally latches the catalog unusable.  Candidate, fragment, and
point attempt counters remain monotone across rejected cells, so the three
caps cannot be evaded by rollback; cap exceptions still prevent publication.
The clipping arithmetic has its own explicit per-catalog limb cap and observed
peak diagnostic; general exact predicates retain their smaller independent
capacity.  Predicate capacity or unrepresentable positive fragments make a
cell unusable rather than silently changing its surface measure.  Focused evidence is
`make -C solvers/cpu immersed_surface_quadrature_test`.

## PR 5.4: adaptive cut-cell volume quadrature

`CutCellVolumeQuadratureCatalog` retains one immutable, cell-id-sorted result
for every Cartesian background cell without retaining surface geometry.  Inside
cells reuse the existing 4x4x4 reference rule, outside cells are certified
empty, and cut cells use a deterministic x-fast octree.  Certified leaves add
scaled positive Gauss rules; unresolved leaves use point classification only
and report lower/upper reference and physical-volume bounds.  Predicate
ambiguity is retained for inspection but makes a rule unusable for assembly.
Rules require an exact grid and canonical-surface binding before assembly;
physical diagnostic estimates sum stored reference weights times the materialized
element Jacobian, while assembly itself continues to own its physical measure.

Focused evidence: `make -C solvers/cpu cut_cell_volume_quadrature_test`.
The target is included in the dependency-free aggregate CPU test.  Surface
quadrature remains explicitly deferred to PR 5.5.

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

## PR 5.3b: Cartesian domain catalog

`CartesianDomainClassification` retains every background cell id and
connectivity while recording x-fast Outside, Inside, or conservative Cut
classification, canonical intersecting surface triangle ids, and ambiguity
diagnostics. Surface contact always wins over point classification. Focused
coverage is `make -C solvers/cpu cartesian_domain_classification_test`.

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

## PR 5.7: Cartesian cubic ghost penalty

`CutCellGhostPenaltyCatalog` is an immutable, exact-domain/volume-bound
catalog for uniform Cartesian degree-three, C2 B-splines.  It selects each
active lower-to-plus neighbor face at most once in x/y/z cell order when at
least one usable positive-volume Cut cell is present.  The static form is
\(s_u=\gamma_u\mu\sum_F h_F^5\int_F[\partial_n^3u]\cdot[\partial_n^3v]\)
and \(s_p=\gamma_p/\mu\sum_F h_F^7\int_F[\partial_n^3p][\partial_n^3q]\),
with positive default gammas 0.01.  There are deliberately no mass or
convection ghost terms.

The face trace uses current Bezier extraction, physical `h_axis^-order`
derivatives, the cubic Bernstein third derivative `{-6,18,-18,6}`, a
same-`+axis` lower-minus-upper jump, and a 4x4 Gauss tangent rule.  Assembly
is deliberately bounded and face-local: `AssembleFace(face_index, ...)`
returns only that face's sorted unique adjacent-cell union (at most 80
background nodes), dense four-field pair tangent, and `-A*x` residual.  There
is no domain-wide dense aggregate API.  `EvaluateFaceJump` exposes the same
bounded trace coefficients for derivative orders zero through three.

Catalog binding is O(1): it records immutable caller-owned domain and volume
object identities plus construction grid/hash/count/options metadata.  Those
objects must remain alive and immutable for the catalog lifetime; validation
compares identities without dereferencing stored pointers.  Diagnostics
precompute nonzero normalized (`mu=1`) trace-jump and coefficient extrema;
selected-face catalogs reject non-representable normalized coefficients before
publication.
A usable positive-volume classified-Cut cell selects faces even when its
volume is certified full; certification only makes it a coverage anchor.
The catalog audits selected-face reachability to a certified-full/Inside
anchor and publishes uncovered cells and components.

The wall adapter now has an explicit catalog-covered overload.  Only a
matching, covered catalog cell may use it; its penalty is
`eta=16*gamma0*mu/h_n`, independent of the cut fraction (the legacy overload
is unchanged at `16*gamma0*mu/(alpha*h_n)`).  Diagnostics state which policy
was used and retain observed alpha bounds.

Focused PETSc coverage is `make -C solvers/cpu ghost_penalty_test PETSC_DIR=…`
and is part of `petsc-test`. It verifies deterministic unique lower-to-plus
face ownership, bounded connectivity, x/y/z trace continuity through order
two and nonzero third jumps, nonzero velocity and pressure blocks, relative
symmetry, Jacobi PSD spectra, constant-pressure nullspace, and central
finite-difference residual linearization. It also checks independent gamma
and viscosity factors, uniform `K_u~s` and `K_p~s^3` scaling, anisotropic
`area/h_n` and `area*h_n` scaling, pure classified-Cut selection semantics,
coverage, caps, nonfinite input, and exact object-identity rejection.

Real retained-rule dyadic sliver evidence currently covers `m=1..4` at octree
depth four for the stated 3x3x3 box, including the selected certified-full
anchor face, nonzero velocity/pressure ghost blocks, covered-wall
`eta*h_n/mu=32`, legacy `32/alpha`, and the diagnostics policy flag. The
requested `m=1..8` at depth eight is not claimed: with the present global
retained 4x4x4 rule catalog it reaches the explicit point cap before
publication (raising it entails tens of millions of retained points).
Likewise, the existing Cartesian domain contract rejects a closed surface
whose bounds enclose the background, so an all-Inside exterior-surface
catalog cannot be formed without changing PR5.3/domain semantics. These are
precise bounded-integration blockers, not substituted evidence.

The remaining Phase 5 closure blockers are the depth-eight `m=1..8`
sliver/global-patch spectra, integration of this pair operator into the
distributed global flow/PSPG assembly, and three-level manufactured Stokes
convergence evidence; none is fabricated by this bounded element/catalog PR.

## PR 5.8: compact cut-volume retention

`CutCellVolumeQuadratureCatalog` now has an explicit storage mode.  Expanded
rules remain the default and preserve the existing point-vector API and
ordering.  Compact mode stores certified material as exact integer boxes on a
`2^max_depth` reference lattice and retains unresolved terminal leaves as a
dyadic key/depth plus a 64-bit inside-sample mask.  It never expands those
records during catalog construction; `ForEachVolumePoint` emits the same
qz/qy/qx 4x4x4 Gauss rule on demand.

The compact builder uses deterministic x-fast octree recursion.  Uniform
children collapse to their parent, while mixed terminal leaves remain sample
records even for all-zero or all-one masks.  Certified boxes then undergo
fixed x/y/z face-adjacent coalescing to a fixed point.  The catalogue exposes
logical emitted-point, final record, monotone record-attempt, rollback-record,
and retained-byte diagnostics separately from the legacy expanded
`output_points` count.  Logical-point, record, and byte caps are checked
before publication.  Bound callers must use
`UsableCompactRule`; the expanded and compact accessors intentionally reject
the opposite storage mode.

`max_points` deliberately retains its legacy expanded-storage meaning.  Compact
rules instead apply the separately configurable, checked `max_logical_points`
cap after block/mask compression (default 64 Mi logical points); node, leaf,
sample, record, and retained-byte caps still bound construction work.  The
compact recursion appends into per-cell shared scratch in Morton child order;
all-full/all-empty rollback releases speculative records, is counted
separately, and cannot bypass the monotone record cap.  Scratch vectors use
deterministic 1.5x planned capacity growth.  `retained_bytes` is the portable,
checked peak accounting model for those planned capacities, including old plus
replacement plans during a reserve.  `std::vector::reserve` may overallocate,
so `max_retained_bytes` deliberately does not claim a strict allocator-byte or
transient-allocation bound.  `observed_retained_bytes` reports only final
published vector capacity for diagnosis; it is not capped and is not an
allocator peak.  Final records and logical emissions are checked before
publication.

Compact validation audits lattice keys, strict block/Morton ordering, sample
prefix intervals, block BVH overlap, every lazily emitted reference point and
weight, and (at bound access) every mapped physical Jacobian and physical
diagnostic sum without materializing an expanded vector.

Focused evidence is build *and* run:
`make -C solvers/cpu cut_cell_volume_quadrature_test && ./solvers/cpu/cut_cell_volume_quadrature_test`.
It covers all 512 degree-0..7 tensor moments for shallow cube/tetra compact
versus expanded rules with operation-scaled tolerances, depth-zero mask-order
bitwise emission, repeat/permutation/translation/scale compact records,
coalescing plus a deterministic shared-helper rollback probe with exact and
cap-minus-one monotone append checks, exact planned-byte-cap boundaries,
transactional work-cap rejection,
dyadic 3^3 slivers through depth eight, and oblique tetrahedral depth-seven/
depth-eight brackets.  It does not claim bitwise equivalence for unresolved
non-polynomial rules beyond the depth-zero mask-order evidence; the existing
deterministic expanded regression suite remains covered.

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
