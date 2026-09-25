# Geometry data contract

This contract gives the centerline/body-fitted IGA, surface/immersed, and
surface/FEM routes the same vocabulary without pretending they use the same
mesh type or already expose one common runtime API.

## Required concepts

- `route` is explicit: `centerline_to_iga_volume`,
  `surface_to_immersed_background`, or `surface_to_fem_volume`. An organ name
  must never choose it implicitly.
- Coordinates are Cartesian. A manifest declares `length_unit`; coupling and
  solver-facing geometry uses metres.
- `reference_geometry.identity_sha256` identifies canonical reference
  geometry, independent of partition ownership or output path.
- `current_geometry` either equals the reference configuration or names a
  separately identified moving configuration and its physical time.
- `stable_ids` describes which IDs key vertices/nodes, elements/cells, and
  interfaces. MPI ownership and file ordering are not identities.
- Every region has an explicit role from `fluid`, `solid`, `shell`, `porous`,
  `network`, or `lumped`. Roles describe physics, not implementation status.
- Every boundary keeps its integer label and names the physical group or
  equivalent lookup. Wall, inlet, outlet, coupling, and constraint semantics
  remain explicit case configuration; they are not inferred from label order.
- The source content hash, applied coordinate scale/transform, preprocessing
  identity, and generated-artifact hash are provenance, not substitutes for
  the canonical reference identity.
- A surface conversion declares whether it exactly preserves the boundary,
  remeshes inside an envelope, or repairs anatomy. Its audit metrics include
  volume, area, volume/area-equivalent radii, triangulation-defined absolute
  curvature, and maximum source-surface distance; these global equivalent
  radii are comparison scalars, not anatomical lumen-radius inference.
- For an envelope-remeshed fTetWild volume, the source surface is retained as
  `source_geometry`, while `reference_geometry.identity_sha256` identifies the
  **output tetrahedral mesh**. Its fixed `current_geometry` equals that output
  reference, not the source surface. `stable_ids` names output MSH tags. This
  prevents a changed boundary from being misreported as the original surface.

## Reference and current configurations

Fixed geometry uses:

```json
"current_geometry": {
  "kind": "same_as_reference",
  "identity_sha256": "..."
}
```

Moving geometry must instead declare a current configuration identity, time in
seconds, predecessor/reference identity, and the mapping or displacement
artifact that produced it. Mesh velocity and fluid velocity are distinct
fields. A solver must reject a current surface, coupling layout, checkpoint, or
traction publication whose reference identity does not match its bound domain.

## Route mappings

| Contract concept | Body-fitted IGA | Immersed | FEM volume |
|---|---|---|---|
| Reference identity | packed/control geometry plus declared transform | `ClosedTriangulatedSurface::CanonicalSha256()` or moving material-map identity | canonical closed-surface SHA-256 |
| Volume discretization | Bézier-extracted spline volume | Cartesian cubic B-spline background plus cut classification | labelled P1 tetrahedra in the initial adapter |
| Stable surface IDs | packed boundary connectivity/labels | canonical surface IDs; material IDs for moving FSI | canonical surface MSH node/element tags |
| Current configuration | fixed in the existing body-fitted route | fixed or explicitly evaluated material motion | same as reference in the initial fixed-wall route |
| Region role | configured 3D fluid/transport domain | configured immersed fluid domain | `fluid` physical volume group |

The FEM manifest emitted by `surface_to_fem_volume.py` contains a
`geometry_contract` object following this definition. For packed body-fitted
IGA data, `iga_inspect DATABASE.ntiga --geometry-manifest FILE.json
--region-role fluid` emits the same top-level contract. Its reference identity
hashes the transform, stable element/node topology, labels, extraction, and
Bezier geometry, but deliberately excludes rank ownership; the local regression
checks that two differently partitioned databases retain one geometry identity
while their artifact hashes differ.

For immersed geometry, `ImmersedGeometryContractJson` serializes the immutable
material/reference identity, canonical surface, evaluated current identity and
time, runtime cut-geometry identity, Cartesian background, stable material and
background IDs, region role, and labelled boundary areas. Its focused test is
available through `make t1-immersed-manifest-test PETSC_DIR=/path/to/petsc`.
The caller still owns where this JSON is persisted alongside a run; this helper
does not change checkpoint or solver output formats.
