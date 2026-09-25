# T1 local surface-to-FEM volume mesh route

The initial external-mesher adapter combines the dependency-free
`solvers/cpu/surface_fem_preflight` utility with
`scripts/surface_to_fem_volume.py`. The C++ stage reads and canonicalizes a
closed triangular VTP or STL surface using the same geometry contract as the
immersed backend. The Python/Gmsh stage adds a first-order tetrahedral volume
and publishes a provenance and quality manifest. It does not solve a fluid
problem and does not repair or smooth the input anatomy.

## Usage

```bash
python3 scripts/surface_to_fem_volume.py \
  examples/vascular_flow/immersed_aneurysm_chain/immersed/surface.vtp \
  /tmp/aneurysm-fem.msh \
  --manifest /tmp/aneurysm-fem.json \
  --target-size-m 0.12
```

Run the local regression with `make t1-fem-mesh-test`.

For large surfaces, the complete C++ preflight may be run once and its
canonical MSH/JSON reused by later Gmsh target-size trials. The reuse path
does not omit a gate: it checks SHA-256 of the current source and canonical
MSH, all preflight options, and the six passing geometry flags before meshing.
The two preflight artifacts must come from the current utility, whose JSON
records `canonical_msh_sha256`; older manifests without that field fail
closed. For example:

```bash
solvers/cpu/surface_fem_preflight input.stl /tmp/validated-surface.msh \
  --manifest /tmp/validated-surface.json --default-boundary-id 1 \
  --length-scale-to-m 0.001 --max-triangles 450000
python3 scripts/surface_to_fem_volume.py input.stl /tmp/volume.msh \
  --manifest /tmp/volume.json --target-size-m 0.005 \
  --default-boundary-id 1 --length-scale-to-m 0.001 --max-triangles 450000 \
  --validated-surface-msh /tmp/validated-surface.msh \
  --validated-surface-manifest /tmp/validated-surface.json
```

The output manifest records `surface_preflight_reused`; source, mesh, option,
or label mismatch is a hard failure, never a reason to skip quality checks.

The adapter requires the Gmsh Python module. The WSL installation used for the
initial gate reports Gmsh 4.8.4 and a GNU GPL license; the repository does not
vendor or redistribute it. PSC module/container availability and deployment
policy remain a separate check before cluster use.

## Fail-closed contract

The VTP route accepts the bounded, uncompressed XML PolyData subset already
supported by `SurfaceReaders`: ASCII, inline base64, and base64-appended arrays,
with triangular `Polys` and a nonnegative integer `boundary_id` per triangle.
ASCII and binary STL are also accepted and receive the configured default
boundary ID. Compressed, multi-piece, polygonal, missing-label, nonfinite, or
invalid-connectivity input is rejected.

`ClosedTriangulatedSurface` performs the shared duplicate-facet, closed-edge,
orientation, vertex-fan, connected-component, self-intersection, positive-volume,
unit scaling, welding, canonical ordering, and canonical-hash checks before
Gmsh sees the surface.

Each input boundary ID becomes a named two-dimensional physical group
`boundary_label_N`; the tetrahedral region is the `fluid` physical group. Gmsh
is instructed to mesh only the empty volume, preserving the validated input
triangles rather than remeshing the anatomical surface.

Before publishing output, the adapter checks:

- the output contains only first-order tetrahedra;
- every tetrahedron has a positive determinant and meets the configured
  minimum corner-scaled Jacobian;
- each boundary label has exactly the original number of triangles;
- the labelled surface facet set exactly equals the tetrahedral boundary;
- summed tetrahedron volume agrees with the closed input surface volume within
  the configured relative tolerance.

The JSON manifest embeds the canonical surface preflight, including input and
canonical SHA-256, label histogram and per-label areas, total area, bounds,
volume- and area-equivalent sphere radii, and triangulation-defined absolute
mean-curvature measures (half the sum of edge length times absolute dihedral),
source/mesh volumes, Gmsh version, target size,
node/element counts, determinant, scaled Jacobian, and thresholds. Generated `.msh` files are case data and remain
untracked according to repository policy.

The embedded `geometry_contract` follows
[`GEOMETRY_DATA_CONTRACT.md`](GEOMETRY_DATA_CONTRACT.md): it declares SI metre
coordinates, reference/current identity, stable ID namespaces, the `fluid`
region role, and explicit boundary groups without guessing wall/port semantics.

## fTetWild route

[fTetWild](https://github.com/wildmeshing/fTetWild) is now an implemented robust
tetrahedralization route, with the current Gmsh adapter retained as the
reference and fallback. This does not change the solver architecture:
fTetWild is a preprocessing executable only; basis functions, quadrature,
residual/Jacobian assembly, boundary conditions, Newton iteration, and QoIs
remain repository-owned native FEM code.

The upstream CLI accepts OFF/OBJ/STL/PLY triangle surfaces and writes MSH or
MESH tetrahedral meshes. The adapter exposes absolute ideal edge length
(`--la`) and relative envelope (`--epsr`), record all quality/optimization
settings, threads, declared version/commit, executable hash, command line,
timings, and peak memory, then
normalize output to the canonical ASCII Gmsh 4.1 interface consumed by the
native solver. Upstream is MPL-2.0 licensed. The repository will invoke a
separately installed executable rather than silently vendoring it.

Because fTetWild may change the boundary triangulation inside its envelope,
labels cannot be copied by facet identity. The adapter projects each output
boundary facet to the validated source surface and requires an
unambiguous source label, distance within the configured envelope, and
consistent orientation. Ambiguous labels, missing/extraneous boundary facets,
nonmanifold output, wrong regions, or a failed volume/Jacobian gate are hard
failures. No automatic hole filling or smoothing is accepted without
being explicitly classified as geometry repair and reporting the resulting
surface deviation.

Both adapters now emit an explicit `geometry_change` record. The Gmsh path
classifies its boundary as exact input preservation and records zero area,
area-radius, curvature, and point-distance change. The fTetWild path classifies
its operation as envelope-constrained boundary remeshing—not repair—and records
volume, surface area, both equivalent radii, discrete-curvature change, and
maximum source-surface distance. Equivalent radii are global audit metrics,
not a claim that a general anatomy is spherical or has one physical lumen
radius.

The comparison card covers straight, bent, bifurcating, and deliberately
imperfect surfaces. It records success/failure class, volume and surface
errors, minimum scaled Jacobian, element count, wall time, and peak RSS for
both fTetWild and Gmsh. Upstream notes that difficult inputs can require more
than 32 GB, so those cases belong on an allocated compute resource rather than
the local WSL mandatory gate.

Run the dependency-independent adapter regression with
`make t1-ftetwild-mesh-test`; it uses a CLI-compatible test double and sends the
normalized result through the native C++ tetrahedral reader. A real local run
was also completed against upstream commit
`d7d99bb4387a07895b9adce058dc7305f6b6e5ab`. It was configured with upstream's
required oneTBB support and GMP from an existing Conda environment; the source
and build remained under `/tmp` and were not vendored.

The first same-input comparison used a 0.1 m by 0.005 m labelled pipe, 0.003 m
target size, and the same preflight surface. fTetWild used a 0.0001 m envelope,
stop energy 12, 40 optimization passes, and two threads. Times include a
separate mesher measurement and the complete adapter pipeline:

| mesher | vertices | tetrahedra | minimum scaled Jacobian | relative volume error | mesher wall (s) | pipeline wall (s) | peak RSS (bytes) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Gmsh 4.8.4 | 322 | 1,176 | 0.001016 | 1.32e-15 | 0.063 | 0.778 | 99,164,160 |
| fTetWild `d7d99bb` | 1,021 | 3,691 | 0.07517 | 0.00332 | 0.892 | 2.536 | 25,174,016 |

The fTetWild boundary stayed within 4.21e-5 m of the source, retained all
wall/inlet/outlet labels, and passed both the native reader and a Gmsh
round-trip. Two further same-input comparisons used the same 40-pass policy:

| case / mesher | vertices | tetrahedra | minimum scaled Jacobian | relative volume error | mesher wall (s) | pipeline wall (s) | peak RSS (bytes) |
|---|---:|---:|---:|---:|---:|---:|---:|
| 90-degree bend / Gmsh | 557 | 1,858 | 0.07841 | 8.73e-16 | 0.097 | 1.922 | 99,090,432 |
| 90-degree bend / fTetWild | 622 | 2,237 | 0.06322 | 0.00133 | 0.316 | 2.660 | 22,245,376 |
| Y bifurcation / Gmsh | 822 | 2,536 | 0.01694 | 4.91e-16 | 0.081 | 3.662 | 99,602,432 |
| Y bifurcation / fTetWild | 732 | 2,337 | 0.10340 | 0.00291 | 0.567 | 5.169 | 29,200,384 |

The fTetWild maximum source distances were 3.47e-5 m for the bend and
3.20e-5 m for the bifurcation, below their configured envelopes. The results
show why selection remains per case: Gmsh had the better minimum quality on
the bend, while fTetWild was substantially better on the straight and
bifurcating fixtures. The regression also sends one shared surface with a
missing triangle through both adapters; both stop in the common closed-surface
preflight, so fTetWild's repair capability cannot silently change anatomy.
The exact local measurements and hashes are retained in
`benchmarks/t1_mesher_comparison.json` and fail-closed checked by
`scripts/validate_t1_mesher_comparison.py`.

An installed executable is invoked as follows:

```bash
python3 scripts/ftetwild_to_fem_volume.py input.vtp output.msh \
  --ftetwild /path/to/FloatTetwild_bin \
  --target-size-m 0.003 --envelope-m 0.0001 \
  --stop-energy 12 --max-optimization-passes 40 --max-threads 2
```

## Remaining T1 scope

The solver-facing surface adapters still create one fluid volume, but a
separate native-contract audit now validates externally prepared multi-region
first-order tetrahedral meshes:

```bash
python3 scripts/audit_multiregion_tet_mesh.py mesh.msh regions.json \
  --output mesh-audit.json --minimum-scaled-jacobian 1e-3
```

Its sidecar explicitly lists every volume physical name and role, every
exterior physical surface and adjacent region, and every interface physical
surface and unordered region pair. The audit requires positive determinants,
the configured scaled-Jacobian floor, exactly one region per tetrahedron,
manifold face incidence, complete exterior-boundary coverage, exact
cross-region interface coverage, no labelled same-region interior faces, and
optional connectedness per region. It emits input hashes, per-region volumes
and element counts, boundary/interface facet counts, and gate status. The
focused regression includes a valid fluid/solid pair and rejects an inverted
tetrahedron, a mislabelled interface, a boundary attached to the wrong region,
and binary Gmsh input. This audit establishes data consistency; it does not
claim the current flow solver can assemble multiple physics regions.

This local slice does not provide automatic capping or surface repair,
multi-physics assembly over multiple material volumes, a common manifest for all
three geometry frontends (the shared contract exists, but persistence remains
route-owned), chamber FEM examples, or PSC installation evidence. Any future
surface repair must be a separate auditable stage using the same area, volume,
equivalent-radius, discrete-curvature, and surface-distance delta vocabulary;
it must never be hidden inside this converter.
