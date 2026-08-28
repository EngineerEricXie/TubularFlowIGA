# Pipeline and Case Preparation

## 1. Generate the control mesh

A 3D case directory starts with a schema-v4 `simulation_config.json`. Its
`geometry` block names an SWC or radius-annotated OBJ, and its `mesh` block
contains smoothing, adaptive centerline, junction-clearance, and quality
controls. The old five-line `mesh_parameter.txt` remains accepted only for
legacy v2/v3 cases; it may not coexist with schema v4.

```json
{
  "schema_version": 4,
  "dimension": "3d",
  "geometry": {
    "kind": "swc_network",
    "file": "skeleton_initial.swc",
    "length_scale_to_m": 0.001
  },
  "mesh": {
    "smoothing": {"iterations": 0, "bifurcation_ratio": 0.2, "noise_ratio": 0.0},
    "centerline": {
      "target_spacing": 1.0,
      "max_spacing_over_diameter": 1.0,
      "max_turn_degrees": 12.0,
      "max_diameter_change_fraction": 0.15,
      "maximum_curvature_radius_product": 0.8
    },
    "junction": {
      "max_spacing_over_diameter": 0.25,
      "upstream_clearance_over_diameter": 1.0,
      "downstream_clearance_over_diameter": 1.5,
      "minimum_angle_degrees": 10.0,
      "maximum_radius_ratio": 8.0,
      "optimization_iterations": 4
    },
    "quality": {
      "minimum_scaled_jacobian": 0.1,
      "check_self_intersection": true,
      "collision_safety_factor": 1.0
    }
  }
}
```

The excerpt omits the unchanged fields, systems, time, and boundary blocks.

Build and test the standalone C++ generator from the repository root:

```bash
make mesh
make mesh-test
./preprocessing/mesh/tubular_mesh pipeline \
  "$CASE_DIR" meshgeneration/template
```

After strict parsing and topology validation, it writes
`skeleton_normalized.swc` and `skeleton.vtp`. OBJ is thereby converted to an
explicitly rooted SWC before smoothing continues. Valid SWC follows the same
normalization path. The original input is never overwritten.

The pipeline then writes `skeleton_smooth.swc`, `mesh_diagnostics.json`,
`skeleton_diagnostics.vtp`, `controlmesh.vtk`, `mesh_quality.json`, and
`initial_velocityfield.txt`. The implementation evaluates its own cubic
B-splines, constructs tube and bifurcation hexahedra, and has no MATLAB, TREES,
Eigen, or VTK-library dependency. A 3D run requires a tree whose nodes have at
most two children after root orientation. It stops before producing a mesh and
reports the offending node ID and child count if this constraint is violated.

The smoothed SWC is intentionally written to eight decimal places and read
back before meshing. This preserves the legacy file-interface behavior at
layer-count boundaries.

MATLAB remains as an optional reference workflow. Install TREES separately,
add both TREES and this repository recursively to the MATLAB path, set
`io_path` near the top of `TreeSmooth.m` and `Hexmesh_main.m`, then run:

```matlab
cd('/path/to/TubularFlowIGA/meshgeneration')
addpath(genpath('/path/to/TREES'))
TreeSmooth
Hexmesh_main
```

Both implementations retain the control-mesh file interfaces. The C++ pipeline
additionally rejects excessive curvature-radius product, diameter transitions,
inadequate bifurcation clearance, unsupported angle/radius ratios, nonpositive
Jacobians, and exterior-surface self-intersections. It uses rotation-minimizing
frames, adaptive arc-length resampling, junction quality optimization, and
corners plus 4x4x4 interior samples. Inspect the diagnostic VTP in ParaView as
a second gate. A positive determinant is mandatory; schema v4 makes `0.1` the
explicit public-case scaled-Jacobian floor.

## 2. Extract the IGA representation

Build with `make spline`. The argument is a directory prefix, so include its
trailing slash:

```bash
OMP_NUM_THREADS=8 ./preprocessing/spline/spline \
  "$CASE_DIR/" --no-legacy-text
```

Set the thread count to the allocated CPU cores. The extractor processes
elements in bounded-memory chunks and formats element records in parallel while
preserving deterministic legacy text output.

It reads `controlmesh.vtk` and writes:

- `bzmeshinfo.txt`: Bezier element connectivity;
- `spline_cache.igacache`: versioned sparse coefficients and Bezier points;
- `bzmesh.vtk`: legacy preprocessing visualization output (its element-local
  points may repeat; production temporal visualization uses the packed
  database and the extraction-signature registry instead);
- `geometry_transform.json`: source origin and normalization scale.

### Coordinate normalization

Before extraction, the spline code subtracts the minimum coordinate on each
axis and divides every coordinate by the smallest domain-axis extent. The
Bezier mesh and packed `.ntiga` geometry therefore use translated, normalized
coordinates rather than the original SWC coordinate units. The extractor writes
the affine map to `geometry_transform.json`; version-5 `.ntiga` also stores it
and the configured source length scale to metres. Transform flow, time,
material, pressure, and transport parameters consistently. The public examples
use internally consistent numerical values and are not presented as
patient-specific SI calibrations.

Omit `--no-legacy-text` to additionally reproduce `cmat.txt` and `bzpt.txt`.
The binary cache records a control-mesh content hash; `iga_pack` rejects stale,
truncated, malformed, or version-incompatible caches.

The production solvers preserve `controlmesh.vtk` as this pipeline interface.
Transient results default to the deduplicated Bézier representation described
in the [visualization output guide](VISUALIZATION.md).

## 3. Partition and pack

CPU rank count, METIS partition suffix, and `mpiexec -np` must agree.

```bash
RANKS=8
mpmetis "$CASE_DIR/bzmeshinfo.txt" "$RANKS"
./solvers/cpu/iga_pack "$CASE_DIR" "$RANKS" "$DATABASE"
./solvers/cpu/iga_inspect "$DATABASE"
```

The packer prefers `spline_cache.igacache` and falls back to legacy text when
the cache is absent. Pass `--legacy-text` after the output path to force that
fallback for regression testing. It creates a binary database with direct
element offsets, sparse extraction rows, adjacency, per-rank touching-element
indices, and six boundary-face labels per element. Version 4 introduced derived
face labels; version 5 additionally stores the source-to-normalized geometry
transform and configured length scale. Readers remain compatible with versions
3 and 4; version 3 face labels are unavailable. Repack
when changing the CPU rank count. CUDA ignores ownership records and may reuse
any valid packed database.

## 4. Validate and solve

Use a workstation only for cases that fit its resources. On a shared cluster,
request a compute allocation through its scheduler. Validate with
`iga_mesh_check` or `iga_cuda mesh-check` before timing either physics
stage. CPU runs use `mpiexec`; CUDA execution is one process on one GPU.

Platform-specific examples:

- [PSC Bridges-2](BRIDGES2.md)

Record the case snapshot, compiler, MPI/PETSc or CUDA versions, rank count or
GPU, scheduler job ID when applicable, solver tolerances, numerical norms,
stage timers, host RSS, and CUDA peak allocation. Compare output fields, not
runtime alone.
