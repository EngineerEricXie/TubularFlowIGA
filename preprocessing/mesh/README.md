# C++ Control-Mesh Generator

This directory replaces the MATLAB/TREES preprocessing stage with a standalone
C++17 program. It reads an SWC or radius-annotated line-OBJ centerline, smooths
and resamples each branch, and writes the hexahedral control mesh consumed by
the existing spline stage. All spline evaluation, frame transport, bifurcation
assembly, element-quality checks, and VTK output are implemented here without a
high-level geometry library. See the repository
[skeleton-format contract](../../docs/SKELETON_FORMATS.md) for the exact input
records and validation rules.

## Build and Test

From the repository root:

```bash
make mesh
make mesh-test
```

Override `CXX` and `CXXFLAGS` when needed. Older GCC installations may require
the default `-lstdc++fs`; remove it through `LDLIBS=` on toolchains where the
filesystem library is built in.

## Run a Case

Create a case directory containing the geometry named by a schema-v4
`simulation_config.json`. Its `geometry` block selects the SWC or radius-annotated
OBJ and its `mesh` block contains all preprocessing controls. Then run from the
repository root:

```bash
./preprocessing/mesh/tubular_mesh pipeline \
  /path/to/case meshgeneration/template
```

The command writes:

- `skeleton_normalized.swc`: strictly validated, explicitly rooted skeleton;
- `skeleton.vtp`: ParaView line data with radius and topology arrays;
- `skeleton_smooth.swc`: smoothed and resampled centerline;
- `mesh_diagnostics.json`: effective limits, per-segment and per-junction metrics,
  warnings, and errors;
- `skeleton_diagnostics.vtp`: the same risk metrics as ParaView point/cell arrays;
- `controlmesh.vtk`: labeled eight-node control elements;
- `cross_section_template.vtk` and `merge_template.vtk`: radius-1 quad template
  previews saved directly alongside the control mesh on every run;
- `mesh_quality.json`: final determinant, scaled-Jacobian, and surface-intersection gates;
- `initial_velocityfield.txt`: branch-aligned initial velocities.

OBJ input is accepted only when it follows the documented radius-annotated
line convention. Surface faces and other OBJ records are rejected. During 3D
preparation, a node with three or more oriented children is rejected before
mesh generation, with its node ID and child count in the error.

Individual stages are also available:

```bash
tubular_mesh smooth INPUT.swc|INPUT.obj mesh_parameter.txt OUTPUT.swc
tubular_mesh generate SMOOTH.swc mesh_parameter.txt TEMPLATE_DIR \
  controlmesh.vtk initial_velocityfield.txt MIN_SCALED_J
tubular_mesh pipeline CASE_DIR TEMPLATE_DIR --allow-preflight-failure
```

`generate` also writes `mesh_quality.json`, `cross_section_template.vtk`, and
`merge_template.vtk` beside `controlmesh.vtk`.
The pipeline override is intended only for explicitly reviewed debug geometry.
It preserves the failed diagnostics and still enforces the final element-quality
and surface-intersection gates.

These two commands retain the strictly parsed five-line legacy format for
regression work. A pipeline case may use either schema v4 or legacy
`mesh_parameter.txt`, never both. Public cases use schema v4.

## Parameters and Assumptions

The schema-v4 `mesh` block separates smoothing, centerline, junction, and
quality controls. Centerline resampling is by approximate arc length and limits
spacing by target length, local diameter, tangent rotation, and fractional
diameter change. Junction controls reserve explicit upstream/downstream
clearance and reject unsupported angles or radius ratios with the SWC node ID.
SWC column 6 remains a radius and is converted internally to diameter.

The topology must be a connected rooted tree. Each nonterminal node must have
one child or exactly two children; higher-order junctions are not yet
supported. Without `mesh.cross_section`, the generator uses the original 201-point tube
and 294-point merge templates under `meshgeneration/template/`.

## Automatic Circular Templates

To generate a circle and its matching bifurcation template, add this optional
entry to the existing schema-v4 `mesh` block in `simulation_config.json`:

```json
"cross_section": {"kind": "circle", "target_size": 0.25}
```

To explicitly use the saved circular and merge templates in
`meshgeneration/template/`, use this entry instead:

```json
"cross_section": {"kind": "default"}
```

Omitting `mesh.cross_section` also selects these defaults. `target_size` is only
accepted with `kind: "circle"`; default mode uses the saved mesh resolution.

`target_size` is an approximate quad edge length on a **radius-1 reference
disk**, in the range `[1/128, 1]`. At a vessel radius `R`, cross-section edges
are approximately `target_size * R` in the geometry's units. Smaller values
create more quads and increase the volume mesh size. Actual lengths vary across
the disk; this is not a uniform or exact edge-length constraint. Axial spacing
is still controlled separately by `mesh.centerline` and `mesh.junction`.

The generator builds a bowed square core with rings extending to the circular
boundary. On the radius-1 disk, side midpoints extend to radius `0.70`, while
corners pull inward to coordinates `(±0.55, ±0.55)`. For normalized grid
coordinates `u,v` in `[-1,1]`, core points are
`x = u*(0.70 - 0.15*v*v)`, `y = v*(0.70 - 0.15*u*u)`.
This increases central spacing and reduces the outer-layer thickness. Resolution
uses `n = 2*ceil(pi/(4*target_size))` divisions per core side and
`ceil(0.30/target_size)` outer layers; intermediate layers interpolate
between the bowed core perimeter and the circle. Edge lengths and quad areas
still vary; the mapping does not impose uniform element sizes.

It folds two half disks and adds a third half disk sharing their
diameter to produce the merge, then derives all three branch connections and
boundary-node mappings from the generated topology. No Cubit installation or
predefined template file is needed in this mode; the required CLI template
path is ignored. The existing determinant, scaled-Jacobian, and intersection
gates also apply to these meshes.

Every pipeline run saves `cross_section_template.vtk` and `merge_template.vtk`
directly in `preprocessing/` when using `scripts/generate_case.sh`. This also
applies when using predefined templates. The manifest lists both paths.

The run exports only the two template VTK previews. It does not create a
`generated_templates/` directory or export template text files and reports.

## Geometry Safety

Preflight reports dimensionless `length/diameter`, diameter gradients,
`curvature*radius`, bifurcation angles/clearance, and broad-phase swept-tube
collision candidates. The curvature-radius safety limit, bifurcation limits,
and diameter-transition limit are hard preflight gates; swept-tube candidates
remain warnings until the exact generated-surface check runs. Element quality
is sampled at corners and a `4 x 4 x 4` Gauss grid. Point updates use
adjacency-aware backtracking; junction interiors undergo configured
quality-improving iterations. Generation fails when a determinant is
nonpositive, the scaled-Jacobian floor is violated, or non-adjacent exterior
triangles intersect. Use `iga_mesh_check` again after Bezier extraction.

The pipeline writes SWC coordinates to eight decimal places and reads that file
back before meshing. This deliberate quantization preserves the MATLAB file
contract and prevents `ceil()` layer counts from changing at roundoff
boundaries. The generated VTK coordinates use double precision so this
topology quantization does not impose a six-decimal absolute geometry scale.

See [the validation report](../../docs/MESH_CPP_VALIDATION.md) for regression
and large-case results.
