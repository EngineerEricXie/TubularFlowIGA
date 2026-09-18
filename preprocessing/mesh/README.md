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

`generate` also writes `mesh_quality.json` beside `controlmesh.vtk`.
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
supported. The 201-point tube and 294-point merge templates under
`meshgeneration/template/` remain part of the stable interface.

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
