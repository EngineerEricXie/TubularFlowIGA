# C++ Control-Mesh Generator

This directory replaces the MATLAB/TREES preprocessing stage with a standalone
C++17 program. It reads an SWC centerline, smooths and resamples each branch,
and writes the hexahedral control mesh consumed by the existing spline stage.
All spline evaluation, frame transport, bifurcation assembly, element-quality
checks, and VTK output are implemented here without a high-level geometry
library.

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

Create a case directory containing `skeleton_initial.swc` and
`mesh_parameter.txt`, then run from the repository root:

```bash
./preprocessing/mesh/tubular_mesh pipeline \
  /path/to/case meshgeneration/template 0.001
```

The final argument is the minimum accepted scaled Jacobian and is optional.
The command writes:

- `skeleton_smooth.swc`: smoothed and resampled centerline;
- `controlmesh.vtk`: labeled eight-node control elements;
- `initial_velocityfield.txt`: branch-aligned initial velocities.

Individual stages are also available:

```bash
tubular_mesh smooth INPUT.swc mesh_parameter.txt OUTPUT.swc
tubular_mesh generate SMOOTH.swc mesh_parameter.txt TEMPLATE_DIR \
  controlmesh.vtk initial_velocityfield.txt MIN_SCALED_J
```

## Parameters and Assumptions

`mesh_parameter.txt` contains five values: smoothing iteration count,
bifurcation smoothing ratio, noise smoothing ratio, target axial segment
length, and bifurcation refinement ratio. SWC column 6 is interpreted as a
radius and converted internally to diameter, matching the MATLAB workflow.

The topology must be a connected rooted tree. Each nonterminal node must have
one child or exactly two children; higher-order junctions are not yet
supported. The 201-point tube and 294-point merge templates under
`meshgeneration/template/` remain part of the stable interface.

## Geometry Safety

Element quality is sampled at corners and a `4 x 4 x 4` Gauss grid. Point
updates use adjacency-aware backtracking and are rejected if they invert a
neighboring element. Generation fails before output when a determinant is
nonpositive or the requested scaled-Jacobian floor is violated. For production
meshes, inspect the VTK geometry and use `iga_mesh_check` again after Bezier
extraction.

The pipeline writes SWC coordinates to eight decimal places and reads that file
back before meshing. This deliberate quantization preserves the MATLAB file
contract and prevents `ceil()` layer counts from changing at roundoff
boundaries.

See [the validation report](../../docs/MESH_CPP_VALIDATION.md) for regression
and large-case results.
