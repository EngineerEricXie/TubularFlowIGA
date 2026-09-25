# Skeleton-to-tetrahedron mesher

`skeleton_to_tet.py` converts a rooted SWC or radius-annotated line OBJ tree
directly into a smooth, labelled, first-order tetrahedral vessel mesh. It uses
cubic Hermite centerlines and overlapping capsules, blends a junction radius
into each child radius with a smoothstep taper, applies bounded surface
smoothing, and fills the resulting closed surface with Gmsh tetrahedra. It does
not use bifurcation templates and accepts two, three, or more children.

```bash
python3 preprocessing/tet/skeleton_to_tet.py \
  examples/vascular_flow/y_bifurcation/skeleton_initial.swc \
  artifacts/y_bifurcation/y \
  --length-scale-to-m 0.001 \
  --surface-target-size-m 0.0001 \
  --volume-target-size-m 0.0002 \
  --radius-transition-fraction 0.5
```

The output prefix receives:

- `.surface.vtp`: ParaView surface with `boundary_id` cell data;
- `.surface.msh`: the same labelled surface as Gmsh 4.1;
- `.volume.msh`: labelled first-order tetrahedra for native FEM;
- `.volume.vtu`: ParaView volume mesh;
- `.json`: geometry controls and mesh quality.

Boundary label `0` is the vessel wall, `1` is the root inlet, and labels `2`
and above are terminal outlets in deterministic tree order. Coordinates and
radii are both multiplied by `--length-scale-to-m`.

`--radius-transition-fraction 0.5` uses the first half of every child edge
leaving a junction to taper from the junction radius to the child radius. A
value of `1` uses the full edge. Non-junction edges interpolate their endpoint
radii over the full edge.

Before capsule construction, degree-2 points closer than
`--minimum-centerline-spacing-m` are reduced along each branch. Root, junction,
and terminal points are always retained. The default spacing is the smallest
input radius, which removes oversampled centerlines without changing sparse
hand-authored trees.

The current OCC construction targets connected, acyclic vessel trees whose
swept branches do not intersect away from declared junctions. Very dense
centerlines can require a larger `--minimum-centerline-spacing-m` to avoid
OpenCASCADE Boolean instability.
