# Radius-1 cross-section creation in Cubit

The reference is a planar, all-quadrilateral disk. Start with Cubit's Circle
scheme; this creates a comparable layout, not a guaranteed reproduction of the
reference's connectivity or the existing template's numbering.

1. Start a fresh Cubit session and play back `create_unit_circle.jou`.
   The journal's `reset` clears the current model. It creates an origin-centered
   circle of radius 1 in the XY plane and meshes it with four-node quads.
2. Adjust `curve all interval 40` to control circumferential resolution. Try 24
   for fewer elements or 48 for more; multiples of eight are a starting point
   for a symmetric center. The Circle scheme's `fraction 0.5` controls the
   interior region's relative size. Inspect the actual mesh counts and quality
   in Cubit rather than assuming a count from the interval setting.
3. Play back `export_unit_circle.py` as a Python script in the same session.
   It exports surface 1 into `unit_circle_template/` under Cubit's working
   directory. It writes coordinate/connectivity text, a VTK preview, and a JSON
   report. For a different surface or output path, call `export_surface` with
   those arguments from Cubit's Python console.

A command-line alternative, from a directory containing both scripts, is:

```bash
coreform_cubit -nogui -batch create_unit_circle.jou export_unit_circle.py
```

Radius 1 is the normalized reference radius. TubularFlowIGA scales this
cross-section to the physical vessel radius while sweeping it along a branch.
The exporter maps Cubit node IDs to contiguous zero-based indices and makes
quadrilaterals counterclockwise in the XY plane. It checks planarity, unit
boundary radius, quad convexity, edge incidence, and Euler characteristic; it
reports interior nodes that the current `radius > 0.95` wall rule would label
incorrectly. These checks are not a proof of bifurcation compatibility or a
replacement for mesh-quality validation in Cubit and the full generator.

## Connecting the candidate to the vessel generator

The candidate is not automatically installed as the active template. The
current C++ loader expects 201 circle points/180 quads and 294 merge points/270
quads. Junction construction uses specific element partitions and reorder
indices. Even a new circle with the same counts can have incompatible ordering.
A new topology requires a matching merge template and updated junction mapping;
only changing the dimension checks is insufficient. Do not overwrite the
existing four active templates until these interfaces have been updated and
validated. Regenerate and validate the control mesh and Bezier representation
after integration.

The Cubit journal could not be executed in the current environment because
Cubit is not available on PATH. The exporter was checked against the existing
201-point/180-quad template using the documented API interface, including
noncontiguous node IDs, reversed quad orientation, and rejection of radius 2.

## Primary references

- [Circular surface creation](https://coreform.com/cubit_help/geometry/geom_creation/bottom_up_creation/surface.htm)
- [Circle meshing scheme](https://coreform.com/cubit_help/mesh_generation/meshing_schemes/traditional/circle.htm)
- [Cubit Python API](https://coreform.com/cubit_help/python/cubit_python_api_method_based.htm)
- [Journal and Python execution](https://coreform.com/cubit_help/environment_control/session_control/execution_command_syntax.htm)
