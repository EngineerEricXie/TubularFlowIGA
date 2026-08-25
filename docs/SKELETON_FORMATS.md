# Skeleton input formats

TubularFlowIGA accepts rooted SWC centerlines and radius-annotated line OBJ
trees. The same input readers feed native 1D networks and the 3D control-mesh
preprocessor.

## SWC

Each non-comment row has:

```text
id type x y z radius parent_id
```

Exactly one node has parent `-1`. IDs must be unique, every other parent must
exist, coordinates and positive radii must be finite, and the parent graph must
be one connected acyclic tree. Rows must contain exactly seven columns.

## Radius-annotated line OBJ

This is a centerline convention, not a surface OBJ. It supports comments plus
only `v` and `l` records:

```text
v x y z radius auxiliary0 auxiliary1
l vertex_index_1 vertex_index_2
```

The fourth value after `v` is the vessel radius. The two auxiliary values are
read and required to be finite but are not used. The supplied
`liver_veins_central.obj` uses zero for both. Line indices are positive,
1-based, vertex-only OBJ indices. A polyline such as `l 1 2 3` is accepted and
creates edges 1–2 and 2–3. Texture/normal index syntax, faces, duplicate edges,
self edges, zero-length edges, cycles, undefined indices, and disconnected
vertices are rejected.

Because OBJ edges are undirected, the default root is the terminal vertex with
the largest radius; ties choose the lowest vertex index. A schema-v3 1D case
can explicitly select a different vertex:

```json
"geometry": {
  "kind": "obj_network",
  "file": "skeleton_initial.obj",
  "length_scale_to_m": 0.001,
  "root_node_id": 1
}
```

`root_node_id` uses the OBJ 1-based vertex index and is valid only for
`obj_network`. If omitted, root inference is deterministic. The original liver
file selects vertex 1 by inference.

## 1D and 3D use

Every successful 1D simulation and 3D preparation writes:

- `skeleton_normalized.swc`, an explicitly rooted canonical SWC;
- `skeleton.vtp`, a ParaView-ready line skeleton.

The VTP contains point arrays `radius`, `diameter`, `node_id`, `parent_id`,
`degree`, and `role`, plus cell arrays `segment_id` and `branch_id`. Numeric
roles are `0` interior, `1` root, `2` junction, and `3` outlet. In ParaView,
apply **Tube**, choose **By Absolute Scalar**, and select `radius`. The original
input is never overwritten. `iga_1d --check` remains validation-only and does
not write these files.

For 1D, `length_scale_to_m` converts both coordinates and radius to SI. The
network solver supports any number of child branches at a junction. Its
`skeleton.vtp` uses those SI coordinates so it overlays the time-series VTP;
the normalized SWC retains source units.

For the 3D preparation pipeline, name the input `skeleton_initial.swc` or
`skeleton_initial.obj`; exactly one must exist. Mesh `seg_length` uses the same
coordinate unit as the skeleton. The resulting spline representation is later
normalized as documented in the [pipeline guide](PIPELINE.md).

The current 3D control-mesh generator supports binary branching: after root
orientation, a node may have at most two children. A 3D preparation stops with
the offending node ID and child count when this is violated. The complete
liver OBJ contains multiway junctions and is therefore valid for 1D but must be
reduced to a binary subtree or otherwise preprocessed before full 3D meshing.
The source-only `liver_vein_obj_segment` examples exercise both supported paths
on vertices 1–12 from that file.
