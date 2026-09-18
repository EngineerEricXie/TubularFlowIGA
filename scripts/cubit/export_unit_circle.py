"""Play back in Cubit after create_unit_circle.jou; exports a candidate only.

Writes to unit_circle_template below Cubit's current working directory.
Does not overwrite the repository's active templates or create a merge template.
"""

from collections import Counter
import json
import math
from pathlib import Path


def export_surface(cubit, surface_id, output_directory):
    selection = "in surface {}".format(surface_id)
    if cubit.parse_cubit_list("tri", selection):
        raise ValueError("The cross-section must contain only quadrilaterals")
    quad_ids = sorted(cubit.parse_cubit_list("quad", selection))
    if not quad_ids:
        raise ValueError("No quadrilaterals found on the selected surface")
    connections = [tuple(cubit.get_expanded_connectivity("quad", q)) for q in quad_ids]
    if any(len(q) != 4 or len(set(q)) != 4 for q in connections):
        raise ValueError("Use four-node quads without higher-order nodes")
    node_ids = sorted({node for q in connections for node in q})
    index = {node: i for i, node in enumerate(node_ids)}
    points = [tuple(cubit.get_nodal_coordinates(n)) for n in node_ids]
    if any(not all(math.isfinite(v) for v in p) or abs(p[2]) > 1e-8 for p in points):
        raise ValueError("All points must be finite and lie in the XY plane, z=0")
    if any(math.hypot(p[0], p[1]) > 1.0 + 1e-6 for p in points):
        raise ValueError("The mesh must fit in an origin-centered unit disk")
    faces = []
    for q in connections:
        face = [index[n] for n in q]
        xy = [points[n] for n in face]
        area2 = sum(xy[i][0]*xy[(i+1)%4][1]-xy[(i+1)%4][0]*xy[i][1] for i in range(4))
        if area2 < 0:
            face = [face[0], face[3], face[2], face[1]]
            xy = [points[n] for n in face]
        for i in range(4):
            a, b, c = xy[i], xy[(i+1)%4], xy[(i+2)%4]
            turn = (b[0]-a[0])*(c[1]-b[1])-(b[1]-a[1])*(c[0]-b[0])
            if turn <= 0:
                raise ValueError("Degenerate, concave, or folded quadrilateral")
        faces.append(face)
    edges = Counter(tuple(sorted((q[i], q[(i+1)%4]))) for q in faces for i in range(4))
    if any(count > 2 for count in edges.values()):
        raise ValueError("Non-manifold mesh edge")
    boundary = sorted({n for edge, count in edges.items() if count == 1 for n in edge})
    if not boundary or any(abs(math.hypot(*points[n][:2])-1.0) > 1e-6 for n in boundary):
        raise ValueError("All exterior boundary nodes must have radius 1 about the origin")
    if len(points)-len(edges)+len(faces) != 1:
        raise ValueError("The mesh must have disk topology")
    boundary_set = set(boundary)
    false_wall = [n for n in range(len(points)) if n not in boundary_set
                  and math.hypot(*points[n][:2]) > 0.95]
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    names = ["template_circle_points90.txt", "template_circle_elements.txt", "unit_circle.vtk", "export_report.json"]
    if any((output / name).exists() for name in names):
        raise ValueError("Candidate output already exists; choose a fresh output directory")
    (output / names[0]).write_text("".join("{:.17g} {:.17g} 0\n".format(*p[:2]) for p in points))
    (output / names[1]).write_text("".join("4 {} {} {} {}\n".format(*q) for q in faces))
    vtk = ["# vtk DataFile Version 2.0", "Cubit unit disk candidate", "ASCII", "DATASET UNSTRUCTURED_GRID",
           "POINTS {} double".format(len(points))]
    vtk.extend("{:.17g} {:.17g} 0".format(*p[:2]) for p in points)
    vtk.append("CELLS {} {}".format(len(faces), 5*len(faces)))
    vtk.extend("4 {} {} {} {}".format(*q) for q in faces)
    vtk.append("CELL_TYPES {}".format(len(faces)))
    vtk.extend("9" for q in faces)
    (output / names[2]).write_text("\n".join(vtk)+"\n")
    report = {"points": len(points), "quads": len(faces), "boundary_nodes": len(boundary),
              "cubit_node_ids_in_export_order": node_ids,
              "interior_nodes_misclassified_by_current_wall_rule": false_wall,
              "matches_current_count_checks": len(points) == 201 and len(faces) == 180,
              "bifurcation_compatibility_verified": False}
    (output / names[3]).write_text(json.dumps(report, indent=2)+"\n")
    print("Candidate exported to {}: {} points, {} quads".format(output.resolve(), len(points), len(faces)))
    print("Matching merge template and junction ordering are required before replacement.")
    if false_wall:
        print("Current radius>0.95 wall rule would mislabel {} interior nodes.".format(len(false_wall)))
    return report


if __name__ == "__main__":
    import cubit
    export_surface(cubit, 1, Path.cwd() / "unit_circle_template")
