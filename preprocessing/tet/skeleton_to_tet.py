#!/usr/bin/env python3
"""Generate a smooth labelled tetrahedral vessel mesh from an SWC or line OBJ skeleton."""

import argparse
from collections import deque
import hashlib
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from preprocessing.tet.skeleton_geometry import (clean_surface_triangles, cross,
	dot, fit_terminal_caps, fit_tree_radii, smooth_surface_groups,
	spline_capsule_paths, sub)


def canonical_tree(points, radii, edges, root, scale):
	adjacency = [[] for _ in points]
	for first, second in edges:
		adjacency[first].append(second)
		adjacency[second].append(first)
	order = []
	parent = {root: None}
	queue = deque([root])
	while queue:
		node = queue.popleft()
		order.append(node)
		for neighbor in sorted(adjacency[node]):
			if neighbor == parent[node]:
				continue
			if neighbor in parent:
				raise ValueError("skeleton contains a cycle")
			parent[neighbor] = node
			queue.append(neighbor)
	if len(order) != len(points):
		raise ValueError("skeleton is disconnected")
	remap = {old: new for new, old in enumerate(order)}
	nodes = [[scale*value for value in points[old]] for old in order]
	node_radii = [scale*radii[old] for old in order]
	segments = []
	children = [[] for _ in order]
	for old in order[1:]:
		start, stop = remap[parent[old]], remap[old]
		segments.append([start, stop, node_radii[stop]])
		children[start].append(stop)
	return {"nodes_m": nodes, "node_radii_m": node_radii,
		"segments": segments, "root": 0,
		"terminals": [index for index, group in enumerate(children) if not group]}


def read_swc(path, scale):
	records = []
	for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
		line = line.strip()
		if not line or line.startswith("#"):
			continue
		fields = line.split()
		if len(fields) != 7:
			raise ValueError(f"{path}:{line_number}: expected seven SWC columns")
		records.append((int(fields[0]), [float(value) for value in fields[2:5]],
			float(fields[5]), int(fields[6])))
	ids = {record[0]: index for index, record in enumerate(records)}
	roots = [index for index, record in enumerate(records) if record[3] == -1]
	if len(roots) != 1:
		raise ValueError("SWC requires exactly one root")
	edges = [(ids[parent], index) for index, (_, _, _, parent) in enumerate(records)
		if parent != -1]
	return canonical_tree([record[1] for record in records],
		[record[2] for record in records], edges, roots[0], scale)


def read_obj(path, scale):
	points, radii, edges = [], [], []
	for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
		line = line.strip()
		if not line or line.startswith("#"):
			continue
		fields = line.split()
		if fields[0] == "v" and len(fields) == 7:
			points.append([float(value) for value in fields[1:4]])
			radii.append(float(fields[4]))
		elif fields[0] == "l" and len(fields) >= 3:
			indices = [int(value)-1 for value in fields[1:]]
			edges.extend(zip(indices[:-1], indices[1:]))
		else:
			raise ValueError(f"{path}:{line_number}: expected radius vertex or line")
	adjacency = [0]*len(points)
	for first, second in edges:
		adjacency[first] += 1
		adjacency[second] += 1
	terminals = [index for index, degree in enumerate(adjacency) if degree == 1]
	root = max(terminals, key=lambda index: (radii[index], -index))
	return canonical_tree(points, radii, list(edges), root, scale)


def read_skeleton(path, scale):
	return read_obj(path, scale) if path.suffix.lower() == ".obj" else read_swc(path, scale)


def simplify_tree(tree, minimum_spacing):
	parents = [None]*len(tree["nodes_m"])
	children = [[] for _ in tree["nodes_m"]]
	for start, stop, _ in tree["segments"]:
		parents[stop] = start
		children[start].append(stop)
	critical = {index for index, group in enumerate(children)
		if index == tree["root"] or len(group) != 1}
	keep = set(critical)
	for start in sorted(critical):
		for first in children[start]:
			section = [start, first]
			while len(children[section[-1]]) == 1:
				section.append(children[section[-1]][0])
			arc = 0.
			for previous, node in zip(section[:-1], section[1:]):
				arc += math.dist(tree["nodes_m"][previous], tree["nodes_m"][node])
				if node == section[-1] or arc >= minimum_spacing:
					keep.add(node)
					arc = 0.
	order = sorted(keep)
	remap = {old: new for new, old in enumerate(order)}
	segments = []
	new_children = [[] for _ in order]
	for old in order[1:]:
		ancestor = parents[old]
		while ancestor not in keep:
			ancestor = parents[ancestor]
		start, stop = remap[ancestor], remap[old]
		segments.append([start, stop, tree["node_radii_m"][old]])
		new_children[start].append(stop)
	return {"nodes_m": [tree["nodes_m"][old] for old in order],
		"node_radii_m": [tree["node_radii_m"][old] for old in order],
		"segments": segments, "root": 0,
		"terminals": [index for index, group in enumerate(new_children) if not group]}


def add_capsule_tree(gmsh, tree, geometry):
	parts = []
	paths = spline_capsule_paths(tree, geometry["subsegments_per_edge"],
		geometry["tangent_scale"])
	radii = tree["node_radii_m"]
	for start, stop, _, points in paths:
		base_radius = min(radii[start], radii[stop])
		for index in range(len(points)-1):
			delta = sub(points[index+1], points[index])
			parts.append((3, gmsh.model.occ.addCylinder(*points[index], *delta,
				base_radius)))
		for point in points[1:-1]:
			parts.append((3, gmsh.model.occ.addSphere(*point, base_radius)))
	incident = [[] for _ in tree["nodes_m"]]
	children = [[] for _ in tree["nodes_m"]]
	for start, stop, _ in tree["segments"]:
		incident[start].append(min(radii[start], radii[stop]))
		incident[stop].append(min(radii[start], radii[stop]))
		children[start].append(stop)
	for index in range(len(tree["nodes_m"])):
		if index not in tree["terminals"] and (index != tree["root"] or len(children[index]) > 1):
			radius = radii[index] if len(children[index]) > 1 else min(incident[index])
			parts.append((3, gmsh.model.occ.addSphere(*tree["nodes_m"][index], radius)))
	fragmented, _ = gmsh.model.occ.fragment([parts[0]], parts[1:])
	return [entity for entity in fragmented if entity[0] == 3]


def tet_quality(points):
	def corner(a, b, c, d):
		u, v, w = sub(b, a), sub(c, a), sub(d, a)
		denominator = math.sqrt(dot(u, u)*dot(v, v)*dot(w, w))
		return dot(u, cross(v, w))/denominator
	a, b, c, d = points
	return corner(a, b, c, d), min(abs(corner(*order)) for order in
		((a, b, c, d), (b, a, d, c), (c, d, a, b), (d, c, b, a)))


def write_surface_vtp(path, points, triangles, labels):
	import vtk
	from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
	import numpy as np
	tags = sorted(points)
	local = {tag: index for index, tag in enumerate(tags)}
	faces = [(label, triangle) for label in labels for triangle in triangles[label]]
	coordinates = vtk.vtkPoints()
	coordinates.SetData(numpy_to_vtk(np.asarray([points[tag] for tag in tags],
		dtype=np.float64), deep=True))
	packed = np.empty((len(faces), 4), dtype=np.int64)
	packed[:, 0] = 3
	packed[:, 1:] = [[local[tag] for tag in triangle] for _, triangle in faces]
	polygons = vtk.vtkCellArray()
	polygons.SetCells(len(faces), numpy_to_vtkIdTypeArray(packed.ravel(), deep=True))
	polydata = vtk.vtkPolyData()
	polydata.SetPoints(coordinates)
	polydata.SetPolys(polygons)
	values = numpy_to_vtk(np.asarray([labels[label] for label, _ in faces],
		dtype=np.int32), deep=True)
	values.SetName("boundary_id")
	polydata.GetCellData().AddArray(values)
	polydata.GetCellData().SetActiveScalars("boundary_id")
	writer = vtk.vtkXMLPolyDataWriter()
	writer.SetFileName(str(path))
	writer.SetInputData(polydata)
	writer.SetDataModeToBinary()
	writer.Write()


def write_volume_vtu(path, coordinates, connectivity):
	import vtk
	from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
	import numpy as np
	tags = sorted(coordinates)
	local = {tag: index for index, tag in enumerate(tags)}
	points = vtk.vtkPoints()
	points.SetData(numpy_to_vtk(np.asarray([coordinates[tag] for tag in tags],
		dtype=np.float64), deep=True))
	packed = np.empty((len(connectivity), 5), dtype=np.int64)
	packed[:, 0] = 4
	packed[:, 1:] = [[local[tag] for tag in cell] for cell in connectivity]
	cells = vtk.vtkCellArray()
	cells.SetCells(len(connectivity), numpy_to_vtkIdTypeArray(packed.ravel(), deep=True))
	grid = vtk.vtkUnstructuredGrid()
	grid.SetPoints(points)
	grid.SetCells(vtk.VTK_TETRA, cells)
	writer = vtk.vtkXMLUnstructuredGridWriter()
	writer.SetFileName(str(path))
	writer.SetInputData(grid)
	writer.SetDataModeToBinary()
	writer.Write()


def generate(tree, output_prefix, controls):
	try:
		import gmsh
	except ImportError as error:
		raise ValueError(f"Gmsh Python module unavailable: {error}") from error
	gmsh.initialize()
	try:
		gmsh.option.setNumber("General.Terminal", 0)
		gmsh.model.add("skeleton-to-tet")
		volumes = add_capsule_tree(gmsh, tree, controls)
		gmsh.model.occ.synchronize()
		oriented = gmsh.model.getBoundary(volumes, combined=True, oriented=True,
			recursive=False)
		signs = {abs(tag): (-1 if tag < 0 else 1) for dimension, tag in oriented
			if dimension == 2}
		endpoints = [tree["root"]]+tree["terminals"]
		labels = {"wall": 0, "inlet": 1,
			**{f"outlet_{index}": index+2 for index in range(len(tree["terminals"]))}}
		surface_groups = {name: [] for name in labels}
		for surface in sorted(signs):
			center = gmsh.model.occ.getCenterOfMass(2, surface)
			distances = [math.dist(center, tree["nodes_m"][node]) for node in endpoints]
			nearest = min(range(len(endpoints)), key=distances.__getitem__)
			radius = tree["node_radii_m"][endpoints[nearest]]
			name = ("inlet" if nearest == 0 else f"outlet_{nearest-1}") \
				if distances[nearest] < 1.e-6*radius else "wall"
			surface_groups[name].append(surface)
		gmsh.option.setNumber("Mesh.MeshSizeMin", controls["surface_target_size_m"])
		gmsh.option.setNumber("Mesh.MeshSizeMax", controls["surface_target_size_m"])
		gmsh.option.setNumber("Mesh.ElementOrder", 1)
		gmsh.model.mesh.generate(2)
		node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
		nodes = {int(tag): tuple(coordinates[3*index:3*index+3])
			for index, tag in enumerate(node_tags)}
		triangles = {name: [] for name in labels}
		for name, surfaces in surface_groups.items():
			for surface in surfaces:
				types, _, blocks = gmsh.model.mesh.getElements(2, surface)
				for element_type, block in zip(types, blocks):
					if element_type != 2:
						raise ValueError("surface mesh is not first-order triangular")
					for start in range(0, len(block), 3):
						face = tuple(int(tag) for tag in block[start:start+3])
						if signs[surface] < 0:
							face = (face[0], face[2], face[1])
						triangles[name].append(face)
		points, triangles, collapsed, dropped = clean_surface_triangles(nodes,
			triangles, controls["join_cleanup_tolerance_m"],
			{"fluid": list(labels)})
		cap_moved, cap_displacement = fit_terminal_caps(points, triangles, tree)
		moved, maximum_displacement = fit_tree_radii(points, triangles, tree,
			controls, ("wall",))
		moved += cap_moved
		maximum_displacement = max(maximum_displacement, cap_displacement)
		if controls["surface_smoothing_iterations"]:
			points = smooth_surface_groups(points, triangles, ("wall",),
				controls["surface_smoothing_iterations"],
				controls["surface_smoothing_pass_band"])
		write_surface_vtp(Path(str(output_prefix)+".surface.vtp"), points,
			triangles, labels)
		gmsh.clear()
		gmsh.model.add("skeleton-to-tet-discrete")
		entities = {}
		for index, name in enumerate(labels, 1):
			gmsh.model.addDiscreteEntity(2, index)
			entities[name] = index
		tags = sorted(points)
		gmsh.model.mesh.addNodes(2, 1, tags,
			[value for tag in tags for value in points[tag]])
		element = 1
		for name, entity in entities.items():
			count = len(triangles[name])
			gmsh.model.mesh.addElementsByType(entity, 2,
				list(range(element, element+count)),
				[tag for triangle in triangles[name] for tag in triangle])
			element += count
			group = gmsh.model.addPhysicalGroup(2, [entity], labels[name]+1)
			gmsh.model.setPhysicalName(2, group, f"boundary_label_{labels[name]}")
		loop = gmsh.model.geo.addSurfaceLoop(list(entities.values()))
		volume = gmsh.model.geo.addVolume([loop])
		gmsh.model.geo.synchronize()
		group = gmsh.model.addPhysicalGroup(3, [volume])
		gmsh.model.setPhysicalName(3, group, "fluid")
		gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
		gmsh.option.setNumber("Mesh.Binary", 0)
		gmsh.write(str(output_prefix)+".surface.msh")
		gmsh.option.setNumber("Mesh.MeshSizeMin", controls["volume_target_size_m"])
		gmsh.option.setNumber("Mesh.MeshSizeMax", controls["volume_target_size_m"])
		gmsh.option.setNumber("Mesh.Algorithm3D", 1)
		gmsh.model.mesh.generate(3)
		types, element_tags, blocks = gmsh.model.mesh.getElements(3)
		if list(types) != [4]:
			raise ValueError("volume mesh is not first-order tetrahedral")
		node_tags, flat, _ = gmsh.model.mesh.getNodes()
		coordinates = {int(tag): tuple(flat[3*index:3*index+3])
			for index, tag in enumerate(node_tags)}
		connectivity = [tuple(int(tag) for tag in blocks[0][start:start+4])
			for start in range(0, len(blocks[0]), 4)]
		minimum_scaled = 1.
		for cell in connectivity:
			signed, scaled = tet_quality([coordinates[tag] for tag in cell])
			if signed <= 0.:
				raise ValueError("tetrahedral mesh has a nonpositive Jacobian")
			minimum_scaled = min(minimum_scaled, scaled)
		if minimum_scaled < controls["minimum_scaled_jacobian"]:
			raise ValueError(f"minimum scaled Jacobian {minimum_scaled:.6g} below gate")
		gmsh.write(str(output_prefix)+".volume.msh")
		write_volume_vtu(Path(str(output_prefix)+".volume.vtu"), coordinates, connectivity)
	finally:
		gmsh.finalize()
	manifest = {"schema_version": 1, "kind": "skeleton_to_tetrahedral_vessel_mesh",
		"source_sha256": hashlib.sha256(controls["input_path"].read_bytes()).hexdigest(),
		"controls": {key: value for key, value in controls.items() if key != "input_path"},
		"skeleton_nodes": len(tree["nodes_m"]), "skeleton_segments": len(tree["segments"]),
		"outlets": len(tree["terminals"]), "collapsed_surface_nodes": collapsed,
		"dropped_surface_triangles": dropped, "radius_fitted_surface_nodes": moved,
		"maximum_radius_fit_displacement_m": maximum_displacement,
		"volume_nodes": len(coordinates), "tetrahedra": len(connectivity),
		"minimum_scaled_jacobian": minimum_scaled,
		"boundary_labels": labels}
	Path(str(output_prefix)+".json").write_text(json.dumps(manifest, indent=2,
		sort_keys=True)+"\n", encoding="utf-8")
	return manifest


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("input_skeleton", type=Path)
	parser.add_argument("output_prefix", type=Path)
	parser.add_argument("--length-scale-to-m", type=float, default=1.)
	parser.add_argument("--surface-target-size-m", type=float, required=True)
	parser.add_argument("--volume-target-size-m", type=float, required=True)
	parser.add_argument("--subsegments-per-edge", type=int, default=2)
	parser.add_argument("--tangent-scale", type=float, default=.25)
	parser.add_argument("--minimum-centerline-spacing-m", type=float)
	parser.add_argument("--radius-transition-fraction", type=float, default=.5)
	parser.add_argument("--join-cleanup-tolerance-m", type=float)
	parser.add_argument("--surface-smoothing-iterations", type=int, default=80)
	parser.add_argument("--surface-smoothing-pass-band", type=float, default=.01)
	parser.add_argument("--minimum-scaled-jacobian", type=float, default=.001)
	args = parser.parse_args()
	path = args.input_skeleton.resolve()
	output_prefix = args.output_prefix.resolve()
	output_prefix.parent.mkdir(parents=True, exist_ok=True)
	tree = read_skeleton(path, args.length_scale_to_m)
	source_nodes = len(tree["nodes_m"])
	minimum_spacing = args.minimum_centerline_spacing_m or min(tree["node_radii_m"])
	tree = simplify_tree(tree, minimum_spacing)
	controls = {"input_path": path,
		"length_scale_to_m": args.length_scale_to_m,
		"source_skeleton_nodes": source_nodes,
		"minimum_centerline_spacing_m": minimum_spacing,
		"surface_target_size_m": args.surface_target_size_m,
		"volume_target_size_m": args.volume_target_size_m,
		"subsegments_per_edge": args.subsegments_per_edge,
		"tangent_scale": args.tangent_scale,
		"radius_transition_fraction": args.radius_transition_fraction,
		"join_cleanup_tolerance_m": args.join_cleanup_tolerance_m or
			.03*min(tree["node_radii_m"]),
		"surface_smoothing_iterations": args.surface_smoothing_iterations,
		"surface_smoothing_pass_band": args.surface_smoothing_pass_band,
		"minimum_scaled_jacobian": args.minimum_scaled_jacobian}
	result = generate(tree, output_prefix, controls)
	print("skeleton_to_tet: PASS nodes={} tetrahedra={} outlets={} minSJ={:.8g}".format(
		result["volume_nodes"], result["tetrahedra"], result["outlets"],
		result["minimum_scaled_jacobian"]))


if __name__ == "__main__":
	try:
		main()
	except (OSError, RuntimeError, ValueError) as error:
		print(f"skeleton_to_tet: ERROR: {error}", file=sys.stderr)
		raise SystemExit(2)
