#!/usr/bin/env python3
"""Run fTetWild and normalize its tetrahedra to the native labelled Gmsh 4.1 contract."""

import argparse
import hashlib
import heapq
import json
import math
from pathlib import Path
import resource
import shutil
import subprocess
import tempfile
import time


def sha256(path):
	return hashlib.sha256(path.read_bytes()).hexdigest()


def vector(a, b):
	return tuple(b[i]-a[i] for i in range(3))


def dot(a, b):
	return sum(a[i]*b[i] for i in range(3))


def cross(a, b):
	return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2],
		a[0]*b[1]-a[1]*b[0])


def norm(a):
	return math.sqrt(dot(a, a))


def determinant(a, b, c, d):
	return dot(vector(a, b), cross(vector(a, c), vector(a, d)))


def scaled_jacobian(a, b, c, d):
	values = []
	for origin, first, second, third in ((a,b,c,d), (b,a,d,c),
			(c,d,a,b), (d,c,b,a)):
		u, v, w = vector(origin, first), vector(origin, second), vector(origin, third)
		denominator = norm(u)*norm(v)*norm(w)
		values.append(abs(dot(u, cross(v, w)))/denominator if denominator else 0.0)
	return min(values)


def boundary_shape_metrics(nodes, labelled, volume):
	edges, area = {}, 0.0
	for face, _ in labelled:
		points = tuple(nodes[tag] for tag in face)
		normal = cross(vector(points[0], points[1]), vector(points[0], points[2]))
		magnitude = norm(normal)
		if magnitude == 0.0:
			raise RuntimeError("degenerate output boundary triangle")
		area += magnitude/2.0
		unit = tuple(value/magnitude for value in normal)
		for first, second in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
			edges.setdefault(tuple(sorted((first, second))), []).append(unit)
	integrated, maximum_dihedral = 0.0, 0.0
	for edge, normals in edges.items():
		if len(normals) != 2:
			raise RuntimeError("output boundary curvature edge is not manifold")
		angle = math.acos(max(-1.0, min(1.0, dot(*normals))))
		integrated += 0.5*math.dist(nodes[edge[0]], nodes[edge[1]])*angle
		maximum_dihedral = max(maximum_dihedral, angle)
	pi = math.pi
	return {"area_m2": area,
		"volume_equivalent_radius_m": (3.0*volume/(4.0*pi))**(1.0/3.0),
		"area_equivalent_radius_m": math.sqrt(area/(4.0*pi)),
		"integrated_absolute_mean_curvature_m": integrated,
		"area_average_absolute_mean_curvature_per_m": integrated/area,
		"maximum_edge_dihedral_rad": maximum_dihedral}


def point_triangle_distance_squared(point, triangle):
	# Closest-point regions from Real-Time Collision Detection, Christer Ericson.
	a, b, c = triangle
	ab, ac, ap = vector(a, b), vector(a, c), vector(a, point)
	d1, d2 = dot(ab, ap), dot(ac, ap)
	if d1 <= 0.0 and d2 <= 0.0:
		return dot(ap, ap)
	bp = vector(b, point)
	d3, d4 = dot(ab, bp), dot(ac, bp)
	if d3 >= 0.0 and d4 <= d3:
		return dot(bp, bp)
	vc = d1*d4-d3*d2
	if vc <= 0.0 and d1 >= 0.0 and d3 <= 0.0:
		v = d1/(d1-d3)
		closest = tuple(a[i]+v*ab[i] for i in range(3))
		return dot(vector(closest, point), vector(closest, point))
	cp = vector(c, point)
	d5, d6 = dot(ab, cp), dot(ac, cp)
	if d6 >= 0.0 and d5 <= d6:
		return dot(cp, cp)
	vb = d5*d2-d1*d6
	if vb <= 0.0 and d2 >= 0.0 and d6 <= 0.0:
		w = d2/(d2-d6)
		closest = tuple(a[i]+w*ac[i] for i in range(3))
		return dot(vector(closest, point), vector(closest, point))
	va = d3*d6-d5*d4
	if va <= 0.0 and d4-d3 >= 0.0 and d5-d6 >= 0.0:
		w = (d4-d3)/((d4-d3)+(d5-d6))
		bc = vector(b, c)
		closest = tuple(b[i]+w*bc[i] for i in range(3))
		return dot(vector(closest, point), vector(closest, point))
	denominator = 1.0/(va+vb+vc)
	v, w = vb*denominator, vc*denominator
	closest = tuple(a[i]+ab[i]*v+ac[i]*w for i in range(3))
	delta = vector(closest, point)
	return dot(delta, delta)


class TriangleTree:
	def __init__(self, triangles):
		self.triangles = triangles
		self.nodes = []
		self.root = self._build(list(range(len(triangles))))

	def _build(self, indices):
		minimum = tuple(min(self.triangles[index][0][axis][coordinate]
			for index in indices for axis in range(3)) for coordinate in range(3))
		maximum = tuple(max(self.triangles[index][0][axis][coordinate]
			for index in indices for axis in range(3)) for coordinate in range(3))
		node = len(self.nodes)
		self.nodes.append(None)
		if len(indices) <= 8:
			self.nodes[node] = (minimum, maximum, indices, None, None)
			return node
		span = [maximum[i]-minimum[i] for i in range(3)]
		axis = max(range(3), key=span.__getitem__)
		indices.sort(key=lambda index: sum(point[axis]
			for point in self.triangles[index][0])/3.0)
		middle = len(indices)//2
		left = self._build(indices[:middle])
		right = self._build(indices[middle:])
		self.nodes[node] = (minimum, maximum, None, left, right)
		return node

	@staticmethod
	def _box_distance_squared(point, minimum, maximum):
		return sum((minimum[i]-point[i])**2 if point[i] < minimum[i]
			else (point[i]-maximum[i])**2 if point[i] > maximum[i] else 0.0
			for i in range(3))

	def nearest(self, point):
		best = math.inf
		best_index = None
		queue = [(0.0, self.root)]
		while queue:
			lower, node_index = heapq.heappop(queue)
			if lower > best:
				break
			minimum, maximum, indices, left, right = self.nodes[node_index]
			if indices is not None:
				for index in indices:
					distance = point_triangle_distance_squared(point,
						self.triangles[index][0])
					if distance < best:
						best, best_index = distance, index
			else:
				for child in (left, right):
					child_min, child_max = self.nodes[child][:2]
					distance = self._box_distance_squared(point, child_min, child_max)
					if distance <= best:
						heapq.heappush(queue, (distance, child))
		return best, self.triangles[best_index]


def read_surface_msh2(path):
	lines = path.read_text(encoding="utf-8").splitlines()
	physical = {}
	nodes = {}
	triangles = []
	index = 0
	while index < len(lines):
		section = lines[index].strip(); index += 1
		if section == "$PhysicalNames":
			count = int(lines[index]); index += 1
			for _ in range(count):
				parts = lines[index].split(maxsplit=2); index += 1
				physical[int(parts[1])] = parts[2].strip().strip('"')
			assert lines[index].strip() == "$EndPhysicalNames"; index += 1
		elif section == "$Nodes":
			count = int(lines[index]); index += 1
			for _ in range(count):
				parts = lines[index].split(); index += 1
				nodes[int(parts[0])] = tuple(map(float, parts[1:4]))
			assert lines[index].strip() == "$EndNodes"; index += 1
		elif section == "$Elements":
			count = int(lines[index]); index += 1
			for _ in range(count):
				parts = list(map(int, lines[index].split())); index += 1
				type_, tags = parts[1], parts[2]
				if type_ == 2:
					name = physical.get(parts[3]) if tags else None
					if not name or not name.startswith("boundary_label_"):
						raise RuntimeError("canonical surface triangle lacks boundary label")
					triangles.append((tuple(parts[3+tags:6+tags]), int(name[15:])))
			assert lines[index].strip() == "$EndElements"; index += 1
	if not nodes or not triangles:
		raise RuntimeError("canonical surface MSH2 is empty")
	return nodes, triangles


def write_obj(path, nodes, triangles):
	ordered = sorted(nodes)
	remap = {tag: index+1 for index, tag in enumerate(ordered)}
	with path.open("w", encoding="utf-8") as output:
		output.write("# canonical preflight surface; labels live in the sidecar manifest\n")
		for tag in ordered:
			output.write("v {:.17g} {:.17g} {:.17g}\n".format(*nodes[tag]))
		for face, _ in triangles:
			output.write("f {} {} {}\n".format(*(remap[tag] for tag in face)))


def read_ftetwild_msh2(path):
	lines = path.read_text(encoding="utf-8").splitlines()
	if len(lines) < 3 or lines[0].strip() != "$MeshFormat" \
			or lines[1].split()[:2] != ["2.2", "0"]:
		raise RuntimeError("fTetWild output must be ASCII Gmsh 2.2 (--no-binary)")
	nodes, tetrahedra = {}, []
	index = 0
	while index < len(lines):
		section = lines[index].strip(); index += 1
		if section == "$Nodes":
			count = int(lines[index]); index += 1
			for _ in range(count):
				parts = lines[index].split(); index += 1
				if len(parts) != 4:
					raise RuntimeError("invalid fTetWild node")
				nodes[int(parts[0])] = tuple(map(float, parts[1:]))
			if lines[index].strip() != "$EndNodes": raise RuntimeError("truncated fTetWild nodes")
			index += 1
		elif section == "$Elements":
			count = int(lines[index]); index += 1
			for _ in range(count):
				parts = list(map(int, lines[index].split())); index += 1
				tags = parts[2]
				if parts[1] != 4:
					raise RuntimeError("fTetWild output contains a non-tetrahedral element")
				cell = tuple(parts[3+tags:7+tags])
				if len(cell) != 4: raise RuntimeError("invalid fTetWild tetrahedron")
				tetrahedra.append(cell)
			if lines[index].strip() != "$EndElements": raise RuntimeError("truncated fTetWild elements")
			index += 1
	if not nodes or not tetrahedra:
		raise RuntimeError("fTetWild output contains no tetrahedra")
	return nodes, tetrahedra


def oriented_boundary(nodes, tetrahedra):
	orientations = ((1,2,3), (0,3,2), (0,1,3), (0,2,1))
	positive = []
	unique = set()
	for cell in tetrahedra:
		if len(set(cell)) != 4 or any(tag not in nodes for tag in cell):
			raise RuntimeError("fTetWild tetrahedron has repeated or missing nodes")
		cell = list(cell)
		value = determinant(*(nodes[tag] for tag in cell))
		if not math.isfinite(value) or value == 0.0:
			raise RuntimeError("fTetWild produced a degenerate tetrahedron")
		if value < 0.0:
			cell[2], cell[3] = cell[3], cell[2]
		cell = tuple(cell)
		key = tuple(sorted(cell))
		if key in unique: raise RuntimeError("fTetWild output contains a duplicate tetrahedron")
		unique.add(key); positive.append(cell)
	faces = {}
	adjacency = [set() for _ in positive]
	for cell_index, cell in enumerate(positive):
		for local in orientations:
			face = tuple(cell[i] for i in local)
			key = tuple(sorted(face))
			if key in faces:
				if faces[key][0] is None:
					raise RuntimeError("fTetWild output has a nonmanifold tetrahedral facet")
				owner = faces[key][1]
				adjacency[owner].add(cell_index); adjacency[cell_index].add(owner)
				faces[key] = (None, owner)
			else:
				faces[key] = (face, cell_index)
	visited, pending = set(), [0]
	while pending:
		cell = pending.pop()
		if cell in visited: continue
		visited.add(cell); pending.extend(adjacency[cell]-visited)
	if len(visited) != len(positive):
		raise RuntimeError("fTetWild output contains disconnected volume regions")
	boundary = [record[0] for record in faces.values() if record[0] is not None]
	if not boundary:
		raise RuntimeError("fTetWild tetrahedral boundary is empty")
	return positive, boundary


def label_boundary(nodes, boundary, source_nodes, source_faces, envelope, ambiguity):
	by_label = {}
	for face, label in source_faces:
		points = tuple(source_nodes[tag] for tag in face)
		normal = cross(vector(points[0], points[1]), vector(points[0], points[2]))
		length = norm(normal)
		if not length: raise RuntimeError("degenerate canonical source triangle")
		by_label.setdefault(label, []).append((points, tuple(value/length for value in normal)))
	trees = {label: TriangleTree(triangles) for label, triangles in by_label.items()}
	labelled, maximum_distance, minimum_normal_dot = [], 0.0, 1.0
	counts = {label: 0 for label in trees}
	for face in boundary:
		points = tuple(nodes[tag] for tag in face)
		centroid = tuple(sum(point[i] for point in points)/3.0 for i in range(3))
		candidates = []
		for label, tree in trees.items():
			distance2, triangle = tree.nearest(centroid)
			candidates.append((math.sqrt(distance2), label, triangle))
		candidates.sort(key=lambda item: (item[0], item[1]))
		if len(candidates) > 1 and candidates[1][0]-candidates[0][0] <= ambiguity:
			raise RuntimeError("ambiguous fTetWild boundary label within mapping tolerance")
		distance, label, source = candidates[0]
		if distance > envelope:
			raise RuntimeError("fTetWild boundary exceeds the configured envelope: "
				f"face={face}, centroid_m={centroid}, source_distance_m={distance:.9g}, "
				f"envelope_m={envelope:.9g}")
		for point in points:
			point_distance = min(math.sqrt(tree.nearest(point)[0]) for tree in trees.values())
			maximum_distance = max(maximum_distance, point_distance)
			if point_distance > envelope:
				raise RuntimeError("fTetWild boundary vertex exceeds the configured envelope: "
				f"vertex={point}, source_distance_m={point_distance:.9g}, "
				f"envelope_m={envelope:.9g}")
		normal = cross(vector(points[0], points[1]), vector(points[0], points[2]))
		length = norm(normal)
		if not length: raise RuntimeError("degenerate fTetWild boundary triangle")
		normal_dot = dot(tuple(value/length for value in normal), source[1])
		minimum_normal_dot = min(minimum_normal_dot, normal_dot)
		if normal_dot <= 0.0:
			raise RuntimeError("fTetWild boundary orientation disagrees with source surface: "
				f"face={face}, centroid_m={centroid}, source_label={label}, "
				f"source_distance_m={distance:.9g}, normal_dot={normal_dot:.9g}")
		labelled.append((face, label)); counts[label] += 1
	if any(count == 0 for count in counts.values()):
		raise RuntimeError("fTetWild output lost an entire source boundary label")
	return labelled, counts, maximum_distance, minimum_normal_dot


def write_msh41(path, nodes, tetrahedra, labelled):
	labels = sorted({label for _, label in labelled})
	physical = {label: index+1 for index, label in enumerate(labels)}
	volume_entity = len(labels)+1
	minimum = tuple(min(point[i] for point in nodes.values()) for i in range(3))
	maximum = tuple(max(point[i] for point in nodes.values()) for i in range(3))
	by_label = {label: [] for label in labels}
	for face, label in labelled: by_label[label].append(face)
	ordered_nodes = sorted(nodes)
	with path.open("w", encoding="utf-8") as output:
		output.write("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$PhysicalNames\n")
		output.write(f"{len(labels)+1}\n")
		for label in labels:
			output.write(f'2 {physical[label]} "boundary_label_{label}"\n')
		output.write(f'3 {volume_entity} "fluid"\n$EndPhysicalNames\n$Entities\n')
		output.write(f"0 0 {len(labels)} 1\n")
		for label in labels:
			entity = physical[label]
			output.write(f"{entity} {' '.join(map(str, minimum+maximum))} 1 {physical[label]} 0\n")
		output.write(f"{volume_entity} {' '.join(map(str, minimum+maximum))} 1 "
			f"{volume_entity} {len(labels)} {' '.join(str(physical[label]) for label in labels)}\n")
		output.write("$EndEntities\n$Nodes\n")
		output.write(f"1 {len(nodes)} {ordered_nodes[0]} {ordered_nodes[-1]}\n")
		output.write(f"3 {volume_entity} 0 {len(nodes)}\n")
		for tag in ordered_nodes: output.write(f"{tag}\n")
		for tag in ordered_nodes: output.write("{:.17g} {:.17g} {:.17g}\n".format(*nodes[tag]))
		total = len(labelled)+len(tetrahedra)
		output.write(f"$EndNodes\n$Elements\n{len(labels)+1} {total} 1 {total}\n")
		element = 1
		for label in labels:
			faces = by_label[label]
			output.write(f"2 {physical[label]} 2 {len(faces)}\n")
			for face in faces:
				output.write(f"{element} {' '.join(map(str, face))}\n"); element += 1
		output.write(f"3 {volume_entity} 4 {len(tetrahedra)}\n")
		for cell in tetrahedra:
			output.write(f"{element} {' '.join(map(str, cell))}\n"); element += 1
		output.write("$EndElements\n")
	return physical, volume_entity


def main():
	start_wall = time.monotonic()
	root = Path(__file__).resolve().parents[1]
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("input_surface")
	parser.add_argument("output_msh")
	parser.add_argument("--manifest")
	parser.add_argument("--ftetwild", default="FloatTetwild_bin")
	parser.add_argument("--surface-preflight",
		default=str(root/"solvers"/"cpu"/"surface_fem_preflight"))
	parser.add_argument("--boundary-array", default="boundary_id")
	parser.add_argument("--length-scale-to-m", type=float, default=1.0)
	parser.add_argument("--weld-tolerance-m", type=float, default=0.0)
	parser.add_argument("--default-boundary-id", type=int, default=0)
	parser.add_argument("--target-size-m", type=float, required=True)
	parser.add_argument("--envelope-m", type=float, required=True)
	parser.add_argument("--stop-energy", type=float, default=10.0)
	parser.add_argument("--max-optimization-passes", type=int, default=80)
	parser.add_argument("--max-threads", type=int, default=1)
	parser.add_argument("--minimum-scaled-jacobian", type=float, default=1.0e-3)
	parser.add_argument("--maximum-volume-relative-error", type=float, default=2.0e-2)
	parser.add_argument("--label-ambiguity-tolerance-m", type=float)
	parser.add_argument("--ftetwild-version", default="unknown")
	args = parser.parse_args()
	for name, value in (("target size", args.target_size_m), ("envelope", args.envelope_m)):
		if not math.isfinite(value) or value <= 0.0: raise RuntimeError(f"{name} must be positive")
	if args.max_threads <= 0: raise RuntimeError("max threads must be positive")
	if args.max_optimization_passes <= 0:
		raise RuntimeError("maximum optimization passes must be positive")
	if not 0.0 < args.minimum_scaled_jacobian <= 1.0:
		raise RuntimeError("minimum scaled Jacobian must be in (0,1]")
	if not 0.0 <= args.maximum_volume_relative_error < 1.0:
		raise RuntimeError("maximum volume relative error must be in [0,1)")
	ambiguity = args.label_ambiguity_tolerance_m
	if ambiguity is None: ambiguity = max(1.0e-12, args.envelope_m*1.0e-6)
	if ambiguity < 0.0 or not math.isfinite(ambiguity):
		raise RuntimeError("label ambiguity tolerance must be finite and nonnegative")
	input_path, output_path = Path(args.input_surface).resolve(), Path(args.output_msh).resolve()
	manifest_path = Path(args.manifest).resolve() if args.manifest else output_path.with_suffix(".json")
	preflight = Path(args.surface_preflight).resolve()
	ftetwild_lookup = shutil.which(args.ftetwild)
	ftetwild = Path(ftetwild_lookup).resolve() if ftetwild_lookup else Path(args.ftetwild).resolve()
	if not preflight.is_file(): raise RuntimeError(f"missing surface preflight: {preflight}")
	if not ftetwild.is_file(): raise RuntimeError(f"missing fTetWild executable: {ftetwild}")
	with tempfile.TemporaryDirectory(prefix="tubularflow-ftetwild-") as directory:
		temporary = Path(directory)
		canonical, preflight_json = temporary/"surface.msh", temporary/"surface.json"
		command = [str(preflight), str(input_path), str(canonical), "--manifest",
			str(preflight_json), "--boundary-array", args.boundary_array,
			"--length-scale-to-m", str(args.length_scale_to_m), "--weld-tolerance-m",
			str(args.weld_tolerance_m), "--default-boundary-id", str(args.default_boundary_id)]
		completed = subprocess.run(command, text=True, capture_output=True)
		if completed.returncode: raise RuntimeError((completed.stderr or completed.stdout).strip())
		surface_data = json.loads(preflight_json.read_text(encoding="utf-8"))
		source_nodes, source_faces = read_surface_msh2(canonical)
		obj, raw = temporary/"surface.obj", temporary/"volume.msh"
		write_obj(obj, source_nodes, source_faces)
		minimum = tuple(min(point[i] for point in source_nodes.values()) for i in range(3))
		maximum = tuple(max(point[i] for point in source_nodes.values()) for i in range(3))
		diagonal = math.dist(minimum, maximum)
		if diagonal <= 0.0: raise RuntimeError("canonical source has a zero bounding-box diagonal")
		command = [str(ftetwild), "--input", str(obj), "--output", str(raw),
			"--la", str(args.target_size_m), "--epsr",
			str(args.envelope_m/diagonal), "--stop-energy", str(args.stop_energy), "--max-threads",
			str(args.max_threads), "--max-its", str(args.max_optimization_passes),
			"--no-binary", "--no-color"]
		before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
		start = time.monotonic()
		completed = subprocess.run(command, text=True, capture_output=True)
		wall = time.monotonic()-start
		peak = max(before, resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)*1024
		if completed.returncode: raise RuntimeError((completed.stderr or completed.stdout).strip())
		if not raw.is_file(): raise RuntimeError("fTetWild did not create its requested output")
		nodes, cells = read_ftetwild_msh2(raw)
		cells, boundary = oriented_boundary(nodes, cells)
		labelled, counts, distance, normal_dot = label_boundary(nodes, boundary,
			source_nodes, source_faces, args.envelope_m, ambiguity)
	qualities = [scaled_jacobian(*(nodes[tag] for tag in cell)) for cell in cells]
	minimum_quality = min(qualities)
	if minimum_quality < args.minimum_scaled_jacobian:
		raise RuntimeError(f"minimum scaled Jacobian {minimum_quality} is below the gate")
	volume = sum(determinant(*(nodes[tag] for tag in cell))/6.0 for cell in cells)
	source_volume = float(surface_data["volume_m3"])
	volume_error = abs(volume-source_volume)/source_volume
	if volume_error > args.maximum_volume_relative_error:
		raise RuntimeError(f"relative volume error {volume_error} exceeds the gate")
	shape = boundary_shape_metrics(nodes, labelled, volume)
	boundary_area = shape["area_m2"]
	source_area = float(surface_data["area_m2"])
	area_error = abs(boundary_area-source_area)/source_area
	source_radii = surface_data["equivalent_radii_m"]
	source_curvature = surface_data["discrete_curvature"]
	geometry_change = {
		"classification": "envelope_constrained_boundary_remesh_not_repair",
		"implicit_repair": False,
		"relative_volume_change": volume_error,
		"relative_surface_area_change": area_error,
		"volume_equivalent_radius_relative_change": abs(
			shape["volume_equivalent_radius_m"]-source_radii["enclosed_volume_sphere"])
			/source_radii["enclosed_volume_sphere"],
		"area_equivalent_radius_relative_change": abs(
			shape["area_equivalent_radius_m"]-source_radii["surface_area_sphere"])
			/source_radii["surface_area_sphere"],
		"integrated_absolute_mean_curvature_relative_change": abs(
			shape["integrated_absolute_mean_curvature_m"]
			-source_curvature["integrated_absolute_mean_curvature_m"])
			/source_curvature["integrated_absolute_mean_curvature_m"],
		"maximum_source_surface_distance_m": distance,
		"output_shape_metrics": shape}
	output_path.parent.mkdir(parents=True, exist_ok=True)
	physical_tags, volume_tag = write_msh41(output_path, nodes, cells, labelled)
	output_sha256 = sha256(output_path)
	peak_rss_bytes = max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
		resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)*1024
	manifest = {
		"schema_version": 1, "route": "surface_to_fem_volume_ftetwild",
		"surface_preflight": surface_data,
		"geometry_contract": {
			"schema_version": 1, "route": "surface_to_fem_volume_ftetwild",
			"coordinate_system": "cartesian", "length_unit": "m",
			"source_geometry": {"identity_sha256": surface_data["canonical_sha256"],
				"source_sha256": surface_data["source"]["sha256"],
				"source_length_scale_to_m": surface_data["applied_length_scale_to_m"]},
			"reference_geometry": {"identity_sha256": output_sha256,
				"derived_from_source_sha256": surface_data["canonical_sha256"]},
			"current_geometry": {"kind": "same_as_reference",
				"identity_sha256": output_sha256},
			"stable_ids": {"volume_nodes": "output_msh_node_tag",
				"volume_elements": "output_msh_element_tag",
				"boundary_triangles": "output_msh_element_tag"},
			"regions": [{"id": "fluid", "role": "fluid", "dimension": 3,
				"physical_group": "fluid", "physical_group_tag": volume_tag}],
			"boundaries": [{"label": label, "physical_group": f"boundary_label_{label}",
				"physical_group_tag": physical_tags[label]}
				for label in sorted(counts)]},
		"mesher": {"name": "fTetWild", "declared_version": args.ftetwild_version,
			"executable_sha256": sha256(ftetwild), "target_size_m": args.target_size_m,
			"envelope_m": args.envelope_m, "stop_energy": args.stop_energy,
			"max_optimization_passes": args.max_optimization_passes,
			"max_threads": args.max_threads, "command": command,
			"wall_s": wall, "child_peak_rss_bytes": peak},
		"resources": {"pipeline_wall_s": time.monotonic()-start_wall,
			"peak_rss_bytes": peak_rss_bytes},
		"geometry_change": geometry_change,
		"volume_mesh": {"format": "gmsh41_ascii", "element_type": "tetrahedron_p1",
			"sha256": output_sha256, "nodes": len(nodes), "elements": len(cells),
			"boundary_triangles": len(labelled),
			"boundary_labels": {str(label): count for label, count in counts.items()},
			"boundary_matches_tetrahedra": True, "volume_m3": volume,
			"volume_relative_error": volume_error,
			"minimum_scaled_jacobian": minimum_quality,
			"boundary_surface_area_m2": boundary_area,
			"boundary_surface_area_relative_error": area_error,
			"maximum_source_surface_distance_m": distance,
			"minimum_source_normal_dot": normal_dot,
			"boundary_remeshed": True,
			"implicit_surface_repair": False},
		"gates": {"required_minimum_scaled_jacobian": args.minimum_scaled_jacobian,
			"maximum_volume_relative_error": args.maximum_volume_relative_error,
			"maximum_surface_distance_m": args.envelope_m,
			"label_ambiguity_tolerance_m": ambiguity, "passed": True}}
	manifest_path.parent.mkdir(parents=True, exist_ok=True)
	manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True)+"\n", encoding="utf-8")
	print(f"ftetwild_to_fem_volume: PASS elements={len(cells)} minSJ={minimum_quality:.8g} "
		f"volume_error={volume_error:.3g} max_distance={distance:.3g}")
	return 0


if __name__ == "__main__":
	try:
		raise SystemExit(main())
	except (AssertionError, KeyError, OSError, RuntimeError, ValueError) as error:
		print(f"ftetwild_to_fem_volume: ERROR: {error}")
		raise SystemExit(2)
