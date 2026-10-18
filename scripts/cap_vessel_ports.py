#!/usr/bin/env python3
"""Cut and cap labelled ports on a closed vessel surface, then remesh it.

Input is a closed triangulated lumen surface, for example the voxel isosurface
written by dicom_seg_to_surface.py, and an explicit list of port planes (a
point, an outward tangent, and a radius per port). The tool

1. optionally smooths the surface with a windowed-sinc filter (explicit, and
   the per-vertex displacement is reported),
2. removes the blind end beyond each port plane, splitting
   the crossing triangles exactly on the plane after snapping nearby vertices,
3. closes every planar cut loop with a cap carrying the port's boundary label,
4. optionally remeshes the closed surface with Gmsh, bounding the element size
   by a target size and a multiple of the local lumen radius (the distance from
   the wall to an SWC centerline), and
5. writes a labelled VTP (cell array ``boundary_id``) plus a JSON manifest.

The output is meant for surface_to_fem_volume.py, which repeats the closed,
manifold, oriented and self-intersection checks in surface_fem_preflight
before meshing the volume. Nothing here infers ports: every cut is listed in
the ports file. All lengths are in the input surface unit unless named _m.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import resource
import sys
import time


def reject(message):
	raise ValueError(message)


def file_sha256(path):
	digest = hashlib.sha256()
	with Path(path).open("rb") as stream:
		for block in iter(lambda: stream.read(1024*1024), b""):
			digest.update(block)
	return digest.hexdigest()


def progress(message):
	print(f"cap_vessel_ports: {message}", file=sys.stderr, flush=True)


def read_surface(path):
	import numpy as np
	import vtk
	from vtk.util.numpy_support import vtk_to_numpy
	suffix = Path(path).suffix.lower()
	if suffix == ".stl":
		reader = vtk.vtkSTLReader()
	elif suffix == ".vtp":
		reader = vtk.vtkXMLPolyDataReader()
	else:
		reject("input surface must be .stl or .vtp")
	reader.SetFileName(str(path))
	reader.Update()
	clean = vtk.vtkCleanPolyData()
	clean.SetInputConnection(reader.GetOutputPort())
	clean.PointMergingOn()
	clean.SetTolerance(0.0)
	clean.Update()
	triangles = vtk.vtkTriangleFilter()
	triangles.SetInputConnection(clean.GetOutputPort())
	triangles.PassVertsOff()
	triangles.PassLinesOff()
	triangles.Update()
	data = triangles.GetOutput()
	if data.GetNumberOfPolys() == 0 or data.GetNumberOfPolys() != data.GetNumberOfCells():
		reject("input surface must contain only triangles")
	points = vtk_to_numpy(data.GetPoints().GetData()).astype(float)
	faces = vtk_to_numpy(data.GetPolys().GetConnectivityArray()).reshape(-1, 3).astype(np.int64)
	return points, faces


def write_vtp(path, points, faces, labels):
	"""Minimal uncompressed ASCII VTP accepted by surface_fem_preflight."""
	lines = ['<?xml version="1.0"?>',
		'<VTKFile type="PolyData" version="1.0" byte_order="LittleEndian" header_type="UInt64">',
		'<PolyData>', f'<Piece NumberOfPoints="{len(points)}" NumberOfPolys="{len(faces)}">',
		'<Points>', '<DataArray type="Float64" NumberOfComponents="3" format="ascii">']
	lines.extend(" ".join(repr(float(value)) for value in point) for point in points)
	lines += ['</DataArray>', '</Points>', '<Polys>',
		'<DataArray type="Int64" Name="connectivity" format="ascii">']
	lines.extend(" ".join(str(int(vertex)) for vertex in face) for face in faces)
	lines += ['</DataArray>', '<DataArray type="Int64" Name="offsets" format="ascii">',
		" ".join(str(3*(index+1)) for index in range(len(faces))), '</DataArray>', '</Polys>',
		'<CellData>', '<DataArray type="Int32" Name="boundary_id" format="ascii">',
		" ".join(str(int(label)) for label in labels), '</DataArray>', '</CellData>',
		'</Piece>', '</PolyData>', '</VTKFile>', '']
	Path(path).write_text("\n".join(lines), encoding="ascii")


def surface_volume_and_area(points, faces):
	import numpy as np
	a, b, c = points[faces[:, 0]], points[faces[:, 1]], points[faces[:, 2]]
	area = 0.5*np.linalg.norm(np.cross(b-a, c-a), axis=1).sum()
	volume = np.einsum("ij,ij->i", a, np.cross(b, c)).sum()/6.0
	return float(volume), float(area)


def triangle_quality(points, faces):
	"""Normalized 4*sqrt(3)*area/sum(edge^2); 1 for equilateral triangles."""
	import numpy as np
	a, b, c = points[faces[:, 0]], points[faces[:, 1]], points[faces[:, 2]]
	area = 0.5*np.linalg.norm(np.cross(b-a, c-a), axis=1)
	squares = (np.sum((b-a)**2, axis=1)+np.sum((c-b)**2, axis=1)+np.sum((a-c)**2, axis=1))
	return 4.0*math.sqrt(3.0)*area/np.maximum(squares, 1e-300)


def smooth(points, faces, iterations, passband):
	import numpy as np
	import vtk
	from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray, vtk_to_numpy
	data = vtk.vtkPolyData()
	vtk_points = vtk.vtkPoints()
	vtk_points.SetData(numpy_to_vtk(np.ascontiguousarray(points), deep=True))
	data.SetPoints(vtk_points)
	cells = vtk.vtkCellArray()
	cells.SetData(numpy_to_vtkIdTypeArray(np.arange(0, 3*len(faces)+1, 3, dtype=np.int64), deep=True),
		numpy_to_vtkIdTypeArray(np.ascontiguousarray(faces.ravel()), deep=True))
	data.SetPolys(cells)
	smoother = vtk.vtkWindowedSincPolyDataFilter()
	smoother.SetInputData(data)
	smoother.SetNumberOfIterations(iterations)
	smoother.SetPassBand(passband)
	smoother.BoundarySmoothingOff()
	smoother.FeatureEdgeSmoothingOff()
	smoother.NonManifoldSmoothingOn()
	smoother.NormalizeCoordinatesOn()
	smoother.Update()
	result = vtk_to_numpy(smoother.GetOutput().GetPoints().GetData()).astype(float)
	if result.shape != points.shape:
		reject("smoothing changed the vertex count")
	return result


def read_ports(path):
	data = json.loads(Path(path).read_text(encoding="utf-8"))
	if not isinstance(data, dict) or data.get("schema_version") != 1:
		reject("ports file must be an object with schema_version 1")
	ports = data.get("ports")
	if not isinstance(ports, list) or not ports:
		reject("ports file needs a nonempty ports list")
	result = []
	for index, port in enumerate(ports):
		if not isinstance(port, dict) or set(port) - {"label", "point", "outward_tangent", "radius",
				"name", "cut_offset_radii"}:
			reject(f"port {index} has unknown or missing fields")
		label = port.get("label")
		point = port.get("point")
		tangent = port.get("outward_tangent")
		radius = port.get("radius")
		if not isinstance(label, int) or isinstance(label, bool) or label < 1:
			reject(f"port {index} label must be a positive integer (0 is the wall)")
		if (not isinstance(point, list) or len(point) != 3 or not isinstance(tangent, list)
				or len(tangent) != 3 or not all(isinstance(v, (int, float)) and math.isfinite(v)
					for v in point+tangent)):
			reject(f"port {index} point and outward_tangent must be finite 3-vectors")
		length = math.sqrt(sum(v*v for v in tangent))
		if not length > 0.0:
			reject(f"port {index} outward_tangent is zero")
		if not isinstance(radius, (int, float)) or not math.isfinite(radius) or radius <= 0.0:
			reject(f"port {index} radius must be finite and positive")
		offset = port.get("cut_offset_radii")
		if offset is not None and (not isinstance(offset, (int, float)) or not math.isfinite(offset)
				or offset < 0.0):
			reject(f"port {index} cut_offset_radii must be finite and nonnegative")
		result.append({"index": index, "name": port.get("name", f"port_{index}"), "label": label,
			"point": [float(v) for v in point], "tangent": [float(v)/length for v in tangent],
			"radius": float(radius), "offset": None if offset is None else float(offset)})
	return result


class SurfaceCutter:
	"""Local plane cuts on a closed, consistently oriented triangle surface."""

	def __init__(self, points, faces):
		self.points = [list(map(float, point)) for point in points]
		self.faces = [list(map(int, face)) for face in faces]
		self.labels = [0]*len(self.faces)
		self.alive = [True]*len(self.faces)

	def edge_faces(self):
		edges = {}
		for index, face in enumerate(self.faces):
			if not self.alive[index]:
				continue
			for corner in range(3):
				key = tuple(sorted((face[corner], face[(corner+1) % 3])))
				edges.setdefault(key, []).append(index)
		return edges

	def cut(self, port, offset_radii, extent_radii, loop_radii, snap_fraction, minimum_cap_fraction):
		import numpy as np
		tangent = np.array(port["tangent"])
		center = np.array(port["point"])-offset_radii*port["radius"]*tangent
		extent = extent_radii*port["radius"]
		points = np.array(self.points)
		alive = [index for index in range(len(self.faces)) if self.alive[index]]
		faces = np.array([self.faces[index] for index in alive])
		signed = (points-center)@tangent
		corner_signed = signed[faces]
		touches_outside = corner_signed.max(axis=1) > 0.0
		centroids = points[faces].mean(axis=1)
		distance = np.linalg.norm(centroids-center, axis=1)
		candidates = np.nonzero(touches_outside & (distance < loop_radii*port["radius"]))[0]
		if len(candidates) == 0:
			return None, "no surface beyond the plane near the port"
		# Seed at the outside triangle closest to the tip direction, then grow
		# through triangles that touch the outside half-space. The growth stops
		# at the plane, so it can only reach the blind end and whatever else is
		# connected to it beyond the plane; a second crossing shows up as a
		# second boundary loop and is rejected below.
		tip = center+port["radius"]*tangent
		seed = candidates[np.argmin(np.linalg.norm(centroids[candidates]-tip, axis=1))]
		neighbor_edges = {}
		for local, face in enumerate(faces):
			for corner in range(3):
				key = (min(face[corner], face[(corner+1) % 3]), max(face[corner], face[(corner+1) % 3]))
				neighbor_edges.setdefault(key, []).append(local)
		region = {int(seed)}
		stack = [int(seed)]
		while stack:
			current = stack.pop()
			face = faces[current]
			for corner in range(3):
				key = (min(face[corner], face[(corner+1) % 3]), max(face[corner], face[(corner+1) % 3]))
				# Only cross edges that have an outside endpoint, so the region
				# stops at the plane instead of leaking along inside triangles.
				if max(signed[key[0]], signed[key[1]]) <= 0.0:
					continue
				for neighbor in neighbor_edges[key]:
					if neighbor in region or not touches_outside[neighbor]:
						continue
					if distance[neighbor] > extent:
						return None, "the beyond-plane region extends too far (plane cuts another branch)"
					region.add(neighbor)
					stack.append(neighbor)
		region_faces = [alive[local] for local in sorted(region)]
		region_vertices = sorted({vertex for index in region_faces for vertex in self.faces[index]})
		# Snap vertices close to the plane onto it, scaled by their shortest edge.
		shortest = {}
		for index in region_faces:
			face = self.faces[index]
			for corner in range(3):
				a, b = face[corner], face[(corner+1) % 3]
				length = float(np.linalg.norm(points[a]-points[b]))
				shortest[a] = min(shortest.get(a, math.inf), length)
				shortest[b] = min(shortest.get(b, math.inf), length)
		snapped = []
		for vertex in region_vertices:
			value = signed[vertex]
			if value != 0.0 and abs(value) < snap_fraction*shortest[vertex]:
				snapped.append((vertex, float(abs(value))))
		trial_points = points.copy()
		trial_signed = signed.copy()
		for vertex, _ in snapped:
			trial_points[vertex] = trial_points[vertex]-trial_signed[vertex]*tangent
			trial_signed[vertex] = 0.0
		# A vertex shared with triangles outside the region must not move,
		# otherwise untouched geometry would change.
		region_set = set(region_faces)
		vertex_faces = {}
		for index in alive:
			for vertex in self.faces[index]:
				vertex_faces.setdefault(vertex, []).append(index)
		for vertex, _ in snapped:
			if any(index not in region_set and max(signed[v] for v in self.faces[index]) > 0.0
					for index in vertex_faces[vertex]):
				return None, "snapping would move a vertex outside the cut region"
		new_points = []
		intersection = {}

		def crossing(a, b):
			key = (min(a, b), max(a, b))
			if key not in intersection:
				da, db = trial_signed[a], trial_signed[b]
				weight = da/(da-db)
				point = trial_points[a]+weight*(trial_points[b]-trial_points[a])
				point = point-((point-center)@tangent)*tangent
				intersection[key] = len(self.points)+len(new_points)
				new_points.append(point.tolist())
			return intersection[key]

		new_faces = []
		removed = []
		for index in region_faces:
			face = self.faces[index]
			sign = [1 if trial_signed[v] > 0.0 else (-1 if trial_signed[v] < 0.0 else 0) for v in face]
			removed.append(index)
			if min(sign) >= 0:
				continue
			# Rotate so the face starts at an inside vertex, keeping orientation.
			kept = []
			for corner in range(3):
				a, b = face[corner], face[(corner+1) % 3]
				sa, sb = sign[corner], sign[(corner+1) % 3]
				if sa <= 0:
					kept.append(a)
				if sa*sb < 0:
					kept.append(crossing(a, b))
			if len(kept) == 3:
				new_faces.append(kept)
			elif len(kept) == 4:
				new_faces.append([kept[0], kept[1], kept[2]])
				new_faces.append([kept[0], kept[2], kept[3]])
			else:
				return None, "unexpected triangle split configuration"
		on_plane = set(intersection.values())
		for vertex in range(len(trial_signed)):
			if trial_signed[vertex] == 0.0:
				on_plane.add(vertex)
		# Boundary edges of the kept surface lying on the plane form the loop.
		counts = {}
		dead = set(removed)
		all_faces = [self.faces[index] for index in alive if index not in dead]+new_faces
		for face in all_faces:
			for corner in range(3):
				a, b = face[corner], face[(corner+1) % 3]
				counts[(a, b)] = counts.get((a, b), 0)+1
		boundary = [(a, b) for (a, b) in counts if (b, a) not in counts]
		touched = {v for index in region_faces for v in self.faces[index]} | on_plane
		loop_edges = [(a, b) for a, b in boundary if a in touched and b in touched]
		if not loop_edges:
			return None, "the cut produced no boundary loop"
		successor = {}
		for a, b in loop_edges:
			if a in successor:
				return None, "the cut loop is not simple"
			successor[a] = b
		start = loop_edges[0][0]
		loop = [start]
		while True:
			following = successor.get(loop[-1])
			if following is None:
				return None, "the cut loop is open"
			if following == start:
				break
			loop.append(following)
			if len(loop) > len(successor):
				return None, "the cut loop does not close"
		if len(loop) != len(successor):
			return None, f"the cut produced {len(successor)} boundary edges but one loop has {len(loop)}"
		if len(loop) < 3:
			return None, "the cut loop has fewer than three vertices"
		all_points = np.vstack([trial_points, np.array(new_points).reshape(-1, 3)])
		loop_distance = np.linalg.norm(all_points[loop]-center, axis=1).max()
		if loop_distance > loop_radii*port["radius"]:
			return None, f"the cut loop reaches {loop_distance/port['radius']:.3g} port radii from the plane center"
		polygon = all_points[loop]-center
		normal_area = 0.5*np.cross(polygon, np.roll(polygon, -1, axis=0)).sum(axis=0)
		equivalent_radius = math.sqrt(abs(float(normal_area@tangent))/math.pi)
		if equivalent_radius < minimum_cap_fraction*port["radius"]:
			return None, (f"the cap radius {equivalent_radius:.3g} is below {minimum_cap_fraction} "
				"port radii")
		return {"equivalent_radius": equivalent_radius, "center": center, "tangent": tangent, "points": trial_points, "new_points": new_points,
			"removed": removed, "new_faces": new_faces, "loop": loop,
			"snapped": snapped, "region_faces": len(region_faces)}, None

	def apply(self, result, label):
		import numpy as np
		points = result["points"]
		for vertex in range(len(self.points)):
			self.points[vertex] = points[vertex].tolist()
		self.points.extend(result["new_points"])
		for index in result["removed"]:
			self.alive[index] = False
		for face in result["new_faces"]:
			self.faces.append(face)
			self.labels.append(0)
			self.alive.append(True)
		# The kept surface traverses the loop with the outside on its left; the
		# cap traverses it backwards so its normal is the outward tangent.
		cap_loop = list(reversed(result["loop"]))
		coordinates = np.array([self.points[v] for v in cap_loop])
		tangent = result["tangent"]
		helper = np.array([1.0, 0.0, 0.0]) if abs(tangent[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
		u = np.cross(tangent, helper)
		u /= np.linalg.norm(u)
		v = np.cross(tangent, u)
		planar = np.stack([(coordinates-result["center"])@u, (coordinates-result["center"])@v], axis=1)
		triangles = ear_clip(planar)
		if triangles is None:
			reject("cap polygon is not simple and cannot be triangulated")
		for a, b, c in triangles:
			self.faces.append([cap_loop[a], cap_loop[b], cap_loop[c]])
			self.labels.append(label)
			self.alive.append(True)
		area = 0.0
		for a, b, c in triangles:
			area += 0.5*abs((planar[b][0]-planar[a][0])*(planar[c][1]-planar[a][1])
				-(planar[c][0]-planar[a][0])*(planar[b][1]-planar[a][1]))
		return area

	def arrays(self):
		import numpy as np
		faces = [face for index, face in enumerate(self.faces) if self.alive[index]]
		labels = [label for index, label in enumerate(self.labels) if self.alive[index]]
		used = sorted({vertex for face in faces for vertex in face})
		remap = {vertex: index for index, vertex in enumerate(used)}
		points = np.array([self.points[vertex] for vertex in used])
		faces = np.array([[remap[v] for v in face] for face in faces], dtype=np.int64)
		return points, faces, np.array(labels, dtype=np.int64)


def ear_clip(polygon):
	"""Triangulate a simple 2D polygon (counter-clockwise after orientation)."""
	count = len(polygon)
	area = 0.0
	for index in range(count):
		x0, y0 = polygon[index]
		x1, y1 = polygon[(index+1) % count]
		area += x0*y1-x1*y0
	order = list(range(count))
	flipped = area < 0.0
	if flipped:
		order.reverse()

	def cross(o, a, b):
		return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])

	for i in range(count):
		a1, a2 = polygon[i], polygon[(i+1) % count]
		for j in range(i+2, count):
			if i == 0 and j == count-1:
				continue
			b1, b2 = polygon[j], polygon[(j+1) % count]
			if (cross(a1, a2, b1)*cross(a1, a2, b2) < 0.0
					and cross(b1, b2, a1)*cross(b1, b2, a2) < 0.0):
				return None
	result = []
	guard = 0
	while len(order) > 3 and guard < 10*count*count:
		guard += 1
		clipped = False
		for position in range(len(order)):
			p, c, n = order[position-1], order[position], order[(position+1) % len(order)]
			if cross(polygon[p], polygon[c], polygon[n]) <= 0.0:
				continue
			inside = False
			for other in order:
				if other in (p, c, n):
					continue
				q = polygon[other]
				if (cross(polygon[p], polygon[c], q) >= 0.0 and cross(polygon[c], polygon[n], q) >= 0.0
						and cross(polygon[n], polygon[p], q) >= 0.0):
					inside = True
					break
			if inside:
				continue
			result.append((p, c, n))
			del order[position]
			clipped = True
			break
		if not clipped:
			return None
	if len(order) == 3:
		result.append(tuple(order))
	if flipped:
		result = [(a, c, b) for a, b, c in result]
	return result


def read_centerline_samples(path, length_scale, spacing):
	"""SWC segments resampled every `spacing`, so the nearest sample distance
	approximates the distance to the centerline polyline within spacing/2."""
	import numpy as np
	nodes, parents = {}, {}
	for line in Path(path).read_text(encoding="utf-8").splitlines():
		line = line.strip()
		if not line or line.startswith("#"):
			continue
		values = line.split()
		if len(values) < 7:
			reject("SWC rows need seven columns")
		identifier = int(values[0])
		nodes[identifier] = np.array([float(values[2]), float(values[3]), float(values[4])])*length_scale
		parents[identifier] = int(values[6])
	if not nodes or any(not np.all(np.isfinite(point)) for point in nodes.values()):
		reject("SWC centerline is empty or nonfinite")
	samples = [point for point in nodes.values()]
	for identifier, parent in parents.items():
		if parent < 0:
			continue
		if parent not in nodes:
			reject(f"SWC node {identifier} has a missing parent {parent}")
		a, b = nodes[parent], nodes[identifier]
		count = int(math.ceil(np.linalg.norm(b-a)/spacing))
		for step in range(1, count):
			samples.append(a+(b-a)*step/count)
	return np.array(samples)


def boundary_loops(faces):
	"""Directed boundary loops of an open, consistently oriented surface."""
	counts = {}
	for face in faces:
		for corner in range(3):
			edge = (int(face[corner]), int(face[(corner+1) % 3]))
			counts[edge] = counts.get(edge, 0)+1
	successor = {}
	for (a, b), count in counts.items():
		if count != 1:
			reject("remeshed wall has a repeated directed edge")
		if (b, a) in counts:
			continue
		if a in successor:
			reject("remeshed wall boundary is not a set of simple loops")
		successor[a] = b
	loops, seen = [], set()
	for start in successor:
		if start in seen:
			continue
		loop = [start]
		seen.add(start)
		while successor[loop[-1]] != start:
			following = successor.get(loop[-1])
			if following is None or following in seen:
				reject("remeshed wall boundary loop is open")
			loop.append(following)
			seen.add(following)
		loops.append(loop)
	return loops


def thin_wall_samples(points, faces, search, factor):
	"""Points where opposite-facing walls are closer than `search`, with the
	size bound factor*distance there. They catch pinches and thin sheets that a
	centerline does not describe; elsewhere the returned arrays are empty."""
	import numpy as np
	from scipy.spatial import cKDTree
	centroids = points[faces].mean(axis=1)
	normals = np.cross(points[faces[:, 1]]-points[faces[:, 0]], points[faces[:, 2]]-points[faces[:, 0]])
	normals /= np.maximum(np.linalg.norm(normals, axis=1), 1e-300)[:, None]
	tree = cKDTree(centroids)
	distance, neighbor = tree.query(centroids, k=min(64, len(centroids)), distance_upper_bound=search)
	sample_points, sample_sizes = [], []
	for index in range(len(centroids)):
		valid = np.isfinite(distance[index])
		candidates = neighbor[index][valid]
		facing = (normals[candidates]@normals[index]) < -0.5
		if facing.any():
			sample_points.append(centroids[index])
			sample_sizes.append(factor*float(distance[index][valid][facing].min()))
	return np.array(sample_points).reshape(-1, 3), np.array(sample_sizes)


def remesh(points, faces, labels, caps, target, minimum, radius_factor, centerline, angle_degrees,
		algorithm=6, thin=None):
	"""Remesh the wall with Gmsh, then mesh every cap as a planar surface on
	the remeshed rim nodes, so boundary labels never depend on feature angles."""
	import numpy as np
	try:
		import gmsh
	except ImportError as error:
		reject(f"Gmsh Python module is unavailable: {error}")
	from scipy.spatial import cKDTree
	tree = cKDTree(centerline) if centerline is not None else None
	thin_tree = cKDTree(thin[0]) if thin is not None and len(thin[0]) else None
	wall = faces[labels == 0]
	used = np.unique(wall)
	local = -np.ones(len(points), dtype=np.int64)
	local[used] = np.arange(len(used))
	gmsh.initialize(sys.argv[:1])
	try:
		gmsh.option.setNumber("General.Terminal", 0)
		gmsh.model.add("wall")
		gmsh.model.addDiscreteEntity(2, 1)
		gmsh.model.mesh.addNodes(2, 1, list(range(1, len(used)+1)), points[used].ravel().tolist())
		gmsh.model.mesh.addElementsByType(1, 2, [], (local[wall]+1).ravel().tolist())
		gmsh.model.mesh.classifySurfaces(angle_degrees*math.pi/180.0, True, True, math.pi)
		gmsh.model.mesh.createGeometry()
		patches = len(gmsh.model.getEntities(2))

		def size(dimension, tag, x, y, z, lc=0.0):  # Gmsh < 4.11 omits lc
			value = target
			if tree is not None:
				# On the wall, the distance to the centerline is the local lumen radius.
				distance, _ = tree.query((x, y, z))
				value = min(value, radius_factor*distance)
			if thin_tree is not None:
				distance, nearest = thin_tree.query((x, y, z))
				if distance <= target:
					value = min(value, thin[1][nearest])
			return max(minimum, value)

		gmsh.model.mesh.setSizeCallback(size)
		for name, value in (("Mesh.MeshSizeExtendFromBoundary", 0), ("Mesh.MeshSizeFromPoints", 0),
				("Mesh.MeshSizeFromCurvature", 0), ("Mesh.Algorithm", algorithm)):
			gmsh.option.setNumber(name, value)
		gmsh.model.mesh.generate(2)
		node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
		index_of = {int(tag): index for index, tag in enumerate(node_tags)}
		new_points = [list(point) for point in np.asarray(coordinates, dtype=float).reshape(-1, 3)]
		new_faces, new_labels = [], []
		for _, surface in gmsh.model.getEntities(2):
			types, _, nodes = gmsh.model.mesh.getElements(2, surface)
			if list(types) != [2]:
				reject("remeshed wall contains non-triangle elements")
			for triangle in np.asarray(nodes[0], dtype=np.int64).reshape(-1, 3):
				new_faces.append([index_of[int(v)] for v in triangle])
				new_labels.append(0)
		wall_faces = np.array(new_faces, dtype=np.int64)
		# Gmsh orients classified patches independently; loops need one orientation.
		closed = orient_consistently(wall_faces)
		loops = boundary_loops(closed)
		if len(loops) != len(caps):
			reject(f"remeshed wall has {len(loops)} boundary loops for {len(caps)} caps")
		coordinates = np.array(new_points)
		assignment = {}
		for index, loop in enumerate(loops):
			ring = coordinates[loop]
			deviation = [float(np.abs((ring-np.array(cap["plane_point"]))@np.array(cap["plane_normal"])).max())
				+float(np.linalg.norm(ring.mean(axis=0)-np.array(cap["plane_point"])))*1e-6 for cap in caps]
			best = int(np.argmin(deviation))
			# Rim nodes sit on the planar cut polyline up to Gmsh's curve
			# placement round-off; the next plane must be far away.
			tolerance = 1e-3*minimum
			second = min((value for k, value in enumerate(deviation) if k != best), default=math.inf)
			if deviation[best] > tolerance or second <= 10.0*tolerance or best in assignment.values():
				order = np.argsort(deviation)[:3]
				reject(f"remeshed rim loop {index} ({len(loop)} nodes) does not lie on exactly one "
					f"port plane: nearest planes {[(int(k), deviation[k]) for k in order]}, "
					f"already assigned {best in assignment.values()}")
			assignment[index] = best
		for index, loop in enumerate(loops):
			cap = caps[assignment[index]]
			normal = np.array(cap["plane_normal"])
			origin = np.array(cap["plane_point"])
			gmsh.model.add(f"cap_{index}")
			ring = [coordinates[v]-((coordinates[v]-origin)@normal)*normal for v in loop]
			spacing = [0.5*(np.linalg.norm(ring[k]-ring[k-1])+np.linalg.norm(ring[(k+1) % len(ring)]-ring[k]))
				for k in range(len(ring))]
			point_tags = [gmsh.model.geo.addPoint(*ring[k], spacing[k]) for k in range(len(ring))]
			line_tags = [gmsh.model.geo.addLine(point_tags[k], point_tags[(k+1) % len(ring)])
				for k in range(len(ring))]
			surface = gmsh.model.geo.addPlaneSurface([gmsh.model.geo.addCurveLoop(line_tags)])
			gmsh.model.geo.synchronize()
			for line in line_tags:
				gmsh.model.mesh.setTransfiniteCurve(line, 2)
			gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
			gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
			gmsh.model.mesh.generate(2)
			cap_index = {}
			for k, point in enumerate(point_tags):
				tags, _, _ = gmsh.model.mesh.getNodes(0, point)
				cap_index[int(tags[0])] = loop[k]
			tags, coordinates_cap, _ = gmsh.model.mesh.getNodes(2, surface)
			for k, tag in enumerate(tags):
				cap_index[int(tag)] = len(new_points)
				new_points.append(list(coordinates_cap[3*k:3*k+3]))
			types, _, nodes = gmsh.model.mesh.getElements(2, surface)
			if list(types) != [2]:
				reject("remeshed cap contains non-triangle elements")
			for triangle in np.asarray(nodes[0], dtype=np.int64).reshape(-1, 3):
				new_faces.append([cap_index[int(v)] for v in triangle])
				new_labels.append(cap["label"])
			for vertex, point in zip(loop, ring):
				new_points[vertex] = list(point)
			coordinates = np.array(new_points)
		version = gmsh.option.getString("General.Version")
	finally:
		gmsh.finalize()
	return (np.array(new_points), np.array(new_faces, dtype=np.int64),
		np.array(new_labels, dtype=np.int64), version, patches)


def orient_consistently(faces):
	"""Flip triangles so neighbors traverse shared edges oppositely."""
	import numpy as np
	edge_faces = {}
	for index, face in enumerate(faces):
		for corner in range(3):
			key = (min(face[corner], face[(corner+1) % 3]), max(face[corner], face[(corner+1) % 3]))
			edge_faces.setdefault(key, []).append(index)
	if any(len(users) > 2 for users in edge_faces.values()):
		reject("remeshed surface has an edge shared by more than two triangles")
	oriented = faces.copy()
	visited = np.zeros(len(faces), dtype=bool)
	for start in range(len(faces)):
		if visited[start]:
			continue
		visited[start] = True
		stack = [start]
		while stack:
			current = stack.pop()
			face = oriented[current]
			for corner in range(3):
				a, b = int(face[corner]), int(face[(corner+1) % 3])
				for neighbor in edge_faces[(min(a, b), max(a, b))]:
					if neighbor == current or visited[neighbor]:
						continue
					other = oriented[neighbor]
					if any(int(other[k]) == a and int(other[(k+1) % 3]) == b for k in range(3)):
						oriented[neighbor] = other[[0, 2, 1]]
					visited[neighbor] = True
					stack.append(neighbor)
	return oriented


def maximum_fold(points, faces):
	"""Largest angle (degrees) between the normals of two edge-adjacent
	triangles, and where it occurs. Near 180 degrees is a folded-back flap,
	which exact intersection predicates do not flag but volume meshers reject."""
	import numpy as np
	normal = np.cross(points[faces[:, 1]]-points[faces[:, 0]], points[faces[:, 2]]-points[faces[:, 0]])
	normal /= np.maximum(np.linalg.norm(normal, axis=1), 1e-300)[:, None]
	edges = {}
	for index, face in enumerate(faces):
		for corner in range(3):
			edges.setdefault((min(face[corner], face[(corner+1) % 3]),
				max(face[corner], face[(corner+1) % 3])), []).append(index)
	worst, where = 0.0, None
	for edge, users in edges.items():
		if len(users) != 2:
			continue
		angle = math.degrees(math.acos(max(-1.0, min(1.0, float(normal[users[0]]@normal[users[1]])))))
		if angle > worst:
			worst, where = angle, points[edge[0]].tolist()
	return worst, where


def orient_outward(points, faces):
	"""Make every triangle agree with its neighbors and the enclosed volume positive."""
	edge_faces = {}
	for face in faces:
		for corner in range(3):
			key = (min(face[corner], face[(corner+1) % 3]), max(face[corner], face[(corner+1) % 3]))
			edge_faces[key] = edge_faces.get(key, 0)+1
	if any(count != 2 for count in edge_faces.values()):
		reject("remeshed surface is not a closed two-manifold")
	oriented = orient_consistently(faces)
	volume, _ = surface_volume_and_area(points, oriented)
	if volume < 0.0:
		oriented = oriented[:, [0, 2, 1]]
	return oriented


def surface_distance(points, reference_points, reference_faces):
	import numpy as np
	import vtk
	from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
	data = vtk.vtkPolyData()
	vtk_points = vtk.vtkPoints()
	vtk_points.SetData(numpy_to_vtk(np.ascontiguousarray(reference_points), deep=True))
	data.SetPoints(vtk_points)
	cells = vtk.vtkCellArray()
	cells.SetData(numpy_to_vtkIdTypeArray(np.arange(0, 3*len(reference_faces)+1, 3, dtype=np.int64), deep=True),
		numpy_to_vtkIdTypeArray(np.ascontiguousarray(reference_faces.ravel()), deep=True))
	data.SetPolys(cells)
	locator = vtk.vtkStaticCellLocator()
	locator.SetDataSet(data)
	locator.BuildLocator()
	closest = [0.0, 0.0, 0.0]
	cell, sub, squared = vtk.reference(0), vtk.reference(0), vtk.reference(0.0)
	result = np.empty(len(points))
	for index, point in enumerate(points):
		locator.FindClosestPoint(point.tolist(), closest, cell, sub, squared)
		result[index] = math.sqrt(float(squared))
	return result


def main():
	start = time.monotonic()
	parser = argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("input_surface", help="closed lumen surface (.stl or .vtp)")
	parser.add_argument("ports", help="ports JSON (schema_version 1)")
	parser.add_argument("output_vtp", help="labelled closed surface for surface_to_fem_volume.py")
	parser.add_argument("--manifest", help="output manifest (default: OUTPUT.json)")
	parser.add_argument("--length-scale-to-m", type=float, default=1.0,
		help="input length unit in metres; only used for reporting")
	parser.add_argument("--smoothing-iterations", type=int, default=0)
	parser.add_argument("--smoothing-passband", type=float, default=0.05)
	parser.add_argument("--max-smoothing-displacement", type=float,
		help="reject if any vertex moves farther than this (input units)")
	parser.add_argument("--cut-offsets-radii", default="1,1.5,2,3,4,6,8",
		help="plane offsets inward from each port point, in port radii, tried in order "
		"(a port's own cut_offset_radii, when given, is tried alone)")
	parser.add_argument("--loop-radii", type=float, default=3.0,
		help="every cut-loop vertex must lie within this many port radii of the plane center")
	parser.add_argument("--extent-radii", type=float, default=40.0,
		help="the removed blind end must lie within this many port radii of the plane center")
	parser.add_argument("--offset-fallback", choices=("largest-cap", "reject"), default="largest-cap",
		help="when no offset passes the cap checks, use the valid cut with the largest cap or fail")
	parser.add_argument("--minimum-cap-fraction", type=float, default=0.5,
		help="reject cuts whose cap equivalent radius is below this fraction of the port radius")
	parser.add_argument("--dome-fraction", type=float, default=0.5,
		help="reject cuts whose cap radius is below this fraction of the cross-section one "
		"port radius deeper (the plane is still in the rounded blind end)")
	parser.add_argument("--snap-fraction", type=float, default=0.3,
		help="snap vertices closer to the plane than this fraction of their shortest edge")
	parser.add_argument("--target-size", type=float,
		help="remesh with this maximum edge size (input units); omit to keep the cut triangulation")
	parser.add_argument("--minimum-size", type=float, help="smallest remesh size (default target/16)")
	parser.add_argument("--centerline-swc",
		help="SWC centerline; the distance from the wall to it (the local radius) bounds the remesh size")
	parser.add_argument("--swc-length-scale", type=float, default=1.0,
		help="SWC unit in input surface units")
	parser.add_argument("--radius-factor", type=float, default=1.0,
		help="remesh size is at most this factor times the local lumen radius")
	parser.add_argument("--thickness-factor", type=float, default=0.5,
		help="remesh size is at most this factor times the local wall-to-wall distance "
		"(0 disables; catches pinches the centerline misses)")
	parser.add_argument("--maximum-fold-degrees", type=float, default=170.0,
		help="reject remeshed surfaces whose adjacent triangle normals differ by more than this")
	parser.add_argument("--feature-angle-degrees", type=float, default=60.0,
		help="dihedral angle at which Gmsh splits the wall into reparametrizable patches")
	args = parser.parse_args()
	input_path = Path(args.input_surface).resolve()
	ports_path = Path(args.ports).resolve()
	output_path = Path(args.output_vtp).resolve()
	manifest_path = Path(args.manifest).resolve() if args.manifest else output_path.with_suffix(".json")
	if output_path.suffix.lower() != ".vtp":
		reject("output must be a .vtp file")
	if not math.isfinite(args.length_scale_to_m) or args.length_scale_to_m <= 0.0:
		reject("length scale must be finite and positive")
	if args.smoothing_iterations < 0 or not 0.0 < args.smoothing_passband <= 2.0:
		reject("smoothing needs nonnegative iterations and a passband in (0,2]")
	if (not args.loop_radii > 1.0 or not args.extent_radii >= args.loop_radii
			or not 0.0 <= args.snap_fraction < 0.5 or not 0.0 <= args.minimum_cap_fraction < 1.0
			or not 0.0 <= args.dome_fraction < 1.0):
		reject("loop radii must exceed 1, extent radii must be at least loop radii, snap fraction "
			"must be in [0,0.5), and the minimum cap fraction in [0,1)")
	try:
		offsets = [float(value) for value in args.cut_offsets_radii.split(",")]
	except ValueError:
		reject("cut offsets must be comma-separated numbers")
	if not offsets or any(not math.isfinite(value) or value < 0.0 for value in offsets):
		reject("cut offsets must be finite and nonnegative")
	if args.target_size is not None and (not math.isfinite(args.target_size) or args.target_size <= 0.0):
		reject("target size must be finite and positive")
	if args.centerline_swc and args.target_size is None:
		reject("a centerline only bounds the remesh size; also pass --target-size")
	if not 0.0 < args.feature_angle_degrees < 90.0 or not args.radius_factor > 0.0:
		reject("feature angle must be in (0,90) degrees and radius factor positive")
	ports = read_ports(ports_path)
	progress(f"reading {input_path.name}")
	points, faces = read_surface(input_path)
	source_volume, source_area = surface_volume_and_area(points, faces)
	if not source_volume > 0.0:
		reject("input surface must be closed and outward oriented (positive enclosed volume)")
	smoothing = {"iterations": args.smoothing_iterations}
	if args.smoothing_iterations:
		progress(f"smoothing {args.smoothing_iterations} iterations")
		smoothed = smooth(points, faces, args.smoothing_iterations, args.smoothing_passband)
		displacement = ((smoothed-points)**2).sum(axis=1)**0.5
		smoothing.update({"passband": args.smoothing_passband,
			"maximum_displacement": float(displacement.max()),
			"median_displacement": float(sorted(displacement)[len(displacement)//2]),
			"p95_displacement": float(sorted(displacement)[int(0.95*(len(displacement)-1))])})
		if (args.max_smoothing_displacement is not None
				and displacement.max() > args.max_smoothing_displacement):
			reject(f"smoothing moved a vertex {displacement.max():.4g} > "
				f"{args.max_smoothing_displacement:.4g}")
		points = smoothed
	reference_points, reference_faces = points.copy(), faces.copy()
	cutter = SurfaceCutter(points, faces)
	port_reports = []
	for port in ports:
		failures, valid, rule = [], [], "first_offset_passing_cap_checks"
		for offset in ([port["offset"]] if port["offset"] is not None else offsets):
			result, failure = cutter.cut(port, offset, args.extent_radii, args.loop_radii,
				args.snap_fraction, 0.0)
			if result is not None and port["offset"] is None:
				valid.append((offset, result))
				# A plane through the rounded or tapering blind end gives a cap much
				# smaller than the port radius or than the tube one radius deeper.
				deeper, _ = cutter.cut(port, offset+1.0, args.extent_radii, args.loop_radii,
					args.snap_fraction, 0.0)
				if result["equivalent_radius"] < args.minimum_cap_fraction*port["radius"]:
					failure = (f"cap radius {result['equivalent_radius']:.3g} is below "
						f"{args.minimum_cap_fraction} port radii")
					result = None
				elif deeper is not None and result["equivalent_radius"] < args.dome_fraction*deeper["equivalent_radius"]:
					failure = (f"cap radius {result['equivalent_radius']:.3g} is below {args.dome_fraction} "
						f"of the radius {deeper['equivalent_radius']:.3g} one port radius deeper (rounded tip)")
					result = None
			if result is not None:
				break
			failures.append(f"offset {offset}: {failure}")
		if result is None and valid and args.offset_fallback == "largest-cap":
			offset, result = max(valid, key=lambda item: item[1]["equivalent_radius"])
			rule = "fallback_largest_valid_cap"
			progress(f"warning: {port['name']} uses the largest valid cap (offset {offset})")
		if result is None:
			reject(f"port {port['name']} (label {port['label']}) cannot be cut: "+"; ".join(failures))
		area = cutter.apply(result, port["label"])
		port_reports.append({"name": port["name"], "label": port["label"], "offset_radii": offset,
			"plane_point": result["center"].tolist(), "plane_normal": result["tangent"].tolist(),
			"loop_vertices": len(result["loop"]), "cap_area": area,
			"equivalent_radius": math.sqrt(area/math.pi), "removed_triangles": result["region_faces"],
			"snapped_vertices": len(result["snapped"]),
			"maximum_snap_distance": max((value for _, value in result["snapped"]), default=0.0),
			"rejected_offsets": failures, "selection": rule if port["offset"] is None else "port_cut_offset"})
		progress(f"capped {port['name']} label {port['label']} offset {offset} "
			f"loop {len(result['loop'])} radius {math.sqrt(area/math.pi):.4g}")
	points, faces, labels = cutter.arrays()
	remesh_report = None
	if args.target_size is not None:
		minimum = args.minimum_size if args.minimum_size is not None else args.target_size/16.0
		centerline = (read_centerline_samples(args.centerline_swc, args.swc_length_scale, minimum/2.0)
			if args.centerline_swc else None)
		progress(f"remeshing to target {args.target_size:g}")
		cut_points, cut_faces = points, faces
		caps = [{"label": report["label"], "plane_point": report["plane_point"],
			"plane_normal": report["plane_normal"]} for report in port_reports]
		# Gmsh's 2D algorithms occasionally fold a tiny wall patch back on
		# itself; retry with the other algorithms before failing closed.
		attempts = []
		thin = thin_wall_samples(cut_points, cut_faces, 2.0*args.target_size, args.thickness_factor) \
			if args.thickness_factor > 0.0 else None
		for algorithm in (6, 5, 1):
			candidate = remesh(cut_points, cut_faces, labels, caps, args.target_size,
				minimum, args.radius_factor, centerline, args.feature_angle_degrees, algorithm, thin)
			oriented = orient_outward(candidate[0], candidate[1])
			fold, where = maximum_fold(candidate[0], oriented)
			attempts.append({"gmsh_algorithm": algorithm, "maximum_fold_degrees": fold})
			if fold <= args.maximum_fold_degrees:
				break
			progress(f"warning: Gmsh algorithm {algorithm} folded the surface ({fold:.1f} deg at {where})")
		else:
			reject(f"every remesh folds the surface (worst {fold:.1f} deg at {where}); "
				"change --target-size or --feature-angle-degrees")
		points, faces, labels, version, patches = candidate
		faces = oriented
		distance = surface_distance(points, cut_points, cut_faces)
		remesh_report = {"mesher": "Gmsh", "version": version, "target_size": args.target_size,
			"minimum_size": minimum, "radius_factor": args.radius_factor,
			"centerline_swc": str(Path(args.centerline_swc).resolve()) if args.centerline_swc else None,
			"centerline_swc_sha256": file_sha256(args.centerline_swc) if args.centerline_swc else None,
			"feature_angle_degrees": args.feature_angle_degrees, "patches": patches,
			"attempts": attempts, "maximum_fold_degrees_gate": args.maximum_fold_degrees,
			"thickness_factor": args.thickness_factor,
			"thin_wall_samples": 0 if thin is None else int(len(thin[0])),
			"minimum_thin_wall_size": None if thin is None or not len(thin[1]) else float(thin[1].min()),
			"maximum_distance_to_cut_surface": float(distance.max())}
	volume, area = surface_volume_and_area(points, faces)
	quality = triangle_quality(points, faces)
	edges = set()
	for face in faces:
		for corner in range(3):
			edges.add((min(face[corner], face[(corner+1) % 3]), max(face[corner], face[(corner+1) % 3])))
	import numpy as np
	edge_array = np.array(sorted(edges))
	lengths = np.linalg.norm(points[edge_array[:, 0]]-points[edge_array[:, 1]], axis=1)
	output_path.parent.mkdir(parents=True, exist_ok=True)
	write_vtp(output_path, points, faces, labels)
	label_counts = {str(int(label)): int((labels == label).sum()) for label in np.unique(labels)}
	manifest = {
		"schema_version": 1,
		"route": "cap_vessel_ports",
		"source": {"path": str(input_path), "sha256": file_sha256(input_path),
			"volume": source_volume, "area": source_area},
		"ports_file": {"path": str(ports_path), "sha256": file_sha256(ports_path)},
		"length_scale_to_m": args.length_scale_to_m,
		"smoothing": smoothing,
		"cut": {"offsets_radii": offsets, "loop_radii": args.loop_radii, "dome_fraction": args.dome_fraction,
			"extent_radii": args.extent_radii, "minimum_cap_fraction": args.minimum_cap_fraction,
			"snap_fraction": args.snap_fraction, "ports": port_reports},
		"remesh": remesh_report,
		"output": {"path": str(output_path), "sha256": file_sha256(output_path),
			"points": int(len(points)), "triangles": int(len(faces)),
			"boundary_array": "boundary_id", "boundary_labels": label_counts,
			"volume": volume, "area": area,
			"relative_volume_change": (volume-source_volume)/source_volume,
			"minimum_triangle_quality": float(quality.min()),
			"maximum_fold_degrees": maximum_fold(points, faces)[0],
			"p1_triangle_quality": float(np.percentile(quality, 1)),
			"minimum_edge_length": float(lengths.min()),
			"median_edge_length": float(np.median(lengths))},
		"geometry_change": {"implicit_repair": False,
			"operations": [name for name, used in (("windowed_sinc_smoothing", args.smoothing_iterations > 0),
				("local_plane_cut_and_cap", True), ("gmsh_remesh", args.target_size is not None)) if used]},
		"resources": {"wall_s": time.monotonic()-start,
			"peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024},
	}
	manifest_path.parent.mkdir(parents=True, exist_ok=True)
	manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True)+"\n", encoding="utf-8")
	print(f"cap_vessel_ports: PASS ports={len(ports)} triangles={len(faces)} "
		f"labels={len(label_counts)} min_quality={quality.min():.4g}")
	return 0


if __name__ == "__main__":
	try:
		raise SystemExit(main())
	except (OSError, ValueError, RuntimeError) as error:
		print(f"cap_vessel_ports: ERROR: {error}", file=sys.stderr)
		raise SystemExit(2)
