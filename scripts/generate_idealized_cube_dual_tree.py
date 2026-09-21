#!/usr/bin/env python3
"""Centerline/radius -> labelled, conforming dual-tree cube mesh (geometry only)."""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile


REGIONS = {"arterial_lumen": 101, "venous_lumen": 102,
	"arterial_wall": 103, "venous_wall": 104, "fixed_darcy_tissue": 105}
SURFACES = {"arterial_fluid_wall": 201, "venous_fluid_wall": 202,
	"arterial_wall_tissue": 203, "venous_wall_tissue": 204,
	"arterial_inlet": 207, "venous_outlet": 208,
	"tissue_exterior": 209, "arterial_wall_root_exterior": 210,
	"venous_wall_root_exterior": 218,
	**{f"arterial_terminal_{index}": 230+index for index in range(4)},
	**{f"venous_terminal_{index}": 240+index for index in range(4)}}
INTERFACES = {
	"arterial_fluid_wall": ("arterial_lumen", "arterial_wall"),
	"venous_fluid_wall": ("venous_lumen", "venous_wall"),
	"arterial_wall_tissue": ("arterial_wall", "fixed_darcy_tissue"),
	"venous_wall_tissue": ("venous_wall", "fixed_darcy_tissue"),
	**{f"arterial_terminal_{index}": ("arterial_lumen", "fixed_darcy_tissue")
		for index in range(4)},
	**{f"venous_terminal_{index}": ("venous_lumen", "fixed_darcy_tissue")
		for index in range(4)}}
REGION_SURFACES = {
	"arterial_lumen": ["arterial_fluid_wall", "arterial_inlet"]+
		[f"arterial_terminal_{index}" for index in range(4)],
	"arterial_wall": ["arterial_fluid_wall", "arterial_wall_tissue",
		"arterial_wall_root_exterior"],
	"venous_lumen": ["venous_fluid_wall", "venous_outlet"]+
		[f"venous_terminal_{index}" for index in range(4)],
	"venous_wall": ["venous_fluid_wall", "venous_wall_tissue",
		"venous_wall_root_exterior"],
	"fixed_darcy_tissue": ["tissue_exterior", "arterial_wall_tissue",
		"venous_wall_tissue"]+[f"arterial_terminal_{index}" for index in range(4)]+
		[f"venous_terminal_{index}" for index in range(4)]}


def require(condition, message):
	if not condition:
		raise ValueError(message)


def dot(a, b):
	return sum(x*y for x, y in zip(a, b))


def sub(a, b):
	return [x-y for x, y in zip(a, b)]


def cross(a, b):
	return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2],
		a[0]*b[1]-a[1]*b[0])


def add(a, b):
	return [x+y for x, y in zip(a, b)]


def scale(a, factor):
	return [factor*x for x in a]


def normalized(a):
	return scale(a, 1./math.sqrt(dot(a, a)))


def spline_capsule_paths(tree, subsegments, tension):
	"""Sample template-free cubic centerlines for every directed tree edge."""
	nodes = tree["nodes_m"]
	parents = [None]*len(nodes)
	children = [[] for _ in nodes]
	for start, stop, _ in tree["segments"]:
		parents[stop] = start
		children[start].append(stop)
	tangents = []
	for index, parent in enumerate(parents):
		incoming = None if parent is None else normalized(sub(nodes[index], nodes[parent]))
		outgoing = [normalized(sub(nodes[child], nodes[index]))
			for child in children[index]]
		if incoming is None:
			tangent = normalized([sum(direction[axis] for direction in outgoing)
				for axis in range(3)])
		elif not outgoing:
			tangent = incoming
		else:
			mean_outgoing = normalized([sum(direction[axis] for direction in outgoing)
				for axis in range(3)])
			tangent = normalized(add(incoming, mean_outgoing))
		tangents.append(tangent)
	paths = []
	for start, stop, radius in tree["segments"]:
		first, last = nodes[start], nodes[stop]
		length = math.dist(first, last)
		first_tangent = scale(tangents[start], tension*length)
		last_tangent = scale(tangents[stop], tension*length)
		points = []
		for sample in range(subsegments+1):
			u = sample/subsegments
			h00 = 2*u**3-3*u**2+1
			h10 = u**3-2*u**2+u
			h01 = -2*u**3+3*u**2
			h11 = u**3-u**2
			points.append([h00*first[axis]+h10*first_tangent[axis]+
				h01*last[axis]+h11*last_tangent[axis] for axis in range(3)])
		paths.append((start, stop, radius, points))
	return paths


def segment_distance(a, b, c, d):
	"""Minimum Euclidean distance between two finite 3D segments."""
	u, v, w = sub(b, a), sub(d, c), sub(a, c)
	aa, bb, cc, dd, ee = dot(u, u), dot(u, v), dot(v, v), dot(u, w), dot(v, w)
	denom = aa*cc-bb*bb
	s = min(1., max(0., (bb*ee-cc*dd)/denom)) if denom > 1e-24 else 0.
	t = min(1., max(0., (aa*ee-bb*dd)/denom)) if denom > 1e-24 else 0.
	for _ in range(3):
		s = min(1., max(0., (bb*t-dd)/aa))
		t = min(1., max(0., (bb*s+ee)/cc))
	return math.sqrt(sum((w[i]+s*u[i]-t*v[i])**2 for i in range(3)))


def tet_quality(points):
	def corner(a, b, c, d):
		u, v, w = sub(b, a), sub(c, a), sub(d, a)
		denom = math.sqrt(dot(u, u)*dot(v, v)*dot(w, w))
		return dot(u, cross(v, w))/denom if denom > 0 else 0.
	a, b, c, d = points
	return (corner(a, b, c, d), min(abs(corner(*order)) for order in
		((a, b, c, d), (b, a, d, c), (c, d, a, b), (d, c, b, a))))


def validate(data):
	require(data.get("schema_version") == 1 and data.get("kind") ==
		"idealized_cube_dual_tree_functional_only" and
		data.get("physiological_validation") is False and data.get("units") == "SI",
		"case must be explicitly SI, versioned and nonphysiological")
	origin, size = data["cube_origin_m"], data["cube_size_m"]
	require(len(origin) == len(size) == 3 and all(math.isfinite(x) for x in origin+size)
		and all(x > 0 for x in size), "cube dimensions must be finite and positive")
	wall = data["wall_thickness_m"]
	clearance = data["minimum_tree_clearance_m"]
	require(math.isfinite(wall) and wall > 0 and math.isfinite(clearance)
		and clearance > 0, "wall and clearance must be positive")
	require(0 < data["minimum_scaled_jacobian"] <= 1 and
		data["target_mesh_size_m"] > 0, "mesh controls are invalid")
	minimum_size = data.get("minimum_mesh_size_m", .3*data["target_mesh_size_m"])
	curvature = data.get("mesh_curvature_divisions", 8)
	vascular = data.get("vascular_geometry", {"kind": "piecewise_cylinders"})
	require(math.isfinite(minimum_size) and 0 < minimum_size <=
		data["target_mesh_size_m"] and isinstance(curvature, int) and
		not isinstance(curvature, bool) and 8 <= curvature <= 128 and
		isinstance(data.get("mesh_size_extend_from_boundary", True), bool) and
		isinstance(data.get("arterial_outer_reverse", True), bool),
		"minimum size, curvature divisions or boundary extension is invalid")
	require(vascular.get("kind") in ("piecewise_cylinders", "spline_capsules"),
		"vascular geometry kind is invalid")
	if vascular["kind"] == "spline_capsules":
		subsegments = vascular.get("subsegments_per_edge", 2)
		tension = vascular.get("tangent_scale", .25)
		cleanup = vascular.get("join_cleanup_tolerance_m", .3*wall)
		smoothing = vascular.get("surface_smoothing_iterations", 0)
		pass_band = vascular.get("surface_smoothing_pass_band", .01)
		require(isinstance(subsegments, int) and not isinstance(subsegments, bool)
			and 2 <= subsegments <= 8 and math.isfinite(tension) and 0 < tension <= .5
			and math.isfinite(cleanup) and 0 < cleanup < .5*wall
			and isinstance(smoothing, int) and not isinstance(smoothing, bool)
			and 0 <= smoothing <= 200 and math.isfinite(pass_band)
			and 0 < pass_band <= 1.,
			"spline capsule controls are invalid")
	for key in ("arterial_tree", "venous_tree"):
		tree = data[key]
		nodes, segments = tree["nodes_m"], tree["segments"]
		require(len(nodes) >= 4 and len(segments) == len(nodes)-1,
			f"{key} must be a connected tree")
		require(all(len(node) == 3 and all(math.isfinite(x) for x in node)
			for node in nodes), f"{key} has invalid node coordinates")
		root, terminals = tree["root"], tree["terminals"]
		require(root == 0 and len(terminals) == 4 and len(set(terminals)) == 4,
			f"{key} requires root 0 and four unique terminals")
		indegree = [0]*len(nodes)
		for segment in segments:
			require(len(segment) == 3, f"{key} segment must be [from,to,radius]")
			i, j, radius = segment
			require(isinstance(i, int) and isinstance(j, int) and
				0 <= i < j < len(nodes) and math.isfinite(radius) and radius > wall,
				f"{key} segment index/radius invalid")
			indegree[j] += 1
			require(dot(sub(nodes[j], nodes[i]), sub(nodes[j], nodes[i])) > 1e-12,
				f"{key} has a zero-length segment")
		require(indegree == [0]+[1]*(len(nodes)-1),
			f"{key} segments do not form one rooted tree")
		outdegree = [0]*len(nodes)
		for i, _, _ in segments:
			outdegree[i] += 1
		require(set(terminals) == {i for i, count in enumerate(outdegree) if count == 0},
			f"{key} terminal list differs from graph leaves")
		for first, (i, j, radius) in enumerate(segments):
			for k, l, other_radius in segments[first+1:]:
				if {i, j} & {k, l}:
					continue
				gap = segment_distance(nodes[i], nodes[j], nodes[k], nodes[l]) \
						-(radius+other_radius+2*wall)
				require(gap > 0., f"{key} unrelated branches overlap")
		for index, node in enumerate(nodes):
			max_radius = max(radius for i, j, radius in segments if i == index or j == index)+wall
			for axis in (1, 2):
				require(origin[axis]+max_radius < node[axis] <
					origin[axis]+size[axis]-max_radius,
					f"{key} branch leaves cube on axis {axis}")
			if index != root:
				require(origin[0]+max_radius < node[0] <
					origin[0]+size[0]-max_radius,
					f"{key} non-root node touches cube end")
		root_x = origin[0] if key == "arterial_tree" else origin[0]+size[0]
		require(abs(nodes[root][0]-root_x) < 1e-12,
			f"{key} root must open on its designated cube face")
	a, v = data["arterial_tree"], data["venous_tree"]
	for ai, aj, ar in a["segments"]:
		for vi, vj, vr in v["segments"]:
			gap = segment_distance(a["nodes_m"][ai], a["nodes_m"][aj],
				v["nodes_m"][vi], v["nodes_m"][vj])-(ar+vr+2*wall)
			require(gap >= clearance,
				f"arterial and venous outer tubes violate clearance: {gap:.6g} m")
	params = data["functional_parameters"]
	for name in ("fluid_density_kg_m3", "fluid_viscosity_pa_s",
		"darcy_mobility_m2_pa_s", "wall_young_modulus_pa", "inlet_speed_m_s"):
		require(math.isfinite(params[name]) and params[name] > 0,
			f"{name} must be finite and positive")
	require(0 < params["wall_poisson_ratio"] < .45 and
		math.isfinite(params["venous_outlet_pressure_pa"]),
		"wall Poisson ratio or outlet pressure invalid")
	return data


def add_tree(gmsh, tree, radius_add, reverse=False, geometry=None):
	nodes = tree["nodes_m"]
	parts = []
	geometry = geometry or {"kind": "piecewise_cylinders"}
	if geometry["kind"] == "spline_capsules":
		paths = spline_capsule_paths(tree, geometry.get("subsegments_per_edge", 2),
			geometry.get("tangent_scale", .25))
		for _, _, radius, points in paths:
			for index in range(len(points)-1):
				delta = sub(points[index+1], points[index])
				parts.append((3, gmsh.model.occ.addCylinder(*points[index], *delta,
					radius+radius_add)))
			for point in points[1:-1]:
				parts.append((3, gmsh.model.occ.addSphere(*point, radius+radius_add)))
	else:
		for i, j, radius in tree["segments"]:
			start, end = (nodes[j], nodes[i]) if reverse else (nodes[i], nodes[j])
			delta = sub(end, start)
			parts.append((3, gmsh.model.occ.addCylinder(*start, *delta,
				radius+radius_add)))
	for index in range(1, len(nodes)):
		if index in tree["terminals"]:
			continue
		radius = max(segment[2] for segment in tree["segments"]
			if index in segment[:2])+radius_add
		parts.append((3, gmsh.model.occ.addSphere(*nodes[index], radius)))
	fused, _ = gmsh.model.occ.fuse([parts[0]], parts[1:])
	require(len(fused) == 1 and fused[0][0] == 3,
		"centerline/radius surface did not fuse into one closed tree")
	return fused[0]


def clean_surface_triangles(nodes, triangles, tolerance):
	"""Collapse short Boolean-seam edges without changing closed-shell topology."""
	parent = {tag: tag for tag in nodes}

	def find(tag):
		while parent[tag] != tag:
			parent[tag] = parent[parent[tag]]
			tag = parent[tag]
		return tag

	def union(first, second):
		first, second = find(first), find(second)
		if first != second:
			parent[max(first, second)] = min(first, second)

	edges = set()
	for group in triangles.values():
		for triangle in group:
			for first, second in ((triangle[0], triangle[1]),
				(triangle[1], triangle[2]), (triangle[2], triangle[0])):
				edge = tuple(sorted((first, second)))
				if edge not in edges:
					edges.add(edge)
					if math.dist(nodes[first], nodes[second]) < tolerance:
						union(first, second)
	for tag in parent:
		parent[tag] = find(tag)
	members = {}
	for tag, representative in parent.items():
		members.setdefault(representative, []).append(tag)
	points = {representative: tuple(sum(nodes[tag][axis] for tag in group)/len(group)
		for axis in range(3)) for representative, group in members.items()}
	cleaned = {}
	dropped = 0
	for name, group in triangles.items():
		unique = set()
		cleaned[name] = []
		for triangle in group:
			mapped = tuple(parent[tag] for tag in triangle)
			key = tuple(sorted(mapped))
			if len(set(mapped)) < 3 or key in unique:
				dropped += 1
				continue
			unique.add(key)
			cleaned[name].append(mapped)
	for name, labels in REGION_SURFACES.items():
		usage = {}
		for label in labels:
			for triangle in cleaned[label]:
				for first, second in ((triangle[0], triangle[1]),
					(triangle[1], triangle[2]), (triangle[2], triangle[0])):
					edge = tuple(sorted((first, second)))
					usage[edge] = usage.get(edge, 0)+1
		require(usage and all(count == 2 for count in usage.values()),
			f"surface cleanup opened the {name} shell")
	return points, cleaned, len(nodes)-len(points), dropped


def smooth_vessel_surfaces(points, triangles, iterations, pass_band):
	try:
		import vtk
	except ImportError as error:
		raise ValueError(f"VTK Python module unavailable: {error}") from error
	for name in ("arterial_fluid_wall", "arterial_wall_tissue",
		"venous_fluid_wall", "venous_wall_tissue"):
		group = triangles[name]
		tags = sorted({tag for triangle in group for tag in triangle})
		local = {tag: index for index, tag in enumerate(tags)}
		edge_use = {}
		for triangle in group:
			for first, second in ((triangle[0], triangle[1]),
				(triangle[1], triangle[2]), (triangle[2], triangle[0])):
				edge = tuple(sorted((first, second)))
				edge_use[edge] = edge_use.get(edge, 0)+1
		boundary = {tag for edge, count in edge_use.items() if count == 1
			for tag in edge}
		vtk_points = vtk.vtkPoints()
		for tag in tags:
			vtk_points.InsertNextPoint(points[tag])
		polygons = vtk.vtkCellArray()
		for triangle in group:
			polygons.InsertNextCell(3)
			for tag in triangle:
				polygons.InsertCellPoint(local[tag])
		polydata = vtk.vtkPolyData()
		polydata.SetPoints(vtk_points)
		polydata.SetPolys(polygons)
		smoother = vtk.vtkWindowedSincPolyDataFilter()
		smoother.SetInputData(polydata)
		smoother.SetNumberOfIterations(iterations)
		smoother.SetPassBand(pass_band)
		smoother.BoundarySmoothingOff()
		smoother.FeatureEdgeSmoothingOff()
		smoother.NonManifoldSmoothingOn()
		smoother.NormalizeCoordinatesOn()
		smoother.Update()
		output = smoother.GetOutput().GetPoints()
		for index, tag in enumerate(tags):
			if tag not in boundary:
				points[tag] = output.GetPoint(index)
	return points


def rebuild_discrete_geometry(gmsh, surface_groups, geometry, wall_thickness):
	"""Replace sliver-prone OCC seams with a conforming discrete surface complex."""
	node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
	nodes = {int(tag): tuple(coordinates[3*index:3*index+3])
		for index, tag in enumerate(node_tags)}
	triangles = {name: [] for name in SURFACES}
	for name, surfaces in surface_groups.items():
		for surface in surfaces:
			types, _, connectivity = gmsh.model.mesh.getElements(2, surface)
			for element_type, flat in zip(types, connectivity):
				require(element_type == 2, "surface mesh is not first-order triangular")
				triangles[name].extend(tuple(int(tag) for tag in flat[start:start+3])
					for start in range(0, len(flat), 3))
	points, triangles, collapsed, dropped = clean_surface_triangles(nodes,
		triangles, geometry.get("join_cleanup_tolerance_m", .3*wall_thickness))
	iterations = geometry.get("surface_smoothing_iterations", 0)
	if iterations:
		points = smooth_vessel_surfaces(points, triangles, iterations,
			geometry.get("surface_smoothing_pass_band", .01))
	gmsh.clear()
	gmsh.model.add("idealized-cube-dual-tree-discrete")
	discrete_surfaces = {}
	for name, entity in SURFACES.items():
		gmsh.model.addDiscreteEntity(2, entity)
		discrete_surfaces[name] = [entity]
	tags = sorted(points)
	gmsh.model.mesh.addNodes(2, next(iter(SURFACES.values())), tags,
		[value for tag in tags for value in points[tag]])
	element = 1
	for name, entity in SURFACES.items():
		count = len(triangles[name])
		gmsh.model.mesh.addElementsByType(entity, 2,
			list(range(element, element+count)),
			[tag for triangle in triangles[name] for tag in triangle])
		element += count
		group = gmsh.model.addPhysicalGroup(2, [entity], SURFACES[name])
		gmsh.model.setPhysicalName(2, group, name)
	regions = {}
	for name, labels in REGION_SURFACES.items():
		loop = gmsh.model.geo.addSurfaceLoop([SURFACES[label] for label in labels])
		volume = gmsh.model.geo.addVolume([loop], REGIONS[name])
		regions[name] = {volume}
	gmsh.model.geo.synchronize()
	for name, tags in regions.items():
		group = gmsh.model.addPhysicalGroup(3, list(tags), REGIONS[name])
		gmsh.model.setPhysicalName(3, group, name)
	return regions, discrete_surfaces, collapsed, dropped


def generate(data, output_prefix):
	try:
		import gmsh
	except ImportError as error:
		raise ValueError(f"Gmsh Python module unavailable: {error}") from error
	gmsh.initialize()
	try:
		gmsh.option.setNumber("General.Terminal", 0)
		gmsh.model.add("idealized-cube-dual-tree")
		origin, size = data["cube_origin_m"], data["cube_size_m"]
		vascular = data.get("vascular_geometry", {"kind": "piecewise_cylinders"})
		cube = (3, gmsh.model.occ.addBox(*origin, *size))
		artery = add_tree(gmsh, data["arterial_tree"], 0., reverse=True,
			geometry=vascular)
		artery_outer = add_tree(gmsh, data["arterial_tree"],
			data["wall_thickness_m"],
			reverse=data.get("arterial_outer_reverse", True), geometry=vascular)
		vein = add_tree(gmsh, data["venous_tree"], 0., geometry=vascular)
		vein_outer = add_tree(gmsh, data["venous_tree"],
			data["wall_thickness_m"], geometry=vascular)
		tissue_shape, _ = gmsh.model.occ.cut([cube],
			[artery_outer, vein_outer], removeTool=False)
		artery_wall_shape, _ = gmsh.model.occ.cut([artery_outer],
			[artery], removeTool=False)
		vein_wall_shape, _ = gmsh.model.occ.cut([vein_outer],
			[vein], removeTool=False)
		require(len(tissue_shape) == len(artery_wall_shape) == len(vein_wall_shape) == 1,
			"Boolean cut did not yield one tissue and two walls")
		parts, mapping = gmsh.model.occ.fragment(tissue_shape,
			[artery, artery_wall_shape[0], vein, vein_wall_shape[0]])
		gmsh.model.occ.synchronize()
		require(all(dim == 3 for dim, _ in parts), "fragment created non-volume part")
		mapped = [{tag for dim, tag in group if dim == 3} for group in mapping]
		require(len(mapped) == 5, "fragment correspondence is incomplete")
		lumen_a = mapped[1]
		wall_a = mapped[2]
		lumen_v = mapped[3]
		wall_v = mapped[4]
		tissue = mapped[0]
		regions = {"arterial_lumen": lumen_a, "venous_lumen": lumen_v,
			"arterial_wall": wall_a, "venous_wall": wall_v,
			"fixed_darcy_tissue": tissue}
		all_parts = {tag for dim, tag in parts if dim == 3}
		require(all(len(tags) == 1 for tags in regions.values()) and
			set().union(*regions.values()) == all_parts and
			sum(map(len, regions.values())) == len(all_parts),
			"Boolean fragment did not yield five disjoint connected regions")
		volume_sum = sum(gmsh.model.occ.getMass(3, next(iter(tags)))
			for tags in regions.values())
		cube_volume = math.prod(size)
		require(abs(volume_sum-cube_volume) <= 1e-5*cube_volume,
			"Boolean regions do not partition the cube volume")
		for name, tags in regions.items():
			group = gmsh.model.addPhysicalGroup(3, list(tags), REGIONS[name])
			gmsh.model.setPhysicalName(3, group, name)
		by_volume = {next(iter(tags)): name for name, tags in regions.items()}
		surface_groups = {name: [] for name in SURFACES}
		for _, surface in gmsh.model.getEntities(2):
			up, _ = gmsh.model.getAdjacencies(2, surface)
			adj = frozenset(by_volume[int(tag)] for tag in up)
			pairs = {
				frozenset(("arterial_lumen", "arterial_wall")): "arterial_fluid_wall",
				frozenset(("venous_lumen", "venous_wall")): "venous_fluid_wall",
				frozenset(("arterial_wall", "fixed_darcy_tissue")): "arterial_wall_tissue",
				frozenset(("venous_wall", "fixed_darcy_tissue")): "venous_wall_tissue"}
			if adj in pairs:
				label = pairs[adj]
			elif adj in (frozenset(("arterial_lumen", "fixed_darcy_tissue")),
				frozenset(("venous_lumen", "fixed_darcy_tissue"))):
				name = "arterial" if "arterial_lumen" in adj else "venous"
				tree = data[name+"_tree"]
				center = gmsh.model.occ.getCenterOfMass(2, surface)
				distances = [math.dist(center, tree["nodes_m"][tip])
					for tip in tree["terminals"]]
				index = min(range(len(distances)), key=distances.__getitem__)
				require(distances[index] < 1e-8,
					f"{name} terminal cap is not centered on a declared leaf")
				label = f"{name}_terminal_{index}"
			elif adj == frozenset(("arterial_lumen",)):
				label = "arterial_inlet"
			elif adj == frozenset(("venous_lumen",)):
				label = "venous_outlet"
			elif adj == frozenset(("fixed_darcy_tissue",)):
				label = "tissue_exterior"
			elif adj == frozenset(("arterial_wall",)):
				label = "arterial_wall_root_exterior"
			elif adj == frozenset(("venous_wall",)):
				label = "venous_wall_root_exterior"
			else:
				raise ValueError(f"unclassified surface {surface} adjacent to {sorted(adj)}")
			surface_groups[label].append(surface)
		for name, surfaces in surface_groups.items():
			require(surfaces, f"missing {name} geometric surface")
			group = gmsh.model.addPhysicalGroup(2, surfaces, SURFACES[name])
			gmsh.model.setPhysicalName(2, group, name)
		for name in ("arterial", "venous"):
			for index in range(4):
				require(len(surface_groups[f"{name}_terminal_{index}"]) == 1,
					f"{name} terminal {index} must have one cap surface")
		gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
		gmsh.option.setNumber("Mesh.Binary", 0)
		gmsh.option.setNumber("Mesh.ElementOrder", 1)
		gmsh.option.setNumber("Mesh.MeshSizeMin", data.get("minimum_mesh_size_m",
			data["target_mesh_size_m"]*.3))
		gmsh.option.setNumber("Mesh.MeshSizeMax", data["target_mesh_size_m"])
		gmsh.option.setNumber("Mesh.MeshSizeFromCurvature",
			data.get("mesh_curvature_divisions", 8))
		gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary",
			int(data.get("mesh_size_extend_from_boundary", True)))
		gmsh.option.setNumber("Mesh.Algorithm3D", 1)
		gmsh.model.mesh.generate(2)
		collapsed_nodes = 0
		dropped_triangles = 0
		if vascular["kind"] == "spline_capsules":
			regions, surface_groups, collapsed_nodes, dropped_triangles = \
				rebuild_discrete_geometry(gmsh, surface_groups, vascular,
					data["wall_thickness_m"])
		with tempfile.TemporaryDirectory(prefix="dual-tree-", dir=output_prefix.parent) as temporary:
			base = Path(temporary)
			surface_file = base/"surface.msh"
			volume_file = base/"volume.msh"
			contract_file = base/"contract.json"
			manifest_file = base/"manifest.json"
			gmsh.write(str(surface_file))
			gmsh.model.mesh.generate(3)
			types, tags, connectivity = gmsh.model.mesh.getElements(3)
			require(len(types) == 1 and types[0] == 4 and len(tags[0]) > 0,
				"volume mesh is not first-order tetrahedral")
			node_tags, flat_coords, _ = gmsh.model.mesh.getNodes()
			coords = {int(tag): flat_coords[3*i:3*i+3] for i, tag in
				enumerate(node_tags)}
			min_scaled = 1.
			for start in range(0, len(connectivity[0]), 4):
				points = [coords[int(tag)] for tag in connectivity[0][start:start+4]]
				signed, scaled = tet_quality(points)
				require(signed > 0., "tetrahedron has nonpositive Jacobian")
				min_scaled = min(min_scaled, scaled)
			require(min_scaled >= data["minimum_scaled_jacobian"],
				f"minimum scaled Jacobian {min_scaled:.6g} below gate")
			region_volumes = {}
			for name, entities in regions.items():
				volume = 0.
				for entity in entities:
					types, _, blocks = gmsh.model.mesh.getElements(3, entity)
					require(len(types) == 1 and types[0] == 4,
						f"{name} is not first-order tetrahedral")
					for start in range(0, len(blocks[0]), 4):
						points = [coords[int(tag)] for tag in blocks[0][start:start+4]]
						volume += dot(sub(points[1], points[0]),
							cross(sub(points[2], points[0]), sub(points[3], points[0])))/6.
				require(volume > 0., f"{name} has nonpositive volume")
				region_volumes[name] = volume
			gmsh.write(str(volume_file))
			contract = {"schema_version": 1, "length_unit": "m",
				"regions": [{"physical_name": name,
					"role": "porous" if name == "fixed_darcy_tissue" else
						("solid" if name.endswith("wall") else "fluid"),
					"require_connected": True} for name in REGIONS],
				"boundaries": [{"physical_name": "tissue_exterior",
					"boundary_id": SURFACES["tissue_exterior"],
					"adjacent_region": "fixed_darcy_tissue",
					"semantic_role": "fixed_exterior_zero_flux_or_explicit_pressure"},
					{"physical_name": "arterial_inlet", "boundary_id": SURFACES["arterial_inlet"],
					"adjacent_region": "arterial_lumen", "semantic_role": "fluid_inlet"},
					{"physical_name": "venous_outlet", "boundary_id": SURFACES["venous_outlet"],
					"adjacent_region": "venous_lumen", "semantic_role": "fluid_outlet"},
					{"physical_name": "arterial_wall_root_exterior",
					"boundary_id": SURFACES["arterial_wall_root_exterior"],
					"adjacent_region": "arterial_wall", "semantic_role": "fixed_wall_support"},
					{"physical_name": "venous_wall_root_exterior",
					"boundary_id": SURFACES["venous_wall_root_exterior"],
					"adjacent_region": "venous_wall", "semantic_role": "fixed_wall_support"}],
				"interfaces": [{"physical_name": name, "regions": list(pair)}
					for name, pair in INTERFACES.items()]}
			contract_file.write_text(json.dumps(contract, indent=2,
				sort_keys=True)+"\n", encoding="utf-8")
			manifest = {"schema_version": 1,
				"kind": "idealized_cube_dual_tree_geometry_only",
				"physiological_validation": False,
				"case_sha256": hashlib.sha256(json.dumps(data, sort_keys=True,
					separators=(",", ":")).encode()).hexdigest(),
				"gmsh_version": gmsh.option.getString("General.Version"),
				"vascular_geometry": vascular,
				"regions": REGIONS, "surfaces": SURFACES,
				"surface_count": {name: len(surfaces)
					for name, surfaces in surface_groups.items()},
				"region_volume_m3": region_volumes,
				"collapsed_surface_nodes": collapsed_nodes,
				"dropped_surface_triangles": dropped_triangles,
				"tetra_count": len(tags[0]),
				"minimum_scaled_jacobian": min_scaled,
				"surface_sha256": hashlib.sha256(surface_file.read_bytes()).hexdigest(),
				"volume_sha256": hashlib.sha256(volume_file.read_bytes()).hexdigest(),
				"contract_sha256": hashlib.sha256(contract_file.read_bytes()).hexdigest()}
			manifest_file.write_text(json.dumps(manifest, indent=2,
				sort_keys=True)+"\n", encoding="utf-8")
			for source, suffix in ((surface_file, ".surface.msh"),
				(volume_file, ".volume.msh"), (contract_file, ".contract.json"),
				(manifest_file, ".json")):
				os.replace(source, Path(str(output_prefix)+suffix))
		return manifest
	finally:
		gmsh.finalize()


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("output_prefix", type=Path)
	parser.add_argument("--case", type=Path,
		default=Path("cases/idealized_cube_dual_tree.json"))
	parser.add_argument("--check-only", action="store_true")
	args = parser.parse_args()
	try:
		data = validate(json.loads(args.case.read_text(encoding="utf-8")))
		if args.check_only:
			print("dual-tree centerline/radius contract: PASS")
			return 0
		prefix = args.output_prefix.resolve()
		require(prefix.parent.is_dir(), "output parent directory must exist")
		require(all(not Path(str(prefix)+suffix).exists() for suffix in
			(".surface.msh", ".volume.msh", ".contract.json", ".json")),
			"output already exists; refusing overwrite")
		result = generate(data, prefix)
		print("dual-tree geometry: PASS tets={} minimum_scaled_jacobian={:.6g}".format(
			result["tetra_count"], result["minimum_scaled_jacobian"]))
		return 0
	except Exception as error:
		print(f"dual-tree geometry rejected: {error}", file=sys.stderr)
		return 2


if __name__ == "__main__":
	sys.exit(main())
