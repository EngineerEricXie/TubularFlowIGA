"""Shared template-free centerline, radius, and surface operations."""

import math


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


def tree_node_radii(tree):
	if "node_radii_m" in tree:
		return tree["node_radii_m"]
	return [max(radius for start, stop, radius in tree["segments"]
		if index in (start, stop)) for index in range(len(tree["nodes_m"]))]


def branch_transition_radius(junction_radius, segment_radius, distance, length):
	u = min(1., distance/length)
	weight = u*u*(3.-2.*u)
	return (1.-weight)*junction_radius+weight*segment_radius


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
		outgoing_sum = [sum(direction[axis] for direction in outgoing)
			for axis in range(3)]
		mean_outgoing = (None if not outgoing else
			(normalized(outgoing_sum) if dot(outgoing_sum, outgoing_sum) > 1.e-24
				else outgoing[0]))
		if incoming is None:
			tangent = mean_outgoing
		elif not outgoing:
			tangent = incoming
		else:
			combined = add(incoming, mean_outgoing)
			tangent = (normalized(combined) if dot(combined, combined) > 1.e-24
				else incoming)
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


def clean_surface_triangles(nodes, triangles, tolerance, shells=None):
	"""Collapse short Boolean-seam edges without changing shell topology."""
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
	if shells:
		for name, labels in shells.items():
			usage = {}
			for label in labels:
				for triangle in cleaned[label]:
					for first, second in ((triangle[0], triangle[1]),
						(triangle[1], triangle[2]), (triangle[2], triangle[0])):
						edge = tuple(sorted((first, second)))
						usage[edge] = usage.get(edge, 0)+1
			if not usage or any(count != 2 for count in usage.values()):
				raise ValueError(f"surface cleanup opened the {name} shell")
	return points, cleaned, len(nodes)-len(points), dropped


def surface_boundary_nodes(triangles):
	edge_use = {}
	for triangle in triangles:
		for first, second in ((triangle[0], triangle[1]),
			(triangle[1], triangle[2]), (triangle[2], triangle[0])):
			edge = tuple(sorted((first, second)))
			edge_use[edge] = edge_use.get(edge, 0)+1
	return {tag for edge, count in edge_use.items() if count == 1 for tag in edge}


def closest_centerline_segment(point, centerline):
	closest = None
	for first, last, arc, length, total, radius, junction_radius, start in centerline:
		direction = sub(last, first)
		position = min(1., max(0., dot(sub(point, first), direction)/
			dot(direction, direction)))
		axis_point = tuple(first[axis]+position*direction[axis]
			for axis in range(3))
		distance = math.dist(point, axis_point)
		candidate = (distance, axis_point, arc+position*length, total,
			radius, junction_radius, start)
		if closest is None or distance < closest[0]:
			closest = candidate
	return closest


def sampled_centerline(tree, geometry):
	paths = spline_capsule_paths(tree, geometry.get("subsegments_per_edge", 2),
		geometry.get("tangent_scale", .25))
	node_radii = tree_node_radii(tree)
	centerline = []
	for start, _, radius, path in paths:
		lengths = [math.dist(path[index], path[index+1])
			for index in range(len(path)-1)]
		total = sum(lengths)
		arc = 0.
		for index, length in enumerate(lengths):
			centerline.append((path[index], path[index+1], arc, length, total,
				radius, node_radii[start], start))
			arc += length
	return centerline


def blend_branch_radii(points, triangles, tree, surface, radius_add, geometry):
	"""Expand child surfaces into a smooth junction-to-segment radius taper."""
	centerline = sampled_centerline(tree, geometry)
	boundary = surface_boundary_nodes(triangles[surface])
	tags = {tag for triangle in triangles[surface] for tag in triangle}
	fraction = geometry["radius_transition_fraction"]
	moved = 0
	maximum_displacement = 0.
	for tag in tags-boundary:
		point = points[tag]
		distance, axis_point, arc, total, radius, junction_radius, _ = \
			closest_centerline_segment(point, centerline)
		transition_length = fraction*total
		if junction_radius > radius and arc < transition_length:
			target = branch_transition_radius(junction_radius+radius_add,
				radius+radius_add, arc, transition_length)
			if target > distance:
				points[tag] = tuple(axis_point[axis]+
					(point[axis]-axis_point[axis])*target/distance for axis in range(3))
				moved += 1
				maximum_displacement = max(maximum_displacement, target-distance)
	return moved, maximum_displacement


def fit_tree_radii(points, triangles, tree, geometry, names=None):
	"""Expand an undersized capsule union to the skeleton's nodal radius field."""
	centerline = sampled_centerline(tree, geometry)
	children = [0]*len(tree["nodes_m"])
	for start, _, _ in tree["segments"]:
		children[start] += 1
	groups = triangles.values() if names is None else (triangles[name] for name in names)
	tags = {tag for group in groups for triangle in group for tag in triangle}
	moved = 0
	maximum_displacement = 0.
	for tag in tags:
		point = points[tag]
		distance, axis_point, arc, total, radius, start_radius, start = \
			closest_centerline_segment(point, centerline)
		transition_length = total
		if children[start] > 1:
			transition_length *= geometry["radius_transition_fraction"]
		target = branch_transition_radius(start_radius, radius, arc, transition_length)
		if distance > 0. and target > distance:
			points[tag] = tuple(axis_point[axis]+
				(point[axis]-axis_point[axis])*target/distance for axis in range(3))
			moved += 1
			maximum_displacement = max(maximum_displacement, target-distance)
	return moved, maximum_displacement


def fit_terminal_caps(points, triangles, tree):
	"""Scale inlet and outlet disks to their endpoint radii."""
	endpoints = [("inlet", tree["root"])]+[
		(f"outlet_{index}", node) for index, node in enumerate(tree["terminals"])]
	moved = 0
	maximum_displacement = 0.
	for name, node in endpoints:
		tags = {tag for triangle in triangles[name] for tag in triangle}
		center = tree["nodes_m"][node]
		current_radius = max(math.dist(points[tag], center) for tag in tags)
		factor = tree_node_radii(tree)[node]/current_radius
		for tag in tags:
			point = points[tag]
			updated = tuple(center[axis]+factor*(point[axis]-center[axis])
				for axis in range(3))
			displacement = math.dist(point, updated)
			if displacement:
				points[tag] = updated
				moved += 1
				maximum_displacement = max(maximum_displacement, displacement)
	return moved, maximum_displacement


def smooth_surface_groups(points, triangles, names, iterations, pass_band):
	try:
		import vtk
	except ImportError as error:
		raise ValueError(f"VTK Python module unavailable: {error}") from error
	for name in names:
		group = triangles[name]
		tags = sorted({tag for triangle in group for tag in triangle})
		local = {tag: index for index, tag in enumerate(tags)}
		boundary = surface_boundary_nodes(group)
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
