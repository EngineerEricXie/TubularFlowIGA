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


def require(condition, message):
	if not condition:
		raise ValueError(message)


def dot(a, b):
	return sum(x*y for x, y in zip(a, b))


def sub(a, b):
	return [x-y for x, y in zip(a, b)]


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
	def cross(a, b):
		return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2],
			a[0]*b[1]-a[1]*b[0])
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
	require(math.isfinite(minimum_size) and 0 < minimum_size <=
		data["target_mesh_size_m"] and isinstance(curvature, int) and
		not isinstance(curvature, bool) and 8 <= curvature <= 128 and
		isinstance(data.get("mesh_size_extend_from_boundary", True), bool) and
		isinstance(data.get("arterial_outer_reverse", True), bool),
		"minimum size, curvature divisions or boundary extension is invalid")
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


def add_tree(gmsh, tree, radius_add, reverse=False):
	nodes = tree["nodes_m"]
	parts = []
	for i, j, radius in tree["segments"]:
		start, end = (nodes[j], nodes[i]) if reverse else (nodes[i], nodes[j])
		delta = sub(end, start)
		parts.append((3, gmsh.model.occ.addCylinder(*start, *delta, radius+radius_add)))
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
		cube = (3, gmsh.model.occ.addBox(*origin, *size))
		artery = add_tree(gmsh, data["arterial_tree"], 0., reverse=True)
		artery_outer = add_tree(gmsh, data["arterial_tree"],
			data["wall_thickness_m"],
			reverse=data.get("arterial_outer_reverse", True))
		vein = add_tree(gmsh, data["venous_tree"], 0.)
		vein_outer = add_tree(gmsh, data["venous_tree"],
			data["wall_thickness_m"])
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
			for name, region in regions.items():
				require(gmsh.model.occ.getMass(3, next(iter(region))) > 0,
					f"{name} has nonpositive volume")
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
				"regions": REGIONS, "surfaces": SURFACES,
				"surface_count": {name: len(surfaces)
					for name, surfaces in surface_groups.items()},
				"region_volume_m3": {name: gmsh.model.occ.getMass(3,
					next(iter(tags))) for name, tags in regions.items()},
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
