#!/usr/bin/env python3
"""Audit first-order multi-region tetrahedra and their boundary/interface facets."""

import argparse
from collections import defaultdict, deque
import hashlib
import json
import math
from pathlib import Path
import sys

from surface_to_fem_volume import tetra_determinant, tetra_minimum_scaled_jacobian


def reject(message):
	raise ValueError(message)


def tokens(line, context):
	value = line.split()
	if not value:
		reject(f"empty Gmsh {context}")
	return value


def parse_msh(path):
	try:
		lines = path.read_text(encoding="ascii").splitlines()
	except (OSError, UnicodeError) as error:
		reject(f"cannot read ASCII Gmsh mesh: {error}")
	index = 0
	physical_names = {}
	entity_physical = {}
	nodes = {}
	elements = []
	mesh_format = False
	while index < len(lines):
		section = lines[index].strip(); index += 1
		if not section:
			continue
		if section == "$MeshFormat":
			if index+1 >= len(lines) or lines[index].strip() != "4.1 0 8" \
					or lines[index+1].strip() != "$EndMeshFormat":
				reject("requires ASCII Gmsh 4.1 with 8-byte reals")
			mesh_format = True; index += 2
		elif section == "$PhysicalNames":
			count = int(lines[index]); index += 1
			for _ in range(count):
				row = lines[index].strip(); index += 1
				parts = row.split(maxsplit=2)
				if len(parts) != 3 or len(parts[2]) < 2 or parts[2][0] != '"' or parts[2][-1] != '"':
					reject("invalid quoted physical name")
				key = (int(parts[0]), int(parts[1]))
				if key in physical_names:
					reject("duplicate physical name tag")
				physical_names[key] = parts[2][1:-1]
			if index >= len(lines) or lines[index].strip() != "$EndPhysicalNames":
				reject("invalid physical-name terminator")
			index += 1
		elif section == "$Entities":
			counts = list(map(int, tokens(lines[index], "entity counts"))); index += 1
			if len(counts) != 4:
				reject("invalid Gmsh entity counts")
			for dimension, count in enumerate(counts):
				for _ in range(count):
					row = tokens(lines[index], "entity"); index += 1
					physical_index = 4 if dimension == 0 else 7
					if len(row) <= physical_index:
						reject("truncated Gmsh entity")
					physical_count = int(row[physical_index])
					if dimension >= 2:
						if physical_count != 1 or len(row) <= physical_index+1:
							reject("every surface and volume entity must have exactly one physical tag")
						entity_physical[(dimension, int(row[0]))] = int(row[physical_index+1])
			if index >= len(lines) or lines[index].strip() != "$EndEntities":
				reject("invalid entity terminator")
			index += 1
		elif section == "$Nodes":
			header = list(map(int, tokens(lines[index], "node header"))); index += 1
			if len(header) != 4:
				reject("invalid node header")
			for _ in range(header[0]):
				block = list(map(int, tokens(lines[index], "node block"))); index += 1
				if len(block) != 4:
					reject("invalid node block")
				dimension, _, parametric, count = block
				tags = []
				for _ in range(count):
					tags.append(int(lines[index])); index += 1
				for tag in tags:
					row = tokens(lines[index], "node coordinates"); index += 1
					if len(row) != 3+(dimension if parametric else 0):
						reject("invalid node coordinate record")
					point = tuple(float(value) for value in row[:3])
					if tag in nodes or not all(math.isfinite(value) for value in point):
						reject("duplicate node tag or nonfinite coordinate")
					nodes[tag] = point
			if len(nodes) != header[1] or index >= len(lines) or lines[index].strip() != "$EndNodes":
				reject("node count mismatch or invalid terminator")
			index += 1
		elif section == "$Elements":
			header = list(map(int, tokens(lines[index], "element header"))); index += 1
			if len(header) != 4:
				reject("invalid element header")
			seen = 0
			for _ in range(header[0]):
				block = list(map(int, tokens(lines[index], "element block"))); index += 1
				if len(block) != 4:
					reject("invalid element block")
				dimension, entity, element_type, count = block
				if element_type not in (2, 4):
					reject(f"unsupported Gmsh element type {element_type}; only P1 triangles/tetrahedra are accepted")
				expected_nodes = 3 if element_type == 2 else 4
				physical_tag = entity_physical.get((dimension, entity))
				if physical_tag is None:
					reject("element entity lacks one physical tag")
				physical_name = physical_names.get((dimension, physical_tag))
				if physical_name is None:
					reject("element physical tag lacks a physical name")
				for _ in range(count):
					row = list(map(int, tokens(lines[index], "element"))); index += 1
					if len(row) != expected_nodes+1:
						reject("invalid element connectivity")
					elements.append((row[0], dimension, element_type, physical_name, tuple(row[1:])))
					seen += 1
			if seen != header[1] or index >= len(lines) or lines[index].strip() != "$EndElements":
				reject("element count mismatch or invalid terminator")
			index += 1
		elif section.startswith("$"):
			reject(f"unsupported Gmsh section {section}")
		else:
			reject("unexpected text outside a Gmsh section")
	if not mesh_format or not nodes or not elements:
		reject("incomplete Gmsh mesh")
	return nodes, elements


def read_contract(path):
	try:
		data = json.loads(path.read_text(encoding="utf-8"))
	except (OSError, json.JSONDecodeError) as error:
		reject(f"cannot read interface contract: {error}")
	if data.get("schema_version") != 1 or data.get("length_unit") != "m":
		reject("contract requires schema_version 1 and metre coordinates")
	regions = data.get("regions")
	boundaries = data.get("boundaries")
	interfaces = data.get("interfaces")
	if not all(isinstance(value, list) for value in (regions, boundaries, interfaces)) or not regions:
		reject("contract regions/boundaries/interfaces must be arrays and regions cannot be empty")
	region_map = {}
	for region in regions:
		if set(region) != {"physical_name", "role", "require_connected"}:
			reject("region entries require physical_name, role, and require_connected")
		name = region["physical_name"]
		if not isinstance(name, str) or not name or name in region_map \
				or region["role"] not in {"fluid", "solid", "porous"} \
				or not isinstance(region["require_connected"], bool):
			reject("invalid or duplicate region entry")
		region_map[name] = region
	boundary_map = {}
	boundary_ids = set()
	for boundary in boundaries:
		if set(boundary) != {"physical_name", "boundary_id", "adjacent_region", "semantic_role"}:
			reject("boundary entries have missing or unknown keys")
		name, region = boundary["physical_name"], boundary["adjacent_region"]
		identifier = boundary["boundary_id"]
		if not isinstance(name, str) or not name or name in boundary_map or region not in region_map \
				or not isinstance(identifier, int) or identifier < 0 or identifier in boundary_ids \
				or not isinstance(boundary["semantic_role"], str) or not boundary["semantic_role"]:
			reject("invalid or duplicate boundary entry")
		boundary_ids.add(identifier); boundary_map[name] = boundary
	interface_map = {}
	for interface in interfaces:
		if set(interface) != {"physical_name", "regions"}:
			reject("interface entries require physical_name and regions")
		name, pair = interface["physical_name"], interface["regions"]
		if not isinstance(name, str) or not name or name in interface_map \
				or not isinstance(pair, list) or len(pair) != 2 or pair[0] == pair[1] \
				or any(region not in region_map for region in pair):
			reject("invalid or duplicate interface entry")
		interface_map[name] = frozenset(pair)
	if set(boundary_map) & set(interface_map):
		reject("a physical surface name cannot be both boundary and interface")
	return data, region_map, boundary_map, interface_map


def audit(nodes, elements, region_map, boundary_map, interface_map, minimum_quality):
	triangles = {}
	physical_triangle_counts = defaultdict(int)
	tetrahedra = []
	region_cells = defaultdict(list)
	minimum_determinant = math.inf
	minimum_scaled = math.inf
	region_volume = defaultdict(float)
	seen_tetrahedra = set()
	for element_id, dimension, element_type, physical_name, connectivity in elements:
		if any(tag not in nodes for tag in connectivity) or len(set(connectivity)) != len(connectivity):
			reject(f"element {element_id} has repeated or missing nodes")
		if element_type == 2:
			if physical_name not in boundary_map and physical_name not in interface_map:
				reject(f"triangle physical group {physical_name!r} is absent from the contract")
			face = tuple(sorted(connectivity))
			if face in triangles:
				reject("duplicate labelled triangle")
			triangles[face] = physical_name
			physical_triangle_counts[physical_name] += 1
		elif element_type == 4:
			if physical_name not in region_map:
				reject(f"tetrahedron physical region {physical_name!r} is absent from the contract")
			key = tuple(sorted(connectivity))
			if key in seen_tetrahedra:
				reject("duplicate tetrahedron")
			seen_tetrahedra.add(key)
			points = [nodes[tag] for tag in connectivity]
			determinant = tetra_determinant(*points)
			quality = tetra_minimum_scaled_jacobian(*points)
			if not math.isfinite(determinant) or determinant <= 0.0:
				reject(f"tetrahedron {element_id} has nonpositive determinant {determinant}")
			if not math.isfinite(quality) or quality < minimum_quality:
				reject(f"tetrahedron {element_id} scaled Jacobian {quality} is below {minimum_quality}")
			minimum_determinant = min(minimum_determinant, determinant)
			minimum_scaled = min(minimum_scaled, quality)
			region_volume[physical_name] += determinant/6.0
			cell_index = len(tetrahedra)
			tetrahedra.append((element_id, physical_name, connectivity))
			region_cells[physical_name].append(cell_index)
	if not tetrahedra or not triangles:
		reject("mesh requires tetrahedra and labelled surface triangles")
	if set(region_cells) != set(region_map):
		reject("every configured region must contain at least one tetrahedron")
	face_uses = defaultdict(list)
	for cell_index, (_, region, cell) in enumerate(tetrahedra):
		for local in ((0,1,2), (0,1,3), (0,2,3), (1,2,3)):
			face_uses[tuple(sorted(cell[index] for index in local))].append((cell_index, region))
	region_adjacency = defaultdict(lambda: defaultdict(set))
	exterior_counts = defaultdict(int)
	interface_counts = defaultdict(int)
	for face, uses in face_uses.items():
		if len(uses) > 2:
			reject("nonmanifold tetrahedral face is used more than twice")
		label = triangles.get(face)
		if len(uses) == 1:
			region = uses[0][1]
			if label is None or label not in boundary_map:
				reject("exterior tetrahedral face lacks exactly one boundary classification")
			if boundary_map[label]["adjacent_region"] != region:
				reject(f"boundary {label!r} is attached to the wrong region")
			exterior_counts[label] += 1
		else:
			left_index, left_region = uses[0]
			right_index, right_region = uses[1]
			if left_region == right_region:
				if label is not None:
					reject("same-region interior face must not carry a boundary/interface triangle")
				region_adjacency[left_region][left_index].add(right_index)
				region_adjacency[left_region][right_index].add(left_index)
			else:
				if label is None or label not in interface_map:
					reject("cross-region face lacks exactly one declared interface triangle")
				if interface_map[label] != frozenset((left_region, right_region)):
					reject(f"interface {label!r} joins the wrong region pair")
				interface_counts[label] += 1
	for face, label in triangles.items():
		if face not in face_uses:
			reject(f"labelled surface triangle in {label!r} is not a tetrahedral facet")
	if set(exterior_counts) != set(boundary_map) or set(interface_counts) != set(interface_map):
		reject("every configured boundary and interface must contain at least one matching facet")
	for region, spec in region_map.items():
		if not spec["require_connected"]:
			continue
		cells = region_cells[region]
		visited = {cells[0]}; queue = deque([cells[0]])
		while queue:
			current = queue.popleft()
			for adjacent in region_adjacency[region][current]:
				if adjacent not in visited:
					visited.add(adjacent); queue.append(adjacent)
		if len(visited) != len(cells):
			reject(f"region {region!r} is disconnected")
	return {
		"nodes": len(nodes), "tetrahedra": len(tetrahedra),
		"labelled_triangles": len(triangles),
		"minimum_determinant_m3": minimum_determinant,
		"minimum_scaled_jacobian": minimum_scaled,
		"regions": {name: {"tetrahedra": len(region_cells[name]),
			"volume_m3": region_volume[name]} for name in sorted(region_map)},
		"boundaries": {name: exterior_counts[name] for name in sorted(boundary_map)},
		"interfaces": {name: interface_counts[name] for name in sorted(interface_map)},
		"exterior_and_interface_partition_exact": True,
		"region_connectivity_checked": True,
	}


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("mesh")
	parser.add_argument("contract")
	parser.add_argument("--output")
	parser.add_argument("--minimum-scaled-jacobian", type=float, default=1.0e-3)
	args = parser.parse_args()
	if not math.isfinite(args.minimum_scaled_jacobian) or not 0.0 < args.minimum_scaled_jacobian <= 1.0:
		reject("minimum scaled Jacobian must be in (0,1]")
	mesh_path = Path(args.mesh).resolve()
	contract_path = Path(args.contract).resolve()
	contract, regions, boundaries, interfaces = read_contract(contract_path)
	nodes, elements = parse_msh(mesh_path)
	result = audit(nodes, elements, regions, boundaries, interfaces,
		args.minimum_scaled_jacobian)
	manifest = {
		"schema_version": 1,
		"kind": "multiregion_tetrahedral_mesh_audit",
		"mesh_sha256": hashlib.sha256(mesh_path.read_bytes()).hexdigest(),
		"contract_sha256": hashlib.sha256(contract_path.read_bytes()).hexdigest(),
		"length_unit": "m",
		"minimum_required_scaled_jacobian": args.minimum_scaled_jacobian,
		"audit": result,
		"gates": {"positive_jacobians": True, "quality": True,
			"boundary_partition": True, "interface_partition": True,
			"region_membership": True, "passed": True},
	}
	if args.output:
		output = Path(args.output).resolve()
		output.parent.mkdir(parents=True, exist_ok=True)
		output.write_text(json.dumps(manifest, indent=2, sort_keys=True)+"\n", encoding="utf-8")
	print("multiregion_tet_mesh_audit: PASS "
		f"regions={len(regions)} tetrahedra={result['tetrahedra']} "
		f"interfaces={sum(result['interfaces'].values())} minSJ={result['minimum_scaled_jacobian']:.8g}")
	return 0


if __name__ == "__main__":
	try:
		raise SystemExit(main())
	except (ValueError, RuntimeError) as error:
		print(f"multiregion_tet_mesh_audit: ERROR: {error}", file=sys.stderr)
		raise SystemExit(2)
