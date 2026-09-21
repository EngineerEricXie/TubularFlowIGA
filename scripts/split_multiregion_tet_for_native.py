#!/usr/bin/env python3
"""Split an audited conforming multi-region mesh into native FEM submeshes.

Each interface receives the same explicit numeric boundary label on both
submeshes. No inlet, outlet, pressure, material, or flow semantics are added.
"""

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import sys

from audit_multiregion_tet_mesh import audit, parse_msh, read_contract


def reject(message):
	raise ValueError(message)


def split(nodes, elements, regions, boundaries, interfaces):
	triangles = {}
	cells = defaultdict(list)
	for _, dimension, _, name, connectivity in elements:
		if dimension == 2:
			triangles[tuple(sorted(connectivity))] = name
		elif dimension == 3:
			cells[name].append(connectivity)
	base = max((item["boundary_id"] for item in boundaries.values()), default=0)
	interface_labels = {name: base+index+1
		for index, name in enumerate(sorted(interfaces))}
	if any(label > 2147483647 for label in interface_labels.values()):
		reject("generated interface boundary label exceeds native int range")
	outputs = {}
	for region in regions:
		if not cells[region]:
			reject(f"region {region!r} has no tetrahedra")
		uses = defaultdict(list)
		for cell in cells[region]:
			a, b, c, d = cell
			for face in ((a, c, b), (a, b, d), (a, d, c), (b, c, d)):
				uses[tuple(sorted(face))].append(face)
		labelled = defaultdict(list)
		for key, faces in uses.items():
			if len(faces) > 2:
				reject(f"region {region!r} has a nonmanifold tetrahedral face")
			if len(faces) == 2:
				continue
			name = triangles.get(key)
			if name is None:
				reject(f"region {region!r} boundary face lacks original label")
			if name in boundaries:
				if boundaries[name]["adjacent_region"] != region:
					reject("original exterior boundary belongs to wrong region")
				label = boundaries[name]["boundary_id"]
			elif name in interfaces:
				if region not in interfaces[name]:
					reject("original interface belongs to wrong region")
				label = interface_labels[name]
			else:
				reject("original face label is neither boundary nor interface")
			labelled[label].append(faces[0])
		used_nodes = {tag for cell in cells[region] for tag in cell}
		outputs[region] = ({tag: nodes[tag] for tag in sorted(used_nodes)},
			cells[region], labelled)
	return outputs, interface_labels


def write_native_mesh(path, nodes, cells, labelled):
	try:
		import gmsh
	except ImportError as error:
		reject(f"Gmsh Python API unavailable: {error}")
	gmsh.initialize()
	try:
		gmsh.model.add("native_single_region_geometry_only")
		gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
		gmsh.option.setNumber("Mesh.Binary", 0)
		gmsh.model.addDiscreteEntity(3, 1)
		gmsh.model.mesh.addNodes(3, 1, list(nodes),
			[value for point in nodes.values() for value in point])
		gmsh.model.mesh.addElementsByType(1, 4, [],
			[node for cell in cells for node in cell])
		gmsh.model.addPhysicalGroup(3, [1], 1)
		gmsh.model.setPhysicalName(3, 1, "fluid")
		for entity, (label, faces) in enumerate(sorted(labelled.items()), 1):
			gmsh.model.addDiscreteEntity(2, entity)
			gmsh.model.mesh.addElementsByType(entity, 2, [],
				[node for face in faces for node in face])
			gmsh.model.addPhysicalGroup(2, [entity], label)
			gmsh.model.setPhysicalName(2, label, f"boundary_label_{label}")
		gmsh.write(str(path))
	finally:
		gmsh.finalize()


def validate_submesh(path, expected_region, expected_cells, expected_labels):
	nodes, elements = parse_msh(path)
	region_map = {"fluid": {"physical_name": "fluid", "role": "fluid",
		"require_connected": False}}
	boundary_map = {f"boundary_label_{label}": {
		"physical_name": f"boundary_label_{label}", "boundary_id": label,
		"adjacent_region": "fluid", "semantic_role": "unclassified_geometry_only"}
		for label in expected_labels}
	result = audit(nodes, elements, region_map, boundary_map, {}, 1.0e-3)
	if result["tetrahedra"] != expected_cells or \
			result["boundaries"] != {f"boundary_label_{label}": len(faces)
				for label, faces in sorted(expected_labels.items())}:
		reject(f"submesh audit changed cells or labelled boundary faces for {expected_region!r}")
	return result


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("multiregion_msh", type=Path)
	parser.add_argument("contract_json", type=Path)
	parser.add_argument("output_directory", type=Path)
	parser.add_argument("--max-tetrahedra", type=int, default=500000)
	args = parser.parse_args()
	try:
		contract, regions, boundaries, interfaces = read_contract(args.contract_json)
		nodes, elements = parse_msh(args.multiregion_msh)
		count = sum(dimension == 3 for _, dimension, _, _, _ in elements)
		if args.max_tetrahedra < 1 or count > args.max_tetrahedra:
			reject(f"input has {count} tetrahedra; allowed 1..{args.max_tetrahedra}")
		original_audit = audit(nodes, elements, regions, boundaries, interfaces, 1.0e-3)
		outputs, interface_labels = split(nodes, elements, regions, boundaries, interfaces)
		folder = args.output_directory.resolve()
		if folder == args.multiregion_msh.resolve().parent or folder == args.contract_json.resolve().parent:
			reject("output directory must be distinct from source directory")
		if folder.exists() and any(folder.iterdir()):
			reject("output directory already contains files")
		folder.mkdir(parents=True, exist_ok=True)
		result = {}
		for region, (sub_nodes, cells, labelled) in outputs.items():
			if not region.replace("_", "").isalnum():
				reject("region name is unsuitable for a safe output filename")
			path = folder/f"{region}.msh"
			write_native_mesh(path, sub_nodes, cells, labelled)
			audit_result = validate_submesh(path, region, len(cells), labelled)
			result[region] = {"file": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
				"audit": audit_result}
		for name, label in interface_labels.items():
			pair = interfaces[name]
			face_sets = []
			for region in pair:
				sub_nodes, _, labelled = outputs[region]
				face_sets.append({tuple(sorted(sub_nodes[tag] for tag in face))
					for face in labelled[label]})
			if face_sets[0] != face_sets[1] or not face_sets[0]:
				reject(f"split interface {name!r} lacks exact matching facets")
		manifest = {"schema_version": 1, "kind": "native_single_region_split_geometry_only",
			"source_mesh_sha256": hashlib.sha256(args.multiregion_msh.read_bytes()).hexdigest(),
			"source_contract_sha256": hashlib.sha256(args.contract_json.read_bytes()).hexdigest(),
			"original_audit": original_audit, "interface_boundary_labels": interface_labels,
			"submeshes": result,
			"limitations": "No vessel inlet/outlet, pressure, tissue material, flow, or perfusion inferred"}
		(folder/"manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True)+"\n",
			encoding="utf-8")
	except (AttributeError, IndexError, KeyError, OSError, ValueError) as error:
		print(f"native multi-region split rejected: {error}", file=sys.stderr)
		return 2
	print(f"native multi-region split: PASS {len(outputs)} regions, "
		f"{len(interface_labels)} exact matching interface labels")
	return 0


if __name__ == "__main__":
	sys.exit(main())
