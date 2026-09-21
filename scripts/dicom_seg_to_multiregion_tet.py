#!/usr/bin/env python3
"""Make a conforming, voxel-faithful tetra ROI from same-grid DICOM SEG labels.

This is a geometry-only route. It does not classify vessel ports, prescribe
flow, smooth anatomy, or assign Darcy material parameters.
"""

import argparse
from collections import defaultdict
import hashlib
from itertools import permutations
import json
import math
from pathlib import Path
import sys

from audit_dicom_seg_region_overlap import audit, face_contact


def reject(message):
	raise ValueError(message)


def parse_roi(value, shape):
	if value == "full":
		return tuple((0, size) for size in shape)
	parts = value.split(",")
	if len(parts) != 3:
		reject("ROI must be z0:z1,y0:y1,x0:x1 or full")
	result = []
	for part, size in zip(parts, shape):
		bounds = part.split(":")
		if len(bounds) != 2:
			reject("ROI must have three half-open index ranges")
		try:
			start, stop = map(int, bounds)
		except ValueError:
			reject("ROI bounds must be integers")
		if not 0 <= start < stop <= size:
			reject("ROI bounds exceed the SEG grid or are empty")
		result.append((start, stop))
	return tuple(result)


def require_manifold_region_boundaries(face_owners):
	"""Reject boundary edge/vertex pinches in each region's voxel solid."""
	region_faces = defaultdict(list)
	for owners in face_owners.values():
		if len(owners) == 1 or owners[0][0] != owners[1][0]:
			for region, face in owners:
				region_faces[region].append(face)
	for region, faces in region_faces.items():
		edges = defaultdict(list)
		vertex_faces = defaultdict(set)
		for face_index, face in enumerate(faces):
			for vertex in face:
				vertex_faces[vertex].add(face_index)
			for left, right in ((face[0], face[1]), (face[1], face[2]),
					(face[2], face[0])):
				edges[tuple(sorted((left, right)))].append(face_index)
		for edge, incident in edges.items():
			if len(incident) != 2:
				reject(f"region {region!r} has a nonmanifold boundary edge")
		fan = defaultdict(lambda: defaultdict(set))
		for (left, right), incident in edges.items():
			for vertex in (left, right):
				fan[vertex][incident[0]].add(incident[1])
				fan[vertex][incident[1]].add(incident[0])
		for vertex, incident in vertex_faces.items():
			visited = {next(iter(incident))}
			pending = list(visited)
			while pending:
				for neighbour in fan[vertex][pending.pop()]:
					if neighbour not in visited:
						visited.add(neighbour)
						pending.append(neighbour)
			if visited != incident:
				reject(f"region {region!r} has a nonmanifold boundary vertex")


def build_mesh(mask_by_region, roi, origin_center_mm, spacing_xyz_mm, max_voxels,
		require_connected=True):
	"""Return shared lattice nodes, positive tets, and exact face owners."""
	import numpy as np
	if len(mask_by_region) < 2:
		reject("at least two region masks are needed")
	shape = next(iter(mask_by_region.values())).shape
	if any(mask.shape != shape or mask.dtype != bool for mask in mask_by_region.values()):
		reject("all region masks must be boolean arrays on one 3D grid")
	occupancy = np.zeros(shape, dtype=np.uint8)
	for index, mask in enumerate(mask_by_region.values(), 1):
		if index > 255 or np.any(occupancy & mask):
			reject("region masks overlap or exceed 255 regions")
		occupancy[mask] = index
	if max_voxels < 1:
		reject("max_voxels must be positive")
	slices = tuple(slice(start, stop) for start, stop in roi)
	selected = occupancy[slices]
	voxel_count = int(np.count_nonzero(selected))
	if voxel_count < 1 or voxel_count > max_voxels:
		reject(f"ROI contains {voxel_count} labelled voxels; allowed 1..{max_voxels}")
	if require_connected:
		from scipy import ndimage
		for index, region in enumerate(mask_by_region, 1):
			if not np.any(selected == index):
				continue
			_, components = ndimage.label(selected == index)
			if components != 1:
				reject(f"region {region!r} has {components} disconnected 6-neighbour ROI components")
	if not all(math.isfinite(value) and value > 0 for value in spacing_xyz_mm):
		reject("invalid voxel spacing")
	if not all(map(math.isfinite, origin_center_mm)):
		reject("invalid voxel origin")
	region_names = tuple(mask_by_region)
	nodes = {}
	vertex_tags = {}
	tets = defaultdict(list)
	face_owners = defaultdict(list)
	perms = tuple(permutations(range(3)))
	for local in np.argwhere(selected):
		k, j, i = (int(local[axis]) + roi[axis][0] for axis in range(3))
		region = region_names[int(occupancy[k, j, i])-1]
		def vertex(delta):
			key = (k+delta[0], j+delta[1], i+delta[2])
			if key not in vertex_tags:
				tag = len(vertex_tags)+1
				vertex_tags[key] = tag
				dx, dy, dz = spacing_xyz_mm
				cx, cy, cz = origin_center_mm
				nodes[tag] = ((cx+(key[2]-.5)*dx)*.001,
					(cy+(key[1]-.5)*dy)*.001,
					(cz+(key[0]-.5)*dz)*.001)
			return vertex_tags[key]
		for perm in perms:
			path = [[0, 0, 0]]
			for axis in perm:
				step = path[-1].copy(); step[axis] = 1; path.append(step)
			cell = [vertex(delta) for delta in path]
			p0, p1, p2, p3 = (nodes[tag] for tag in cell)
			ab = [p1[axis]-p0[axis] for axis in range(3)]
			ac = [p2[axis]-p0[axis] for axis in range(3)]
			ad = [p3[axis]-p0[axis] for axis in range(3)]
			det = (ab[0]*(ac[1]*ad[2]-ac[2]*ad[1])
				-ab[1]*(ac[0]*ad[2]-ac[2]*ad[0])
				+ab[2]*(ac[0]*ad[1]-ac[1]*ad[0]))
			if det < 0:
				cell[1], cell[2] = cell[2], cell[1]
			elif det == 0:
				reject("degenerate voxel tetrahedron")
			tets[region].append(tuple(cell))
			a, b, c, d = cell
			for face in ((a, c, b), (a, b, d), (a, d, c), (b, c, d)):
				face_owners[tuple(sorted(face))].append((region, face))
	tag_grid = {tag: grid for grid, tag in vertex_tags.items()}
	interfaces = defaultdict(list)
	boundaries = defaultdict(list)
	for key, owners in face_owners.items():
		if len(owners) > 2:
			reject("nonmanifold tetrahedron face")
		if len(owners) == 2:
			left, right = owners
			if left[0] != right[0]:
				pair = tuple(sorted((left[0], right[0])))
				if "tissue_candidate" not in pair:
					reject("direct vessel-to-vessel contact needs an explicit contract")
				vessel_side = left if left[0] != "tissue_candidate" else right
				interfaces[pair].append(vessel_side[1])
		else:
			region, face = owners[0]
			grid_vertices = [tag_grid[node] for node in face]
			is_cut = any(all(vertex[axis] == bound for vertex in grid_vertices)
				for axis, limits in enumerate(roi) for bound in limits)
			boundaries[(region, "roi_cut" if is_cut else "unclassified_exterior")].append(face)
	if not interfaces:
		reject("ROI contains no tissue-vessel tetrahedral interface")
	require_manifold_region_boundaries(face_owners)
	return nodes, tets, boundaries, interfaces, voxel_count


def write_gmsh(path, nodes, tets, boundaries, interfaces):
	try:
		import gmsh
	except ImportError as error:
		reject(f"Gmsh Python API unavailable: {error}")
	gmsh.initialize()
	try:
		gmsh.model.add("seg_multiregion_geometry_only")
		gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
		gmsh.option.setNumber("Mesh.Binary", 0)
		regions = list(tets)
		for index, region in enumerate(regions, 1):
			gmsh.model.addDiscreteEntity(3, index)
			if index == 1:
				gmsh.model.mesh.addNodes(3, index, list(nodes),
					[value for point in nodes.values() for value in point])
			gmsh.model.mesh.addElementsByType(index, 4, [],
				[node for cell in tets[region] for node in cell])
			gmsh.model.addPhysicalGroup(3, [index], 100+index)
			gmsh.model.setPhysicalName(3, 100+index, region)
		surfaces = {}
		for name, faces in [*((f"{region}_{kind}", faces)
			for (region, kind), faces in sorted(boundaries.items())),
			*(("interface_"+"_".join(pair), faces)
			for pair, faces in sorted(interfaces.items()))]:
			entity = len(surfaces)+1
			gmsh.model.addDiscreteEntity(2, entity)
			gmsh.model.mesh.addElementsByType(entity, 2, [],
				[node for face in faces for node in face])
			gmsh.model.addPhysicalGroup(2, [entity], entity)
			gmsh.model.setPhysicalName(2, entity, name)
			surfaces[name] = entity
		gmsh.write(str(path))
	finally:
		gmsh.finalize()
	return surfaces


def make_contract(tets, boundaries, interfaces, surfaces):
	return {"schema_version": 1, "length_unit": "m",
		"regions": [{"physical_name": name,
			"role": "porous" if name == "tissue_candidate" else "fluid",
			"require_connected": True} for name in tets],
		"boundaries": [{"physical_name": f"{region}_{kind}",
			"boundary_id": surfaces[f"{region}_{kind}"],
			"adjacent_region": region,
			"semantic_role": kind+"_not_a_flow_port"}
			for region, kind in boundaries],
		"interfaces": [{"physical_name": "interface_"+"_".join(pair),
			"regions": list(pair)} for pair in interfaces]}


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("seg_dicom", type=Path)
	parser.add_argument("tissue_segment_number", type=int)
	parser.add_argument("vessel_segment_numbers", type=int, nargs="+")
	parser.add_argument("output_msh", type=Path)
	parser.add_argument("--roi", required=True, help="half-open z0:z1,y0:y1,x0:x1 or full")
	parser.add_argument("--max-voxels", type=int, default=50000)
	args = parser.parse_args()
	try:
		result, masks, origin = audit(args.seg_dicom, args.tissue_segment_number,
			args.vessel_segment_numbers, return_masks=True)
		if not result["vessel_masks_mutually_exclusive"]:
			reject("vessel labels overlap; cannot create an unambiguous partition")
		roi = parse_roi(args.roi, masks[args.tissue_segment_number].shape)
		import numpy as np
		vessel_union = np.zeros_like(masks[args.tissue_segment_number])
		for number in args.vessel_segment_numbers:
			vessel_union |= masks[number]
		regions = {"tissue_candidate": masks[args.tissue_segment_number] & ~vessel_union}
		for vessel in result["vessels"]:
			regions[f"seg_{vessel['segment_number']}_{vessel['segment_label'].lower()}"] = \
				masks[vessel["segment_number"]]
		nodes, tets, boundaries, interfaces, count = build_mesh(regions, roi,
			origin, result["voxel_spacing_mm"], args.max_voxels)
		output = args.output_msh.resolve()
		if output == args.seg_dicom.resolve():
			reject("mesh output must not overwrite source SEG")
		contract_path = output.with_suffix(".contract.json")
		manifest_path = output.with_suffix(".manifest.json")
		if any(path.exists() for path in (output, contract_path, manifest_path)):
			reject("mesh, contract, or manifest output already exists")
		output.parent.mkdir(parents=True, exist_ok=True)
		surfaces = write_gmsh(output, nodes, tets, boundaries, interfaces)
		contract = make_contract(tets, boundaries, interfaces, surfaces)
		contract_path.write_text(json.dumps(contract, indent=2, sort_keys=True)+"\n",
			encoding="utf-8")
		from audit_multiregion_tet_mesh import audit as audit_mesh, parse_msh, read_contract
		_, region_map, boundary_map, interface_map = read_contract(contract_path)
		audit_result = audit_mesh(*parse_msh(output), region_map, boundary_map,
			interface_map, 1.0e-3)
		voxel_volume_m3 = math.prod(result["voxel_spacing_mm"])*1.0e-9
		for region, cells in tets.items():
			if len(cells) % 6:
				reject("region tetra count is not six per source voxel")
			expected_volume = len(cells)/6*voxel_volume_m3
			actual_volume = audit_result["regions"][region]["volume_m3"]
			if abs(actual_volume-expected_volume) > 1e-10*expected_volume:
				reject(f"region {region!r} tetra volume differs from source voxels")
		roi_slices = tuple(slice(start, stop) for start, stop in roi)
		for pair, faces in interfaces.items():
			vessel = next(name for name in pair if name != "tissue_candidate")
			contact = face_contact(regions["tissue_candidate"][roi_slices],
				regions[vessel][roi_slices], result["voxel_spacing_mm"])
			if len(faces) != 2*contact["total_contact_faces"]:
				reject("tetra interface differs from source voxel-face contact")
		manifest = {"schema_version": 1, "kind": "voxel_faithful_multiregion_roi_geometry_only",
			"source_sha256": result["source_sha256"], "mesh_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
			"roi_zyx_half_open": roi, "labelled_voxels": count,
			"nodes": len(nodes), "tetrahedra": sum(map(len, tets.values())),
			"interfaces": {"_".join(pair): len(faces) for pair, faces in interfaces.items()},
			"boundary_labels": {f"{region}_{kind}": len(faces)
				for (region, kind), faces in boundaries.items()},
			"geometry_audit": audit_result,
			"limitations": "ROI geometry only: no vessel ports, physiological BC, material, flow, or perfusion claim"}
		manifest_path.write_text(
			json.dumps(manifest, indent=2, sort_keys=True)+"\n", encoding="utf-8")
	except (AttributeError, IndexError, KeyError, OSError, ValueError) as error:
		print(f"DICOM SEG multiregion tetra rejected: {error}", file=sys.stderr)
		return 2
	print(f"DICOM SEG multiregion tetra: geometry-only ROI, {count} voxels, "
		f"{sum(map(len, tets.values()))} tetrahedra")
	return 0


if __name__ == "__main__":
	sys.exit(main())
