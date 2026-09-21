#!/usr/bin/env python3
"""Inventory vessel voxel faces outside tissue; never infer flow ports or BCs."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

try:
	from scripts.audit_dicom_seg_region_overlap import audit, reject
except ModuleNotFoundError:
	from audit_dicom_seg_region_overlap import audit, reject


def count_exposed_faces(components, tissue_only, vessel_union, spacing_xyz_mm):
	"""Per 6-connected vessel component: tissue, background and SEG-grid faces."""
	import numpy as np
	if components.ndim != 3 or tissue_only.shape != components.shape \
			or vessel_union.shape != components.shape \
			or not np.issubdtype(components.dtype, np.integer) \
			or np.any(components < 0) \
			or np.any((components > 0) & ~vessel_union) \
			or np.any(tissue_only & vessel_union):
		reject("exposed-face arrays must be disjoint, aligned 3D masks")
	dx, dy, dz = spacing_xyz_mm
	if not all(math.isfinite(value) and value > 0 for value in (dx, dy, dz)):
		reject("exposed-face spacing must be finite and positive")
	count = int(components.max())
	if set(np.unique(components)) != set(range(count+1)):
		reject("exposed-face component IDs must be contiguous from 1")
	background = ~tissue_only & ~vessel_union
	areas = (dx*dy, dx*dz, dy*dz)
	counts = {role: np.zeros((count+1, 3), dtype=np.int64)
		for role in ("tissue", "background", "other_vessel", "grid_extent")}
	for axis in range(3):
		lower = [slice(None)]*3
		upper = [slice(None)]*3
		lower[axis] = slice(None, -1)
		upper[axis] = slice(1, None)
		lo, hi = tuple(lower), tuple(upper)
		for role, neighbour in (("tissue", tissue_only),
				("background", background),
				("other_vessel", vessel_union & (components == 0))):
			counts[role][:, axis] += np.bincount(
				components[lo][(components[lo] > 0) & neighbour[hi]],
				minlength=count+1)
			counts[role][:, axis] += np.bincount(
				components[hi][(components[hi] > 0) & neighbour[lo]],
				minlength=count+1)
		first = [slice(None)]*3
		last = [slice(None)]*3
		first[axis] = 0
		last[axis] = -1
		for boundary in (tuple(first), tuple(last)):
			counts["grid_extent"][:, axis] += np.bincount(
				components[boundary][components[boundary] > 0],
				minlength=count+1)
	sizes = np.bincount(components.ravel(), minlength=count+1)
	internal = np.zeros(count+1, dtype=np.int64)
	for axis in range(3):
		lower = [slice(None)]*3
		upper = [slice(None)]*3
		lower[axis] = slice(None, -1)
		upper[axis] = slice(1, None)
		lo, hi = components[tuple(lower)], components[tuple(upper)]
		internal += np.bincount(lo[(lo > 0) & (lo == hi)], minlength=count+1)
	partition = sum(values.sum(axis=1) for values in counts.values())
	if not np.array_equal(partition[1:], (6*sizes-2*internal)[1:]):
		reject("vessel surface faces do not partition into audited categories")
	return [{"component_id": component_id,
		**{role: {"faces_normal_z_y_x": counts[role][component_id].tolist(),
			"total_faces": int(counts[role][component_id].sum()),
			"area_mm2": float(np.dot(counts[role][component_id], areas))}
			for role in counts}}
		for component_id in range(1, count+1)]


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("seg_dicom", type=Path)
	parser.add_argument("tissue_segment_number", type=int)
	parser.add_argument("vessel_segment_numbers", type=int, nargs="+")
	parser.add_argument("--manifest", type=Path, required=True)
	args = parser.parse_args()
	try:
		if args.manifest.resolve() == args.seg_dicom.resolve() \
				or args.manifest.exists():
			reject("manifest already exists or overlaps source SEG")
		base, masks, _ = audit(args.seg_dicom, args.tissue_segment_number,
			args.vessel_segment_numbers, return_masks=True)
		if not base["vessel_masks_mutually_exclusive"]:
			reject("overlapping vessel masks make exposed-face ownership ambiguous")
		import numpy as np
		from scipy import ndimage
		vessel_union = np.logical_or.reduce([masks[number]
			for number in args.vessel_segment_numbers])
		tissue_only = masks[args.tissue_segment_number] & ~vessel_union
		vessels = []
		for item in base["vessels"]:
			number = item["segment_number"]
			components, count = ndimage.label(masks[number])
			faces = count_exposed_faces(components, tissue_only, vessel_union,
				base["voxel_spacing_mm"])
			if count != item["components_6_neighbour"] or any(
				face["tissue"]["total_faces"] != recorded["tissue_contact_faces"]
				for face, recorded in zip(faces,
					item["component_tissue_contact"]["components"])):
				reject("exposed-face and component-contact audits disagree")
			vessels.append({"segment_number": number,
				"segment_label": item["segment_label"], "components": faces})
		result = {"schema_version": 1,
			"kind": "seg_vessel_face_inventory_not_ports_or_boundary_conditions",
			"source_seg_sha256": base["source_sha256"],
			"auditor_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
			"tissue_segment_number": args.tissue_segment_number,
			"vessel_segment_numbers": args.vessel_segment_numbers,
			"source_audit_sha256": hashlib.sha256(json.dumps(base,
				sort_keys=True).encode()).hexdigest(),
			"voxel_spacing_mm": base["voxel_spacing_mm"],
			"vessels": vessels,
			"limitations": "Background and grid-extent faces are not automatically ports; no flow direction, BC, material, or physiology inferred"}
		args.manifest.parent.mkdir(parents=True, exist_ok=True)
		with args.manifest.open("x", encoding="utf-8") as stream:
			stream.write(json.dumps(result, indent=2, sort_keys=True)+"\n")
	except (OSError, ValueError) as error:
		print(f"SEG exposed vessel faces rejected: {error}", file=sys.stderr)
		return 2
	print(f"SEG exposed vessel faces: PASS {args.manifest}")
	return 0


if __name__ == "__main__":
	sys.exit(main())
