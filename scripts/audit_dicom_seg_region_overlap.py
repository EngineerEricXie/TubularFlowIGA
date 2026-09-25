#!/usr/bin/env python3
"""Audit same-SEG voxel overlap; never infer a vascular port or tissue interface."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys


def reject(message):
	raise ValueError(message)


def face_contact(tissue, vessel, spacing_xyz_mm):
	"""Count only shared voxel faces, never diagonal contact or inferred ports."""
	import numpy as np
	if tissue.shape != vessel.shape or tissue.ndim != 3:
		reject("face-contact masks must share one 3D grid")
	if np.any(tissue & vessel):
		reject("face-contact regions must be disjoint")
	dx, dy, dz = spacing_xyz_mm
	if not all(math.isfinite(value) and value > 0 for value in (dx, dy, dz)):
		reject("face-contact spacing must be positive and finite")
	counts = []
	for axis in range(3):
		lower = [slice(None)] * 3
		upper = [slice(None)] * 3
		lower[axis] = slice(None, -1)
		upper[axis] = slice(1, None)
		first = tuple(lower)
		second = tuple(upper)
		counts.append(int(np.count_nonzero(tissue[first] & vessel[second]))
			+ int(np.count_nonzero(vessel[first] & tissue[second])))
	areas = [dx * dy, dx * dz, dy * dz]
	return {"faces_normal_z_y_x": counts,
		"face_area_mm2_normal_z_y_x": areas,
		"total_contact_faces": sum(counts),
		"total_contact_area_mm2": sum(count * area
			for count, area in zip(counts, areas))}


def component_face_contacts(tissue, components, count):
	"""Return shared-face counts per 6-neighbour vessel component ID."""
	import numpy as np
	if tissue.shape != components.shape or tissue.ndim != 3 or count < 0:
		reject("component-contact arrays must share one 3D grid")
	if np.any(tissue & (components > 0)) or np.any(components < 0) \
			or np.any(components > count):
		reject("component-contact labels overlap tissue or are out of range")
	result = np.zeros(count+1, dtype=np.int64)
	for axis in range(3):
		lower = [slice(None)] * 3
		upper = [slice(None)] * 3
		lower[axis] = slice(None, -1)
		upper[axis] = slice(1, None)
		first, second = tuple(lower), tuple(upper)
		left = components[first]
		right = components[second]
		selected = tissue[first] & (right > 0)
		result += np.bincount(right[selected], minlength=count+1)
		selected = tissue[second] & (left > 0)
		result += np.bincount(left[selected], minlength=count+1)
	return result


def component_bounding_boxes_zyx(components, count):
	"""Return half-open voxel-index boxes; these do not classify flow ports."""
	import numpy as np
	from scipy import ndimage
	if components.ndim != 3 or count < 0 or np.any(components < 0) \
			or np.any(components > count):
		reject("component bounding-box labels are invalid")
	boxes = ndimage.find_objects(components, max_label=count)
	if len(boxes) != count or any(box is None for box in boxes):
		reject("component IDs are missing from bounding-box inventory")
	return [{"component_id": index+1,
		"bbox_zyx_half_open": [[axis.start, axis.stop] for axis in box],
		"touches_seg_grid_extent": any(axis.start == 0
			or axis.stop == components.shape[dimension]
			for dimension, axis in enumerate(box))}
		for index, box in enumerate(boxes)]


def parse_roi_zyx(value, shape):
	if value == "full":
		return [[0, extent] for extent in shape]
	parts = value.split(",")
	if len(parts) != 3:
		reject("ROI must be half-open z0:z1,y0:y1,x0:x1")
	try:
		roi = [[int(bound) for bound in part.split(":")] for part in parts]
	except ValueError:
		reject("ROI bounds must be integers")
	if any(len(bounds) != 2 or bounds[0] < 0 or bounds[0] >= bounds[1]
			or bounds[1] > shape[axis] for axis, bounds in enumerate(roi)):
		reject("ROI is empty or exceeds the SEG grid")
	return roi


def roi_component_membership(mask, roi):
	"""Count selected voxels by full-grid 6-neighbour component ID."""
	import numpy as np
	from scipy import ndimage
	if mask.ndim != 3 or mask.dtype != np.bool_:
		reject("ROI component mask must be a 3D boolean array")
	labels, count = ndimage.label(mask)
	slices = tuple(slice(*bounds) for bounds in roi)
	ids, counts = np.unique(labels[slices][labels[slices] > 0],
		return_counts=True)
	if any(component_id < 1 or component_id > count for component_id in ids):
		reject("ROI component ID is out of range")
	return [{"component_id": int(component_id), "roi_voxels": int(voxels)}
		for component_id, voxels in zip(ids, counts)]


def audit(path, tissue_number, vessel_numbers, return_masks=False):
	try:
		import numpy as np
		import pydicom
		from scipy import ndimage
	except ImportError as error:
		reject(f"DICOM SEG overlap dependencies unavailable: {error}")
	if tissue_number in vessel_numbers or len(set(vessel_numbers)) != len(vessel_numbers):
		reject("tissue and vessel segment numbers must be distinct")
	data = pydicom.dcmread(str(path))
	if data.Modality != "SEG" or data.SegmentationType != "BINARY":
		reject("expected a binary DICOM SEG")
	segments = {int(item.SegmentNumber): str(item.SegmentLabel)
		for item in data.SegmentSequence}
	selected = [tissue_number, *vessel_numbers]
	if any(number not in segments for number in selected):
		reject("requested segment number is absent")
	if len(data.ReferencedSeriesSequence) != 1:
		reject("expected exactly one referenced image series")
	shared = data.SharedFunctionalGroupsSequence[0]
	orientation = tuple(map(float,
		shared.PlaneOrientationSequence[0].ImageOrientationPatient))
	if len(orientation) != 6 or max(abs(value-expected) for value, expected in
			zip(orientation, (1, 0, 0, 0, 1, 0))) > 1e-6:
		reject("only unflipped axial LPS orientation is supported")
	spacing = tuple(map(float, shared.PixelMeasuresSequence[0].PixelSpacing))
	if len(spacing) != 2 or not all(math.isfinite(x) and x > 0 for x in spacing):
		reject("invalid in-plane spacing")
	frames = data.PerFrameFunctionalGroupsSequence
	if len(frames) != int(data.NumberOfFrames):
		reject("SEG frame count mismatch")
	indices = {number: {} for number in selected}
	xy = None
	for index, frame in enumerate(frames):
		number = int(frame.SegmentIdentificationSequence[0].ReferencedSegmentNumber)
		if number not in indices:
			continue
		position = tuple(map(float, frame.PlanePositionSequence[0].ImagePositionPatient))
		if len(position) != 3 or not all(map(math.isfinite, position)):
			reject("nonfinite or incomplete frame position")
		if xy is None:
			xy = position[:2]
		elif max(abs(position[axis]-xy[axis]) for axis in (0, 1)) > 1e-4:
			reject("selected segments do not share one in-plane grid")
		z = round(position[2], 4)
		if z in indices[number] or abs(z-position[2]) > 1e-4:
			reject("duplicate or imprecise segment z position")
		indices[number][z] = index
	if any(not indices[number] for number in selected):
		reject("a selected segment has no frames")
	z_values = sorted({z for by_z in indices.values() for z in by_z})
	if len(z_values) < 2:
		reject("at least two distinct z positions are required")
	z_step = (z_values[-1]-z_values[0])/(len(z_values)-1)
	if z_step <= 0 or any(abs(z-(z_values[0]+i*z_step)) > 1e-4
			for i, z in enumerate(z_values)):
		reject("selected segments are not on one regular z grid")
	array = np.asarray(data.pixel_array)
	if array.shape != (len(frames), int(data.Rows), int(data.Columns)):
		reject("unexpected SEG pixel array shape")
	if not np.isin(array, [0, 1]).all():
		reject("binary SEG contains nonbinary pixel values")
	shape = (len(z_values), int(data.Rows), int(data.Columns))
	masks = {}
	for number in selected:
		mask = np.zeros(shape, dtype=bool)
		for z_index, z in enumerate(z_values):
			if z in indices[number]:
				mask[z_index] = array[indices[number][z]]
		if not mask.any():
			reject("selected segment has no nonzero voxels")
		masks[number] = mask
	tissue = masks[tissue_number]
	results = []
	vessel_union = np.zeros(shape, dtype=bool)
	vessel_overlap = False
	for number in vessel_numbers:
		mask = masks[number]
		vessel_overlap |= bool(np.any(vessel_union & mask))
		vessel_union |= mask
	tissue_only = tissue & ~vessel_union
	if not tissue_only.any():
		reject("tissue minus vessel union has no voxels")
	for number in vessel_numbers:
		mask = masks[number]
		components, count = ndimage.label(mask)
		sizes = np.bincount(components.ravel())[1:]
		boxes = component_bounding_boxes_zyx(components, count)
		inside = int(np.count_nonzero(mask & tissue))
		outside = int(np.count_nonzero(mask & ~tissue))
		component_summary = None
		if not vessel_overlap:
			contact_counts = component_face_contacts(tissue_only, components, count)
			if int(contact_counts.sum()) != face_contact(tissue_only, mask,
					(spacing[1], spacing[0], z_step))["total_contact_faces"]:
				reject("component and whole-mask tissue contacts disagree")
			component_summary = {
				"components_touching_tissue_by_face": int(np.count_nonzero(contact_counts[1:])),
				"components_without_tissue_face_contact": int(np.count_nonzero(contact_counts[1:] == 0)),
				"voxels_without_tissue_face_contact": int(sizes[contact_counts[1:] == 0].sum()),
				"largest_component_tissue_contact_faces": int(contact_counts[int(np.argmax(sizes))+1]),
				"components": [{**boxes[index], "voxels": int(size),
					"tissue_contact_faces": int(contact_counts[index+1])}
					for index, size in enumerate(sizes)]}
		results.append({"segment_number": number, "segment_label": segments[number],
			"nonzero_voxels": int(mask.sum()), "inside_tissue_mask_voxels": inside,
			"outside_tissue_mask_voxels": outside,
			"components_6_neighbour": int(count),
			"largest_component_voxels": int(sizes.max()),
			"component_tissue_contact": component_summary})
	pairwise_overlap = [{"segments": [left, right],
		"overlap_voxels": int(np.count_nonzero(masks[left] & masks[right]))}
		for offset, left in enumerate(vessel_numbers)
		for right in vessel_numbers[offset+1:]]
	contacts = None if vessel_overlap else [
		{"segment_number": number, "segment_label": segments[number],
			**face_contact(tissue_only, masks[number],
				(spacing[1], spacing[0], z_step))}
		for number in vessel_numbers]
	result = {"schema_version": 2, "kind": "same_seg_voxel_overlap_and_face_contact_not_mesh_or_flow",
		"source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
		"study_instance_uid": str(data.StudyInstanceUID),
		"seg_series_instance_uid": str(data.SeriesInstanceUID),
		"referenced_series_instance_uid": str(
			data.ReferencedSeriesSequence[0].SeriesInstanceUID),
		"coordinate_system": "DICOM LPS", "length_unit": "mm",
		"voxel_spacing_mm": [spacing[1], spacing[0], z_step],
		"z_positions": len(z_values), "tissue": {"segment_number": tissue_number,
			"segment_label": segments[tissue_number], "nonzero_voxels": int(tissue.sum())},
		"vessels": results, "vessel_pairwise_overlap": pairwise_overlap,
		"vessel_masks_mutually_exclusive": not vessel_overlap,
		"tissue_minus_vessel_union_voxels": int(np.count_nonzero(tissue_only)),
		"tissue_minus_vessel_union_components_6_neighbour": None,
		"tissue_minus_vessel_union_role": "voxel set difference only; not a meshed or validated Darcy region",
		"tissue_vessel_voxel_face_contacts": contacts,
		"face_contact_role": "candidate grid-face adjacency of disjoint voxel sets only; no smooth/conforming tetra interface, port, BC, or physiology",
		"interpretation": "Mask overlap and grid-face adjacency only; no port, conforming tetra interface, BC, or physiology inferred"}
	tissue_components, tissue_count = ndimage.label(tissue_only)
	tissue_sizes = np.bincount(tissue_components.ravel())[1:]
	result["tissue_minus_vessel_union_components_6_neighbour"] = {
		"count": int(tissue_count),
		"largest_component_voxels": int(tissue_sizes.max()),
		"voxels_outside_largest_component": int(tissue_sizes.sum()-tissue_sizes.max()),
		"components": [{**box, "voxels": int(tissue_sizes[index])}
			for index, box in enumerate(component_bounding_boxes_zyx(
				tissue_components, tissue_count))]}
	if return_masks:
		return result, masks, (xy[0], xy[1], z_values[0])
	return result


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("seg_dicom", type=Path)
	parser.add_argument("tissue_segment_number", type=int)
	parser.add_argument("vessel_segment_numbers", type=int, nargs="+")
	parser.add_argument("--manifest", type=Path, required=True)
	parser.add_argument("--roi", help="optional half-open z0:z1,y0:y1,x0:x1 membership audit")
	args = parser.parse_args()
	try:
		if args.manifest.resolve() == args.seg_dicom.resolve():
			reject("manifest must not overwrite source SEG")
		if args.manifest.exists():
			reject("manifest already exists; choose an immutable output path")
		if args.roi is None:
			result = audit(args.seg_dicom, args.tissue_segment_number,
				args.vessel_segment_numbers)
		else:
			result, masks, _ = audit(args.seg_dicom, args.tissue_segment_number,
				args.vessel_segment_numbers, return_masks=True)
			roi = parse_roi_zyx(args.roi, masks[args.tissue_segment_number].shape)
			vessel_union = masks[args.vessel_segment_numbers[0]].copy()
			for number in args.vessel_segment_numbers[1:]:
				vessel_union |= masks[number]
			tissue_only = masks[args.tissue_segment_number] & ~vessel_union
			result["roi_component_membership"] = {
				"roi_zyx_half_open": roi,
				"tissue_candidate": roi_component_membership(tissue_only, roi),
				"vessels": {str(number): roi_component_membership(masks[number], roi)
					for number in args.vessel_segment_numbers},
				"role": "full-grid component IDs intersecting ROI only; no port, BC, or component-selection policy"}
		args.manifest.parent.mkdir(parents=True, exist_ok=True)
		with args.manifest.open("x", encoding="utf-8") as stream:
			stream.write(json.dumps(result, indent=2, sort_keys=True)+"\n")
	except (AttributeError, IndexError, KeyError, OSError, ValueError) as error:
		print(f"DICOM SEG overlap rejected: {error}", file=sys.stderr)
		return 2
	print(f"DICOM SEG overlap: PASS tissue={args.tissue_segment_number} "
		f"vessels={args.vessel_segment_numbers} (voxel audit only)")
	return 0


if __name__ == "__main__":
	sys.exit(main())
