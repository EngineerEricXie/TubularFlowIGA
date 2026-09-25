#!/usr/bin/env python3
"""Extract one DICOM SEG label as a provenance-tracked voxel-isosurface.

This intentionally supports only regular axis-aligned SEG grids. It does not
infer vascular ports, smooth anatomy, or build volume mesh. Optional largest-
component selection is explicit and records the excluded voxel count.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys


def reject(message):
	raise ValueError(message)


def extract(seg_path, segment_number, output_path, manifest_path, keep_largest_component=False,
		smooth_iterations=0, smooth_passband=None, max_displacement_mm=None,
		max_volume_relative_error=None):
	try:
		import numpy as np
		import pydicom
		from scipy import ndimage
		import vtk
		from vtk.util.numpy_support import numpy_to_vtk
	except ImportError as error:
		reject(f"DICOM SEG surface dependencies unavailable: {error}")

	data = pydicom.dcmread(str(seg_path))
	if data.Modality != "SEG" or data.SegmentationType != "BINARY":
		reject("expected a binary DICOM SEG object")
	segments = {int(item.SegmentNumber): item.SegmentLabel for item in data.SegmentSequence}
	if segment_number not in segments:
		reject(f"segment {segment_number} is missing")
	if len(data.ReferencedSeriesSequence) != 1:
		reject("expected exactly one referenced source series")
	spacing = [float(value) for value in data.SharedFunctionalGroupsSequence[0]
		.PixelMeasuresSequence[0].PixelSpacing]
	if len(spacing) != 2 or min(spacing) <= 0 or not all(np.isfinite(spacing)):
		reject("invalid pixel spacing")
	frames = data.PerFrameFunctionalGroupsSequence
	if len(frames) != int(data.NumberOfFrames):
		reject("frame count mismatch")
	xy = None
	indexed = {}
	for index, frame in enumerate(frames):
		if int(frame.SegmentIdentificationSequence[0].ReferencedSegmentNumber) != segment_number:
			continue
		position = tuple(float(value) for value in frame.PlanePositionSequence[0].ImagePositionPatient)
		if len(position) != 3 or not all(np.isfinite(position)):
			reject("invalid frame position")
		if xy is None:
			xy = position[:2]
		if max(abs(position[axis]-xy[axis]) for axis in (0, 1)) > 1e-4:
			reject("inconsistent in-plane frame position")
		if position[2] in indexed:
			reject("duplicate z position for selected segment")
		indexed[position[2]] = index
	if not indexed:
		reject("selected segment has no frames")
	if len(indexed) < 2:
		reject("at least two z positions are required")
	z_values = sorted(indexed)
	z_spacing = (z_values[-1]-z_values[0])/(len(z_values)-1)
	if z_spacing <= 0 or any(abs(z_values[i]-(z_values[0]+i*z_spacing)) > 1e-4
			for i in range(len(z_values))):
		reject("selected segment is not on a regular z grid")
	# The current voxel-index-to-LPS transform is deliberately narrow. Reject
	# oblique and flipped planes until their affine transform is implemented.
	orientation = tuple(float(value) for value in
		data.SharedFunctionalGroupsSequence[0].PlaneOrientationSequence[0].ImageOrientationPatient)
	if max(abs(orientation[i]-expected) for i, expected in enumerate((1, 0, 0, 0, 1, 0))) > 1e-6:
		reject("only unflipped axial LPS orientation is supported")
	array = np.asarray(data.pixel_array)
	if array.shape != (len(frames), int(data.Rows), int(data.Columns)):
		reject("unexpected pixel array shape")
	mask = np.stack([array[indexed[z]] for z in z_values]).astype(bool)
	if not mask.any():
		reject("selected segment is empty")
	component_ids, component_count = ndimage.label(mask)
	component_sizes = np.bincount(component_ids.ravel())[1:]
	excluded_voxels = 0
	if keep_largest_component:
		largest_id = int(np.argmax(component_sizes))+1
		excluded_voxels = int(mask.sum()-component_sizes[largest_id-1])
		mask = component_ids == largest_id
	# Padding closes masks that touch the outermost z slice; no anatomy is
	# fabricated inside the source SEG. Components are only removed by request.
	padded = np.pad(mask, 1, constant_values=False).astype(np.uint8)
	image = vtk.vtkImageData()
	image.SetDimensions(int(padded.shape[2]), int(padded.shape[1]), int(padded.shape[0]))
	image.SetSpacing(spacing[1], spacing[0], z_spacing)
	image.SetOrigin(xy[0]-spacing[1], xy[1]-spacing[0], z_values[0]-z_spacing)
	image.GetPointData().SetScalars(numpy_to_vtk(padded.ravel(order="C"), deep=True,
		array_type=vtk.VTK_UNSIGNED_CHAR))
	contour = vtk.vtkFlyingEdges3D()
	contour.SetInputData(image)
	contour.SetValue(0, 0.5)
	contour.Update()
	surface = contour.GetOutput()
	if surface.GetNumberOfPoints() == 0 or surface.GetNumberOfPolys() == 0:
		reject("isosurface is empty")
	point_displacement_mm = {"maximum": 0.0, "median": 0.0, "p95": 0.0}
	if smooth_iterations:
		if smooth_iterations < 0 or smooth_passband is None or not 0 < smooth_passband < 1:
			reject("smoothing needs positive iterations and passband in (0,1)")
		if max_displacement_mm is None or not math.isfinite(max_displacement_mm) \
				or max_displacement_mm <= 0:
			reject("smoothing needs a positive maximum point displacement in mm")
		original_points = np.array(vtk.util.numpy_support.vtk_to_numpy(
			surface.GetPoints().GetData()), dtype=np.float64, copy=True)
		smoother = vtk.vtkWindowedSincPolyDataFilter()
		smoother.SetInputData(surface)
		smoother.SetNumberOfIterations(smooth_iterations)
		smoother.SetPassBand(smooth_passband)
		smoother.NormalizeCoordinatesOn()
		smoother.BoundarySmoothingOff()
		smoother.FeatureEdgeSmoothingOff()
		smoother.Update()
		surface = smoother.GetOutput()
		if surface.GetNumberOfPoints() != len(original_points):
			reject("smoothing changed surface vertex count")
		new_points = np.asarray(vtk.util.numpy_support.vtk_to_numpy(
			surface.GetPoints().GetData()), dtype=np.float64)
		displacement = np.linalg.norm(new_points-original_points, axis=1)
		point_displacement_mm = {"maximum": float(displacement.max()),
			"median": float(np.median(displacement)),
			"p95": float(np.percentile(displacement, 95))}
		if not math.isfinite(point_displacement_mm["maximum"]) \
				or point_displacement_mm["maximum"] > max_displacement_mm:
			reject("smoothed surface exceeds maximum allowed point displacement")
	elif smooth_passband is not None or max_displacement_mm is not None:
		reject("smoothing parameters were supplied without positive iterations")
	if max_volume_relative_error is None or not math.isfinite(max_volume_relative_error) \
			or not 0 <= max_volume_relative_error < 1:
		reject("an explicit --max-volume-relative-error in [0,1) is required")
	voxel_volume_mm3 = float(mask.sum()*spacing[0]*spacing[1]*z_spacing)
	mass = vtk.vtkMassProperties()
	mass.SetInputData(surface)
	mass.Update()
	surface_volume_mm3 = float(mass.GetVolume())
	volume_relative_error = abs(surface_volume_mm3-voxel_volume_mm3)/voxel_volume_mm3
	if not math.isfinite(volume_relative_error) or volume_relative_error > max_volume_relative_error:
		reject("surface volume differs from source voxel volume beyond allowed relative error")
	surface.GetPointData().Initialize()
	surface.GetCellData().Initialize()
	output_path.parent.mkdir(parents=True, exist_ok=True)
	if output_path.suffix.lower() != ".stl":
		reject("output must have .stl extension")
	writer = vtk.vtkSTLWriter()
	writer.SetFileTypeToBinary()
	writer.SetFileName(str(output_path))
	writer.SetInputData(surface)
	if writer.Write() != 1:
		reject("could not write surface STL")
	manifest = {
		"route": "dicom_seg_to_surface", "segment_number": segment_number,
		"segment_label": segments[segment_number], "source_sha256":
		hashlib.sha256(seg_path.read_bytes()).hexdigest(),
		"surface_sha256": hashlib.sha256(output_path.read_bytes()).hexdigest(),
		"study_instance_uid": str(data.StudyInstanceUID),
		"seg_series_instance_uid": str(data.SeriesInstanceUID),
		"referenced_series_instance_uid": str(data.ReferencedSeriesSequence[0].SeriesInstanceUID),
		"source_coordinate_system": "DICOM LPS", "source_length_unit": "mm",
		"voxel_spacing_mm": [spacing[1], spacing[0], z_spacing],
		"z_range_mm": [z_values[0], z_values[-1]],
		"nonzero_voxels": int(mask.sum()),
		"source_nonzero_voxels": int(mask.sum())+excluded_voxels,
		"source_components_6_neighbour": component_count,
		"source_component_sizes_voxels": sorted(map(int, component_sizes), reverse=True),
		"keep_largest_component": bool(keep_largest_component),
		"excluded_voxels": excluded_voxels,
		"voxel_volume_mm3": voxel_volume_mm3,
		"surface_volume_mm3": surface_volume_mm3,
		"volume_relative_error": volume_relative_error,
		"max_volume_relative_error": max_volume_relative_error,
		"smoothing": {"iterations": smooth_iterations, "passband": smooth_passband,
			"max_displacement_mm": max_displacement_mm,
			"point_displacement_mm": point_displacement_mm},
		"surface_points": surface.GetNumberOfPoints(),
		"surface_triangles": surface.GetNumberOfPolys(),
		"geometry_operation": "padded 0.5 voxel isosurface; "
			+ ("explicit largest-component selection" if keep_largest_component
				else "no component removal")
			+ ("; explicit windowed-sinc smoothing" if smooth_iterations else "; no smoothing"),
		"boundary_semantics": "unassigned; original segment number only",
		"surface_format": "binary STL",
		"stl_default_boundary_id": segment_number,
	}
	manifest_path.parent.mkdir(parents=True, exist_ok=True)
	manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True)+"\n", encoding="utf-8")
	print(json.dumps(manifest, indent=2, sort_keys=True))


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("seg_dicom", type=Path)
	parser.add_argument("segment_number", type=int)
	parser.add_argument("output_stl", type=Path)
	parser.add_argument("--manifest", type=Path)
	parser.add_argument("--keep-largest-component", action="store_true",
		help="explicitly discard all but the largest 6-neighbour voxel component")
	parser.add_argument("--smooth-iterations", type=int, default=0)
	parser.add_argument("--smooth-passband", type=float)
	parser.add_argument("--max-displacement-mm", type=float)
	parser.add_argument("--max-volume-relative-error", type=float, required=True)
	args = parser.parse_args()
	manifest = args.manifest or args.output_stl.with_suffix(".json")
	try:
		extract(args.seg_dicom, args.segment_number, args.output_stl, manifest,
			args.keep_largest_component, args.smooth_iterations, args.smooth_passband,
			args.max_displacement_mm, args.max_volume_relative_error)
	except (AttributeError, KeyError, IndexError, ValueError) as error:
		print(f"DICOM SEG surface rejected: {error}", file=sys.stderr)
		return 2
	return 0


if __name__ == "__main__":
	sys.exit(main())
