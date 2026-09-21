#!/usr/bin/env python3
"""Export artificial dual-tree pipe and wall surfaces as ParaView PVD state series."""

import argparse
import json
import math
import os
from pathlib import Path
import tempfile
import xml.etree.ElementTree as ET

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy


def read_grid(path):
	reader = vtk.vtkXMLPUnstructuredGridReader()
	reader.SetFileName(str(path))
	reader.Update()
	grid = reader.GetOutput()
	if not grid or not grid.GetNumberOfPoints() or not grid.GetNumberOfCells():
		raise ValueError(f"empty or unreadable field: {path}")
	return grid


def reference_and_displacement(grid, name):
	data = grid.GetPointData()
	for field in ("GlobalPointIds", "reference_position_m", "displacement_m"):
		if data.GetArray(field) is None:
			raise ValueError(f"{name} lacks {field}")
	ref = np.asarray(vtk_to_numpy(data.GetArray("reference_position_m")))
	disp = np.asarray(vtk_to_numpy(data.GetArray("displacement_m")))
	current = np.asarray(vtk_to_numpy(grid.GetPoints().GetData()))
	if ref.shape != disp.shape or ref.shape != current.shape or ref.shape[1] != 3:
		raise ValueError(f"{name} coordinates or displacement have invalid shape")
	if not all(np.isfinite(array).all() for array in (ref, disp, current)):
		raise ValueError(f"{name} geometry is nonfinite")
	if np.max(np.abs(current-(ref+disp))) > 1e-12:
		raise ValueError(f"{name} points disagree with reference + displacement")
	return ref, disp


def surface_at_state(grid, scale, initial, name):
	ref, disp = reference_and_displacement(grid, name)
	if not initial and np.linalg.norm(disp, axis=1).max() <= 0:
		raise ValueError(f"{name} final state does not move")
	state = vtk.vtkUnstructuredGrid()
	state.DeepCopy(grid)
	warped = np.ascontiguousarray(ref if initial else ref+scale*disp,
		dtype=np.float64)
	if not initial:
		# A display-only warp must not turn an otherwise valid pipe inside out.
		# This corner-tetra check is necessary, though not a full curved-face
		# self-intersection or high-order Jacobian certificate.
		tetra_volumes(grid, warped, name+f" display x{scale:g}")
	points = vtk.vtkPoints()
	points.SetData(numpy_to_vtk(warped, deep=True))
	state.SetPoints(points)
	if initial:
		zero = numpy_to_vtk(np.zeros_like(disp), deep=True)
		zero.SetName("displacement_m")
		state.GetPointData().RemoveArray("displacement_m")
		state.GetPointData().AddArray(zero)
	filter_surface = vtk.vtkDataSetSurfaceFilter()
	filter_surface.SetInputData(state)
	filter_surface.Update()
	surface = vtk.vtkPolyData()
	surface.DeepCopy(filter_surface.GetOutput())
	if surface.GetNumberOfPolys() <= 0:
		raise ValueError(f"{name} surface extraction failed")
	# High-order tetra faces create interpolated surface points. Their sampled
	# point-data fields, not input point indices, are the authoritative mapping.
	surface_ref = np.asarray(vtk_to_numpy(surface.GetPointData().GetArray(
		"reference_position_m")))
	surface_disp = np.asarray(vtk_to_numpy(surface.GetPointData().GetArray(
		"displacement_m")))
	if surface_ref.shape != surface_disp.shape or surface_ref.shape[1] != 3 or \
		not np.isfinite(surface_ref).all() or not np.isfinite(surface_disp).all():
		raise ValueError(f"{name} interpolated surface fields are invalid")
	double_points = vtk.vtkPoints()
	double_points.SetData(numpy_to_vtk(np.ascontiguousarray(
		surface_ref+scale*surface_disp,
		dtype=np.float64), deep=True))
	surface.SetPoints(double_points)
	metadata = vtk.vtkDoubleArray()
	metadata.SetName("display_displacement_scale")
	metadata.InsertNextValue(scale)
	surface.GetFieldData().AddArray(metadata)
	return surface


def write_surface(path, surface):
	writer = vtk.vtkXMLPolyDataWriter()
	writer.SetFileName(str(path))
	writer.SetInputData(surface)
	writer.SetDataModeToBinary()
	if writer.Write() != 1:
		raise ValueError(f"could not write {path}")


def write_volume(path, grid):
	writer = vtk.vtkXMLUnstructuredGridWriter()
	writer.SetFileName(str(path))
	writer.SetInputData(grid)
	writer.SetDataModeToBinary()
	if writer.Write() != 1:
		raise ValueError(f"could not write {path}")


def add_cell_array(grid, name, values):
	array = numpy_to_vtk(np.ascontiguousarray(values, dtype=np.float64), deep=True)
	array.SetName(name)
	grid.GetCellData().AddArray(array)


def tetra_volumes(grid, positions, name):
	if not all(grid.GetCellType(i) in (10, 24) for i in range(grid.GetNumberOfCells())):
		raise ValueError(f"{name} contains non-tetrahedral cells")
	connectivity = np.array([[grid.GetCell(i).GetPointId(j) for j in range(4)]
		for i in range(grid.GetNumberOfCells())], dtype=np.int64)
	p = positions[connectivity]
	volumes = np.einsum("ij,ij->i",
		np.cross(p[:, 1]-p[:, 0], p[:, 2]-p[:, 0]),
		p[:, 3]-p[:, 0]) / 6.
	if not np.isfinite(volumes).all() or np.any(volumes <= 0):
		raise ValueError(f"{name} has a nonpositive tetrahedron Jacobian")
	return volumes


def deformation_volume(grid, name):
	ref, _ = reference_and_displacement(grid, name)
	current = vtk_to_numpy(grid.GetPoints().GetData())
	initial = tetra_volumes(grid, ref, name+" reference")
	final = tetra_volumes(grid, current, name+" final")
	result = vtk.vtkUnstructuredGrid()
	result.DeepCopy(grid)
	add_cell_array(result, "reference_cell_volume_m3", initial)
	add_cell_array(result, "deformed_cell_volume_m3", final)
	add_cell_array(result, "deformation_jacobian", final/initial)
	add_cell_array(result, "relative_cell_volume_change", final/initial-1.)
	return result, {"reference_volume_m3": float(initial.sum()),
		"deformed_volume_m3": float(final.sum()),
		"volume_change_m3": float((final-initial).sum()),
		"relative_volume_change_percent": float(100*((final-initial).sum()/initial.sum())),
		"minimum_deformation_jacobian": float(np.min(final/initial)),
		"maximum_deformation_jacobian": float(np.max(final/initial))}


def darcy_volume(grid):
	field = grid.GetCellData().GetArray("darcy_rt0_centroid_flux_m_s")
	if field is None:
		raise ValueError("Darcy result lacks conservative RT0 centroid velocity")
	velocity = np.asarray(vtk_to_numpy(field))
	if velocity.shape != (grid.GetNumberOfCells(), 3) or \
		not np.isfinite(velocity).all():
		raise ValueError("Darcy velocity is invalid")
	result = vtk.vtkUnstructuredGrid()
	result.DeepCopy(grid)
	add_cell_array(result, "darcy_velocity_m_s", velocity)
	speed = np.linalg.norm(velocity, axis=1)
	add_cell_array(result, "darcy_speed_m_s", speed)
	return result, velocity, speed


def darcy_arrows(grid, velocity, speed):
	centers = vtk.vtkPolyData()
	points = vtk.vtkPoints()
	vectors = vtk.vtkDoubleArray()
	vectors.SetName("darcy_velocity_m_s")
	vectors.SetNumberOfComponents(3)
	scalars = vtk.vtkDoubleArray()
	scalars.SetName("darcy_speed_m_s")
	vertices = vtk.vtkCellArray()
	step = max(1, grid.GetNumberOfCells() // 180)
	for index in range(0, grid.GetNumberOfCells(), step):
		cell = grid.GetCell(index)
		xyz = np.mean([grid.GetPoint(cell.GetPointId(j)) for j in range(4)], axis=0)
		point = points.InsertNextPoint(*xyz)
		vertices.InsertNextCell(1)
		vertices.InsertCellPoint(point)
		vectors.InsertNextTuple(velocity[index])
		scalars.InsertNextValue(speed[index])
	centers.SetPoints(points)
	centers.SetVerts(vertices)
	centers.GetPointData().AddArray(vectors)
	centers.GetPointData().SetActiveVectors("darcy_velocity_m_s")
	centers.GetPointData().AddArray(scalars)
	centers.GetPointData().SetActiveScalars("darcy_speed_m_s")
	arrow = vtk.vtkArrowSource()
	glyph = vtk.vtkGlyph3D()
	glyph.SetSourceConnection(arrow.GetOutputPort())
	glyph.SetInputData(centers)
	glyph.OrientOn()
	glyph.SetVectorModeToUseVector()
	glyph.SetScaleModeToDataScalingOff()
	glyph.SetScaleFactor(.001)  # Display-only 1 mm arrows; speed remains in SI.
	glyph.Update()
	result = vtk.vtkPolyData()
	result.DeepCopy(glyph.GetOutput())
	if result.GetNumberOfPoints() <= 0:
		raise ValueError("Darcy velocity arrow generation failed")
	return result


def write_pvd(path, names):
	root = ET.Element("VTKFile", type="Collection", version="0.1",
		byte_order="LittleEndian")
	collection = ET.SubElement(root, "Collection")
	for index, name in enumerate(names):
		ET.SubElement(collection, "DataSet", timestep=str(index), group="",
			part="0", file=name)
	ET.indent(root, space="  ")
	ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("run_directory", type=Path)
	parser.add_argument("--output", type=Path)
	parser.add_argument("--magnification", type=float, default=1000.)
	args = parser.parse_args()
	if not math.isfinite(args.magnification) or args.magnification <= 1:
		raise ValueError("magnification must be finite and greater than one")
	run = args.run_directory.resolve()
	summary = json.loads((run/"summary.json").read_text(encoding="utf-8"))
	if summary.get("kind") != "idealized_cube_dual_tree_quasi_steady_fsi_functional_only" or \
		not summary.get("vessel_wall_fsi") or summary.get("physiological_validation"):
		raise ValueError("run is not a classified artificial dual-tree FSI result")
	output = args.output.resolve() if args.output else run/"paraview"
	if output.exists() or not output.parent.is_dir():
		raise ValueError("output exists or its parent is missing; refusing overwrite")
	with tempfile.TemporaryDirectory(prefix="paraview-dual-tree-", dir=output.parent) as staging:
		stage = Path(staging)
		volume_report = {"schema_version": 1,
			"kind": "idealized_cube_linear_tetra_deformation_audit",
			"physiological_validation": False,
			"case_sha256": summary["case_sha256"],
			"solver_sha256": summary["solver_sha256"]}
		for vessel in ("arterial", "venous"):
			initial = read_grid(run/"fields"/vessel/"snapshot.pvtu")
			final = read_grid(run/"fields"/(vessel+"_fsi")/"snapshot.pvtu")
			wall = read_grid(run/"fields"/(vessel+"_wall")/"snapshot.pvtu")
			for name, grid in (("lumen", final), ("wall", wall)):
				deformed, report = deformation_volume(grid, vessel+"_"+name)
				volume_report[vessel+"_"+name] = report
				write_volume(stage/f"{vessel}_{name}_deformation.vtu", deformed)
			for name, first, second in (("pipe", initial, final),
				("wall", wall, wall)):
				for suffix, scale in (("actual", 1.),
					(f"x{args.magnification:g}", args.magnification)):
					stem = f"{vessel}_{name}_{suffix}"
					files = [f"{stem}_state{state}.vtp" for state in (0, 1)]
					for state, file in enumerate(files):
						grid = first if state == 0 else second
						write_surface(stage/file, surface_at_state(grid, scale,
							state == 0, stem))
					write_pvd(stage/(stem+".pvd"), files)
		arterial_delta = volume_report["arterial_lumen"]["volume_change_m3"] + \
			volume_report["arterial_wall"]["volume_change_m3"]
		venous_delta = volume_report["venous_lumen"]["volume_change_m3"] + \
			volume_report["venous_wall"]["volume_change_m3"]
		if max(abs(arterial_delta), abs(venous_delta)) > 1e-10 * max(
			volume_report["arterial_lumen"]["reference_volume_m3"],
			volume_report["venous_lumen"]["reference_volume_m3"]):
			raise ValueError("fixed outer wall and lumen volumes do not balance")
		volume_report["arterial_lumen_plus_wall_change_m3"] = arterial_delta
		volume_report["venous_lumen_plus_wall_change_m3"] = venous_delta
		tissue_initial = read_grid(run/"fields"/"tissue"/"snapshot.pvtu")
		tissue = read_grid(run/"fields"/"tissue_fsi"/"snapshot.pvtu")
		initial_ids = vtk_to_numpy(tissue_initial.GetPointData().GetArray("GlobalPointIds"))
		final_ids = vtk_to_numpy(tissue.GetPointData().GetArray("GlobalPointIds"))
		initial_coordinates = vtk_to_numpy(tissue_initial.GetPoints().GetData())
		final_coordinates = vtk_to_numpy(tissue.GetPoints().GetData())
		if not np.array_equal(initial_ids, final_ids) or \
			not np.array_equal(initial_coordinates, final_coordinates):
			raise ValueError("fixed Darcy tissue geometry changed")
		volume_report["fixed_tissue_volume_m3"] = float(tetra_volumes(tissue,
			final_coordinates, "fixed Darcy tissue").sum())
		tissue_grid, velocity, speed = darcy_volume(tissue)
		write_volume(stage/"tissue_darcy_velocity.vtu", tissue_grid)
		write_surface(stage/"tissue_darcy_direction_arrows.vtp",
			darcy_arrows(tissue_grid, velocity, speed))
		volume_report["tissue_darcy_speed_m_s"] = {
			"minimum": float(np.min(speed)), "median": float(np.median(speed)),
			"maximum": float(np.max(speed))}
		(stage/"volume_audit.json").write_text(json.dumps(volume_report, indent=2,
			sort_keys=True)+"\n", encoding="utf-8")
		(stage/"README.txt").write_text(
			"Artificial cube functional case, NOT anatomy or physiology.\n"
			"Open a .pvd file in ParaView and press Play; steps 0/1 are states, "
			"not physical seconds.\n"
			"*_actual uses the solved displacement in metres.\n"
			f"*_x{args.magnification:g} scales displacement for visibility; "
			"the physical field displacement_m remains unscaled.\n"
			"pipe = fluid lumen boundary; wall = solid vessel-wall exterior/interior.\n"
			"tissue_darcy_velocity.vtu contains the FIXED tetra mesh and cell-centroid "
			"conservative RT0 velocity (m/s). Use Clip or Slice to see inside.\n"
			"tissue_darcy_direction_arrows.vtp is a sparse display-only sample: "
			"arrows have fixed 1 mm display length, not speed-proportional length.\n"
			"*_deformation.vtu includes physical cell-volume ratio/Jacobian; "
			"volume_audit.json reports reference/final volumes.\n"
			"Tissue remains fixed; hydraulic pressure continuity is not coupled.\n",
			encoding="utf-8")
		os.replace(stage, output)
	print(f"ParaView surface, Darcy and volume fields: PASS {output}")


if __name__ == "__main__":
	main()
