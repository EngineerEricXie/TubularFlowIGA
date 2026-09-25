#!/usr/bin/env python3
"""Extract checked, global-ID-ordered frozen FSI fields for native transport."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


def sha256(path):
	digest = hashlib.sha256()
	with path.open("rb") as stream:
		for block in iter(lambda: stream.read(1024*1024), b""):
			digest.update(block)
	return digest.hexdigest()


def read_grid(path):
	reader = vtk.vtkXMLPUnstructuredGridReader()
	reader.SetFileName(str(path))
	reader.Update()
	grid = reader.GetOutput()
	if not grid or not grid.GetNumberOfCells():
		raise ValueError(f"empty FSI field: {path}")
	return grid


def ordered_points(grid, fields):
	data = grid.GetPointData()
	ids_array = data.GetArray("GlobalPointIds")
	if ids_array is None:
		raise ValueError("FSI point global IDs are absent")
	ids = np.asarray(vtk_to_numpy(ids_array), dtype=np.int64)
	if ids.size == 0 or ids.min() != 0:
		raise ValueError("FSI point global IDs do not start at zero")
	values = []
	for name in fields:
		field = data.GetArray(name)
		if field is None:
			raise ValueError(f"FSI point field is absent: {name}")
		array = np.asarray(vtk_to_numpy(field), dtype=np.float64)
		if array.shape != (ids.size, 3) or not np.isfinite(array).all():
			raise ValueError(f"FSI point field is invalid: {name}")
		values.append(array)
	order = np.argsort(ids, kind="stable")
	sorted_ids = ids[order]
	unique, first = np.unique(sorted_ids, return_index=True)
	if not np.array_equal(unique, np.arange(unique.size, dtype=np.int64)):
		raise ValueError("FSI point global IDs are incomplete")
	result = []
	for name, array in zip(fields, values):
		sorted_values = array[order]
		if np.max(np.abs(sorted_values - sorted_values[first][np.searchsorted(
			unique, sorted_ids)])) > 1e-12:
			raise ValueError(f"FSI shared point field differs: {name}")
		result.append(sorted_values[first])
	return unique, result


def fluid_rows(path):
	grid = read_grid(path)
	ids, (velocity, reference, displacement) = ordered_points(grid,
		("velocity_m_s", "reference_position_m", "displacement_m"))
	current = np.asarray(vtk_to_numpy(grid.GetPoints().GetData()), dtype=np.float64)
	point_ids = np.asarray(vtk_to_numpy(
		grid.GetPointData().GetArray("GlobalPointIds")), dtype=np.int64)
	if current.shape != (point_ids.size, 3) or not np.isfinite(current).all():
		raise ValueError("FSI fluid coordinates are invalid")
	if np.max(np.abs(current-(reference+displacement)[point_ids])) > 1e-12:
		raise ValueError("FSI fluid coordinates disagree with displacement")
	return np.column_stack((ids, velocity, reference, reference+displacement))


def tissue_rows(path):
	grid = read_grid(path)
	data = grid.GetCellData()
	ids_array = data.GetArray("GlobalCellIds")
	flow_array = data.GetArray("conservative_outward_face_flow_m3_s")
	if ids_array is None or flow_array is None:
		raise ValueError("FSI Darcy cell IDs or conservative RT0 flows are absent")
	ids = np.asarray(vtk_to_numpy(ids_array), dtype=np.int64)
	flows = np.asarray(vtk_to_numpy(flow_array), dtype=np.float64)
	if ids.size != grid.GetNumberOfCells() or flows.shape != (ids.size, 4) or \
		(ids < 0).any() or len(np.unique(ids)) != ids.size or \
		not np.isfinite(flows).all():
		raise ValueError("FSI Darcy cell IDs or RT0 flows are invalid")
	order = np.argsort(ids)
	return np.column_stack((ids[order], flows[order]))


def write_rows(path, rows):
	formats = ["%d"] + ["%.17g"]*(rows.shape[1]-1)
	with path.open("w", encoding="utf-8") as stream:
		stream.write(f"{rows.shape[0]}\n")
		np.savetxt(stream, rows, fmt=formats)


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("flow_run", type=Path)
	parser.add_argument("split_directory", type=Path)
	parser.add_argument("output_directory", type=Path)
	args = parser.parse_args()
	run, split, output = (path.resolve() for path in
		(args.flow_run, args.split_directory, args.output_directory))
	summary_path = run/"summary.json"
	summary = json.loads(summary_path.read_text(encoding="utf-8"))
	if summary.get("kind") != "idealized_cube_dual_tree_quasi_steady_fsi_functional_only" \
		or summary.get("physiological_validation") is not False or \
		sha256(split/"manifest.json") != summary.get("split_manifest_sha256"):
		raise ValueError("flow classification or independently rebuilt split mesh differs")
	if output.exists() or not output.parent.is_dir():
		raise ValueError("oxygen field output exists or parent is absent")
	artery = fluid_rows(run/"fields"/"arterial_fsi"/"snapshot.pvtu")
	vein = fluid_rows(run/"fields"/"venous_fsi"/"snapshot.pvtu")
	tissue = tissue_rows(run/"fields"/"tissue_fsi"/"snapshot.pvtu")
	output.mkdir()
	files = {"arterial_fsi_state": artery, "venous_fsi_state": vein,
		"tissue_rt0_face_flow": tissue}
	for name, rows in files.items():
		write_rows(output/(name+".txt"), rows)
	(output/"source_manifest.json").write_text(json.dumps({
		"schema_version": 1,
		"kind": "idealized_cube_frozen_fsi_oxygen_input_functional_only",
		"physiological_validation": False,
		"flow_summary_sha256": sha256(summary_path),
		"case_sha256": summary["case_sha256"],
		"split_manifest_sha256": summary["split_manifest_sha256"],
		"field_sha256": {name: sha256(output/(name+".txt")) for name in files},
		"counts": {name: int(rows.shape[0]) for name, rows in files.items()},
	}, indent=2, sort_keys=True)+"\n", encoding="utf-8")
	print(f"frozen FSI oxygen fields: PASS {output}")


if __name__ == "__main__":
	main()
