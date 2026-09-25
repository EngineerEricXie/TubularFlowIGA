#!/usr/bin/env python3
"""Compare two independently reconstructed cube FSI runs by global VTK IDs."""

import argparse
import json
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


FIELDS = ("arterial", "tissue", "venous", "arterial_wall", "venous_wall",
	"arterial_fsi", "tissue_fsi", "venous_fsi")


def canonical(data, identifiers, field):
	ids = vtk_to_numpy(data.GetArray(identifiers)).astype(np.int64)
	if len(ids) == 0:
		raise ValueError(f"empty {field} identifiers")
	unique, first, counts = np.unique(ids, return_index=True, return_counts=True)
	arrays = {}
	for index in range(data.GetNumberOfArrays()):
		name = data.GetArrayName(index)
		if name == identifiers:
			continue
		values = np.asarray(vtk_to_numpy(data.GetArray(index)))
		if not np.isfinite(values).all():
			raise ValueError(f"nonfinite {field}/{name}")
		ordered = values[np.argsort(ids)]
		start = np.r_[0, np.cumsum(counts)[:-1]]
		selected = ordered[start]
		if not np.allclose(ordered, np.repeat(selected, counts, axis=0),
			rtol=1e-10, atol=1e-12):
			raise ValueError(f"duplicate {field}/{name} global IDs disagree")
		arrays[name] = selected
	return unique, arrays


def read_field(run, name):
	reader = vtk.vtkXMLPUnstructuredGridReader()
	reader.SetFileName(str(run / "fields" / name / "snapshot.pvtu"))
	reader.Update()
	grid = reader.GetOutput()
	if grid.GetNumberOfPoints() == 0 or grid.GetNumberOfCells() == 0:
		raise ValueError(f"empty {name} output in {run}")
	points = canonical(grid.GetPointData(), "GlobalPointIds", name + "/point")
	cells = canonical(grid.GetCellData(), "GlobalCellIds", name + "/cell")
	return points, cells


def compare_array(one, two, label):
	if one.shape != two.shape:
		raise ValueError(f"{label} shape mismatch")
	maximum = float(np.max(np.abs(one-two)))
	scale = max(float(np.max(np.abs(one))), float(np.max(np.abs(two))), 1e-12)
	if maximum > 1e-8*scale + 1e-12:
		raise ValueError(f"{label} mismatch: maximum={maximum:.6g}, scale={scale:.6g}")
	return maximum


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("first_run", type=Path)
	parser.add_argument("second_run", type=Path)
	args = parser.parse_args()
	paths = [args.first_run.resolve(), args.second_run.resolve()]
	summaries = [json.loads((path / "summary.json").read_text(encoding="utf-8"))
		for path in paths]
	for key in ("case_sha256", "solver_sha256", "geometry_manifest_sha256",
		"split_manifest_sha256", "kind", "fsi_model"):
		if summaries[0][key] != summaries[1][key]:
			raise ValueError(f"{key} differs")
	if not all(s["vessel_wall_fsi"] and not s["physiological_validation"]
		for s in summaries):
		raise ValueError("case classification differs")
	for key in summaries[0]["flow"]:
		compare_array(np.array([summaries[0]["flow"][key]]),
			np.array([summaries[1]["flow"][key]]), "flow/" + key)
	for ledger in ("terminals", "fsi_terminals"):
		if summaries[0][ledger].keys() != summaries[1][ledger].keys():
			raise ValueError(f"{ledger} IDs differ")
		for index, first in summaries[0][ledger].items():
			second = summaries[1][ledger][index]
			if first.keys() != second.keys():
				raise ValueError(f"{ledger}/{index} fields differ")
			for field in first:
				compare_array(np.array([first[field]]), np.array([second[field]]),
					f"{ledger}/{index}/{field}")
	count = 0
	maximum = 0.
	for name in FIELDS:
		first, second = [read_field(path, name) for path in paths]
		for dimension in range(2):
			ids_a, arrays_a = first[dimension]
			ids_b, arrays_b = second[dimension]
			if not np.array_equal(ids_a, ids_b) or arrays_a.keys() != arrays_b.keys():
				raise ValueError(f"{name} point/cell ownership differs")
			for field in arrays_a:
				maximum = max(maximum, compare_array(arrays_a[field],
					arrays_b[field], f"{name}/{field}"))
				count += 1
	print(f"dual-tree field comparison: PASS ranks={summaries[0]['ranks']}/"
		f"{summaries[1]['ranks']} arrays={count} max_abs_difference={maximum:.6g}")


if __name__ == "__main__":
	main()
