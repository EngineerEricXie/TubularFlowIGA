#!/usr/bin/env python3
"""Compare independently published artificial oxygen time-series fields by global ID."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


def field(path):
	reader = vtk.vtkXMLPUnstructuredGridReader()
	reader.SetFileName(str(path))
	reader.Update()
	grid = reader.GetOutput()
	if not grid or not grid.GetNumberOfCells():
		raise ValueError(f"oxygen PVTU is empty: {path}")
	point_data = grid.GetPointData()
	ids_array = point_data.GetArray("GlobalPointIds")
	concentration_array = point_data.GetArray("concentration_mol_m3")
	if ids_array is None or concentration_array is None:
		raise ValueError(f"oxygen PVTU field schema differs: {path}")
	ids = np.asarray(vtk_to_numpy(ids_array), dtype=np.int64)
	concentration = np.asarray(vtk_to_numpy(concentration_array), dtype=np.float64)
	order = np.argsort(ids, kind="stable")
	unique, first = np.unique(ids[order], return_index=True)
	if not np.array_equal(unique, np.arange(unique.size, dtype=np.int64)) or \
		concentration.shape != ids.shape or not np.isfinite(concentration).all():
		raise ValueError(f"oxygen PVTU global IDs or concentration differ: {path}")
	sorted_values = concentration[order]
	if np.max(np.abs(sorted_values-sorted_values[first][np.searchsorted(
		unique, ids[order])])) > 1e-10:
		raise ValueError(f"oxygen shared point values differ: {path}")
	return unique, sorted_values[first], grid.GetNumberOfCells()


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("first", type=Path)
	parser.add_argument("second", type=Path)
	args = parser.parse_args()
	roots = (args.first.resolve(), args.second.resolve())
	provenance = [json.loads((root/"provenance.json").read_text(encoding="utf-8"))
		for root in roots]
	for key in ("kind", "case_sha256", "flow_summary_sha256",
		"split_manifest_sha256", "frozen_field_manifest_sha256", "solver_sha256"):
		if provenance[0].get(key) != provenance[1].get(key):
			raise ValueError(f"oxygen rerun provenance differs: {key}")
	if provenance[0]["ranks"] == provenance[1]["ranks"]:
		raise ValueError("oxygen comparison requires different MPI ranks")
	summary = [json.loads((root/"transport"/"summary.json").read_text(
		encoding="utf-8")) for root in roots]
	steps = summary[0]["accepted_steps"]
	if steps != summary[1]["accepted_steps"] or steps < 2:
		raise ValueError("oxygen accepted step counts differ")
	max_difference = 0.
	arrays = 0
	for region in ("artery", "tissue", "vein"):
		for step in range(steps+1):
			values = [field(root/"transport"/f"{region}_step_{step}"/
				"snapshot.pvtu") for root in roots]
			if not np.array_equal(values[0][0], values[1][0]) or \
				values[0][2] != values[1][2]:
				raise ValueError(f"oxygen {region} step {step} mesh differs")
			max_difference = max(max_difference, float(np.max(np.abs(
				values[0][1]-values[1][1]))))
			arrays += 1
	if max_difference > 1e-9:
		raise ValueError(f"oxygen cross-rank concentration difference {max_difference:g}")
	ledgers = []
	for root in roots:
		with (root/"transport"/"ledger.csv").open(newline="", encoding="utf-8") as stream:
			ledgers.append(list(csv.DictReader(stream)))
	if len(ledgers[0]) != steps or len(ledgers[1]) != steps:
		raise ValueError("oxygen ledger step counts differ")
	for first, second in zip(*ledgers):
		if first.keys() != second.keys() or first["step"] != second["step"]:
			raise ValueError("oxygen ledger columns or steps differ")
		for key in first:
			if key == "step":
				continue
			a, b = float(first[key]), float(second[key])
			if key.endswith("_reason"):
				if a <= 0 or b <= 0:
					raise ValueError(f"oxygen {key} did not converge")
				continue
			if not np.isfinite(a) or not np.isfinite(b) or \
				abs(a-b) > 1e-9*max(1e-10, abs(a), abs(b)):
				raise ValueError(f"oxygen ledger differs at {key}")
	print(f"artificial oxygen cross-rank comparison: PASS "
		f"ranks={provenance[0]['ranks']}/{provenance[1]['ranks']} "
		f"arrays={arrays} max_abs_difference_mol_m3={max_difference:.12g}")


if __name__ == "__main__":
	main()
