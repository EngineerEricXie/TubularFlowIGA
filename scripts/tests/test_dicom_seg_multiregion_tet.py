#!/usr/bin/env python3
"""Synthetic geometry-only checks for the shared-lattice SEG tetra route."""

import json
from pathlib import Path
import sys
import tempfile

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from audit_multiregion_tet_mesh import audit, parse_msh, read_contract
from dicom_seg_to_multiregion_tet import build_mesh, make_contract, write_gmsh


def require_rejection(operation, needle):
	try:
		operation()
	except ValueError as error:
		if needle in str(error):
			return
		raise RuntimeError(f"wrong rejection: {error}") from error
	raise RuntimeError(f"expected rejection containing {needle!r}")


def main():
	tissue = np.zeros((1, 1, 2), dtype=bool)
	vessel = np.zeros_like(tissue)
	tissue[0, 0, 0] = True
	vessel[0, 0, 1] = True
	regions = {"tissue_candidate": tissue, "seg_4_portal": vessel}
	roi = ((0, 1), (0, 1), (0, 2))
	nodes, tets, boundaries, interfaces, count = build_mesh(
		regions, roi, (0, 0, 0), (1, 1, 1), 2)
	if count != 2 or len(nodes) != 12 or sum(map(len, tets.values())) != 12 \
			or list(map(len, interfaces.values())) != [2]:
		raise RuntimeError("two-voxel shared-lattice geometry changed")
	with tempfile.TemporaryDirectory(prefix="tubularflow-seg-multiregion-") as directory:
		path = Path(directory)/"two-voxel.msh"
		surfaces = write_gmsh(path, nodes, tets, boundaries, interfaces)
		contract_path = path.with_suffix(".json")
		contract_path.write_text(json.dumps(make_contract(tets, boundaries,
			interfaces, surfaces)), encoding="utf-8")
		_, region_map, boundary_map, interface_map = read_contract(contract_path)
		result = audit(*parse_msh(path), region_map, boundary_map, interface_map, 1e-3)
		if result["tetrahedra"] != 12 or result["interfaces"] != {
			"interface_seg_4_portal_tissue_candidate": 2} \
			or result["minimum_determinant_m3"] <= 0 \
			or result["minimum_scaled_jacobian"] < 1e-3:
			raise RuntimeError("two-voxel multiregion audit failed")
		second = Path(directory)/"repeat.msh"
		write_gmsh(second, nodes, tets, boundaries, interfaces)
		if path.read_bytes() != second.read_bytes():
			raise RuntimeError("two-voxel mesh output is not deterministic")
	require_rejection(lambda: build_mesh(regions, roi, (0, 0, 0), (1, 1, 1), 1),
		"allowed 1..1")
	overlapped = vessel.copy(); overlapped[0, 0, 0] = True
	require_rejection(lambda: build_mesh({"tissue_candidate": tissue,
		"seg_4_portal": overlapped}, roi, (0, 0, 0), (1, 1, 1), 3),
		"overlap")
	disconnected = vessel.copy(); disconnected[0, 0, 1] = False
	require_rejection(lambda: build_mesh({"tissue_candidate": tissue,
		"seg_4_portal": disconnected}, roi, (0, 0, 0), (1, 1, 1), 2),
		"no tissue-vessel")
	pinched_tissue = np.zeros((2, 2, 3), dtype=bool)
	pinched_vessel = np.zeros_like(pinched_tissue)
	pinched_tissue[0, 0, 0] = True
	pinched_tissue[1, 1, 1] = True
	pinched_vessel[0, 0, 1] = True
	require_rejection(lambda: build_mesh({"tissue_candidate": pinched_tissue,
		"seg_4_portal": pinched_vessel}, ((0, 2), (0, 2), (0, 3)),
		(0, 0, 0), (1, 1, 1), 3), "disconnected 6-neighbour")
	require_rejection(lambda: build_mesh({"tissue_candidate": pinched_tissue,
		"seg_4_portal": pinched_vessel}, ((0, 2), (0, 2), (0, 3)),
		(0, 0, 0), (1, 1, 1), 3, require_connected=False), "nonmanifold boundary")
	print("dicom_seg_multiregion_tet_test: PASS")


if __name__ == "__main__":
	main()
