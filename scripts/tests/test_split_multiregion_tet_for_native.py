#!/usr/bin/env python3
"""End-to-end synthetic multi-region split and native reader/interface gates."""

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"scripts"))
from dicom_seg_to_multiregion_tet import build_mesh, make_contract, write_gmsh


def run(command):
	return subprocess.run(command, cwd=ROOT, text=True,
		stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def main():
	tissue = np.zeros((1, 1, 2), dtype=bool)
	vessel = np.zeros_like(tissue)
	tissue[0, 0, 0] = True
	vessel[0, 0, 1] = True
	with tempfile.TemporaryDirectory(prefix="tubularflow-native-split-") as directory:
		folder = Path(directory)
		mesh = folder/"source.msh"
		nodes, cells, boundaries, interfaces, _ = build_mesh(
			{"tissue_candidate": tissue, "seg_4_portal": vessel},
			((0, 1), (0, 1), (0, 2)), (0, 0, 0), (1, 1, 1), 2)
		surfaces = write_gmsh(mesh, nodes, cells, boundaries, interfaces)
		contract = folder/"source.contract.json"
		contract.write_text(json.dumps(make_contract(cells, boundaries,
			interfaces, surfaces)), encoding="utf-8")
		script = ROOT/"scripts"/"split_multiregion_tet_for_native.py"
		first = folder/"first"
		result = run([sys.executable, str(script), str(mesh), str(contract), str(first)])
		if result.returncode:
			raise RuntimeError(result.stdout+result.stderr)
		manifest = json.loads((first/"manifest.json").read_text(encoding="utf-8"))
		label = manifest["interface_boundary_labels"][
			"interface_seg_4_portal_tissue_candidate"]
		if label != 3 or set(manifest["submeshes"]) != {
			"tissue_candidate", "seg_4_portal"}:
			raise RuntimeError("split manifest changed")
		for region in ("tissue_candidate", "seg_4_portal"):
			check = run([str(ROOT/"solvers"/"cpu"/"native_tet_mesh_inspect"),
				str(first/f"{region}.msh")])
			if check.returncode or "tetrahedra=6" not in check.stdout:
				raise RuntimeError(check.stdout+check.stderr)
		matching = run([str(ROOT/"solvers"/"cpu"/"native_tet_matching_mesh_inspect"),
			str(first/"seg_4_portal.msh"), str(first/"tissue_candidate.msh"), str(label)])
		if matching.returncode or "facets=2" not in matching.stdout:
			raise RuntimeError(matching.stdout+matching.stderr)
		wrong = run([str(ROOT/"solvers"/"cpu"/"native_tet_matching_mesh_inspect"),
			str(first/"seg_4_portal.msh"), str(first/"tissue_candidate.msh"), "999"])
		if wrong.returncode != 2 or "label is absent" not in wrong.stderr:
			raise RuntimeError("wrong interface label was not rejected")
		second = folder/"second"
		result = run([sys.executable, str(script), str(mesh), str(contract), str(second)])
		if result.returncode:
			raise RuntimeError(result.stdout+result.stderr)
		for region in ("tissue_candidate", "seg_4_portal"):
			if hashlib.sha256((first/f"{region}.msh").read_bytes()).digest() != \
					hashlib.sha256((second/f"{region}.msh").read_bytes()).digest():
				raise RuntimeError("split mesh is not deterministic")
		again = run([sys.executable, str(script), str(mesh), str(contract), str(first)])
		if again.returncode != 2 or "already contains files" not in again.stderr:
			raise RuntimeError("existing outputs were not protected")
		overlimit = run([sys.executable, str(script), str(mesh), str(contract),
			str(folder/"overlimit"), "--max-tetrahedra", "1"])
		if overlimit.returncode != 2 or "allowed 1..1" not in overlimit.stderr:
			raise RuntimeError("oversized mesh was not rejected")
	print("split_multiregion_tet_for_native_test: PASS")


if __name__ == "__main__":
	main()
