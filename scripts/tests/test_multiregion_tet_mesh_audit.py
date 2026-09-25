#!/usr/bin/env python3
"""Positive and fail-closed tests for the multi-region tetrahedral audit."""

import json
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[2]
TOOL = ROOT/"scripts"/"audit_multiregion_tet_mesh.py"


def mesh():
	return """$MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
5
2 1 "outer_fluid"
2 2 "outer_solid"
2 3 "fluid_solid_interface"
3 11 "fluid_domain"
3 12 "solid_domain"
$EndPhysicalNames
$Entities
0 0 7 2
1 0 0 0 1 1 1 1 1 0
2 0 0 0 1 1 1 1 1 0
3 0 0 0 1 1 1 1 1 0
4 0 0 -1 1 1 0 1 2 0
5 0 0 -1 1 1 0 1 2 0
6 0 0 -1 1 1 0 1 2 0
7 0 0 0 1 1 0 1 3 0
1 0 0 0 1 1 1 1 11 0
2 0 0 -1 1 1 0 1 12 0
$EndEntities
$Nodes
1 5 1 5
3 1 0 5
1
2
3
4
5
0 0 0
1 0 0
0 1 0
0 0 1
0 0 -1
$EndNodes
$Elements
9 9 1 9
2 1 2 1
1 1 2 4
2 2 2 1
2 1 4 3
2 3 2 1
3 2 3 4
2 4 2 1
4 1 3 5
2 5 2 1
5 1 5 2
2 6 2 1
6 3 2 5
2 7 2 1
7 1 2 3
3 1 4 1
8 1 2 3 4
3 2 4 1
9 1 3 2 5
$EndElements
"""


def contract():
	return {
		"schema_version": 1,
		"length_unit": "m",
		"regions": [
			{"physical_name": "fluid_domain", "role": "fluid", "require_connected": True},
			{"physical_name": "solid_domain", "role": "solid", "require_connected": True},
		],
		"boundaries": [
			{"physical_name": "outer_fluid", "boundary_id": 10,
				"adjacent_region": "fluid_domain", "semantic_role": "explicit_fluid_boundary"},
			{"physical_name": "outer_solid", "boundary_id": 20,
				"adjacent_region": "solid_domain", "semantic_role": "explicit_solid_boundary"},
		],
		"interfaces": [
			{"physical_name": "fluid_solid_interface",
				"regions": ["fluid_domain", "solid_domain"]},
		],
	}


def run(mesh_path, contract_path, output=None):
	command = [sys.executable, str(TOOL), str(mesh_path), str(contract_path)]
	if output is not None:
		command += ["--output", str(output)]
	return subprocess.run(command, cwd=ROOT, text=True,
		stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def require_failure(result, needle):
	if result.returncode != 2 or needle not in result.stderr:
		raise RuntimeError(f"expected failure containing {needle!r}:\n{result.stdout}{result.stderr}")


def main():
	with tempfile.TemporaryDirectory(prefix="tubularflow-multiregion-") as directory:
		directory = Path(directory)
		mesh_path = directory/"two-region.msh"
		contract_path = directory/"contract.json"
		manifest_path = directory/"audit.json"
		mesh_path.write_text(mesh(), encoding="ascii")
		contract_path.write_text(json.dumps(contract()), encoding="utf-8")
		result = run(mesh_path, contract_path, manifest_path)
		if result.returncode != 0:
			raise RuntimeError(result.stdout+result.stderr)
		data = json.loads(manifest_path.read_text(encoding="utf-8"))
		if not data["gates"]["passed"] or data["audit"]["tetrahedra"] != 2 \
				or data["audit"]["interfaces"] != {"fluid_solid_interface": 1} \
				or data["audit"]["boundaries"] != {"outer_fluid": 3, "outer_solid": 3} \
				or data["audit"]["regions"]["fluid_domain"]["volume_m3"] != 1.0/6.0 \
				or data["audit"]["regions"]["solid_domain"]["volume_m3"] != 1.0/6.0:
			raise RuntimeError("healthy two-region audit manifest changed")

		wrong_interface = directory/"wrong-interface.msh"
		wrong_interface.write_text(mesh().replace(
			"7 0 0 0 1 1 0 1 3 0", "7 0 0 0 1 1 0 1 1 0"), encoding="ascii")
		require_failure(run(wrong_interface, contract_path),
			"cross-region face lacks exactly one declared interface")

		inverted = directory/"inverted.msh"
		inverted.write_text(mesh().replace("9 1 3 2 5", "9 1 2 3 5"), encoding="ascii")
		require_failure(run(inverted, contract_path), "nonpositive determinant")

		wrong_boundary_contract = contract()
		wrong_boundary_contract["boundaries"][0]["adjacent_region"] = "solid_domain"
		wrong_contract_path = directory/"wrong-contract.json"
		wrong_contract_path.write_text(json.dumps(wrong_boundary_contract), encoding="utf-8")
		require_failure(run(mesh_path, wrong_contract_path), "attached to the wrong region")

		binary = directory/"binary.msh"
		binary.write_text(mesh().replace("4.1 0 8", "4.1 1 8"), encoding="ascii")
		require_failure(run(binary, contract_path), "requires ASCII Gmsh 4.1")
	print("multiregion_tet_mesh_audit_test: PASS")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
