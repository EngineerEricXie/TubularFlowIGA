#!/usr/bin/env python3
"""Check shared-entry routing for editable and directly imported solver cases."""

import json
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[2]
SOLVER = ROOT/"scripts/solver.py"


def run(configuration, *options):
	completed = subprocess.run([sys.executable, str(SOLVER), str(configuration),
		"--dry-run", *map(str, options)], cwd=ROOT, text=True,
		stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	assert completed.returncode == 0, completed.stdout+completed.stderr
	return completed.stdout


def main():
	count = 0
	for configuration in sorted((ROOT/"examples/solver").glob("*.json")):
		data = json.loads(configuration.read_text(encoding="utf-8"))
		if data.get("schema_version") == "solver-v1":
			run(configuration)
			count += 1
	flow_darcy = run(ROOT/"examples/solver/fem_flow_darcy_species_gpu.json")
	assert "native_tet_flow_cuda" in flow_darcy
	assert "native_tet_darcy_cuda" in flow_darcy and "--source-map" in flow_darcy
	assert "native_tet_species_cuda" in flow_darcy

	assert "iga_0d" in run(ROOT/"examples/zero_d/steady_resistive_straight/simulation_config.json")
	assert "iga_1d" in run(ROOT/"examples/one_d/rigid_straight/simulation_config.json")
	iga_case = ROOT/"examples/vascular_flow/straight_tube/simulation_config.json"
	assert "iga_navier_stokes" in run(iga_case,
		"--database", "/tmp/direct-import.ntiga")
	assert "iga_cuda navier-stokes" in run(iga_case,
		"--database", "/tmp/direct-import.ntiga", "--device", "GPU")
	assert "iga_solve" in run(iga_case,
		"--database", "/tmp/direct-import.ntiga", "--physics", "species")
	assert "iga_cuda transport" in run(iga_case,
		"--database", "/tmp/direct-import.ntiga", "--device", "GPU",
		"--physics", "species")
	assert "run_native_tet_workflow.py" in run(
		ROOT/"examples/solver/fem_hydraulic_workflow.json")
	assert "native_tet_flow_cuda" in run(
		ROOT/"examples/solver/fem_flow_gpu_case.json",
		"--backend", "native_tet_flow_cuda")
	assert "native_tet_darcy" in run(
		ROOT/"examples/solver/fem_darcy_case.json",
		"--backend", "native_tet_darcy")
	assert "native_tet_species_transport" in run(
		ROOT/"examples/solver/fem_species_case.json",
		"--backend", "native_tet_species_transport")

	with tempfile.TemporaryDirectory(prefix="solver-direct-import-") as directory:
		directory = Path(directory)
		fsi = directory/"fsi.json"
		fsi.write_text(json.dumps({"schema_version": 1,
			"fluid_mesh_file": "fluid.msh", "solid_mesh_file": "wall.msh",
			"time": {"dt_s": 0.01, "steps": 2}}))
		assert "native_tet_fsi_steps" in run(fsi)

		chain = directory/"chain.json"
		chain.write_text(json.dumps({"schema_version": 1,
			"meshes": {"artery": "a.msh", "tissue": "t.msh", "vein": "v.msh"}}))
		assert "native_tet_flow_darcy_chain" in run(chain)

	print(f"solver_entry_test: PASS editable_configs={count} direct_import_routes=12")


if __name__ == "__main__":
	main()
