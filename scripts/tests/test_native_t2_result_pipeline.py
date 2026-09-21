#!/usr/bin/env python3
"""Test functional/validation separation in the native T2 result pipeline."""

import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from t2_contract import exact_contract_sha256, mesh_contract_sha256


FAKE_SOLVER = r'''#!/usr/bin/env python3
import json
import os
from pathlib import Path
import sys

mesh, contract_path = Path(sys.argv[1]), Path(sys.argv[2])
args = sys.argv[3:]
def value(option):
	return args[args.index(option)+1]
contract = json.loads(contract_path.read_text())
validation = value("-native_tet_validation") == "true"
index = int(mesh.stem.split("-")[-1])-1
flow = contract["boundary_conditions"]["prescribed_flow_m3_s"]
result = {
	"schema_version": 1,
	"result_classification": "physical_validation_candidate" if validation else "functional_smoke",
	"validation_gates_enforced": validation,
	"nonlinear_absolute_tolerance": 1e-14 if validation else 1e-8,
	"backend": "native_cpp_petsc_tetrahedral_fem",
	"backend_version": "test-double",
	"velocity_space": "continuous tetrahedral Lagrange P2 vector",
	"pressure_space": "continuous tetrahedral Lagrange P1 scalar",
	"stabilization": "none (inf-sup-stable Taylor-Hood)",
	"nonlinear_convergence_reason": "converged_test_double",
	"linear_convergence_reason": "converged_test_double",
	"linear_solver": "test-double",
	"mpi_ranks": int(os.environ.get("FAKE_MPI_RANKS", "1")),
	"inlet_outward_flow_m3_s": -flow,
	"outlet_outward_flow_m3_s": flow,
	"pressure_drop_pa": contract["analytic_reference"]["pressure_drop_pa"],
	"velocity_relative_l2": 0.016/(4**index),
	"wall_shear_relative_l2": 0.04/(2**index),
	"outlet_backflow_area_fraction": 0.0,
	"relative_mass_imbalance": 0.0,
	"surface_volume_divergence_relative_error": 0.0,
	"assembly_wall_s": 0.1,
	"solve_wall_s": 0.2,
	"peak_rss_bytes": 1024,
}
Path(value("-native_tet_level_result")).write_text(json.dumps(result))
'''

FAKE_MPIEXEC = r'''#!/usr/bin/env python3
import os
import subprocess
import sys

index = sys.argv.index("-np")
environment = os.environ.copy()
environment["FAKE_MPI_RANKS"] = sys.argv[index+1]
raise SystemExit(subprocess.run(sys.argv[index+2:], env=environment).returncode)
'''

SLEEPER = r'''#!/usr/bin/env python3
import time
time.sleep(10)
'''


def run(command, expected=0):
	completed = subprocess.run(command, text=True, capture_output=True)
	assert completed.returncode == expected, completed.stdout+completed.stderr
	return completed


def main():
	root = Path(__file__).resolve().parents[2]
	contract_path = root/"benchmarks"/"t2_fixed_flow_contract.json"
	contract_bytes = contract_path.read_bytes()
	contract = json.loads(contract_bytes)
	runner = root/"scripts"/"run_native_t2_level.py"
	collector = root/"scripts"/"collect_native_t2_results.py"
	validator = root/"scripts"/"validate_t2_fixed_flow_result.py"
	with tempfile.TemporaryDirectory(prefix="native-t2-pipeline-") as directory:
		temporary = Path(directory)
		fake = temporary/"fake-solver.py"
		fake.write_text(FAKE_SOLVER, encoding="utf-8")
		fake.chmod(0o755)
		fake_mpiexec = temporary/"fake-mpiexec.py"
		fake_mpiexec.write_text(FAKE_MPIEXEC, encoding="utf-8")
		fake_mpiexec.chmod(0o755)
		sleeper = temporary/"sleeper.py"
		sleeper.write_text(SLEEPER, encoding="utf-8")
		sleeper.chmod(0o755)
		levels = []
		for index, target in enumerate(
				contract["discretization_contract"]["mesh_levels_target_size_m"], 1):
			mesh = temporary/f"mesh-{index}.msh"
			mesh.write_text(f"synthetic mesh {index}\n", encoding="utf-8")
			mesh_hash = hashlib.sha256(mesh.read_bytes()).hexdigest()
			manifest = temporary/f"mesh-{index}.json"
			manifest.write_text(json.dumps({"volume_mesh": {"sha256": mesh_hash}}),
				encoding="utf-8")
			levels.append({"name": f"level-{index}", "target_size_m": target,
				"volume_mesh": mesh.name, "volume_manifest": manifest.name,
				"mesh_sha256": mesh_hash, "geometry_error": {
					"maximum_boundary_distance_m": target**2,
					"relative_volume_error": target,
					"boundary_label_counts_match": True}})
		series = temporary/"mesh-series.json"
		series.write_text(json.dumps({"schema_version": 1,
			"kind": "t2_mesh_preparation_not_solver_evidence",
			"case_id": contract["case_id"],
			"contract_sha256": exact_contract_sha256(contract_bytes),
			"mesh_contract_sha256": mesh_contract_sha256(contract), "levels": levels}),
			encoding="utf-8")
		mpi_output = temporary/"functional-mpi.json"
		run(["python3", str(runner), str(series), "1", str(mpi_output),
			"--contract", str(contract_path), "--executable", str(fake),
			"--mode", "functional", "--mpi-ranks", "2", "--mpiexec", str(fake_mpiexec),
			"--mpi-launcher-option=--test-launcher-option"])
		mpi_data = json.loads(mpi_output.read_text())
		assert mpi_data["requested_mpi_ranks"] == 2
		assert mpi_data["mpi_launcher_options"] == ["--test-launcher-option"]
		timed_out = run(["python3", str(runner), str(series), "1",
			str(temporary/"must-not-exist.json"), "--contract", str(contract_path),
			"--executable", str(sleeper), "--timeout-s", "0.05"], expected=2)
		assert "exceeded 0.05 s; no result was promoted" in timed_out.stdout
		assert not (temporary/"must-not-exist.json").exists()
		for mode in ("functional", "validation"):
			outputs = []
			for index in range(1, 4):
				output = temporary/f"{mode}-{index}.json"
				run(["python3", str(runner), str(series), str(index), str(output),
					"--contract", str(contract_path), "--executable", str(fake),
					"--mode", mode])
				outputs.append(output)
			combined = temporary/f"{mode}.json"
			run(["python3", str(collector), str(combined), *map(str, outputs),
				"--contract", str(contract_path)])
			validated = subprocess.run(["python3", str(validator), str(combined),
				"--contract", str(contract_path)], text=True, capture_output=True)
			if mode == "functional":
				assert validated.returncode != 0
				assert "functional smoke evidence" in validated.stderr
			else:
				assert validated.returncode == 0, validated.stderr
	print("native_t2_result_pipeline_test: PASS")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
