#!/usr/bin/env python3
"""Exercise healthy and fail-closed T2 result validation."""

import copy
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile


def main():
	root = Path(__file__).resolve().parents[2]
	contract_path = root/"benchmarks"/"t2_fixed_flow_contract.json"
	contract_bytes = contract_path.read_bytes()
	contract = json.loads(contract_bytes)
	validator = root/"scripts"/"validate_t2_fixed_flow_result.py"
	pressure = contract["analytic_reference"]["pressure_drop_pa"]
	levels = []
	for index, target in enumerate(contract["discretization_contract"]["mesh_levels_target_size_m"]):
		level = {field: 0.0 for field in contract["required_result_fields"]}
		level.update({
			"result_classification": "physical_validation_candidate",
			"validation_gates_enforced": True,
			"target_size_m": target,
			"backend": "test_backend", "backend_version": "1",
			"source_commit": "1"*40, "source_tree_state": "dirty_with_source_manifest",
			"source_files_sha256": "9"*64, "input_sha256": "2"*64,
			"mesh_sha256": f"{index+3}"*64,
			"velocity_space": "P2", "pressure_space": "P1", "stabilization": "none",
			"nonlinear_convergence_reason": "converged_exact_fixture",
			"linear_convergence_reason": "converged_exact_fixture",
			"linear_solver": "synthetic_exact",
			"inlet_outward_flow_m3_s": -1e-6, "outlet_outward_flow_m3_s": 1e-6,
			"pressure_drop_pa": pressure, "velocity_relative_l2": 0.016/(4**index),
			"wall_shear_relative_l2": 0.04/(2**index), "relative_mass_imbalance": 0.0,
			"outlet_backflow_area_fraction": 0.0,
			"surface_volume_divergence_relative_error": 1e-8,
			"geometry_error": {"maximum_boundary_distance_m": target**2,
				"relative_volume_error": target, "boundary_label_counts_match": True},
			"assembly_wall_s": 1.0, "solve_wall_s": 2.0, "peak_rss_bytes": 1024})
		levels.append(level)
	result = {"schema_version": 1, "case_id": contract["case_id"],
		"result_classification": "physical_validation",
		"contract_sha256": hashlib.sha256(contract_bytes).hexdigest(), "levels": levels,
		"minimum_observed_velocity_convergence_order": 2.0}
	with tempfile.TemporaryDirectory(prefix="tubularflow-t2-result-") as directory:
		temporary = Path(directory)
		def run(value, name):
			path = temporary/name
			path.write_text(json.dumps(value), encoding="utf-8")
			return subprocess.run(["python3", str(validator), str(path), "--contract",
				str(contract_path)], text=True, capture_output=True)
		healthy = run(result, "healthy.json")
		assert healthy.returncode == 0, healthy.stderr
		mutations = []
		bad_hash = copy.deepcopy(result); bad_hash["contract_sha256"] = "0"*64; mutations.append(bad_hash)
		bad_sign = copy.deepcopy(result); bad_sign["levels"][0]["inlet_outward_flow_m3_s"] = 1e-6; mutations.append(bad_sign)
		bad_mass = copy.deepcopy(result); bad_mass["levels"][0]["outlet_outward_flow_m3_s"] = .9e-6; mutations.append(bad_mass)
		bad_order = copy.deepcopy(result); bad_order["levels"][2]["velocity_relative_l2"] = .01; mutations.append(bad_order)
		bad_labels = copy.deepcopy(result); bad_labels["levels"][1]["geometry_error"]["boundary_label_counts_match"] = False; mutations.append(bad_labels)
		functional = copy.deepcopy(result); functional["result_classification"] = "functional_smoke"; mutations.append(functional)
		ungated = copy.deepcopy(result); ungated["levels"][2]["validation_gates_enforced"] = False; mutations.append(ungated)
		for index, mutation in enumerate(mutations):
			failed = run(mutation, f"invalid-{index}.json")
			assert failed.returncode != 0, f"invalid result mutation {index} passed"
			assert "T2 fixed-flow result error:" in failed.stderr
	print("t2_fixed_flow_result_test: PASS")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
