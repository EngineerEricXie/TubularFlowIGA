#!/usr/bin/env python3
"""Negative regressions for the frozen T2 fixed-flow validation card."""

import copy
import json
from pathlib import Path
import subprocess
import tempfile


def main():
	root = Path(__file__).resolve().parents[2]
	validator = root/"scripts"/"validate_t2_fixed_flow_contract.py"
	baseline = json.loads((root/"benchmarks"/"t2_fixed_flow_contract.json").read_text(
		encoding="utf-8"))
	with tempfile.TemporaryDirectory(prefix="tubularflow-t2-contract-") as directory:
		temporary = Path(directory)
		healthy = temporary/"healthy.json"
		healthy.write_text(json.dumps(baseline), encoding="utf-8")
		result = subprocess.run(["python3", str(validator), "--contract", str(healthy)],
			text=True, capture_output=True)
		assert result.returncode == 0, result.stderr
		mutations = []
		wrong_pressure = copy.deepcopy(baseline)
		wrong_pressure["analytic_reference"]["pressure_drop_pa"] *= 2.0
		mutations.append(wrong_pressure)
		wrong_sign = copy.deepcopy(baseline)
		wrong_sign["boundary_conditions"]["expected_inlet_outward_flow_m3_s"] *= -1.0
		mutations.append(wrong_sign)
		loose_mass = copy.deepcopy(baseline)
		loose_mass["qoi_gates"]["maximum_relative_mass_imbalance"] = 1e-4
		mutations.append(loose_mass)
		missing_evidence = copy.deepcopy(baseline)
		missing_evidence["required_result_fields"].remove("geometry_error")
		mutations.append(missing_evidence)
		external_backend = copy.deepcopy(baseline)
		external_backend["backend_decision"]["architecture"] = "external_adapter"
		external_backend["backend_decision"]["external_discretization_frameworks_allowed"] = True
		mutations.append(external_backend)
		for index, mutation in enumerate(mutations):
			path = temporary/f"invalid-{index}.json"
			path.write_text(json.dumps(mutation), encoding="utf-8")
			result = subprocess.run(["python3", str(validator), "--contract", str(path)],
				text=True, capture_output=True)
			assert result.returncode != 0, f"invalid mutation {index} passed"
			assert "T2 fixed-flow contract error:" in result.stderr
	print("t2_fixed_flow_contract_test: PASS")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
