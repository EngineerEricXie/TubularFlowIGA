#!/usr/bin/env python3
"""Collect three native T2 level results without promoting smoke data to validation."""

import argparse
import hashlib
import json
import math
from pathlib import Path


def main():
	root = Path(__file__).resolve().parents[1]
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("output")
	parser.add_argument("levels", nargs=3)
	parser.add_argument("--contract", default=str(root/"benchmarks"/"t2_fixed_flow_contract.json"))
	args = parser.parse_args()
	contract_path = Path(args.contract).resolve()
	contract_bytes = contract_path.read_bytes()
	contract = json.loads(contract_bytes)
	levels = [json.loads(Path(path).read_text(encoding="utf-8")) for path in args.levels]
	targets = contract["discretization_contract"]["mesh_levels_target_size_m"]
	for index, (level, target) in enumerate(zip(levels, targets), 1):
		if level.get("mesh_level") != index or level.get("target_size_m") != target:
			raise RuntimeError(f"level result {index} does not match frozen ordering")
		if level.get("contract_sha256") != hashlib.sha256(contract_bytes).hexdigest():
			raise RuntimeError(f"level result {index} contract SHA-256 is stale")
	errors = [float(level["velocity_relative_l2"]) for level in levels]
	if any(not math.isfinite(error) or error <= 0.0 for error in errors):
		raise RuntimeError("velocity errors must be finite and positive")
	orders = [math.log(errors[index]/errors[index+1])
		/math.log(targets[index]/targets[index+1]) for index in range(2)]
	validated = all(level.get("result_classification") == "physical_validation_candidate"
		and level.get("validation_gates_enforced") is True for level in levels)
	result = {
		"schema_version": 1,
		"case_id": contract["case_id"],
		"contract_sha256": hashlib.sha256(contract_bytes).hexdigest(),
		"result_classification": "physical_validation" if validated else "functional_smoke",
		"levels": levels,
		"minimum_observed_velocity_convergence_order": min(orders),
	}
	output = Path(args.output).resolve()
	output.parent.mkdir(parents=True, exist_ok=True)
	output.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n", encoding="utf-8")
	print(f"collect_native_t2_results: PASS classification={result['result_classification']} "
		f"minimum_order={min(orders):.6g}")
	return 0


if __name__ == "__main__":
	try:
		raise SystemExit(main())
	except (KeyError, OSError, RuntimeError, TypeError, ValueError) as error:
		print(f"collect_native_t2_results: ERROR: {error}")
		raise SystemExit(2)
