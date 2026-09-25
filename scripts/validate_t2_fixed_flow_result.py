#!/usr/bin/env python3
"""Validate a three-level T2 solver result against the frozen contract."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys


HASH = re.compile(r"[0-9a-f]{64}")
COMMIT = re.compile(r"(?:[0-9a-f]{40}|[0-9a-f]{64})")


def fail(message):
	print(f"T2 fixed-flow result error: {message}", file=sys.stderr)
	return 1


def finite_number(value, name, nonnegative=False):
	if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
		raise ValueError(f"{name} must be a finite number")
	if nonnegative and value < 0.0:
		raise ValueError(f"{name} must be nonnegative")
	return float(value)


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("result")
	parser.add_argument("--contract", default="benchmarks/t2_fixed_flow_contract.json")
	args = parser.parse_args()
	try:
		contract_path = Path(args.contract)
		contract_bytes = contract_path.read_bytes()
		contract = json.loads(contract_bytes)
		result = json.loads(Path(args.result).read_text(encoding="utf-8"))
	except (OSError, json.JSONDecodeError) as error:
		return fail(f"cannot read input: {error}")
	try:
		if result.get("schema_version") != 1 or result.get("case_id") != contract["case_id"]:
			raise ValueError("schema or case id does not match the contract")
		if result.get("result_classification") != "physical_validation":
			raise ValueError("result is functional smoke evidence, not physical validation")
		expected_contract_hash = hashlib.sha256(contract_bytes).hexdigest()
		if result.get("contract_sha256") != expected_contract_hash:
			raise ValueError("contract SHA-256 does not match the exact validated bytes")
		levels = result.get("levels")
		targets = contract["discretization_contract"]["mesh_levels_target_size_m"]
		if not isinstance(levels, list) or len(levels) != len(targets):
			raise ValueError("result must contain every frozen mesh level exactly once")
		gates = contract["qoi_gates"]
		flow = contract["boundary_conditions"]["prescribed_flow_m3_s"]
		reference_pressure = contract["analytic_reference"]["pressure_drop_pa"]
		velocity_errors = []
		backend_identity = None
		for index, (level, target) in enumerate(zip(levels, targets)):
			name = f"levels[{index}]"
			if not isinstance(level, dict):
				raise ValueError(f"{name} must be an object")
			if level.get("result_classification") != "physical_validation_candidate" \
					or level.get("validation_gates_enforced") is not True:
				raise ValueError(f"{name} was not produced with validation gates enabled")
			for field in contract["required_result_fields"]:
				if field not in level:
					raise ValueError(f"{name} is missing required field {field!r}")
			if finite_number(level.get("target_size_m"), f"{name}.target_size_m") != target:
				raise ValueError(f"{name} target size differs from the frozen series")
			identity = tuple(level.get(field) for field in
				("backend", "backend_version", "velocity_space", "pressure_space", "stabilization"))
			if not all(isinstance(value, str) and value.strip() for value in identity):
				raise ValueError(f"{name} backend and discretization fields must be nonempty")
			if backend_identity is None:
				backend_identity = identity
			elif identity != backend_identity:
				raise ValueError("backend or discrete formulation changed across mesh levels")
			if not isinstance(level["source_commit"], str) \
					or not COMMIT.fullmatch(level["source_commit"]):
				raise ValueError(f"{name}.source_commit must be a full lowercase Git object id")
			if level["source_tree_state"] not in ("clean", "dirty_with_source_manifest"):
				raise ValueError(f"{name}.source_tree_state is invalid")
			for field in ("source_files_sha256", "input_sha256", "mesh_sha256"):
				if not isinstance(level[field], str) or not HASH.fullmatch(level[field]):
					raise ValueError(f"{name}.{field} must be a lowercase SHA-256")
			for field in ("nonlinear_convergence_reason", "linear_convergence_reason"):
				if not isinstance(level[field], str) or not level[field].startswith("converged"):
					raise ValueError(f"{name}.{field} is not a successful convergence reason")
			if not isinstance(level["linear_solver"], str) or not level["linear_solver"].strip():
				raise ValueError(f"{name}.linear_solver must be explicit")
			inlet = finite_number(level["inlet_outward_flow_m3_s"], f"{name}.inlet flow")
			outlet = finite_number(level["outlet_outward_flow_m3_s"], f"{name}.outlet flow")
			if inlet >= 0.0 or outlet <= 0.0:
				raise ValueError(f"{name} violates outward-positive port signs")
			mass = abs(inlet+outlet)/max(abs(inlet), abs(outlet), flow)
			reported_mass = finite_number(level["relative_mass_imbalance"], f"{name}.mass", True)
			if not math.isclose(mass, reported_mass, rel_tol=1e-10, abs_tol=1e-15):
				raise ValueError(f"{name} reported mass imbalance is inconsistent with port flows")
			pressure = finite_number(level["pressure_drop_pa"], f"{name}.pressure drop")
			pressure_error = abs(pressure-reference_pressure)/reference_pressure
			velocity_error = finite_number(level["velocity_relative_l2"], f"{name}.velocity error", True)
			wall_error = finite_number(level["wall_shear_relative_l2"], f"{name}.WSS error", True)
			backflow = finite_number(level["outlet_backflow_area_fraction"],
				f"{name}.outlet backflow", True)
			divergence = finite_number(level["surface_volume_divergence_relative_error"],
				f"{name}.divergence error", True)
			if reported_mass > gates["maximum_relative_mass_imbalance"] \
					or divergence > gates["maximum_surface_volume_divergence_relative_error"] \
					or backflow > contract["boundary_conditions"][
						"maximum_outlet_backflow_area_fraction"]:
				raise ValueError(f"{name} fails a conservation or backflow gate")
			if index == len(levels)-1 and (velocity_error > gates["maximum_velocity_relative_l2"]
					or pressure_error > gates["maximum_pressure_drop_relative_error"]
					or wall_error > gates["maximum_wall_shear_relative_l2"]):
				raise ValueError(f"{name} fails a finest-level absolute QoI gate")
			geometry = level["geometry_error"]
			if not isinstance(geometry, dict) or geometry.get("boundary_label_counts_match") is not True:
				raise ValueError(f"{name}.geometry_error must confirm boundary labels")
			finite_number(geometry.get("maximum_boundary_distance_m"),
				f"{name}.maximum boundary distance", True)
			finite_number(geometry.get("relative_volume_error"), f"{name}.volume error", True)
			for field in ("assembly_wall_s", "solve_wall_s", "peak_rss_bytes"):
				finite_number(level[field], f"{name}.{field}", True)
			velocity_errors.append(velocity_error)
		orders = []
		for index in range(len(levels)-1):
			if velocity_errors[index+1] <= 0.0 or velocity_errors[index] <= 0.0:
				raise ValueError("velocity convergence order requires positive errors")
			orders.append(math.log(velocity_errors[index]/velocity_errors[index+1])
				/math.log(targets[index]/targets[index+1]))
		observed = finite_number(result.get("minimum_observed_velocity_convergence_order"),
			"minimum observed convergence order")
		if not math.isclose(observed, min(orders), rel_tol=1e-10, abs_tol=1e-12):
			raise ValueError("reported convergence order is inconsistent with level errors")
		if observed < gates["minimum_observed_velocity_convergence_order"]:
			raise ValueError("velocity convergence order is below the frozen gate")
	except (KeyError, TypeError, ValueError) as error:
		return fail(str(error))
	print(f"T2 fixed-flow result: PASS ({len(levels)} levels, minimum order {observed:.6g})")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
