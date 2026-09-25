#!/usr/bin/env python3
"""Validate the frozen native tetrahedral ALE card."""

import json
from pathlib import Path


def require(condition, message):
	if not condition:
		raise RuntimeError(message)


def main():
	root = Path(__file__).resolve().parents[1]
	data = json.loads((root/"benchmarks"/"t4_native_ale_contract.json").read_text())
	require(data.get("schema_version") == 1, "unsupported T4 schema")
	require(data.get("contract_status") == "frozen_before_runtime_implementation",
		"T4 contract status changed")
	require(data.get("implementation_status") ==
		"native_coarse_mesh_deforming_runtime_verified_chamber_comparison_pending",
		"T4 implementation status is inaccurate")
	require(data.get("backend") == "native_cpp_petsc_tetrahedral_ALE",
		"T4 backend must remain native")
	require(data.get("external_discretization_or_assembly_framework") is False,
		"external ALE assembly cannot count as native completion")
	kinematics = data["kinematics"]
	require(kinematics.get("current_coordinate") == "x(X,t)=X+d(X,t)",
		"reference-to-current map is missing")
	require(kinematics.get("mesh_velocity") == "w=(d_next-d_committed)/dt",
		"mesh-velocity time difference changed")
	require(kinematics.get("fluid_convective_velocity") == "u-w",
		"ALE relative convection is missing")
	require(kinematics.get("state_lifecycle") == ["committed", "trial", "accepted_or_rollback"],
		"ALE lifecycle is not transactional")
	gates = data["quality_gates"]
	require(gates.get("positive_current_tetrahedron_determinant") is True,
		"positive Jacobian gate is missing")
	require(0 < gates.get("minimum_current_to_reference_determinant_ratio", 0) < 1,
		"determinant-ratio gate is invalid")
	require(gates.get("quality_failure_action") == "rollback_then_reduce_step_or_stop",
		"quality failure must roll back or stop")
	require(gates.get("implicit_remeshing") is False,
		"remeshing cannot be implicit")
	boundary = data["boundary_policy"]
	require(boundary.get("every_surface_label_requires_an_explicit_rule") is True,
		"every ALE surface label must have an explicit rule")
	require(boundary.get("fixed_port") == "zero_mesh_displacement_and_velocity",
		"fixed-port mesh motion changed")
	require("fluid_velocity_is_a_separate_boundary_condition" in boundary.get("moving_wall", ""),
		"mesh velocity must not silently become fluid velocity")
	cards = set(data.get("verification_cards", []))
	for required in ("zero_motion_reduces_exactly_to_fixed_mesh_operator",
			"discrete_geometric_conservation_uniform_field",
			"inverted_trial_mesh_rolls_back_without_history_leak",
			"temporal_convergence_on_prescribed_motion"):
		require(required in cards, f"missing ALE verification card {required}")
	print("T4 native ALE contract: PASS")
	return 0


if __name__ == "__main__":
	try:
		raise SystemExit(main())
	except (KeyError, OSError, RuntimeError, TypeError, ValueError) as error:
		print(f"T4 native ALE contract: ERROR: {error}")
		raise SystemExit(2)
