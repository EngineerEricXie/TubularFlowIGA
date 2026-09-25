#!/usr/bin/env python3
"""Validate the frozen native Kirchhoff-Love shell implementation card."""

import json
from pathlib import Path


def require(condition, message):
	if not condition:
		raise RuntimeError(message)


def main():
	root = Path(__file__).resolve().parents[1]
	path = root/"benchmarks"/"t3_native_iga_shell_contract.json"
	data = json.loads(path.read_text(encoding="utf-8"))
	require(data.get("schema_version") == 1, "unsupported schema")
	require(data.get("contract_status") == "frozen", "T3 contract is not frozen")
	require(data.get("implementation_status") == "verified_single_patch_reference_runtime",
		"T3 implementation status is stale")
	require(data.get("backend") == "native_cpp_iga_kirchhoff_love_shell",
		"first shell backend must remain native")
	require(data.get("external_discretization_or_assembly_framework") is False,
		"external assembly cannot satisfy native shell completion")
	theory = data["theory"]
	require(theory.get("name") == "Kirchhoff-Love thin shell", "shell theory changed")
	require(theory.get("transverse_shear") == "neglected", "KL shear assumption is missing")
	require(0 < theory.get("thickness_ratio_recommended_maximum", 0) <= 0.05,
		"thin-shell applicability bound is missing")
	require(theory.get("within_patch_continuity") == "C1_or_higher_across_element_boundaries",
		"KL continuity requirement changed")
	required = set(data["input_contract"]["required"])
	for name in ("tensor_product_B_spline_or_NURBS_midsurface",
			"degrees_p_and_q_at_least_2", "positive_thickness",
			"Young_modulus_and_Poisson_ratio", "explicit_boundary_constraints_and_loads"):
		require(name in required, f"missing required input {name}")
	cards = {card["id"]: card for card in data.get("verification_cards", [])}
	for name in ("rigid-motion-invariance", "flat-membrane-patch",
			"pure-bending-patch", "consistent-tangent", "spatial-convergence"):
		require(name in cards, f"missing verification card {name}")
	require(cards["spatial-convergence"].get("minimum_levels", 0) >= 3,
		"spatial convergence requires at least three levels")
	policy = data["classification_policy"]
	require(policy.get("rendered_deformation_only") == "demonstration_visualization",
		"visualization cannot be numerical verification")
	require(policy.get("external_solver_comparison") == "cross_validation_not_native_completion",
		"external solver cannot be relabelled native completion")
	print("T3 native IGA shell contract: PASS")
	return 0


if __name__ == "__main__":
	try:
		raise SystemExit(main())
	except (KeyError, OSError, RuntimeError, TypeError, ValueError) as error:
		print(f"T3 native IGA shell contract: ERROR: {error}")
		raise SystemExit(2)
