#!/usr/bin/env python3
"""Fail closed on the checked-in native T2 functional-smoke evidence."""

import argparse
import hashlib
import json
import math
from pathlib import Path

from run_native_t2_level import SOURCE_FILES, combined_sha256


def require(condition, message):
	if not condition:
		raise RuntimeError(message)


def main():
	root = Path(__file__).resolve().parents[1]
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("evidence", nargs="?",
		default=str(root/"benchmarks"/"t2_native_functional_evidence.json"))
	parser.add_argument("--skip-source-hash", action="store_true",
		help="test-only: validate a mutated copy after source binding was already checked")
	args = parser.parse_args()
	evidence = json.loads(Path(args.evidence).read_text(encoding="utf-8"))
	require(evidence.get("schema_version") == 1, "unsupported schema version")
	require(evidence.get("kind") == "native_fem_functional_smoke_not_physical_validation",
		"evidence must remain explicitly non-validation")
	require(evidence.get("backend") == "native_cpp_petsc_tetrahedral_fem",
		"unexpected backend")
	require(evidence.get("external_fem_framework") is False,
		"external FEM framework cannot be native evidence")
	policy = evidence.get("solver_policy", {})
	require(policy.get("result_classification") == "functional_smoke",
		"functional evidence classification changed")
	require(policy.get("validation_gates_enforced") is False,
		"functional evidence cannot claim validation gates")
	contract = root/"benchmarks"/"t2_fixed_flow_contract.json"
	require(hashlib.sha256(contract.read_bytes()).hexdigest() == evidence.get("contract_sha256"),
		"contract hash is stale")
	if not args.skip_source_hash:
		sources = [(name, root/name) for name in SOURCE_FILES]
		current_hash=combined_sha256(sources)
		require(current_hash == evidence.get("current_compatible_source_files_sha256"),
			"native FEM compatible source hash is stale")
		require(evidence.get("current_source_compatibility") ==
			"affine body acceleration defaults identically to zero; zero-force element and prior temporal regressions pass",
			"native FEM source compatibility is not explicitly bounded")
	require(isinstance(evidence.get("executable_sha256"), str)
		and len(evidence["executable_sha256"]) == 64,
		"native FEM executable SHA-256 is missing")
	levels = evidence.get("passed_levels")
	require(isinstance(levels, list) and [x.get("mesh_level") for x in levels] == [1, 2],
		"exactly coarse and medium functional levels must be recorded")
	for level in levels:
		for key in ("target_size_m", "assembly_wall_s", "solve_wall_s", "total_wall_s",
				"velocity_relative_l2", "wall_shear_relative_l2",
				"relative_mass_imbalance", "surface_volume_divergence_relative_error"):
			value = level.get(key)
			require(isinstance(value, (int, float)) and math.isfinite(value) and value >= 0,
				f"invalid {key} at level {level.get('mesh_level')}")
		require(level.get("total_wall_s") >= level.get("assembly_wall_s")
			+level.get("solve_wall_s"), "total time is smaller than phase times")
		require(level.get("outlet_backflow_area_fraction") == 0,
			"functional observation has outlet backflow")
	require(levels[1]["target_size_m"] < levels[0]["target_size_m"],
		"mesh size did not decrease")
	require(levels[1]["tetrahedra"] > levels[0]["tetrahedra"],
		"tetrahedron count did not increase")
	require(levels[1]["velocity_relative_l2"] < levels[0]["velocity_relative_l2"],
		"velocity error did not decrease")
	require(levels[1]["wall_shear_relative_l2"] < levels[0]["wall_shear_relative_l2"],
		"wall-shear error did not decrease")
	failures = evidence.get("failed_attempts")
	require(isinstance(failures, list) and {x.get("mpi_ranks") for x in failures} == {8, 16},
		"8-rank and 16-rank fine timeout observations are required")
	for fine in failures:
		require(fine.get("mesh_level") == 3 and fine.get("failure") == "wall_time_timeout"
			and fine.get("timeout_s") == 600 and fine.get("result_promoted") is False,
			"fine timeout must remain failed and unpromoted")
	for experiment in evidence.get("rejected_solver_experiments", []):
		require(experiment.get("accepted") is False,
			"a rejected solver experiment was silently promoted")
	print("validate_t2_native_functional_evidence: PASS")
	return 0


if __name__ == "__main__":
	try:
		raise SystemExit(main())
	except (KeyError, OSError, RuntimeError, TypeError, ValueError) as error:
		print(f"validate_t2_native_functional_evidence: ERROR: {error}")
		raise SystemExit(2)
