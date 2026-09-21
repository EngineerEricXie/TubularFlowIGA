#!/usr/bin/env python3
"""Validate the checked-in T1 same-input mesher comparison evidence."""

import argparse
import json
import math
from pathlib import Path
import re
import sys


HASH = re.compile(r"[0-9a-f]{64}")


def fail(message):
	print(f"T1 mesher comparison error: {message}", file=sys.stderr)
	return 1


def positive(value, name, allow_zero=False):
	if isinstance(value, bool) or not isinstance(value, (int,float)) \
			or not math.isfinite(value) or value < 0.0 or (not allow_zero and value == 0.0):
		raise ValueError(f"{name} must be {'nonnegative' if allow_zero else 'positive'}")
	return float(value)


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("report", nargs="?", default="benchmarks/t1_mesher_comparison.json")
	args = parser.parse_args()
	try:
		data = json.loads(Path(args.report).read_text(encoding="utf-8"))
		if data.get("schema_version") != 1 \
				or data.get("kind") != "local_wsl_mesher_comparison_evidence":
			raise ValueError("schema or evidence kind changed")
		policy = data["policy"]
		minimum_quality = positive(policy["required_minimum_scaled_jacobian"], "quality gate")
		maximum_volume = positive(policy["maximum_relative_volume_error"], "volume gate")
		cases = data["valid_cases"]
		if [case.get("id") for case in cases] != ["straight_pipe", "bent_pipe_90_degree",
				"y_bifurcation"]:
			raise ValueError("valid comparison case set changed")
		for case in cases:
			name = case["id"]
			for field in ("source_sha256", "canonical_surface_sha256"):
				if not isinstance(case.get(field), str) or not HASH.fullmatch(case[field]):
					raise ValueError(f"{name}.{field} is not a SHA-256")
			envelope = positive(case["ftetwild_envelope_m"], f"{name}.envelope")
			results = case.get("results")
			if not isinstance(results, dict) or set(results) != {"gmsh","ftetwild"}:
				raise ValueError(f"{name} must compare exactly Gmsh and fTetWild")
			label_sets = []
			for mesher, result in results.items():
				context = f"{name}.{mesher}"
				positive(result["nodes"], context+".nodes")
				positive(result["tetrahedra"], context+".tetrahedra")
				if positive(result["minimum_scaled_jacobian"], context+".quality") < minimum_quality:
					raise ValueError(f"{context} fails the quality gate")
				if positive(result["relative_volume_error"], context+".volume", True) > maximum_volume:
					raise ValueError(f"{context} fails the volume gate")
				distance = positive(result["maximum_source_surface_distance_m"],
					context+".surface distance", True)
				if mesher == "ftetwild" and distance > envelope:
					raise ValueError(f"{context} exceeds its envelope")
				for field in ("mesher_wall_s","pipeline_wall_s","peak_rss_bytes"):
					positive(result[field], context+"."+field)
				if not HASH.fullmatch(result.get("mesh_sha256", "")):
					raise ValueError(f"{context}.mesh_sha256 is invalid")
				if result.get("native_reader_passed") is not True:
					raise ValueError(f"{context} did not pass the native reader")
				labels = result.get("boundary_labels")
				if not isinstance(labels, dict) or any(not str(key).isdigit()
						or not isinstance(value,int) or value <= 0 for key,value in labels.items()):
					raise ValueError(f"{context} boundary labels are invalid")
				label_sets.append(set(labels))
			if label_sets[0] != label_sets[1]:
				raise ValueError(f"{name} meshers did not preserve the same labels")
		invalid = data.get("invalid_cases")
		if not isinstance(invalid,list) or len(invalid) != 1 \
				or invalid[0].get("gmsh_rejected") is not True \
				or invalid[0].get("ftetwild_rejected") is not True \
				or invalid[0].get("expected_failure_class") != "shared_surface_preflight_not_closed":
			raise ValueError("shared imperfect-surface rejection evidence is missing")
	except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError) as error:
		return fail(str(error))
	print(f"T1 mesher comparison: PASS ({len(cases)} valid, {len(invalid)} invalid cases)")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
