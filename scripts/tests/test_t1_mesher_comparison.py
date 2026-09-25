#!/usr/bin/env python3
"""Check healthy and fail-closed T1 mesher comparison validation."""

import copy
import json
from pathlib import Path
import subprocess
import tempfile


def main():
	root = Path(__file__).resolve().parents[2]
	report = root/"benchmarks"/"t1_mesher_comparison.json"
	validator = root/"scripts"/"validate_t1_mesher_comparison.py"
	healthy = subprocess.run(["python3",str(validator),str(report)], text=True,
		capture_output=True)
	assert healthy.returncode == 0, healthy.stdout+healthy.stderr
	data = json.loads(report.read_text(encoding="utf-8"))
	mutations = []
	bad_quality = copy.deepcopy(data)
	bad_quality["valid_cases"][0]["results"]["gmsh"]["minimum_scaled_jacobian"] = 0.0
	mutations.append(bad_quality)
	bad_labels = copy.deepcopy(data)
	del bad_labels["valid_cases"][2]["results"]["ftetwild"]["boundary_labels"]["3"]
	mutations.append(bad_labels)
	bad_rejection = copy.deepcopy(data)
	bad_rejection["invalid_cases"][0]["ftetwild_rejected"] = False
	mutations.append(bad_rejection)
	with tempfile.TemporaryDirectory(prefix="t1-mesher-comparison-") as directory:
		for index, mutation in enumerate(mutations):
			path = Path(directory)/f"invalid-{index}.json"
			path.write_text(json.dumps(mutation), encoding="utf-8")
			failed = subprocess.run(["python3",str(validator),str(path)], text=True,
				capture_output=True)
			assert failed.returncode != 0, f"invalid comparison mutation {index} passed"
	print("t1_mesher_comparison_test: PASS")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
