#!/usr/bin/env python3
"""Mutation checks for native T2 functional evidence classification."""

import json
from pathlib import Path
import subprocess
import tempfile


def run(command, expected):
	completed = subprocess.run(command, text=True, capture_output=True)
	assert completed.returncode == expected, completed.stdout+completed.stderr
	return completed


def main():
	root = Path(__file__).resolve().parents[2]
	validator = root/"scripts"/"validate_t2_native_functional_evidence.py"
	evidence_path = root/"benchmarks"/"t2_native_functional_evidence.json"
	run(["python3", str(validator), str(evidence_path)], 0)
	baseline = json.loads(evidence_path.read_text(encoding="utf-8"))
	with tempfile.TemporaryDirectory(prefix="t2-functional-evidence-") as directory:
		for name, mutate, expected_text in (
			("promoted", lambda x: x["failed_attempts"][0].update(result_promoted=True),
				"fine timeout must remain failed and unpromoted"),
			("validation", lambda x: x["solver_policy"].update(
				result_classification="physical_validation"),
				"functional evidence classification changed"),
		):
			candidate = json.loads(json.dumps(baseline))
			mutate(candidate)
			path = Path(directory)/f"{name}.json"
			path.write_text(json.dumps(candidate), encoding="utf-8")
			failed = run(["python3", str(validator), str(path), "--skip-source-hash"], 2)
			assert expected_text in failed.stdout
	print("t2_native_functional_evidence_test: PASS")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
