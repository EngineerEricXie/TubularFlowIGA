"""Fail-closed checks for the versioned liver ROI functional case."""

import json
from pathlib import Path
import tempfile
import unittest

from scripts.run_liver_roi_functional_case import validate_case


CASE = Path(__file__).resolve().parents[2]/"cases/liver_roi_functional.json"


class LiverRoiFunctionalCaseTest(unittest.TestCase):
	def test_missing_source_data_is_not_a_valid_case_run(self):
		with tempfile.TemporaryDirectory() as directory:
			with self.assertRaisesRegex(ValueError, "source_seg file is missing"):
				validate_case(CASE, Path(directory))

	def test_physiological_claim_is_rejected_before_data_access(self):
		case = json.loads(CASE.read_text(encoding="utf-8"))
		case["physiological_validation"] = True
		with tempfile.TemporaryDirectory() as directory:
			path = Path(directory)/"case.json"
			path.write_text(json.dumps(case), encoding="utf-8")
			with self.assertRaisesRegex(ValueError, "classification is invalid"):
				validate_case(path, Path(directory))

	def test_traversal_is_rejected_before_data_access(self):
		case = json.loads(CASE.read_text(encoding="utf-8"))
		case["files_relative_to_data_dir"]["source_seg"]["path"] = "../seg.dcm"
		with tempfile.TemporaryDirectory() as directory:
			path = Path(directory)/"case.json"
			path.write_text(json.dumps(case), encoding="utf-8")
			with self.assertRaisesRegex(ValueError, "path escapes"):
				validate_case(path, Path(directory))


if __name__ == "__main__":
	unittest.main()
