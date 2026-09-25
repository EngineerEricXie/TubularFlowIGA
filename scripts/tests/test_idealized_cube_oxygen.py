"""Artificial three-domain tracer input must state its restricted assumptions."""

import copy
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_idealized_cube_dual_tree_oxygen import validate


CASE = Path(__file__).resolve().parents[2]/"cases/idealized_cube_dual_tree_oxygen.json"


class ArtificialOxygenCaseTest(unittest.TestCase):
	def setUp(self):
		self.case = json.loads(CASE.read_text(encoding="utf-8"))

	def test_valid_artificial_contract(self):
		self.assertTrue(validate(self.case, CASE).is_file())
		self.assertLess(self.case["vascular_diffusivity_m2_s"],
			self.case["tissue_diffusivity_m2_s"])

	def test_rejects_claim_of_physiology(self):
		case = copy.deepcopy(self.case)
		case["physiological_validation"] = True
		with self.assertRaisesRegex(ValueError, "functional-only"):
			validate(case, CASE)

	def test_rejects_unmodelled_tissue_consumption(self):
		case = copy.deepcopy(self.case)
		case["tissue_consumption_mol_m3_s"] = 1e-3
		with self.assertRaisesRegex(ValueError, "functional-only"):
			validate(case, CASE)

	def test_rejects_invalid_time_or_breakthrough(self):
		for key, value in (("time_step_s", 0), ("maximum_steps", 1),
			("venous_breakthrough_fraction", 1.0)):
			with self.subTest(key=key):
				case = copy.deepcopy(self.case)
				case[key] = value
				with self.assertRaisesRegex(ValueError, "contract"):
					validate(case, CASE)

	def test_rejects_unsplit_diffusivity(self):
		case = copy.deepcopy(self.case)
		case["vascular_diffusivity_m2_s"] = case["tissue_diffusivity_m2_s"]
		with self.assertRaisesRegex(ValueError, "contract"):
			validate(case, CASE)


if __name__ == "__main__":
	unittest.main()
