"""Fail-closed contract tests for the artificial centerline/radius case."""

import copy
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from generate_idealized_cube_dual_tree import validate


CASE = Path(__file__).resolve().parents[2] / "cases/idealized_cube_dual_tree.json"
REFINED_CASE = Path(__file__).resolve().parents[2] / \
	"cases/idealized_cube_dual_tree_refined.json"


class IdealizedCubeDualTreeCaseTest(unittest.TestCase):
	def setUp(self):
		self.case = json.loads(CASE.read_text(encoding="utf-8"))

	def test_valid_centerlines_and_radii(self):
		self.assertEqual(validate(copy.deepcopy(self.case)), self.case)

	def test_refined_thin_wall_and_curvature_contract(self):
		refined = json.loads(REFINED_CASE.read_text(encoding="utf-8"))
		self.assertEqual(validate(refined), refined)
		self.assertLess(refined["wall_thickness_m"],
			self.case["wall_thickness_m"])
		self.assertGreater(refined["mesh_curvature_divisions"], 8)

	def test_reject_invalid_refined_mesh_controls(self):
		refined = json.loads(REFINED_CASE.read_text(encoding="utf-8"))
		refined["minimum_mesh_size_m"] = refined["target_mesh_size_m"] * 2
		with self.assertRaisesRegex(ValueError, "mesh controls|minimum size"):
			validate(refined)

	def test_reject_overlapping_trees(self):
		case = copy.deepcopy(self.case)
		for node in case["venous_tree"]["nodes_m"]:
			node[2] -= 0.008
		with self.assertRaisesRegex(ValueError, "clearance"):
			validate(case)

	def test_reject_wrong_units(self):
		case = copy.deepcopy(self.case)
		case["units"] = "mm"
		with self.assertRaisesRegex(ValueError, "SI"):
			validate(case)

	def test_reject_branch_outside_cube(self):
		case = copy.deepcopy(self.case)
		case["arterial_tree"]["nodes_m"][4][1] = -0.02
		with self.assertRaisesRegex(ValueError, "leaves cube"):
			validate(case)

	def test_reject_nonpositive_wall_stiffness(self):
		case = copy.deepcopy(self.case)
		case["functional_parameters"]["wall_young_modulus_pa"] = 0
		with self.assertRaisesRegex(ValueError, "wall_young_modulus_pa"):
			validate(case)

	def test_reject_duplicate_terminal(self):
		case = copy.deepcopy(self.case)
		case["arterial_tree"]["terminals"][3] = case["arterial_tree"]["terminals"][0]
		with self.assertRaisesRegex(ValueError, "unique terminals"):
			validate(case)


if __name__ == "__main__":
	unittest.main()
