"""Fail-closed contract tests for the artificial centerline/radius case."""

import copy
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from generate_idealized_cube_dual_tree import (branch_transition_radius,
	spline_capsule_paths, tree_node_radii, validate)


CASE = Path(__file__).resolve().parents[2] / "cases/idealized_cube_dual_tree.json"
REFINED_CASE = Path(__file__).resolve().parents[2] / \
	"cases/idealized_cube_dual_tree_refined.json"
SMOOTH_CASE = Path(__file__).resolve().parents[2] / \
	"cases/idealized_cube_dual_tree_smooth.json"


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

	def test_template_free_spline_capsule_contract(self):
		smooth = json.loads(SMOOTH_CASE.read_text(encoding="utf-8"))
		self.assertEqual(validate(smooth), smooth)
		self.assertEqual(smooth["vascular_geometry"]["kind"], "spline_capsules")
		self.assertEqual(smooth["vascular_geometry"]["radius_transition_fraction"], .5)
		self.assertEqual(smooth["vascular_geometry"]["surface_smoothing_iterations"], 80)

	def test_spline_capsules_accept_a_trifurcation(self):
		tree = {
			"nodes_m": [[0., 0., 0.], [1., 0., 0.], [2., -1., 0.],
				[2., 0., 0.], [2., 1., 0.], [3., .5, -.5], [3., 1.5, .5]],
			"segments": [[0, 1, 1.], [1, 2, .8], [1, 3, .8], [1, 4, .8],
				[4, 5, .6], [4, 6, .6]],
			"root": 0, "terminals": [2, 3, 5, 6]}
		paths = spline_capsule_paths(tree, 2, .25)
		self.assertEqual(len(paths), len(tree["segments"]))
		self.assertEqual(sum(start == 1 for start, _, _, _ in paths), 3)
		self.assertTrue(all(len(points) == 3 for _, _, _, points in paths))
		radii = tree_node_radii(tree)
		self.assertEqual(radii[1], 1.)
		self.assertEqual(radii[4], .8)
		self.assertEqual(branch_transition_radius(1., .8, 0., 2.), 1.)
		self.assertEqual(branch_transition_radius(1., .8, 1., 2.), .9)
		self.assertEqual(branch_transition_radius(1., .8, 2., 2.), .8)

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
