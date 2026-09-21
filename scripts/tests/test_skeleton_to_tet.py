"""Tests for the template-free skeleton-to-tetrahedron input path."""

from pathlib import Path
import tempfile
import unittest

from preprocessing.tet.skeleton_geometry import (branch_transition_radius,
	fit_terminal_caps, spline_capsule_paths)
from preprocessing.tet.skeleton_to_tet import read_skeleton, simplify_tree


class SkeletonToTetTest(unittest.TestCase):
	def test_swc_trifurcation_and_radius_transition(self):
		with tempfile.TemporaryDirectory() as directory:
			path = Path(directory)/"tree.swc"
			path.write_text("""1 2 0 0 0 1.2 -1
2 2 1 0 0 1.2 1
3 2 2 -1 0 0.8 2
4 2 2 0 0 0.7 2
5 2 2 1 0 0.6 2
""", encoding="utf-8")
			tree = read_skeleton(path, .001)
		self.assertEqual(tree["root"], 0)
		self.assertEqual(tree["terminals"], [2, 3, 4])
		self.assertEqual(len(spline_capsule_paths(tree, 2, .25)), 4)
		self.assertAlmostEqual(tree["node_radii_m"][1], .0012)
		self.assertAlmostEqual(branch_transition_radius(.0012, .0008, .5, 1.), .001)

	def test_centerline_simplification_preserves_topology_nodes(self):
		tree = {"nodes_m": [[0., 0., 0.], [.1, 0., 0.], [.2, 0., 0.],
			[1., 1., 0.], [1., -1., 0.]],
			"node_radii_m": [1., 1., 1., .5, .5],
			"segments": [[0, 1, 1.], [1, 2, 1.], [2, 3, .5], [2, 4, .5]],
			"root": 0, "terminals": [3, 4]}
		simplified = simplify_tree(tree, .15)
		self.assertEqual(simplified["nodes_m"],
			[[0., 0., 0.], [.2, 0., 0.], [1., 1., 0.], [1., -1., 0.]])
		self.assertEqual(simplified["terminals"], [2, 3])

	def test_obj_root_is_largest_terminal_radius(self):
		with tempfile.TemporaryDirectory() as directory:
			path = Path(directory)/"tree.obj"
			path.write_text("""v 0 0 0 0.4 0 0
v 1 0 0 0.5 0 0
v 2 1 0 0.3 0 0
v 2 -1 0 0.2 0 0
l 1 2
l 2 3
l 2 4
""", encoding="utf-8")
			tree = read_skeleton(path, 1.)
		self.assertEqual(tree["nodes_m"][0], [0., 0., 0.])
		self.assertEqual(tree["terminals"], [2, 3])

	def test_terminal_caps_follow_endpoint_radii(self):
		points = {1: (0., 0., 0.), 2: (0., 1., 0.), 3: (0., 0., 1.),
			4: (2., 0., 0.), 5: (2., .5, 0.), 6: (2., 0., .5)}
		triangles = {"inlet": [(1, 2, 3)], "outlet_0": [(4, 5, 6)]}
		tree = {"nodes_m": [[0., 0., 0.], [2., 0., 0.]],
			"node_radii_m": [2., 1.], "root": 0, "terminals": [1]}
		fit_terminal_caps(points, triangles, tree)
		self.assertAlmostEqual(points[2][1], 2.)
		self.assertAlmostEqual(points[5][1], 1.)


if __name__ == "__main__":
	unittest.main()
