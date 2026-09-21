"""Artificial grids only; do not interpret exposed faces as anatomical ports."""

import unittest

import numpy as np

from scripts.audit_seg_exposed_vessel_faces import count_exposed_faces


class ExposedFacesTest(unittest.TestCase):
	def test_tissue_background_and_grid_are_distinct(self):
		components = np.zeros((3, 3, 3), dtype=np.int32)
		components[0, 1, 1] = 1
		tissue = np.zeros_like(components, dtype=bool)
		tissue[1, 1, 1] = True
		vessel = components > 0
		item, = count_exposed_faces(components, tissue, vessel, (1, 2, 3))
		self.assertEqual(item["tissue"]["total_faces"], 1)
		self.assertEqual(item["background"]["total_faces"], 4)
		self.assertEqual(item["other_vessel"]["total_faces"], 0)
		self.assertEqual(item["grid_extent"]["total_faces"], 1)
		self.assertEqual(item["tissue"]["area_mm2"], 2)
		self.assertEqual(item["grid_extent"]["area_mm2"], 2)

	def test_invalid_overlap_and_component_ids_rejected(self):
		components = np.zeros((2, 2, 2), dtype=np.int32)
		components[0, 0, 0] = 2
		vessel = components > 0
		with self.assertRaisesRegex(ValueError, "contiguous"):
			count_exposed_faces(components, np.zeros_like(vessel), vessel, (1, 1, 1))
		components[0, 0, 0] = 1
		with self.assertRaisesRegex(ValueError, "disjoint"):
			count_exposed_faces(components, vessel, vessel, (1, 1, 1))

	def test_other_vessel_is_not_background_or_tissue(self):
		components = np.zeros((2, 2, 2), dtype=np.int32)
		components[0, 0, 0] = 1
		vessel_union = components > 0
		vessel_union[0, 0, 1] = True
		item, = count_exposed_faces(components,
			np.zeros_like(vessel_union), vessel_union, (1, 1, 1))
		self.assertEqual(item["other_vessel"]["total_faces"], 1)
		self.assertEqual(item["background"]["total_faces"], 2)
		self.assertEqual(item["grid_extent"]["total_faces"], 3)


if __name__ == "__main__":
	unittest.main()
