#!/usr/bin/env python3
"""Small synthetic checks for voxel-face contact arithmetic, not anatomy."""

import unittest

import numpy as np

from scripts.audit_dicom_seg_region_overlap import (
	component_bounding_boxes_zyx, component_face_contacts, face_contact,
	parse_roi_zyx, roi_component_membership)


class FaceContactTest(unittest.TestCase):
	def test_one_face_per_axis_with_physical_area(self):
		for axis, expected_area in ((0, 1 * 2), (1, 1 * 3), (2, 2 * 3)):
			with self.subTest(axis=axis):
				tissue = np.zeros((3, 3, 3), dtype=bool)
				vessel = np.zeros_like(tissue)
				center = [1, 1, 1]
				neighbour = center.copy()
				neighbour[axis] += 1
				tissue[tuple(center)] = True
				vessel[tuple(neighbour)] = True
				result = face_contact(tissue, vessel, (1, 2, 3))
				self.assertEqual(result["total_contact_faces"], 1)
				self.assertEqual(result["total_contact_area_mm2"], expected_area)
				self.assertEqual(result["faces_normal_z_y_x"][axis], 1)

	def test_diagonal_is_not_face_contact(self):
		tissue = np.zeros((3, 3, 3), dtype=bool)
		vessel = np.zeros_like(tissue)
		tissue[1, 1, 1] = True
		vessel[2, 2, 1] = True
		self.assertEqual(face_contact(tissue, vessel, (1, 2, 3))["total_contact_faces"], 0)

	def test_overlap_and_bad_grid_rejected(self):
		tissue = np.zeros((2, 2, 2), dtype=bool)
		vessel = np.zeros_like(tissue)
		tissue[0, 0, 0] = True
		vessel[0, 0, 0] = True
		with self.assertRaisesRegex(ValueError, "disjoint"):
			face_contact(tissue, vessel, (1, 2, 3))
		vessel[0, 0, 0] = False
		with self.assertRaisesRegex(ValueError, "positive"):
			face_contact(tissue, vessel, (1, 0, 3))
		with self.assertRaisesRegex(ValueError, "share"):
			face_contact(tissue, vessel[0], (1, 2, 3))

	def test_component_face_contacts_distinguish_disconnected_island(self):
		tissue = np.zeros((3, 3, 5), dtype=bool)
		components = np.zeros_like(tissue, dtype=np.int32)
		tissue[1, 1, 1] = True
		components[1, 1, 2] = 1
		components[1, 1, 4] = 2
		self.assertEqual(component_face_contacts(tissue, components, 2).tolist(),
			[0, 1, 0])
		with self.assertRaisesRegex(ValueError, "overlap"):
			component_face_contacts(components > 0, components, 2)

	def test_component_boxes_are_half_open_and_report_image_edge(self):
		labels = np.zeros((4, 5, 6), dtype=np.int32)
		labels[1:3, 2:4, 2:4] = 1
		labels[0, 0, 5] = 2
		self.assertEqual(component_bounding_boxes_zyx(labels, 2), [
			{"component_id": 1,
				"bbox_zyx_half_open": [[1, 3], [2, 4], [2, 4]],
				"touches_seg_grid_extent": False},
			{"component_id": 2,
				"bbox_zyx_half_open": [[0, 1], [0, 1], [5, 6]],
				"touches_seg_grid_extent": True}])
		with self.assertRaisesRegex(ValueError, "missing"):
			component_bounding_boxes_zyx(labels, 3)

	def test_roi_membership_uses_full_grid_component_ids(self):
		mask = np.zeros((4, 5, 6), dtype=bool)
		mask[1:3, 1, 1] = True
		mask[2, 3, 4] = True
		roi = parse_roi_zyx("2:3,0:5,0:6", mask.shape)
		self.assertEqual(roi_component_membership(mask, roi), [
			{"component_id": 1, "roi_voxels": 1},
			{"component_id": 2, "roi_voxels": 1}])
		self.assertEqual(roi_component_membership(mask,
			parse_roi_zyx("0:1,0:5,0:6", mask.shape)), [])
		with self.assertRaisesRegex(ValueError, "exceeds"):
			parse_roi_zyx("0:5,0:5,0:6", mask.shape)


if __name__ == "__main__":
	unittest.main()
