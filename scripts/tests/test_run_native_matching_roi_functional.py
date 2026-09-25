#!/usr/bin/env python3
"""No-MPI parser and rejection checks for the reproducible ROI runner."""

import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_native_matching_roi_functional import parse_result


class FunctionalResultTest(unittest.TestCase):
	def test_complete_converged_line(self):
		line = ("native_tet_matching_solved_roi_smoke: PASS facets=158 "
			"fluid_inlet_m3_s=-1.5e-8 interface_source_m3_s=1.5e-8 "
			"darcy_outward_m3_s=1.5e-8 fluid_newton_iterations=4 "
			"fluid_linear_reason=4 darcy_linear_reason=3 "
			"darcy_max_cell_defect_m3_s=1e-21\n")
		self.assertEqual(parse_result(line)["facets"], 158)

	def test_missing_duplicate_nonfinite_and_unconverged_rejected(self):
		valid = ("native_tet_matching_solved_roi_smoke: PASS facets=2 "
			"fluid_inlet_m3_s=-1 interface_source_m3_s=1 "
			"darcy_outward_m3_s=1 fluid_newton_iterations=2 "
			"fluid_linear_reason=4 darcy_linear_reason=3 "
			"darcy_max_cell_defect_m3_s=0\n")
		for malformed in (valid.replace(" facets=2", ""),
				valid.replace("facets=2", "facets=2 facets=2"),
				valid.replace("facets=2", "facets=nan"),
				valid.replace("darcy_linear_reason=3", "darcy_linear_reason=-1"),
				valid+valid):
			with self.subTest(malformed=malformed):
				with self.assertRaises(ValueError):
					parse_result(malformed)


if __name__ == "__main__":
	unittest.main()
