"""Double-precision surface state tests for the ParaView dual-tree export."""

from pathlib import Path
import sys
import unittest

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from export_idealized_cube_paraview import deformation_volume, surface_at_state


class ParaViewSurfaceTest(unittest.TestCase):
	def test_actual_and_magnified_coordinates(self):
		ref = np.array([[0., 0., 0.], [.01, 0., 0.],
			[0., .01, 0.], [0., 0., .01]], dtype=np.float64)
		disp = np.array([[0., 0., 0.], [2e-8, 0., 0.],
			[0., -1e-8, 0.], [0., 0., 1e-8]], dtype=np.float64)
		grid = vtk.vtkUnstructuredGrid()
		points = vtk.vtkPoints()
		points.SetData(numpy_to_vtk(ref + disp, deep=True))
		grid.SetPoints(points)
		cell = vtk.vtkTetra()
		for index in range(4):
			cell.GetPointIds().SetId(index, index)
		grid.InsertNextCell(cell.GetCellType(), cell.GetPointIds())
		for name, values in (("reference_position_m", ref),
			("displacement_m", disp)):
			array = numpy_to_vtk(values, deep=True)
			array.SetName(name)
			grid.GetPointData().AddArray(array)
		identifiers = numpy_to_vtk(np.arange(4, dtype=np.int64), deep=True)
		identifiers.SetName("GlobalPointIds")
		grid.GetPointData().AddArray(identifiers)
		for initial, scale in ((True, 1.), (False, 1.), (False, 30000.)):
			surface = surface_at_state(grid, scale, initial, "synthetic")
			position = vtk_to_numpy(surface.GetPoints().GetData())
			surface_ref = vtk_to_numpy(surface.GetPointData().GetArray(
				"reference_position_m"))
			surface_disp = vtk_to_numpy(surface.GetPointData().GetArray(
				"displacement_m"))
			self.assertEqual(position.dtype, np.float64)
			np.testing.assert_allclose(position,
				surface_ref + scale*surface_disp, atol=1e-15, rtol=0)
			if initial:
				self.assertEqual(float(np.linalg.norm(surface_disp, axis=1).max()), 0.)
			else:
				self.assertGreater(float(np.linalg.norm(surface_disp, axis=1).max()), 0.)
		deformed, report = deformation_volume(grid, "synthetic")
		self.assertGreater(report["volume_change_m3"], 0.)
		self.assertGreater(report["minimum_deformation_jacobian"], 1.)
		self.assertIsNotNone(deformed.GetCellData().GetArray(
			"relative_cell_volume_change"))

	def test_rejects_display_warp_with_inverted_tetrahedron(self):
		ref = np.array([[0., 0., 0.], [.01, 0., 0.],
			[0., .01, 0.], [0., 0., .01]], dtype=np.float64)
		disp = np.zeros_like(ref)
		disp[1, 0] = -2e-5
		grid = vtk.vtkUnstructuredGrid()
		points = vtk.vtkPoints()
		points.SetData(numpy_to_vtk(ref+disp, deep=True))
		grid.SetPoints(points)
		cell = vtk.vtkTetra()
		for index in range(4):
			cell.GetPointIds().SetId(index, index)
		grid.InsertNextCell(cell.GetCellType(), cell.GetPointIds())
		for name, values in (("reference_position_m", ref),
			("displacement_m", disp)):
			array = numpy_to_vtk(values, deep=True)
			array.SetName(name)
			grid.GetPointData().AddArray(array)
		identifiers = numpy_to_vtk(np.arange(4, dtype=np.int64), deep=True)
		identifiers.SetName("GlobalPointIds")
		grid.GetPointData().AddArray(identifiers)
		with self.assertRaisesRegex(ValueError, "nonpositive tetrahedron Jacobian"):
			surface_at_state(grid, 1000., False, "inverted")


if __name__ == "__main__":
	unittest.main()
