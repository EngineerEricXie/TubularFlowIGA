#!/usr/bin/env python3
"""Exercise the surface-label GUI state changes without a display server."""

import json
from pathlib import Path
import subprocess
import sys
import tempfile

import vtk

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from scripts.label_surface_bc_gui import SurfaceLabeler, boundary_loops, read_surface


def tetra_surface():
	points = vtk.vtkPoints()
	for point in ((0., 0., 0.), (1., 0., 0.), (0., 1., 0.), (0., 0., 1.)):
		points.InsertNextPoint(*point)
	polys = vtk.vtkCellArray()
	for nodes in ((0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3)):
		cell = vtk.vtkTriangle()
		for local, node in enumerate(nodes):
			cell.GetPointIds().SetId(local, node)
		polys.InsertNextCell(cell)
	surface = vtk.vtkPolyData()
	surface.SetPoints(points)
	surface.SetPolys(polys)
	labels = vtk.vtkIntArray()
	labels.SetName("boundary_id")
	labels.SetNumberOfTuples(4)
	labels.Fill(0)
	surface.GetCellData().AddArray(labels)
	surface.GetCellData().SetScalars(labels)
	return surface


def main():
	with tempfile.TemporaryDirectory(prefix="surface-label-gui-") as directory:
		output = Path(directory)/"labelled.vtp"
		labeler = SurfaceLabeler(tetra_surface(), output, interactive=False)

		labeler.select_cells([0])
		preview = labeler.surface.GetCellData().GetArray("selection_preview")
		assert tuple(preview.GetTuple3(0)) == (255., 215., 0.)
		labeler.assign(1)
		labels = labeler.surface.GetCellData().GetArray("boundary_id")
		assert labels.GetValue(0) == 1

		labeler.select_cells([1, 2])
		labeler.select_cells([2, 3], "toggle")
		assert labeler.selected == {1, 3}
		labeler.assign(2)
		labeler.undo()
		labels = labeler.surface.GetCellData().GetArray("boundary_id")
		assert [labels.GetValue(i) for i in range(4)] == [1, 0, 0, 0]

		labeler.select_cells([3])
		labeler.cut_selected()
		assert labeler.surface.GetNumberOfCells() == 3
		assert [len(loop) for loop in boundary_loops(labeler.surface)] == [3]
		labeler.cap(0, 3)
		assert labeler.surface.GetNumberOfCells() == 6
		assert boundary_loops(labeler.surface) == []
		labeler.save()

		saved = read_surface(output)
		assert saved.GetNumberOfCells() == 6
		metadata = json.loads(output.with_suffix(".labels.json").read_text())
		assert metadata["wall_label"] == 0
		assert metadata["inflow_label"] == 1
		assert metadata["outflow_labels"] == [3]
		assert sum(metadata["boundary_triangles_by_label"].values()) == 6
		preflight = subprocess.run([str(ROOT/"solvers/cpu/surface_fem_preflight"),
			str(output), str(Path(directory)/"surface.msh")], text=True,
			capture_output=True)
		assert preflight.returncode == 0, preflight.stdout+preflight.stderr
	print("label_surface_bc_gui_test: PASS")


if __name__ == "__main__":
	main()
