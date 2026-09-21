#!/usr/bin/env python3
"""Interactively label a triangulated vessel surface for the FEM mesh route."""

import argparse
import json
from pathlib import Path
import tkinter as tk
from tkinter import messagebox, simpledialog

import vtk


def triangles(polydata):
	return [tuple(polydata.GetCell(i).GetPointId(j) for j in range(3))
		for i in range(polydata.GetNumberOfCells())]


def boundary_loops(polydata):
	edges = {}
	for cell in triangles(polydata):
		for a, b in ((cell[0], cell[1]), (cell[1], cell[2]),
				(cell[2], cell[0])):
			key = tuple(sorted((a, b)))
			edges.setdefault(key, []).append((a, b))
	oriented = {items[0][0]: items[0][1] for items in edges.values()
		if len(items) == 1}
	loops = []
	while oriented:
		start = next(iter(oriented))
		loop = [start]
		point = start
		while point in oriented:
			point = oriented.pop(point)
			if point == start:
				break
			loop.append(point)
		loops.append(loop)
	return loops


def replace_geometry(polydata, cells, labels, extra_points=()):
	all_points = [polydata.GetPoint(index)
		for index in range(polydata.GetNumberOfPoints())] + list(extra_points)
	used = sorted({node for cell in cells for node in cell})
	local_map = {node: index for index, node in enumerate(used)}
	points = vtk.vtkPoints()
	points.SetDataTypeToDouble()
	for node in used:
		points.InsertNextPoint(*all_points[node])
	polys = vtk.vtkCellArray()
	for cell in cells:
		triangle = vtk.vtkTriangle()
		for local, node in enumerate(cell):
			triangle.GetPointIds().SetId(local, local_map[node])
		polys.InsertNextCell(triangle)
	result = vtk.vtkPolyData()
	result.SetPoints(points)
	result.SetPolys(polys)
	boundary = vtk.vtkIntArray()
	boundary.SetName("boundary_id")
	for label in labels:
		boundary.InsertNextValue(label)
	result.GetCellData().AddArray(boundary)
	result.GetCellData().SetScalars(boundary)
	return result


def remove_selected(polydata, selected):
	labels = polydata.GetCellData().GetArray("boundary_id")
	cells = triangles(polydata)
	keep = [index for index in range(len(cells)) if index not in selected]
	return replace_geometry(polydata, [cells[index] for index in keep],
		[labels.GetValue(index) for index in keep])


def cap_loop(polydata, loop, label):
	points = [polydata.GetPoint(node) for node in loop]
	center = tuple(sum(point[axis] for point in points)/len(points)
		for axis in range(3))
	centroid = polydata.GetNumberOfPoints()
	cells = triangles(polydata)
	labels = polydata.GetCellData().GetArray("boundary_id")
	new_cells = [(loop[(index+1) % len(loop)], loop[index], centroid)
		for index in range(len(loop))]
	return replace_geometry(polydata, cells+new_cells,
		[labels.GetValue(index) for index in range(len(cells))]
		+[label]*len(new_cells), [center])


def read_surface(path):
	readers = {".vtp": vtk.vtkXMLPolyDataReader,
		".stl": vtk.vtkSTLReader, ".ply": vtk.vtkPLYReader,
		".obj": vtk.vtkOBJReader}
	reader = readers[path.suffix.lower()]()
	reader.SetFileName(str(path))
	reader.Update()
	triangulate = vtk.vtkTriangleFilter()
	triangulate.SetInputConnection(reader.GetOutputPort())
	triangulate.Update()
	polydata = vtk.vtkPolyData()
	polydata.DeepCopy(triangulate.GetOutput())
	boundary = polydata.GetCellData().GetArray("boundary_id")
	if boundary is None:
		boundary = vtk.vtkIntArray()
		boundary.SetName("boundary_id")
		boundary.SetNumberOfTuples(polydata.GetNumberOfCells())
		boundary.Fill(0)
		polydata.GetCellData().AddArray(boundary)
	polydata.GetCellData().SetScalars(boundary)
	return polydata


def write_surface(path, polydata):
	points = " ".join(f"{value:.17g}" for index in range(polydata.GetNumberOfPoints())
		for value in polydata.GetPoint(index))
	cells = triangles(polydata)
	connectivity = " ".join(str(node) for cell in cells for node in cell)
	offsets = " ".join(str(3*(index+1)) for index in range(len(cells)))
	array = polydata.GetCellData().GetArray("boundary_id")
	labels = " ".join(str(array.GetValue(index)) for index in range(len(cells)))
	path.write_text(f'''<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">
  <PolyData>
    <Piece NumberOfPoints="{polydata.GetNumberOfPoints()}" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="{len(cells)}">
      <Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">{points}</DataArray></Points>
      <Polys>
        <DataArray type="Int64" Name="connectivity" format="ascii">{connectivity}</DataArray>
        <DataArray type="Int64" Name="offsets" format="ascii">{offsets}</DataArray>
      </Polys>
      <CellData Scalars="boundary_id">
        <DataArray type="Int32" Name="boundary_id" format="ascii">{labels}</DataArray>
      </CellData>
    </Piece>
  </PolyData>
</VTKFile>
''')


class SurfaceLabeler(vtk.vtkInteractorStyleTrackballCamera):
	def __init__(self, surface, output, interactive=True):
		super().__init__()
		self.surface = surface
		self.output = output
		self.interactive = interactive
		self.selected = set()
		self.history = []
		self.outlet_label = 2
		self.drag_start = None
		self.box = False
		if not interactive:
			self.dialog = None
			self.renderer = None
			self.mapper = None
			self.actor = None
			self.window = None
			self.interactor = None
			self.refresh()
			return
		self.dialog = tk.Tk()
		self.dialog.withdraw()
		self.renderer = vtk.vtkRenderer()
		self.renderer.SetBackground(0.08, 0.10, 0.15)
		self.mapper = vtk.vtkPolyDataMapper()
		self.actor = vtk.vtkActor()
		self.actor.SetMapper(self.mapper)
		self.renderer.AddActor(self.actor)
		self.actor.GetProperty().EdgeVisibilityOn()
		self.actor.GetProperty().SetEdgeColor(0.15, 0.17, 0.2)
		self.window = vtk.vtkRenderWindow()
		self.window.AddRenderer(self.renderer)
		self.window.SetSize(1100, 800)
		self.interactor = vtk.vtkRenderWindowInteractor()
		self.interactor.SetRenderWindow(self.window)
		self.interactor.SetInteractorStyle(self)
		self.AddObserver("LeftButtonPressEvent", self.press)
		self.AddObserver("LeftButtonReleaseEvent", self.release)
		self.AddObserver("KeyPressEvent", self.key)
		self.refresh()
		self.renderer.ResetCamera()

	def refresh(self):
		labels = self.surface.GetCellData().GetArray("boundary_id")
		colors = vtk.vtkUnsignedCharArray()
		colors.SetName("selection_preview")
		colors.SetNumberOfComponents(3)
		for index in range(self.surface.GetNumberOfCells()):
			label = labels.GetValue(index)
			if index in self.selected:
				color = (255, 215, 0)
			elif label == 0:
				color = (150, 158, 170)
			elif label == 1:
				color = (230, 72, 64)
			else:
				color = ((41*label) % 180+50, (97*label) % 180+50,
					(173*label) % 180+50)
			colors.InsertNextTuple3(*color)
		self.surface.GetCellData().AddArray(colors)
		self.surface.GetCellData().SetActiveScalars("selection_preview")
		if not self.interactive:
			return
		self.mapper.SetInputData(self.surface)
		self.mapper.SetScalarModeToUseCellData()
		self.mapper.SetColorModeToDirectScalars()
		self.mapper.ScalarVisibilityOn()
		self.window.SetWindowName(
			f"Surface BC labels | selected {len(self.selected)} | outlet {self.outlet_label} "
			"| i inlet, o outlet, w wall, n next outlet, u undo, x cut, c cap, s save")
		self.window.Render()

	def snapshot(self):
		copy = vtk.vtkPolyData()
		copy.DeepCopy(self.surface)
		self.history.append(copy)

	def press(self, caller, event):
		self.drag_start = self.interactor.GetEventPosition()
		self.box = bool(self.interactor.GetShiftKey())
		if not self.box:
			self.OnLeftButtonDown()

	def select_cells(self, indices, mode="replace"):
		indices = set(indices)
		if mode == "replace":
			self.selected = indices
		elif mode == "add":
			self.selected.update(indices)
		else:
			self.selected.symmetric_difference_update(indices)
		self.refresh()

	def release(self, caller, event):
		end = self.interactor.GetEventPosition()
		if self.box:
			x0, x1 = sorted((self.drag_start[0], end[0]))
			y0, y1 = sorted((self.drag_start[1], end[1]))
			selector = vtk.vtkHardwareSelector()
			selector.SetRenderer(self.renderer)
			selector.SetFieldAssociation(vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS)
			selector.SetArea(x0, y0, x1, y1)
			selection = selector.Select()
			indices = []
			if selection is not None:
				for node_index in range(selection.GetNumberOfNodes()):
					node = selection.GetNode(node_index)
					ids = node.GetSelectionList()
					for item in range(ids.GetNumberOfTuples()):
						indices.append(int(ids.GetValue(item)))
			self.select_cells(indices,
				"add" if self.interactor.GetControlKey() else "replace")
		else:
			self.OnLeftButtonUp()
			if abs(end[0]-self.drag_start[0])+abs(end[1]-self.drag_start[1]) < 5:
				picker = vtk.vtkCellPicker()
				picker.Pick(end[0], end[1], 0., self.renderer)
				index = picker.GetCellId()
				if index >= 0:
					self.select_cells([index],
						"toggle" if self.interactor.GetControlKey() else "replace")
		self.drag_start = None

	def assign(self, label):
		if not self.selected:
			return
		self.snapshot()
		values = self.surface.GetCellData().GetArray("boundary_id")
		for index in self.selected:
			values.SetValue(index, label)
		self.selected.clear()
		self.refresh()

	def undo(self):
		if not self.history:
			return
		self.surface = self.history.pop()
		self.selected.clear()
		self.refresh()

	def cut_selected(self):
		if not self.selected:
			return
		self.snapshot()
		self.surface = remove_selected(self.surface, self.selected)
		self.selected.clear()
		self.refresh()

	def cap(self, index, label=None):
		loops = boundary_loops(self.surface)
		self.snapshot()
		self.surface = cap_loop(self.surface, loops[index],
			self.outlet_label if label is None else label)
		self.selected.clear()
		self.refresh()

	def choose_cap(self):
		loops = boundary_loops(self.surface)
		if not loops:
			messagebox.showinfo("Surface BC labels", "No open boundary loops found")
			return
		index = simpledialog.askinteger("Cap opening",
			"Open loops (number of vertices): "
			+", ".join(f"{i}: {len(loop)}" for i, loop in enumerate(loops))
			+f"\nChoose a loop to cap with outlet label {self.outlet_label}:",
			minvalue=0, maxvalue=len(loops)-1)
		if index is None:
			return
		if messagebox.askyesno("Confirm end cap",
				f"Cap loop {index} with boundary label {self.outlet_label}?"):
			self.cap(index)

	def save(self):
		self.output.parent.mkdir(parents=True, exist_ok=True)
		write_surface(self.output, self.surface)
		labels = self.surface.GetCellData().GetArray("boundary_id")
		counts = {}
		for index in range(self.surface.GetNumberOfCells()):
			label = str(labels.GetValue(index))
			counts[label] = counts.get(label, 0)+1
		self.output.with_suffix(".labels.json").write_text(json.dumps({
			"surface": self.output.name, "boundary_array": "boundary_id",
			"wall_label": 0, "inflow_label": 1,
			"outflow_labels": sorted(int(label) for label in counts if int(label) >= 2),
			"boundary_triangles_by_label": counts}, indent=2)+"\n")
		if self.interactive:
			messagebox.showinfo("Surface BC labels", f"Saved {self.output}")

	def key(self, caller, event):
		key = self.interactor.GetKeySym().lower()
		if key == "i": self.assign(1)
		elif key == "o": self.assign(self.outlet_label)
		elif key == "w": self.assign(0)
		elif key == "n":
			self.outlet_label += 1
			self.refresh()
		elif key == "u": self.undo()
		elif key == "x" and self.selected:
			if messagebox.askyesno("Confirm cut",
					f"Cut {len(self.selected)} selected triangles? Use C to cap the opening."):
				self.cut_selected()
		elif key == "c": self.choose_cap()
		elif key == "s": self.save()

	def start(self):
		self.interactor.Initialize()
		self.window.Render()
		self.interactor.Start()


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("input_surface", type=Path)
	parser.add_argument("output_surface", type=Path)
	args = parser.parse_args()
	SurfaceLabeler(read_surface(args.input_surface), args.output_surface).start()


if __name__ == "__main__":
	main()
