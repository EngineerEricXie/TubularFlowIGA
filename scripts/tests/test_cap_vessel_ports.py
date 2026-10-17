"""Tests for the explicit port cut-and-cap route (scripts/cap_vessel_ports.py)."""

from importlib.util import find_spec
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
TOOL = ROOT/"scripts"/"cap_vessel_ports.py"
PREFLIGHT = ROOT/"solvers"/"cpu"/"surface_fem_preflight"
SEGMENTS = (((0., 0., 0.), (0., 0., 10.), 1.5), ((0., 0., 10.), (5., 0., 18.), 1.0),
	((0., 0., 10.), (-5., 0., 18.), 1.0))


def write_y_surface(path):
	"""Closed blind-ended Y lumen: marching cubes of a union of capsules."""
	import numpy as np
	import vtk
	from vtk.util.numpy_support import numpy_to_vtk
	spacing = 0.25
	lower, upper = np.array([-8., -3., -3.]), np.array([8., 3., 21.])
	shape = np.ceil((upper-lower)/spacing).astype(int)+1
	axes = [lower[axis]+spacing*np.arange(shape[axis]) for axis in range(3)]
	x, y, z = np.meshgrid(*axes, indexing="ij")
	points = np.stack([x, y, z], axis=-1)
	field = np.full(x.shape, np.inf)
	for a, b, radius in SEGMENTS:
		a, b = np.array(a), np.array(b)
		t = np.clip(((points-a)@(b-a))/((b-a)@(b-a)), 0., 1.)
		field = np.minimum(field, np.linalg.norm(points-(a+t[..., None]*(b-a)), axis=-1)-radius)
	image = vtk.vtkImageData()
	image.SetDimensions(*map(int, shape))
	image.SetSpacing(spacing, spacing, spacing)
	image.SetOrigin(*lower)
	array = numpy_to_vtk(np.ascontiguousarray(field.transpose(2, 1, 0).ravel()), deep=True)
	array.SetName("distance")
	image.GetPointData().SetScalars(array)
	contour = vtk.vtkMarchingCubes()
	contour.SetInputData(image)
	contour.SetValue(0, 0.)
	contour.Update()
	reverse = vtk.vtkReverseSense()
	reverse.SetInputConnection(contour.GetOutputPort())
	reverse.ReverseCellsOn()
	reverse.Update()
	writer = vtk.vtkSTLWriter()
	writer.SetFileName(str(path))
	writer.SetInputConnection(reverse.GetOutputPort())
	writer.Write()


def ports_file(path, extra=()):
	branch = 1/math.sqrt(5.**2+8.**2)
	ports = [{"name": "inlet", "label": 1, "point": [0., 0., 0.], "outward_tangent": [0., 0., -1.],
		"radius": 1.5},
		{"name": "right", "label": 2, "point": [5., 0., 18.], "outward_tangent": [5*branch, 0., 8*branch],
		"radius": 1.0},
		{"name": "left", "label": 3, "point": [-5., 0., 18.], "outward_tangent": [-5*branch, 0., 8*branch],
		"radius": 1.0}, *extra]
	path.write_text(json.dumps({"schema_version": 1, "ports": ports}))
	return path


def read_vtp(path):
	import numpy as np
	import vtk
	from vtk.util.numpy_support import vtk_to_numpy
	reader = vtk.vtkXMLPolyDataReader()
	reader.SetFileName(str(path))
	reader.Update()
	data = reader.GetOutput()
	points = vtk_to_numpy(data.GetPoints().GetData()).astype(float)
	faces = vtk_to_numpy(data.GetPolys().GetConnectivityArray()).reshape(-1, 3)
	labels = vtk_to_numpy(data.GetCellData().GetArray("boundary_id"))
	return points, faces, labels, np


@unittest.skipUnless(find_spec("vtk") and find_spec("numpy"), "VTK and NumPy are required")
class CapVesselPortsTest(unittest.TestCase):
	def run_tool(self, *arguments, expect_success=True):
		result = subprocess.run([sys.executable, str(TOOL), *map(str, arguments)],
			text=True, capture_output=True, timeout=900, check=False)
		if expect_success and result.returncode:
			self.fail(result.stdout+result.stderr)
		if not expect_success and result.returncode == 0:
			self.fail("cap_vessel_ports unexpectedly accepted")
		return result

	def check_surface(self, path, minimum_quality):
		points, faces, labels, np = read_vtp(path)
		self.assertEqual(set(np.unique(labels).tolist()), {0, 1, 2, 3})
		directed = {}
		for face in faces:
			for corner in range(3):
				edge = (int(face[corner]), int(face[(corner+1) % 3]))
				self.assertNotIn(edge, directed, "surface is not consistently oriented")
				directed[edge] = True
		self.assertTrue(all((b, a) in directed for a, b in directed), "surface is not closed")
		a, b, c = points[faces[:, 0]], points[faces[:, 1]], points[faces[:, 2]]
		self.assertGreater(np.einsum("ij,ij->i", a, np.cross(b, c)).sum(), 0.)
		area = 0.5*np.linalg.norm(np.cross(b-a, c-a), axis=1)
		quality = 4*math.sqrt(3)*area/((b-a)**2+(c-b)**2+(a-c)**2).sum(axis=1)
		self.assertGreater(quality.min(), minimum_quality)
		for label, radius in ((1, 1.5), (2, 1.0), (3, 1.0)):
			cap = faces[labels == label]
			normal = np.cross(points[cap[0, 1]]-points[cap[0, 0]], points[cap[0, 2]]-points[cap[0, 0]])
			normal /= np.linalg.norm(normal)
			offsets = (points[np.unique(cap)]-points[cap[0, 0]])@normal
			self.assertLess(np.abs(offsets).max(), 1e-9, f"cap {label} is not planar")
			# Marching cubes slightly shrinks the lumen; the cap must still be a full section.
			self.assertGreater(area[labels == label].sum(), 0.8*math.pi*radius**2)
			self.assertLess(area[labels == label].sum(), 1.1*math.pi*radius**2)
		if PREFLIGHT.is_file():
			with tempfile.TemporaryDirectory() as directory:
				checked = subprocess.run([str(PREFLIGHT), str(path), str(Path(directory)/"surface.msh")],
					text=True, capture_output=True, timeout=300, check=False)
				self.assertEqual(checked.returncode, 0, checked.stdout+checked.stderr)

	def test_cut_cap_and_reject(self):
		with tempfile.TemporaryDirectory() as directory:
			directory = Path(directory)
			surface = directory/"y.stl"
			write_y_surface(surface)
			output = directory/"capped.vtp"
			self.run_tool(surface, ports_file(directory/"ports.json"), output, "--smoothing-iterations", 10)
			manifest = json.loads(output.with_suffix(".json").read_text())
			self.assertEqual([port["label"] for port in manifest["cut"]["ports"]], [1, 2, 3])
			self.assertEqual(manifest["source"]["sha256"], manifest["source"]["sha256"].lower())
			self.check_surface(output, 1e-6)
			far = {"name": "far", "label": 4, "point": [50., 50., 50.], "outward_tangent": [1., 0., 0.],
				"radius": 1.0}
			rejected = self.run_tool(surface, ports_file(directory/"bad.json", (far,)),
				directory/"bad.vtp", expect_success=False)
			self.assertIn("port far (label 4) cannot be cut", rejected.stderr)
			self.assertFalse((directory/"bad.vtp").exists())

	@unittest.skipUnless(find_spec("gmsh") and find_spec("scipy"), "Gmsh and SciPy are required")
	def test_remeshed_surface_has_exact_labels(self):
		with tempfile.TemporaryDirectory() as directory:
			directory = Path(directory)
			surface = directory/"y.stl"
			write_y_surface(surface)
			swc = directory/"y.swc"
			swc.write_text("1 3 0 0 0 1.5 -1\n2 3 0 0 10 1.5 1\n3 3 5 0 18 1.0 2\n4 3 -5 0 18 1.0 2\n")
			output = directory/"remeshed.vtp"
			self.run_tool(surface, ports_file(directory/"ports.json"), output, "--smoothing-iterations", 10,
				"--target-size", 0.5, "--centerline-swc", swc)
			manifest = json.loads(output.with_suffix(".json").read_text())
			self.assertEqual(manifest["remesh"]["mesher"], "Gmsh")
			self.check_surface(output, 0.2)


if __name__ == "__main__":
	unittest.main()
