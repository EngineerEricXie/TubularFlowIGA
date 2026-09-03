#!/usr/bin/env python3

import tempfile
import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from vtk_centerline_to_obj import write_radius_annotated_obj


class ScalableValues:
    def __init__(self, values):
        self.values = values

    def __rmul__(self, scale):
        if self.values and isinstance(self.values[0], (list, tuple)):
            return [[scale * component for component in value] for value in self.values]
        return [scale * value for value in self.values]


class FakeMesh:
    def __init__(self, radii=(0.5, 0.25)):
        self.points = ScalableValues(((1.0, 2.0, 3.0), (4.0, 5.0, 6.0)))
        self.point_data = {"radius": ScalableValues(radii)}
        self.n_points = 2
        self.lines = [2, 0, 1]


class VtkCenterlineToObjTest(unittest.TestCase):
    def test_scale_applies_to_coordinates_and_radii(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "centerline.obj"
            write_radius_annotated_obj(
                FakeMesh(), output, scale=2.0, radius_array="radius"
            )
            text = output.read_text(encoding="utf-8")
            self.assertIn("v 2 4 6 1 0 0", text)
            self.assertIn("v 8 10 12 0.5 0 0", text)
            self.assertIn("l 1 2", text)

    def test_missing_radius_array_is_rejected(self):
        mesh = FakeMesh()
        mesh.point_data = {}
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "was not found"):
                write_radius_annotated_obj(
                    mesh,
                    Path(directory) / "centerline.obj",
                    scale=1.0,
                    radius_array="radius",
                )

    def test_nonpositive_radius_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "must be finite and positive"):
                write_radius_annotated_obj(
                    FakeMesh((0.5, 0.0)),
                    Path(directory) / "centerline.obj",
                    scale=1.0,
                    radius_array="radius",
                )

    def test_missing_lines_are_rejected_without_creating_output(self):
        mesh = FakeMesh()
        mesh.lines = []
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "centerline.obj"
            with self.assertRaisesRegex(ValueError, "no centerline polylines"):
                write_radius_annotated_obj(
                    mesh, output, scale=1.0, radius_array="radius"
                )
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
