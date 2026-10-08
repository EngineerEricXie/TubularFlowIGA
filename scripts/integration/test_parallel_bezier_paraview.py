#!/usr/bin/env pvpython
"""Validate actual cubic Bezier cells produced by the MPI partition builder."""
import json
import sys
from pathlib import Path
from vtkmodules.vtkCommonCore import reference
from vtkmodules.vtkIOXML import vtkXMLPUnstructuredGridReader

root = Path(sys.argv[1])
checked = 0
maximum_velocity_coordinate_error = 0.0
for group in ('world', 'self', 'pair'):
    for layout in range(3):
        reader = vtkXMLPUnstructuredGridReader()
        reader.SetFileName(str(root / group / f'layout{layout}' / 'snapshot.pvtu'))
        reader.Update()
        grid = reader.GetOutput()
        assert grid.GetNumberOfCells() == 2
        degrees = grid.GetCellData().GetHigherOrderDegrees()
        assert degrees is not None
        ids = grid.GetPointData().GetGlobalIds()
        seen = {}
        for point in range(grid.GetNumberOfPoints()):
            identity = ids.GetValue(point)
            xyz = grid.GetPoint(point)
            velocity = grid.GetPointData().GetArray('velocity').GetTuple(point)
            scalar = grid.GetPointData().GetArray('scalar').GetValue(point)
            error = max(abs(a-b) for a, b in zip(velocity, xyz))
            maximum_velocity_coordinate_error = max(maximum_velocity_coordinate_error, error)
            assert error < 1e-14
            assert abs(scalar - (xyz[0] + 2*xyz[1] + 3*xyz[2])) < 1e-14
            data = (xyz, velocity, scalar)
            if identity in seen:
                assert seen[identity] == data
            seen[identity] = data
        assert len(seen) == 112
        cell_ids = grid.GetCellData().GetGlobalIds()
        assert sorted(cell_ids.GetValue(i) for i in range(2)) == [0, 1]
        for cell_index in range(2):
            assert degrees.GetTuple(cell_index) == (3., 3., 3.)
            cell = grid.GetCell(cell_index)
            assert cell.GetCellType() == 79 and cell.GetNumberOfPoints() == 64
            location, weights = [0., 0., 0.], [0.] * 64
            cell.EvaluateLocation(reference(0), [.2, .4, .6], location, weights)
            expected = (cell_ids.GetValue(cell_index) + .2, .4, .6)
            assert max(abs(a-b) for a, b in zip(location, expected)) < 1e-13
            scalar = sum(weights[i] * grid.GetPointData().GetArray('scalar').GetValue(cell.GetPointId(i))
                         for i in range(64))
            assert abs(scalar - (expected[0] + 2*expected[1] + 3*expected[2])) < 1e-13
        checked += 1
print(json.dumps({'status': 'passed', 'snapshots': checked, 'unique_points': 112,
                  'bezier_cells': 2, 'parametric_interpolation': 'passed',
                  'maximum_velocity_coordinate_error': maximum_velocity_coordinate_error}))
