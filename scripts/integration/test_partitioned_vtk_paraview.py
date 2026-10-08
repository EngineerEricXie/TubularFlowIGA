#!/usr/bin/env pvpython
"""Read partitioned output with ParaView's VTK reader and verify physical IDs."""
import json
import sys
from pathlib import Path
from paraview import servermanager
from paraview.simple import PVDReader
from vtkmodules.vtkIOXML import vtkXMLPUnstructuredGridReader, vtkXMLUnstructuredGridReader

root = Path(sys.argv[1])
base = 9007199254740993
for step, time in enumerate((.25, .5)):
    for rank in range(3):
        reader = vtkXMLUnstructuredGridReader()
        reader.SetFileName(str(root / f'step{step}.rank{rank}.vtu'))
        reader.Update()
        grid = reader.GetOutput()
        assert grid.GetNumberOfCells() == (0 if rank == 2 else 1)
        assert grid.GetNumberOfPoints() == (0 if rank == 2 else 4)
        assert grid.GetFieldData().GetArray('TimeValue').GetValue(0) == time
    reader = vtkXMLPUnstructuredGridReader()
    reader.SetFileName(str(root / f'step{step}.pvtu'))
    reader.Update()
    grid = reader.GetOutput()
    assert grid.GetNumberOfCells() == 2
    assert grid.GetNumberOfPoints() == 8  # shared IDs do not imply reader merging
    point_ids = grid.GetPointData().GetGlobalIds()
    cell_ids = grid.GetCellData().GetGlobalIds()
    assert point_ids is not None and cell_ids is not None
    assert sorted(cell_ids.GetValue(i) for i in range(2)) == [100, 101]
    for cell in range(2):
        identity = cell_ids.GetValue(cell)
        assert grid.GetCellType(cell) == 10
        assert grid.GetCellData().GetArray('owner').GetValue(cell) == identity - 100
        connectivity = grid.GetCell(cell).GetPointIds()
        ids = {point_ids.GetValue(connectivity.GetId(i)) for i in range(4)}
        assert ids == {base, base + 1, base + 2, base + 3 + identity - 100}
    seen = {}
    for point in range(8):
        identity = point_ids.GetValue(point)  # avoid floating-point GetTuple()
        xyz = grid.GetPoint(point)
        velocity = grid.GetPointData().GetArray('velocity').GetTuple(point)
        temperature = grid.GetPointData().GetArray('temperature & tracer').GetValue(point)
        assert velocity == tuple(x + time for x in xyz)
        assert temperature == sum(xyz) + time
        value = (xyz, velocity, temperature)
        if identity in seen:
            assert seen[identity] == value
        seen[identity] = value
    assert set(seen) == set(range(base, base + 5))
series = PVDReader(FileName=str(root / 'series.pvd'))
assert list(series.TimestepValues) == [.25, .5]
for time in series.TimestepValues:
    series.UpdatePipeline(time)
    result = servermanager.Fetch(series)
    assert result.GetNumberOfCells() == 2
    assert result.GetPointData().GetArray('velocity').GetTuple(0) == (time, time, time)
print(json.dumps({'status': 'passed', 'steps': 2, 'cells': 2, 'unique_point_ids': 5,
                  'empty_piece': True, 'paraview': str(servermanager.vtkSMProxyManager.GetVersionMajor()) + '.' + str(servermanager.vtkSMProxyManager.GetVersionMinor())}))
