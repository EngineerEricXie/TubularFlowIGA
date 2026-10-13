#!/usr/bin/env pvpython
"""Read the real MPI writer's independent communicator snapshots."""
from pathlib import Path
import json
import sys
from vtkmodules.vtkIOXML import vtkXMLPUnstructuredGridReader

root = Path(sys.argv[1])
count = 0
for group, ranks, seed in [('world', 3, 10), ('self', 1, 20), ('pair', 2, 30)]:
    expected = 1 if ranks == 1 else ranks - 1
    indices = sorted((root / group).glob('*/snapshot.pvtu'))
    assert len(indices) == (3 if ranks == 1 else 5)
    for path in indices:
        reader = vtkXMLPUnstructuredGridReader()
        reader.SetFileName(str(path))
        reader.Update()
        grid = reader.GetOutput()
        assert grid.GetNumberOfCells() == expected
        assert grid.GetNumberOfPoints() == expected
        ids = grid.GetPointData().GetGlobalIds()
        cell_ids = grid.GetCellData().GetGlobalIds()
        values = grid.GetPointData().GetArray('value')
        for point in range(expected):
            rank = ids.GetValue(point)
            assert rank == cell_ids.GetValue(point)
            assert grid.GetPoint(point) == (rank, seed, 0)
            assert values.GetValue(point) == seed + rank
        assert sorted(ids.GetValue(i) for i in range(expected)) == list(range(expected))
        count += 1
    assert not (root / group / 'write' / 'snapshot.pvtu').exists()
print(json.dumps({'status': 'passed', 'snapshots': count, 'groups': [3, 1, 2]}))
