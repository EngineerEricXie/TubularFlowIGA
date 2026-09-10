#!/usr/bin/env pvpython
"""Compare CPU flow PVD time frames against the serial VTKHDF reference."""
import json
import sys
from pathlib import Path
import xml.etree.ElementTree as ET
from paraview.simple import OpenDataFile, Delete
from paraview import servermanager
from vtkmodules.vtkIOHDF import vtkHDFReader
from vtkmodules.vtkCommonDataModel import vtkBezierHexahedron

root = Path(sys.argv[1]).resolve()
checked = 0
maximum_error = 0.0
for name, ranks in [('parallel-1', 1), ('parallel-2', 2), ('final-1', 1),
                    ('sparse-2', 2), ('stopped-2', 2), ('restarted-2', 2)]:
    path = root / name / 'result/flow.pvd'
    times = [float(node.attrib['timestep']) for node in ET.parse(path).iter('DataSet')]
    series = OpenDataFile(str(path))
    assert list(series.TimestepValues) == times
    reference = vtkHDFReader()
    reference.SetFileName(str(root / f'reference-{ranks}/result/flow.vtkhdf'))
    for time in times:
        reference.SetStep(round(time / .01))
        reference.Update()
        expected = reference.GetOutput()
        series.UpdatePipeline(time)
        actual = servermanager.Fetch(series)
        cell_ids = expected.GetCellData().GetArray('element_id')
        # The serial mesh deduplicates signatures. Its smallest occurrence is
        # the partition-independent physical ID used by the parallel writer.
        occurrences = {}
        for c in range(expected.GetNumberOfCells()):
            cell = expected.GetCell(c)
            for tensor in range(64):
                local = vtkBezierHexahedron.PointIndexFromIJK(tensor % 4, tensor // 4 % 4, tensor // 16, [3, 3, 3])
                point = cell.GetPointId(local)
                identity = 64 * cell_ids.GetValue(c) + tensor
                occurrences[point] = min(identity, occurrences.get(point, identity))
        lookup = {identity: point for point, identity in occurrences.items()}
        assert len(lookup) == expected.GetNumberOfPoints()
        actual_cells = actual.GetCellData().GetGlobalIds()
        assert sorted(actual_cells.GetValue(c) for c in range(actual.GetNumberOfCells())) == sorted(
            cell_ids.GetValue(c) for c in range(expected.GetNumberOfCells()))
        ids = actual.GetPointData().GetGlobalIds()
        seen = {}
        for point in range(actual.GetNumberOfPoints()):
            identity = ids.GetValue(point)
            target = lookup[identity]
            values = tuple(actual.GetPoint(point))
            baseline = tuple(expected.GetPoint(target))
            for field in ('velocity', 'pressure'):
                values += actual.GetPointData().GetArray(field).GetTuple(point)
                baseline += expected.GetPointData().GetArray(field).GetTuple(target)
            assert all(a == b for a, b in zip(values, baseline)), (name, time, identity, values, baseline)
            maximum_error = max(maximum_error, max(abs(a-b) for a, b in zip(values, baseline)))
            if identity in seen:
                assert seen[identity] == values
            seen[identity] = values
        assert set(seen) == set(lookup)
        checked += 1
    Delete(series)
print(json.dumps(dict(status='passed', frames=checked, maximum_absolute_error=maximum_error,
                      fields=['coordinates', 'velocity', 'pressure'], comparison='exact tuple equality')))
