#!/usr/bin/env pvpython
"""Read transport PVTU time series and compare fields across rank/frequency layouts."""

import argparse
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

from paraview import servermanager
from paraview.simple import Delete, OpenDataFile


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('root', type=Path)
parser.add_argument('--relative-tolerance', type=float, default=1e-6)
parser.add_argument('--absolute-tolerance', type=float, default=1e-12)
parser.add_argument('--output', type=Path)
args = parser.parse_args()
root = args.root.resolve()
fields = ('three_red', 'three_blue')


def read_case(name):
    path = root/name/'result/field.pvd'
    times = [float(node.attrib['timestep']) for node in ET.parse(path).iter('DataSet')]
    reader = OpenDataFile(str(path))
    assert list(reader.TimestepValues) == times
    frames = {}
    for time in times:
        reader.UpdatePipeline(time)
        data = servermanager.Fetch(reader)
        point_ids = data.GetPointData().GetGlobalIds()
        cell_ids = data.GetCellData().GetGlobalIds()
        assert point_ids is not None and cell_ids is not None
        cells = sorted(cell_ids.GetValue(index) for index in range(data.GetNumberOfCells()))
        values = {}
        for point in range(data.GetNumberOfPoints()):
            identity = point_ids.GetValue(point)
            row = tuple(data.GetPoint(point))
            for field in fields:
                array = data.GetPointData().GetArray(field)
                assert array is not None
                row += tuple(array.GetTuple(point))
            if identity in values:
                assert values[identity] == row
            values[identity] = row
        frames[time] = {'cells': cells, 'points': values}
    Delete(reader)
    return frames


cases = {name: read_case(name) for name in (
    'transport-ranks-1-frequency-1', 'transport-ranks-1-frequency-2',
    'transport-ranks-2-frequency-1', 'transport-ranks-2-frequency-2')}
reference = cases['transport-ranks-1-frequency-1']
comparisons = []
maximum_absolute = 0.0
for name, frames in cases.items():
    assert set(frames).issubset(reference)
    for time, frame in frames.items():
        expected = reference[time]
        assert frame['cells'] == expected['cells']
        assert set(frame['points']) == set(expected['points'])
        defect2 = {field: 0.0 for field in fields}
        reference2 = {field: 0.0 for field in fields}
        for identity, row in frame['points'].items():
            baseline = expected['points'][identity]
            assert row[:3] == baseline[:3]
            for offset, field in enumerate(fields, 3):
                defect = row[offset]-baseline[offset]
                defect2[field] += defect*defect
                reference2[field] += baseline[offset]*baseline[offset]
                maximum_absolute = max(maximum_absolute, abs(defect))
        for field in fields:
            absolute = math.sqrt(defect2[field])
            norm = math.sqrt(reference2[field])
            relative = absolute/norm if norm else None
            assert relative <= args.relative_tolerance if relative is not None else absolute <= args.absolute_tolerance
            comparisons.append({'case': name, 'time': time, 'field': field,
                                'absolute_l2': absolute, 'relative_l2': relative})

report = {'status': 'passed', 'frames': sum(map(len, cases.values())),
          'maximum_absolute_error': maximum_absolute, 'fields': list(fields),
          'comparisons': comparisons, 'relative_tolerance': args.relative_tolerance,
          'absolute_tolerance': args.absolute_tolerance}
text = json.dumps(report)
if args.output:
    args.output.write_text(json.dumps(report, indent=2)+'\n')
print(text)
