#!/usr/bin/env python3
"""Read every accepted frame with VTK, including multiblock and unit metadata."""
import argparse
import csv
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


def read(reader_type, path):
    reader = reader_type()
    errors = []
    reader.AddObserver('ErrorEvent', lambda *_: errors.append(str(path)))
    reader.SetFileName(str(path))
    reader.Update()
    assert not errors and reader.GetErrorCode() == 0, f'VTK reader failed: {path}'
    return reader.GetOutput()


def verify(run):
    manifest = json.loads((run / 'run.json').read_text())
    assert manifest['status'] == 'passed'
    with (run / 'history.csv').open() as stream:
        history = list(csv.DictReader(stream))
    entries = ET.parse(run / 'fsi.pvd').getroot().findall('./Collection/DataSet')
    assert len(entries) == len(history) == manifest['steps']
    frames = []
    for entry, row in zip(entries, history):
        time = float(row['time_s'])
        assert float(entry.attrib['timestep']) == time
        path = run / entry.attrib['file']
        multi = read(vtk.vtkXMLMultiBlockDataReader, path)
        assert multi.GetNumberOfBlocks() == 2
        counts = {}
        for index, name in enumerate(('fluid', 'membrane')):
            assert multi.GetMetaData(index).Get(vtk.vtkCompositeDataSet.NAME()) == name
            block = multi.GetBlock(index)
            grid = read(vtk.vtkXMLUnstructuredGridReader, path.parent / f'{name}.vtu')
            assert grid.GetNumberOfPoints() > 0 and grid.GetNumberOfCells() > 0
            assert block.GetNumberOfPoints() == grid.GetNumberOfPoints()
            assert block.GetNumberOfCells() == grid.GetNumberOfCells()
            for data in (grid, block):
                assert data.GetFieldData().GetArray('time_s').GetTuple1(0) == time
                units = data.GetFieldData().GetAbstractArray('field_units')
                assert units is not None and units.GetNumberOfTuples() == 1
                assert 'velocity_m_per_s' in units.GetValue(0)
                assert np.isfinite(vtk_to_numpy(data.GetPoints().GetData())).all()
                fields = data.GetPointData()
                assert fields.GetNumberOfArrays() > 0
                for field_index in range(fields.GetNumberOfArrays()):
                    field = fields.GetArray(field_index)
                    assert field.GetNumberOfTuples() == data.GetNumberOfPoints()
                    assert np.isfinite(vtk_to_numpy(field)).all()
            counts[name] = {'points': grid.GetNumberOfPoints(), 'cells': grid.GetNumberOfCells()}
        frames.append({'time_s': time, 'datasets': counts})
    return {'status': 'passed', 'vtk_version': vtk.vtkVersion.GetVTKVersion(), 'frames': frames}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.run)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
