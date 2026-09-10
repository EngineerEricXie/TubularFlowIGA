#!/usr/bin/env pvpython
"""Compare CPU flow PVD time frames against the serial VTKHDF reference."""
import argparse
import ctypes
import ctypes.util
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET
from paraview.simple import OpenDataFile, Delete
from paraview import servermanager
from vtkmodules.vtkIOHDF import vtkHDFReader
from vtkmodules.vtkCommonDataModel import vtkBezierHexahedron

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('root', type=Path)
parser.add_argument('--parallel-ranks', type=int, nargs='+',
                    help='compare only parallel-R against reference-R for these rank counts')
parser.add_argument('--relative-tolerance', type=float, default=0.)
parser.add_argument('--absolute-tolerance', type=float, default=0.)
args = parser.parse_args()
assert all(math.isfinite(v) and v >= 0 for v in (args.relative_tolerance, args.absolute_tolerance))
root = args.root.resolve()
checked = 0
maximum_error = 0.0
comparisons = []
cases = ([(f'parallel-{r}', r) for r in args.parallel_ranks] if args.parallel_ranks else
         [('parallel-1', 1), ('parallel-2', 2), ('final-1', 1),
          ('sparse-2', 2), ('stopped-2', 2), ('restarted-2', 2)])
for name, ranks in cases:
    path = root / name / 'result/flow.pvd'
    times = [float(node.attrib['timestep']) for node in ET.parse(path).iter('DataSet')]
    series = OpenDataFile(str(path))
    assert list(series.TimestepValues) == times
    reference = vtkHDFReader()
    reference.SetFileName(str(root / f'reference-{ranks}/result/flow.vtkhdf'))
    reference.UpdateInformation()
    # ParaView 5.13 treats a one-step HDF as static (GetTimeValue returns 0).
    # Read the actual Steps/Values dataset rather than inventing a start time.
    hdf = ctypes.CDLL(ctypes.util.find_library('hdf5_serial'))
    high = ctypes.CDLL(ctypes.util.find_library('hdf5_serial_hl'))
    hdf.H5Fopen.argtypes = [ctypes.c_char_p, ctypes.c_uint, ctypes.c_longlong]
    hdf.H5Fopen.restype = ctypes.c_longlong
    hdf.H5Fclose.argtypes = [ctypes.c_longlong]
    high.H5LTget_dataset_ndims.argtypes = [ctypes.c_longlong, ctypes.c_char_p, ctypes.POINTER(ctypes.c_int)]
    high.H5LTget_dataset_info.argtypes = [ctypes.c_longlong, ctypes.c_char_p, ctypes.POINTER(ctypes.c_ulonglong), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_size_t)]
    high.H5LTread_dataset_double.argtypes = [ctypes.c_longlong, ctypes.c_char_p, ctypes.POINTER(ctypes.c_double)]
    handle = hdf.H5Fopen(reference.GetFileName().encode(), 0, 0)
    assert handle >= 0
    try:
        dimensions = (ctypes.c_ulonglong * 32)()
        kind, size = ctypes.c_int(), ctypes.c_size_t()
        dataset = b'/VTKHDF/Steps/Values'
        ndims = ctypes.c_int()
        assert high.H5LTget_dataset_ndims(handle, dataset, ctypes.byref(ndims)) == 0 and ndims.value == 1
        assert high.H5LTget_dataset_info(handle, dataset, dimensions, ctypes.byref(kind), ctypes.byref(size)) == 0
        assert dimensions[0] == reference.GetNumberOfSteps()
        values = (ctypes.c_double * dimensions[0])()
        assert high.H5LTread_dataset_double(handle, dataset, values) == 0
        reference_times = {value: step for step, value in enumerate(values)}
        assert len(reference_times) == len(values)
    finally:
        assert hdf.H5Fclose(handle) == 0
    for time in times:
        assert time in reference_times, (time, reference_times)
        reference.SetStep(reference_times[time])
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
        squared = {field: [0., 0.] for field in ('velocity', 'pressure')}
        for point in range(actual.GetNumberOfPoints()):
            identity = ids.GetValue(point)
            target = lookup[identity]
            values = tuple(actual.GetPoint(point))
            baseline = tuple(expected.GetPoint(target))
            for field in ('velocity', 'pressure'):
                values += actual.GetPointData().GetArray(field).GetTuple(point)
                baseline += expected.GetPointData().GetArray(field).GetTuple(target)
            assert values[:3] == baseline[:3], (name, time, identity, 'coordinates')
            if identity not in seen:
                for field, span in [('velocity', slice(3, 6)), ('pressure', slice(6, 7))]:
                    squared[field][0] += sum((a-b)**2 for a, b in zip(values[span], baseline[span]))
                    squared[field][1] += sum(b*b for b in baseline[span])
            maximum_error = max(maximum_error, max(abs(a-b) for a, b in zip(values, baseline)))
            if identity in seen:
                assert seen[identity] == values
            seen[identity] = values
        assert set(seen) == set(lookup)
        for field, (defect2, reference2) in squared.items():
            defect, norm = math.sqrt(defect2), math.sqrt(reference2)
            relative = defect / norm if norm else None
            assert (relative <= args.relative_tolerance if relative is not None else defect <= args.absolute_tolerance), (name, time, field, defect, relative)
            comparisons.append(dict(case=name, time=time, field=field, absolute_l2=defect, relative_l2=relative))
        checked += 1
    Delete(series)
print(json.dumps(dict(status='passed', frames=checked, maximum_absolute_error=maximum_error,
                      fields=['coordinates', 'velocity', 'pressure'], comparison='unique point field L2; exact coordinates and shared tuples', comparisons=comparisons,
                      relative_tolerance=args.relative_tolerance, zero_reference_absolute_tolerance=args.absolute_tolerance)))
