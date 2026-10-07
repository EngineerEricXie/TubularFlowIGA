#!/usr/bin/env python3
"""Read all registered HDF/VTU pairs using ParaView's actual VTK readers.

Run with pvpython --no-mpi. This verifies the file reader, not the GUI.
"""
import json
import sys
from pathlib import Path

import numpy as np
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkIOHDF import vtkHDFReader
from vtkmodules.vtkIOXML import vtkXMLPUnstructuredGridReader, vtkXMLUnstructuredGridReader


def normalized(grid):
    points = grid.GetPointData()
    cells = grid.GetCellData()
    point_ids = vtk_to_numpy(points.GetArray("GlobalPointIds"))
    cell_ids = vtk_to_numpy(cells.GetArray("GlobalCellIds"))
    assert point_ids.dtype == np.int64 and cell_ids.dtype == np.int64
    assert len(set(cell_ids)) == len(cell_ids)
    unique, indices = np.unique(point_ids, return_index=True)
    inverse = np.searchsorted(unique, point_ids)
    cell_order = np.argsort(cell_ids)
    result = {"point_ids": unique, "cell_ids": cell_ids[cell_order]}
    xyz = vtk_to_numpy(grid.GetPoints().GetData())
    result["coordinates"] = xyz[indices]
    np.testing.assert_array_equal(xyz, result["coordinates"][inverse])
    for prefix, data, order in (("point", points, indices), ("cell", cells, cell_order)):
        for index in range(data.GetNumberOfArrays()):
            array = data.GetArray(index)
            values = vtk_to_numpy(array)
            assert np.isfinite(values).all()
            result[prefix+":"+array.GetName()] = values[order]
            if prefix == "point":
                np.testing.assert_array_equal(values, values[indices][inverse])
    offsets = vtk_to_numpy(grid.GetCells().GetOffsetsArray())
    connectivity = vtk_to_numpy(grid.GetCells().GetConnectivityArray())
    rows = [tuple(point_ids[connectivity[offsets[i]:offsets[i+1]]]) for i in cell_order]
    result["connectivity"] = np.asarray(rows, dtype=np.int64)
    result["types"] = vtk_to_numpy(grid.GetCellTypesArray())[cell_order]
    return result


def main():
    manifest = Path(sys.argv[1])
    records = json.loads(manifest.read_text())
    checks = []
    for record in records["comparisons"]:
        reader = vtkHDFReader()
        reader.SetFileName(record["hdf"])
        reader.UpdateInformation()
        assert reader.GetNumberOfSteps() == len(record["snapshots"])
        for step, snapshot in enumerate(record["snapshots"]):
            reader.SetStep(step)
            reader.Update()
            assert reader.GetErrorCode() == 0
            grid = reader.GetOutput()
            assert abs(reader.GetTimeValue()-snapshot["time_s"]) < 1e-12
            actual = normalized(grid)
            assert len(actual["point_ids"]) == grid.GetNumberOfPoints()
            legacy = (vtkXMLPUnstructuredGridReader() if snapshot["file"].endswith(".pvtu")
                      else vtkXMLUnstructuredGridReader())
            legacy.SetFileName(snapshot["file"])
            legacy.Update()
            expected = normalized(legacy.GetOutput())
            assert actual.keys() == expected.keys()
            maximum = 0.
            for name in actual:
                if actual[name].dtype.kind in "iu":
                    np.testing.assert_array_equal(actual[name], expected[name])
                else:
                    np.testing.assert_allclose(actual[name], expected[name], rtol=1e-10, atol=1e-11)
                    maximum = max(maximum, float(np.max(np.abs(actual[name]-expected[name]))))
            checks.append({"hdf": record["hdf"], "step": step,
                           "maximum_absolute_difference": maximum})
    result = {"passed": True, "snapshots": checks, "reader": "ParaView VTK HDF reader"}
    manifest.with_name("reader_verification.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"passed": True, "snapshots": len(checks)}))


if __name__ == "__main__":
    main()
