#!/usr/bin/env python3
"""Read heartbeat demonstration frames with native VTK, without formal acceptance."""
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from verify_vtk import read

root=Path(sys.argv[1])
manifest=json.loads((root/'run.json').read_text())
assert manifest['status']=='visualization_only'
frames=[]
for entry in ET.parse(root/'fsi.pvd').findall('.//DataSet'):
    path=root/entry.attrib['file'];time=float(entry.attrib['timestep'])
    blocks=read(vtk.vtkXMLMultiBlockDataReader,path)
    assert blocks.GetNumberOfBlocks()==2
    counts={}
    for name in ('fluid','membrane','membrane-display'):
        grid=read(vtk.vtkXMLUnstructuredGridReader,path.parent/f'{name}.vtu')
        assert grid.GetNumberOfPoints()>0 and grid.GetNumberOfCells()>0
        assert grid.GetFieldData().GetArray('time_s').GetTuple1(0)==time
        assert np.isfinite(vtk_to_numpy(grid.GetPoints().GetData())).all()
        data=grid.GetPointData()
        for i in range(data.GetNumberOfArrays()):
            array=data.GetArray(i)
            assert array.GetNumberOfTuples()==grid.GetNumberOfPoints()
            assert np.isfinite(vtk_to_numpy(array)).all()
        counts[name]=grid.GetNumberOfPoints()
        if name=='membrane-display':
            assert grid.GetFieldData().GetArray('display_displacement_scale').GetTuple1(0)==manifest['display_scale']
    frames.append({'time_s':time,'points':counts})
assert len(frames)==manifest['steps']
report={'status':'visualization_vtk_verified','formal_acceptance':False,'vtk_version':vtk.vtkVersion.GetVTKVersion(),'frames':frames}
(root/'heartbeat-vtk-verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
