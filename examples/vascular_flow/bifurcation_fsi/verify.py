#!/usr/bin/env python3
"""Independent checks of saved, accepted Y-vessel FSI datasets."""
import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np


def dataset(path):
    root=ET.parse(path).getroot()
    saved_time=float(root.find(".//FieldData/DataArray[@Name='time_s']").text)
    assert np.isfinite(saved_time)
    piece=root.find('.//Piece')
    n=int(piece.attrib['NumberOfPoints'])
    points=np.fromstring(piece.find('./Points/DataArray').text,sep=' ').reshape(n,3)
    fields={}
    for a in piece.findall('./PointData/DataArray'):
        fields[a.attrib['Name']]=np.fromstring(a.text,sep=' ').reshape(n,int(a.attrib.get('NumberOfComponents',1)))
    assert np.isfinite(points).all() and all(np.isfinite(v).all() for v in fields.values())
    return points,fields,saved_time


def verify(run):
    manifest=json.loads((run/'run.json').read_text());assert manifest['status']=='passed'
    history=np.atleast_1d(np.genfromtxt(run/'history.csv',delimiter=',',names=True))
    coupling=np.atleast_1d(np.genfromtxt(run/'coupling.csv',delimiter=',',names=True))
    assert len(history)==manifest['steps'] and np.array_equal(history['step'],np.arange(1,len(history)+1))
    assert all(np.isfinite(history[n]).all() for n in history.dtype.names)
    assert all(np.isfinite(coupling[n]).all() for n in coupling.dtype.names)
    assert np.allclose(np.diff(np.r_[0.,history['time_s']]),manifest['dt_s'],rtol=1e-13,atol=1e-15)
    assert set(coupling['step'])==set(history['step'])
    assert np.all(history['inlet_m3_s']<0) and np.all(history['lower_m3_s']>0) and np.all(history['upper_m3_s']>0)
    assert np.all((history['max_displacement_m']>0)&(history['max_displacement_m']<1e-4))
    assert np.all(history['mass_defect']<.03) and np.all(history['wall_leakage']<.03) and np.all(history['continuity']<1e-8)
    scale=history['normalization_scale_m3_s']
    expected_scale=np.maximum.reduce([
        np.full_like(scale,1e-6),
        np.abs(history['backward_euler_volume_rate_m3_s']),
        np.abs(history['total_fluid_surface_outward_flow_m3_s']),
        np.abs(history['total_material_surface_outward_flow_m3_s']),
    ])
    assert np.array_equal(scale,expected_scale)
    assert np.all(history['continuity']>=0)
    assert np.allclose(history['mass_defect'],np.abs(history['moving_mass_defect_m3_s'])/scale,rtol=1e-12,atol=1e-15)
    assert np.allclose(history['wall_leakage'],np.abs(history['wall_relative_leakage_m3_s'])/scale,rtol=1e-12,atol=1e-15)
    reconstructed=history['backward_euler_volume_rate_m3_s']+history['total_fluid_surface_outward_flow_m3_s']-history['total_material_surface_outward_flow_m3_s']
    assert np.allclose(history['moving_mass_defect_m3_s'],reconstructed,rtol=1e-12,atol=1e-18)
    collection=ET.parse(run/'fsi.pvd').getroot().find('Collection')
    assert len(collection)==len(history)
    for entry,row in zip(collection,history):
        assert float(entry.attrib['timestep'])==row['time_s']
        assert entry.attrib['file']==f"step-{int(row['step'])}/fsi.vtm"
        frame=run/entry.attrib['file'];assert frame.is_file()
        blocks=ET.parse(frame).getroot().find('vtkMultiBlockDataSet')
        assert {b.attrib['name'] for b in blocks}=={'fluid','membrane'}
        for block in blocks:assert (frame.parent/block.attrib['file']).is_file()
    previous=None;reference=None;normals=None;maximum_speed=0.;maximum_velocity_error=0.
    for row in history:
        i=int(row['step']);h=coupling[coupling['step']==i]
        assert len(h)==int(row['iterations']) and h[-1]['rms_m']<=h[-1]['threshold_m']
        points,w,saved_time=dataset(run/f'step-{i}/membrane.vtu')
        assert saved_time==row['time_s']
        d=w['displacement_m'];v=w['velocity_m_per_s'];ref=w['reference_coordinates_m']
        assert np.allclose(points,ref+d,rtol=1e-13,atol=1e-15)
        if previous is None:previous=np.zeros_like(d);reference=ref.copy();normals=w['reference_normal'].copy()
        assert np.array_equal(ref,reference) and np.array_equal(w['reference_normal'],normals)
        error=float(np.max(np.abs(v-(d-previous)/manifest['dt_s'])));maximum_velocity_error=max(maximum_velocity_error,error)
        assert error<1e-10
        assert np.isclose(np.linalg.norm(d,axis=1).max(),row['max_displacement_m'],rtol=1e-10,atol=1e-15)
        assert np.linalg.norm(w['fluid_traction_on_structure_pa'])>0
        previous=d.copy()
        _,f,saved_time=dataset(run/f'step-{i}/fluid.vtu')
        assert saved_time==row['time_s']
        assert np.allclose(np.linalg.norm(f['velocity_m_per_s'],axis=1),f['speed_m_per_s'][:,0],rtol=1e-12,atol=1e-14)
        maximum_speed=max(maximum_speed,float(f['speed_m_per_s'].max()))
    assert maximum_speed>0
    return {'status':'passed','steps':len(history),'maximum_speed_m_s':maximum_speed,'maximum_wall_displacement_m':float(history['max_displacement_m'].max()),'maximum_wall_velocity_consistency_error_m_s':maximum_velocity_error,'maximum_moving_mass_defect':float(history['mass_defect'].max()),'maximum_wall_leakage':float(history['wall_leakage'].max()),'maximum_continuity_defect':float(history['continuity'].max()),'final_inlet_m3_s':float(history[-1]['inlet_m3_s']),'final_lower_m3_s':float(history[-1]['lower_m3_s']),'final_upper_m3_s':float(history[-1]['upper_m3_s'])}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();result=verify(a.run);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
