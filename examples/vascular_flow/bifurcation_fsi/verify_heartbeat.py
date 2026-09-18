#!/usr/bin/env python3
"""Check actual solver frames and expansion/contraction, without granting acceptance."""
import json
import sys
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np


def verify(root):
    manifest=json.loads((root/'run.json').read_text())
    assert manifest['status']=='visualization_only'
    h=np.atleast_1d(np.genfromtxt(root/'history.csv',delimiter=',',names=True))
    assert len(h)==manifest['steps']
    assert np.isfinite(np.column_stack([h[k] for k in h.dtype.names])).all()
    assert list(h['step'])==list(range(1,len(h)+1))
    assert np.allclose(h['inlet_pressure_pa'],manifest['inlet_pressure_pa']+manifest['pulse_amplitude_pa']*.5*(1-np.cos(2*np.pi*h['time_s']/manifest['period_s'])),rtol=1e-12,atol=1e-12)
    for name in ('fsi.pvd','wall-display.pvd'):
        frames=ET.parse(root/name).findall('.//DataSet')
        assert len(frames)==len(h)
        assert np.allclose([float(x.attrib['timestep']) for x in frames],h['time_s'])
        assert all((root/x.attrib['file']).is_file() for x in frames)
    coupling=np.atleast_1d(np.genfromtxt(root/'coupling.csv',delimiter=',',names=True))
    assert np.isfinite(np.column_stack([coupling[k] for k in coupling.dtype.names])).all()
    for row in h:
        rounds=coupling[coupling['step']==row['step']]
        assert len(rounds)==int(row['iterations']) and len(rounds)>0
        assert rounds[-1]['rms_m']<=rounds[-1]['threshold_m']
    means=[]
    previous=None
    for step in range(1,len(h)+1):
        directory=root/f'step-{step}'
        fields={}
        for name in ('fluid','membrane','membrane-display'):
            tree=ET.parse(directory/f'{name}.vtu')
            data={a.attrib['Name']:np.fromstring(a.text,sep=' ') for a in tree.findall('.//PointData/DataArray')}
            assert all(np.isfinite(v).all() for v in data.values())
            xyz=np.fromstring(tree.find('.//Points/DataArray').text,sep=' ').reshape(-1,3)
            assert np.isfinite(xyz).all()
            if name.startswith('membrane'):
                ref=data['reference_coordinates_m'].reshape(-1,3)
                disp=data['displacement_m'].reshape(-1,3)
                scale=manifest['display_scale'] if name=='membrane-display' else 1.
                assert np.allclose(xyz,ref+scale*disp,rtol=1e-12,atol=1e-14)
            fields[name]=data
        actual=fields['membrane']
        displacement=actual['displacement_m'].reshape(-1,3)
        if previous is None:
            previous=np.zeros_like(displacement)
        velocity=actual['velocity_m_per_s'].reshape(-1,3)
        assert np.allclose(velocity,(displacement-previous)/manifest['dt_s'],rtol=1e-8,atol=1e-12)
        previous=displacement
        means.append(float(actual['normal_displacement_m'].mean()))
    delta=np.diff(means)
    assert delta.max()>1e-10 and delta.min()<-1e-10, 'Wall did not both expand and contract'
    report={'status':'visualization_verified','formal_acceptance':False,'frames':len(h),'duration_s':float(h['time_s'][-1]),'mean_normal_displacement_m':means,'expansion_observed':True,'contraction_observed':True,'conservation_failed_steps':int(np.sum(h['conservation_passed']==0)),'display_scale':manifest['display_scale']}
    (root/'heartbeat-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    verify(Path(sys.argv[1]))
