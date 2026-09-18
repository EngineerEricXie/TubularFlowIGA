#!/usr/bin/env python3
"""Render accepted FSI VTU samples; animation never invents intermediate states."""
import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import imageio.v2 as imageio


def read(path):
    piece=ET.parse(path).getroot().find('.//Piece')
    points=np.fromstring(piece.find('./Points/DataArray').text,sep=' ').reshape(-1,3)
    fields={a.attrib['Name']:np.fromstring(a.text,sep=' ').reshape(-1,int(a.attrib.get('NumberOfComponents',1))) for a in piece.findall('./PointData/DataArray')}
    cells=np.fromstring(piece.find('./Cells/DataArray[@Name="connectivity"]').text,sep=' ',dtype=int)
    return points,fields,cells


def section_mask(x,y):
    trunk=np.maximum(np.abs(y)-.002,x-.0003)
    def union(a,b,k=.0007):
        h=np.maximum(k-np.abs(a-b),0)/k
        return np.minimum(a,b)-h*h*k*.25
    sdf=trunk
    for sign in [-1,1]:
        axial=(x+.001)*.81923192+sign*y*.57346234
        radial=-(x+.001)*.57346234+sign*y*.81923192
        sdf=union(sdf,np.maximum(np.abs(radial)-.00165,-axial))
    return np.maximum(sdf,np.maximum(-.009-x,x-.009))<-.00012


def render(run,output,warp,preview_step=None):
    history=np.atleast_1d(np.genfromtxt(run/'history.csv',delimiter=',',names=True))
    if preview_step is None:
        manifest=json.loads((run/'run.json').read_text())
        if manifest.get('status')!='passed': raise ValueError('Only passed solver runs may be rendered as a showcase')
        steps=[int(x) for x in history['step']]
        if steps!=list(range(1,manifest['steps']+1)): raise ValueError('Incomplete accepted history')
    else:
        history=history[history['step']==preview_step]
        if len(history)!=1: raise ValueError('Preview requires one existing accepted step')
        row=history[0]
        coupling=np.atleast_1d(np.genfromtxt(run/'coupling.csv',delimiter=',',names=True))
        coupling=coupling[coupling['step']==preview_step]
        valid=all(np.isfinite(row[name]) for name in history.dtype.names)
        valid=valid and len(coupling)==int(row['iterations']) and len(coupling)>0
        valid=valid and np.isfinite(coupling['rms_m']).all() and np.isfinite(coupling['threshold_m']).all()
        valid=valid and coupling[-1]['rms_m']<=coupling[-1]['threshold_m']
        valid=valid and 0<row['max_displacement_m']<1e-4 and 0<=row['mass_defect']<.03
        valid=valid and 0<=row['wall_leakage']<.03 and 0<=row['continuity']<1e-8
        valid=valid and row['inlet_m3_s']<0 and row['lower_m3_s']>0 and row['upper_m3_s']>0
        if not valid: raise ValueError('Preview step did not satisfy numerical acceptance gates')
        steps=[preview_step]
    frames=[];vmax=0.;pmin=0.;pmax=0.;dmax=0.
    for i in steps:
        _,f,_=read(run/f'step-{i}/fluid.vtu');_,w,_=read(run/f'step-{i}/membrane.vtu')
        vmax=max(vmax,float(f['speed_m_per_s'].max()));pmin=min(pmin,float(f['pressure_pa'].min()));pmax=max(pmax,float(f['pressure_pa'].max()));dmax=max(dmax,float(np.linalg.norm(w['displacement_m'],axis=1).max()*1e6))
    output.mkdir(parents=True,exist_ok=True)
    x=np.linspace(-.009,.009,240);y=np.linspace(-.009,.009,230);xx,yy=np.meshgrid(x,y)
    query=np.column_stack([xx.ravel(),yy.ravel(),np.zeros(xx.size)]);mask=section_mask(xx,yy)
    plt.style.use('dark_background')
    for i,row in zip(steps,history):
        points,f,_=read(run/f'step-{i}/fluid.vtu');wall,w,cells=read(run/f'step-{i}/membrane.vtu')
        tree=cKDTree(points);distance,index=tree.query(query,k=8)
        weights=1/np.maximum(distance,1e-10)**2;weights/=weights.sum(axis=1)[:,None]
        def interpolate(values):
            out=(values[index]*weights[:,:,None]).sum(axis=1)
            return out.reshape(xx.shape+(out.shape[-1],))
        uv=interpolate(f['velocity_m_per_s']);speed=np.linalg.norm(uv,axis=2);pressure=interpolate(f['pressure_pa'])[:,:,0]
        fig=plt.figure(figsize=(16,8),facecolor='#101925');gs=fig.add_gridspec(2,3,height_ratios=[4,1.35])
        a=fig.add_subplot(gs[0,0]);b=fig.add_subplot(gs[0,1]);c=fig.add_subplot(gs[0,2],projection='3d')
        for ax,data,title,cmap,limits in [(a,speed,'Blood velocity · m/s','turbo',(0,vmax)),(b,pressure,'Fluid pressure · Pa','coolwarm',(pmin,pmax))]:
            m=ax.pcolormesh(xx*1e3,yy*1e3,np.ma.masked_where(~mask,data),cmap=cmap,vmin=limits[0],vmax=limits[1],shading='auto');fig.colorbar(m,ax=ax,shrink=.65,pad=.02)
            ax.set(title=title,xlabel='x · mm',ylabel='y · mm',aspect='equal');ax.set_facecolor('#101925')
        a.streamplot(x*1e3,y*1e3,np.ma.masked_where(~mask,uv[:,:,0]),np.ma.masked_where(~mask,uv[:,:,1]),density=.8,color='white',linewidth=.6,arrowsize=.7)
        displacement=w['displacement_m'];reference=wall-displacement;warped=reference+warp*displacement;tris=cells.reshape(-1,3)
        scalar=np.linalg.norm(displacement,axis=1)*1e6;norm=Normalize(0,max(dmax,1e-12));colors=plt.get_cmap('plasma')(norm(scalar[tris].mean(axis=1)))
        c.add_collection3d(Poly3DCollection(warped[tris]*1e3,facecolors=colors,edgecolors='none',alpha=1.))
        c.set(xlim=(-10,10),ylim=(-9,9),zlim=(-4,4),title=f'Wall displacement · µm\nDeformation displayed ×{warp:g}')
        c.set_box_aspect((20,18,8));c.view_init(42,-105);c.set_axis_off();fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='plasma'),ax=c,shrink=.65,pad=.01)
        q=fig.add_subplot(gs[1,:2]);q.plot(history['time_s']*1000,-history['inlet_m3_s']*1e6,label='Inlet',marker='o' if preview_step is not None else None);q.plot(history['time_s']*1000,history['lower_m3_s']*1e6,label='Lower outlet',marker='o' if preview_step is not None else None);q.plot(history['time_s']*1000,history['upper_m3_s']*1e6,'--',label='Upper outlet',marker='o' if preview_step is not None else None);q.axvline(row['time_s']*1000,color='white',alpha=.5);q.set(xlabel='Simulation time · ms',ylabel='Flow · mL/s');q.legend(ncol=3,fontsize=9);q.grid(alpha=.15)
        note=fig.add_subplot(gs[1,2]);note.set_axis_off();note.text(0,.9,f"Accepted step {i} · {row['time_s']*1000:.0f} ms\nStrong coupling: {int(row['iterations'])} iterations\nMax wall motion: {row['max_displacement_m']*1e6:.3g} µm\nMoving continuity: {row['continuity']:.2e}",va='top',fontsize=11)
        title='TubularFlowIGA  |  Bifurcating vessel · two-way FSI'
        if preview_step is not None: title=f'PREVIEW · accepted step {preview_step} only · full run incomplete'
        fig.suptitle(title,fontsize=22,y=.98)
        fig.text(.5,.015,'Blood-like Newtonian fluid · small-displacement pre-tensioned membrane · startup under constant pressure\nCenter-plane fields interpolated from solver samples; wall deformation scale explicitly amplified',ha='center',fontsize=10,color='#b5c5d5')
        fig.subplots_adjust(top=.88,bottom=.14,hspace=.32,wspace=.3)
        frame=output/f'frame-{i:03d}.png';fig.savefig(frame,dpi=125,facecolor=fig.get_facecolor());plt.close(fig);frames.append(imageio.imread(frame))
    if preview_step is None:
        imageio.mimsave(output/'bifurcation-fsi.gif',frames,duration=[450]*(len(frames)-1)+[1500],loop=0)
        imageio.imwrite(output/'bifurcation-fsi.png',frames[-1])
    else:
        imageio.imwrite(output/f'preview-step-{preview_step}.png',frames[-1])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('output',type=Path);p.add_argument('--warp',type=float,default=100.);p.add_argument('--preview-step',type=int,help='Render one numerically accepted step with an explicit incomplete-run label');a=p.parse_args();render(a.run,a.output,a.warp,a.preview_step)
