#!/usr/bin/env python3
"""Generate an idealized, closed Y lumen in metres (wall=7, ports=1/2/3)."""
import argparse
from pathlib import Path
import numpy as np
from skimage.measure import marching_cubes


def generate(spacing):
    low=np.array([-.012,-.010,-.004]); high=np.array([.012,.010,.004])
    axes=[np.linspace(a,b,int(np.ceil((b-a)/spacing))+1) for a,b in zip(low,high)]
    x,y,z=np.meshgrid(*axes,indexing='ij')
    # A parent cylinder and two oblique daughters joined with a smooth minimum.
    trunk=np.sqrt(y*y+z*z)-.002
    trunk=np.maximum(trunk,x-.0003)
    branches=[]
    for sign in [-1,1]:
        axial=(x+.001)*.81923192+sign*y*.57346234
        radial=-(x+.001)*.57346234+sign*y*.81923192
        branches.append(np.maximum(np.sqrt(radial**2+z*z)-.00165,-axial))
    def union(a,b,k=.0007):
        h=np.maximum(k-np.abs(a-b),0)/k
        return np.minimum(a,b)-h*h*k*.25
    sdf=union(union(trunk,branches[0]),branches[1])
    sdf=np.maximum(sdf,np.maximum(-.009-x,x-.009))
    step=np.array([a[1]-a[0] for a in axes])
    vertices,faces,_,_=marching_cubes(sdf,0,spacing=step,allow_degenerate=False)
    vertices+=low
    signed=np.einsum('ij,ij->i',vertices[faces[:,0]],np.cross(vertices[faces[:,1]],vertices[faces[:,2]])).sum()/6
    if signed<0: faces=faces[:,[0,2,1]]
    labels=np.full(len(faces),7,dtype=int)
    for i,f in enumerate(faces):
        if np.max(np.abs(vertices[f,0]+.009))<1e-8: labels[i]=1
        if np.max(np.abs(vertices[f,0]-.009))<1e-8: labels[i]=2 if vertices[f,1].mean()<0 else 3
    assert set(labels)=={1,2,3,7}
    edges={}
    for f in faces:
        for a,b in zip(f,np.roll(f,-1)):
            key=tuple(sorted((int(a),int(b))));edges.setdefault(key,[]).append((int(a),int(b)))
    assert all(len(v)==2 and v[0]==v[1][::-1] for v in edges.values()),'not a closed oriented manifold'
    return vertices,faces,labels

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);p.add_argument('--spacing',type=float,default=.0008);a=p.parse_args()
    v,f,l=generate(a.spacing);a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('x') as out:
        out.write(f'{len(v)} {len(f)}\n')
        np.savetxt(out,v,fmt='%.17g');np.savetxt(out,np.column_stack([f,l]),fmt='%d')
    print('vertices',len(v),'triangles',len(f),'labels',dict(zip(*np.unique(l,return_counts=True))))
