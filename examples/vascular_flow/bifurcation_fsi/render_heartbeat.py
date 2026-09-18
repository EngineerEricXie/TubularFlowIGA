#!/usr/bin/env python3
"""Render only actual heartbeat solver snapshots, visibly labelled as a demonstration."""
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import imageio.v2 as imageio

root=Path(sys.argv[1])
verified=json.loads((root/'heartbeat-verification.json').read_text())
assert verified['status']=='visualization_verified'
h=np.atleast_1d(np.genfromtxt(root/'history.csv',delimiter=',',names=True))
frames=[]
maximum=float(h['max_displacement_m'].max()*1e6)
for row in h:
    tree=ET.parse(root/f"step-{int(row['step'])}/membrane-display.vtu")
    points=np.fromstring(tree.find('.//Points/DataArray').text,sep=' ').reshape(-1,3)*1000
    tris=np.fromstring(tree.find('.//Cells/DataArray[@Name="connectivity"]').text,sep=' ',dtype=int).reshape(-1,3)
    d=np.fromstring(tree.find('.//PointData/DataArray[@Name="normal_displacement_m"]').text,sep=' ')*1e6
    fig=plt.figure(figsize=(10,7),facecolor='white')
    ax=fig.add_subplot(211,projection='3d')
    norm=plt.Normalize(-maximum,maximum)
    ax.add_collection3d(Poly3DCollection(points[tris],facecolors=plt.get_cmap('coolwarm')(norm(d[tris].mean(axis=1))),edgecolors='none'))
    ax.set(xlim=(-10,10),ylim=(-9,9),zlim=(-4,4))
    ax.set_box_aspect((20,18,8));ax.view_init(65,-90);ax.set_axis_off()
    fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='coolwarm'),ax=ax,shrink=.7,label='Actual normal displacement (µm)')
    ax.set_title(f"Two-way FSI heartbeat · t = {row['time_s']:.2f} s\nWall displacement displayed ×{verified['display_scale']:g}")
    plot=fig.add_subplot(212)
    plot.plot(h['time_s'],verified['mean_normal_displacement_m'],color='darkred',label='Mean wall displacement (m)')
    plot.set(xlabel='Time (s)',ylabel='Mean wall displacement (m)')
    pressure=plot.twinx();pressure.plot(h['time_s'],h['inlet_pressure_pa'],color='navy',linestyle='--');pressure.set_ylabel('Inlet pressure (Pa)',color='navy')
    plot.axvline(row['time_s'],color='gray');plot.grid(alpha=.2)
    fig.text(.5,.01,'VISUALIZATION ONLY · coarse quadrature · conservation gate bypass enabled',ha='center',fontsize=10)
    fig.tight_layout(rect=(0,.04,1,1))
    fig.canvas.draw();frames.append(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy());plt.close(fig)
imageio.mimsave(root/'heartbeat.gif',frames,duration=100,loop=0)
imageio.imwrite(root/'heartbeat.png',frames[len(frames)//2-1])
