#!/usr/bin/env python3
"""Generate a watertight labelled Y-pipe comparison surface with Gmsh OCC."""

import argparse
import math
from pathlib import Path


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("output_vtp")
	parser.add_argument("--radius-m", type=float, default=0.005)
	parser.add_argument("--trunk-length-m", type=float, default=0.05)
	parser.add_argument("--branch-x-m", type=float, default=0.04)
	parser.add_argument("--branch-y-m", type=float, default=0.03)
	parser.add_argument("--target-size-m", type=float, default=0.003)
	args = parser.parse_args()
	values = (args.radius_m, args.trunk_length_m, args.branch_x_m,
		args.branch_y_m, args.target_size_m)
	if any(not math.isfinite(value) or value <= 0.0 for value in values):
		raise ValueError("Y-pipe dimensions and target size must be finite and positive")
	if args.trunk_length_m <= 2.0*args.radius_m:
		raise ValueError("Y-pipe trunk must be longer than its diameter")
	try:
		import gmsh
	except ImportError as error:
		raise ValueError(f"Gmsh Python module is unavailable: {error}")
	gmsh.initialize()
	try:
		gmsh.option.setNumber("General.Terminal", 0)
		gmsh.model.add("labelled-y-pipe-surface")
		trunk = gmsh.model.occ.addCylinder(0,0,0,args.trunk_length_m,0,0,args.radius_m)
		start = args.trunk_length_m-2.0*args.radius_m
		upper = gmsh.model.occ.addCylinder(start,0,0,args.branch_x_m,args.branch_y_m,0,
			args.radius_m)
		lower = gmsh.model.occ.addCylinder(start,0,0,args.branch_x_m,-args.branch_y_m,0,
			args.radius_m)
		volumes, _ = gmsh.model.occ.fuse([(3,trunk)], [(3,upper),(3,lower)],
			removeObject=True, removeTool=True)
		gmsh.model.occ.synchronize()
		volume_tags = [tag for dimension, tag in volumes if dimension == 3]
		if len(volume_tags) != 1:
			raise RuntimeError("Y-pipe Boolean union did not produce one fluid volume")
		volume = volume_tags[0]
		oriented = gmsh.model.getBoundary([(3,volume)], oriented=True, recursive=False)
		surface_sign = {abs(tag): (-1 if tag < 0 else 1)
			for dimension, tag in oriented if dimension == 2}
		endpoints = ((0.0,0.0,0.0),
			(start+args.branch_x_m,args.branch_y_m,0.0),
			(start+args.branch_x_m,-args.branch_y_m,0.0))
		labels = {}
		for dimension, tag in gmsh.model.getEntities(2):
			center = gmsh.model.occ.getCenterOfMass(dimension, tag)
			distances = [math.dist(center, endpoint) for endpoint in endpoints]
			nearest = min(range(3), key=distances.__getitem__)
			labels[tag] = (1,2,3)[nearest] if distances[nearest] < 0.75*args.radius_m else 0
		if set(labels.values()) != {0,1,2,3}:
			raise RuntimeError("Y-pipe could not identify wall/inlet/two outlet surfaces")
		gmsh.option.setNumber("Mesh.MeshSizeMin", args.target_size_m)
		gmsh.option.setNumber("Mesh.MeshSizeMax", args.target_size_m)
		gmsh.option.setNumber("Mesh.ElementOrder", 1)
		gmsh.model.mesh.generate(2)
		node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
		nodes = {int(tag): tuple(coordinates[3*index:3*index+3])
			for index, tag in enumerate(node_tags)}
		triangles, triangle_labels = [], []
		for tag, label in labels.items():
			types, _, connectivity = gmsh.model.mesh.getElements(2, tag)
			for type_, values in zip(types, connectivity):
				if type_ != 2: raise RuntimeError("Y-pipe surface contains a non-triangle element")
				values = list(map(int, values))
				for index in range(0, len(values), 3):
					face = tuple(values[index:index+3])
					if surface_sign.get(tag, 1) < 0: face = (face[0],face[2],face[1])
					triangles.append(face); triangle_labels.append(label)
		if not triangles: raise RuntimeError("Y-pipe surface meshing produced no triangles")
		ordered = sorted(nodes); remap = {tag: index for index, tag in enumerate(ordered)}
		points = [nodes[tag] for tag in ordered]
		faces = [tuple(remap[tag] for tag in face) for face in triangles]
	finally:
		gmsh.finalize()
	coordinates = " ".join(f"{value:.17g}" for point in points for value in point)
	connectivity = " ".join(str(node) for face in faces for node in face)
	offsets = " ".join(str(3*(index+1)) for index in range(len(faces)))
	boundary_ids = " ".join(map(str, triangle_labels))
	contents = f'''<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian"><PolyData>
<Piece NumberOfPoints="{len(points)}" NumberOfPolys="{len(faces)}">
<Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">{coordinates}</DataArray></Points>
<Polys><DataArray type="Int64" Name="connectivity" format="ascii">{connectivity}</DataArray>
<DataArray type="Int64" Name="offsets" format="ascii">{offsets}</DataArray></Polys>
<CellData Scalars="boundary_id"><DataArray type="Int32" Name="boundary_id" format="ascii">{boundary_ids}</DataArray></CellData>
</Piece></PolyData></VTKFile>
'''
	output = Path(args.output_vtp); output.parent.mkdir(parents=True, exist_ok=True)
	output.write_text(contents, encoding="utf-8")
	counts = {label: triangle_labels.count(label) for label in sorted(set(triangle_labels))}
	print(f"y_pipe_surface: PASS points={len(points)} triangles={len(faces)} labels={counts}")
	return 0


if __name__ == "__main__":
	try:
		raise SystemExit(main())
	except (OSError, RuntimeError, ValueError) as error:
		print(f"y_pipe_surface: ERROR: {error}")
		raise SystemExit(2)
