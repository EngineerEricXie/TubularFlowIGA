#!/usr/bin/env python3
"""Generate a closed labelled constant-radius toroidal pipe segment."""

import argparse
import math
from pathlib import Path


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("output_vtp")
	parser.add_argument("--bend-radius-m", type=float, required=True)
	parser.add_argument("--tube-radius-m", type=float, required=True)
	parser.add_argument("--bend-angle-deg", type=float, required=True)
	parser.add_argument("--circumferential-segments", type=int, required=True)
	parser.add_argument("--axial-segments", type=int, required=True)
	parser.add_argument("--radial-segments", type=int, default=1)
	args = parser.parse_args()
	angle = math.radians(args.bend_angle_deg)
	if not math.isfinite(args.bend_radius_m) or args.bend_radius_m <= 0.0:
		raise ValueError("bend radius must be finite and positive")
	if not math.isfinite(args.tube_radius_m) or not 0.0 < args.tube_radius_m < args.bend_radius_m:
		raise ValueError("tube radius must be positive and smaller than bend radius")
	if not math.isfinite(angle) or not 0.0 < angle < math.pi:
		raise ValueError("bend angle must be in (0,180) degrees")
	if args.circumferential_segments < 8 or args.axial_segments < 2 \
			or args.radial_segments < 1:
		raise ValueError("bent pipe requires at least 8 circumferential and 2 axial segments")
	n, m = args.circumferential_segments, args.axial_segments
	points = []
	def frame(theta):
		center = (args.bend_radius_m*math.sin(theta),
			args.bend_radius_m*(1.0-math.cos(theta)), 0.0)
		normal = (-math.sin(theta), math.cos(theta), 0.0)
		binormal = (0.0, 0.0, 1.0)
		return center, normal, binormal
	def ring_point(theta, radius, phi):
		center, normal, binormal = frame(theta)
		return tuple(center[i]+radius*(math.cos(phi)*normal[i]
			+math.sin(phi)*binormal[i]) for i in range(3))
	for axial in range(m+1):
		theta = angle*axial/m
		for circumferential in range(n):
			phi = 2.0*math.pi*circumferential/n
			points.append(ring_point(theta, args.tube_radius_m, phi))
	triangles, labels = [], []
	for axial in range(m):
		for circumferential in range(n):
			next_circumferential = (circumferential+1) % n
			a = axial*n+circumferential
			b = (axial+1)*n+circumferential
			c = (axial+1)*n+next_circumferential
			d = axial*n+next_circumferential
			triangles.extend(((a,c,b), (a,d,c))); labels.extend((0,0))
	def cap(theta, outer, label, positive):
		center_index = len(points); points.append(frame(theta)[0])
		rings = []
		for radial in range(1, args.radial_segments):
			ring = []
			for circumferential in range(n):
				phi = 2.0*math.pi*circumferential/n
				ring.append(len(points))
				points.append(ring_point(theta,
					args.tube_radius_m*radial/args.radial_segments, phi))
			rings.append(ring)
		rings.append(outer)
		first = rings[0]
		for circumferential in range(n):
			next_circumferential = (circumferential+1) % n
			triangle = (center_index, first[circumferential], first[next_circumferential])
			triangles.append(triangle if positive else (triangle[0],triangle[2],triangle[1]))
			labels.append(label)
		for inner, outer_ring in zip(rings, rings[1:]):
			for circumferential in range(n):
				next_circumferential = (circumferential+1) % n
				a, b = inner[circumferential], outer_ring[circumferential]
				c, d = outer_ring[next_circumferential], inner[next_circumferential]
				triangles.extend(((a,b,c),(a,c,d)) if positive else ((a,c,b),(a,d,c)))
				labels.extend((label,label))
	cap(0.0, list(range(n)), 1, False)
	cap(angle, list(range(m*n,(m+1)*n)), 2, True)
	coordinates = " ".join(f"{value:.17g}" for point in points for value in point)
	connectivity = " ".join(str(node) for face in triangles for node in face)
	offsets = " ".join(str(3*(index+1)) for index in range(len(triangles)))
	boundary_ids = " ".join(map(str, labels))
	contents = f'''<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian"><PolyData>
<Piece NumberOfPoints="{len(points)}" NumberOfPolys="{len(triangles)}">
<Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">{coordinates}</DataArray></Points>
<Polys><DataArray type="Int64" Name="connectivity" format="ascii">{connectivity}</DataArray>
<DataArray type="Int64" Name="offsets" format="ascii">{offsets}</DataArray></Polys>
<CellData Scalars="boundary_id"><DataArray type="Int32" Name="boundary_id" format="ascii">{boundary_ids}</DataArray></CellData>
</Piece></PolyData></VTKFile>
'''
	output = Path(args.output_vtp); output.parent.mkdir(parents=True, exist_ok=True)
	output.write_text(contents, encoding="utf-8")
	polygon_area = 0.5*n*args.tube_radius_m**2*math.sin(2.0*math.pi/n)
	print(f"bent_pipe_surface: PASS points={len(points)} triangles={len(triangles)} "
		f"polygon_volume_m3={polygon_area*args.bend_radius_m*angle:.17g}")
	return 0


if __name__ == "__main__":
	try:
		raise SystemExit(main())
	except ValueError as error:
		print(f"bent_pipe_surface: ERROR: {error}")
		raise SystemExit(2)
