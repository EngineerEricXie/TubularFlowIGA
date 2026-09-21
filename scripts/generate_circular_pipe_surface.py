#!/usr/bin/env python3
"""Generate a closed labelled triangulated circular-pipe VTP surface."""

import argparse
import math
from pathlib import Path


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("output_vtp")
	parser.add_argument("--length-m", type=float, required=True)
	parser.add_argument("--radius-m", type=float, required=True)
	parser.add_argument("--circumferential-segments", type=int, required=True)
	parser.add_argument("--axial-segments", type=int, required=True)
	parser.add_argument("--radial-segments", type=int, default=1)
	args = parser.parse_args()
	if not math.isfinite(args.length_m) or args.length_m <= 0.0:
		raise ValueError("length must be finite and positive")
	if not math.isfinite(args.radius_m) or args.radius_m <= 0.0:
		raise ValueError("radius must be finite and positive")
	if args.circumferential_segments < 8 or args.axial_segments < 1 \
			or args.radial_segments < 1:
		raise ValueError("pipe requires at least 8 circumferential and positive axial/radial segments")
	n = args.circumferential_segments
	m = args.axial_segments
	points = []
	for axial in range(m+1):
		x = args.length_m*axial/m
		for circumferential in range(n):
			theta = 2.0*math.pi*circumferential/n
			points.append((x, args.radius_m*math.cos(theta),
				args.radius_m*math.sin(theta)))
	triangles = []
	labels = []
	for axial in range(m):
		for circumferential in range(n):
			next_circumferential = (circumferential+1) % n
			a = axial*n+circumferential
			b = (axial+1)*n+circumferential
			c = (axial+1)*n+next_circumferential
			d = axial*n+next_circumferential
			triangles.extend(((a, c, b), (a, d, c)))
			labels.extend((0, 0))
	def cap(x, outer, label, positive):
		center = len(points)
		points.append((x, 0.0, 0.0))
		rings = []
		for radial in range(1, args.radial_segments):
			ring = []
			radius = args.radius_m*radial/args.radial_segments
			for circumferential in range(n):
				theta = 2.0*math.pi*circumferential/n
				ring.append(len(points))
				points.append((x, radius*math.cos(theta), radius*math.sin(theta)))
			rings.append(ring)
		rings.append(outer)
		first = rings[0]
		for circumferential in range(n):
			next_circumferential = (circumferential+1) % n
			triangles.append((center, first[circumferential], first[next_circumferential])
				if positive else (center, first[next_circumferential], first[circumferential]))
			labels.append(label)
		for inner, outer_ring in zip(rings, rings[1:]):
			for circumferential in range(n):
				next_circumferential = (circumferential+1) % n
				a = inner[circumferential]
				b = outer_ring[circumferential]
				c = outer_ring[next_circumferential]
				d = inner[next_circumferential]
				triangles.extend(((a, b, c), (a, c, d)) if positive
					else ((a, c, b), (a, d, c)))
				labels.extend((label, label))
	cap(0.0, list(range(n)), 1, False)
	cap(args.length_m, list(range(m*n, (m+1)*n)), 2, True)
	connectivity = " ".join(str(node) for triangle in triangles for node in triangle)
	offsets = " ".join(str(3*(index+1)) for index in range(len(triangles)))
	coordinates = " ".join(f"{value:.17g}" for point in points for value in point)
	boundary_ids = " ".join(str(value) for value in labels)
	contents = f'''<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">
  <PolyData>
    <Piece NumberOfPoints="{len(points)}" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="{len(triangles)}">
      <Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">{coordinates}</DataArray></Points>
      <Polys>
        <DataArray type="Int64" Name="connectivity" format="ascii">{connectivity}</DataArray>
        <DataArray type="Int64" Name="offsets" format="ascii">{offsets}</DataArray>
      </Polys>
      <CellData Scalars="boundary_id">
        <DataArray type="Int32" Name="boundary_id" format="ascii">{boundary_ids}</DataArray>
      </CellData>
    </Piece>
  </PolyData>
</VTKFile>
'''
	output = Path(args.output_vtp)
	output.parent.mkdir(parents=True, exist_ok=True)
	output.write_text(contents, encoding="utf-8")
	polygon_area = 0.5*n*args.radius_m**2*math.sin(2.0*math.pi/n)
	sagitta = args.radius_m*(1.0-math.cos(math.pi/n))
	print(f"circular_pipe_surface: PASS points={len(points)} triangles={len(triangles)} "
		f"polygon_volume_m3={polygon_area*args.length_m:.17g} "
		f"maximum_boundary_distance_m={sagitta:.17g}")
	return 0


if __name__ == "__main__":
	try:
		raise SystemExit(main())
	except ValueError as error:
		print(f"circular_pipe_surface: ERROR: {error}")
		raise SystemExit(2)
