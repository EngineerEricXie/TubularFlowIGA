#!/usr/bin/env python3
"""Write a structured tetrahedral box mesh (Gmsh 4.1) for solver benchmarks.

Each of the N x N x N hexahedra of a cube with side LENGTH is split into six
tetrahedra around its main diagonal, so the mesh has 6*N^3 tetrahedra and
(N+1)^3 vertices. Boundary labels: 1 on x=0, 2 on x=LENGTH, 3 elsewhere.
"""

import argparse
from collections import defaultdict
from pathlib import Path

from generate_fsi_channel_fixture import oriented, write_msh


def box_mesh(divisions, length):
	count = divisions+1
	def node(i, j, k):
		return (i*count+j)*count+k
	points = [(length*i/divisions, length*j/divisions, length*k/divisions)
		for i in range(count) for j in range(count) for k in range(count)]
	cells = []
	for i in range(divisions):
		for j in range(divisions):
			for k in range(divisions):
				a, b, c, d = node(i, j, k), node(i+1, j, k), node(i, j+1, k), node(i+1, j+1, k)
				e, f, g, h = node(i, j, k+1), node(i+1, j, k+1), node(i, j+1, k+1), node(i+1, j+1, k+1)
				for vertices in ((a, b, d, h), (a, d, c, h), (a, c, g, h),
						(a, g, e, h), (a, e, f, h), (a, f, b, h)):
					cells.append(oriented(points, vertices))
	uses = defaultdict(list)
	for cell in cells:
		for opposite in range(4):
			face = tuple(cell[local] for local in range(4) if local != opposite)
			uses[tuple(sorted(face))].append(face)
	triangles = []
	for occurrences in uses.values():
		if len(occurrences) != 1:
			continue
		face = occurrences[0]
		x = [points[index][0] for index in face]
		label = 1 if max(x) == 0 else (2 if min(x) == length else 3)
		triangles.append((face, label))
	return points, cells, triangles


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("divisions", type=int, help="hexahedra per side (N)")
	parser.add_argument("output", type=Path, help="output .msh file")
	parser.add_argument("--length", type=float, default=0.01, help="cube side in metres")
	args = parser.parse_args()
	if args.divisions < 1:
		parser.error("divisions must be positive")
	points, cells, triangles = box_mesh(args.divisions, args.length)
	write_msh(args.output, points, cells, triangles)
	print(f"{args.output}: {len(points)} vertices, {len(cells)} tetrahedra, {len(triangles)} boundary triangles")


if __name__ == "__main__":
	main()
