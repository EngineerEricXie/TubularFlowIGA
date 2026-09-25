#!/usr/bin/env python3
"""Render measured dual-tree geometry and wall displacement; display scale is explicit."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


def wall_points(output, name):
	reader = vtk.vtkXMLPUnstructuredGridReader()
	reader.SetFileName(str(output / "fields" / name / "snapshot.pvtu"))
	reader.Update()
	grid = reader.GetOutput()
	if grid.GetNumberOfPoints() == 0:
		raise ValueError(f"empty {name} wall field")
	position = vtk_to_numpy(grid.GetPointData().GetArray("reference_position_m"))
	displacement = vtk_to_numpy(grid.GetPointData().GetArray("displacement_m"))
	if position.shape != displacement.shape or position.shape[1] != 3:
		raise ValueError(f"invalid {name} wall field")
	return position, displacement


def cube_edges(axis, origin, size):
	for i in (0, 1):
		for j in (0, 1):
			for direction in range(3):
				p = origin.copy()
				q = origin.copy()
				fixed = [k for k in range(3) if k != direction]
				p[fixed[0]] += i * size[fixed[0]]
				p[fixed[1]] += j * size[fixed[1]]
				q[:] = p
				q[direction] += size[direction]
				axis.plot(*np.array([p, q]).T, color="0.75", lw=0.6)


def plot_tree(axis, tree, color):
	nodes = np.asarray(tree["nodes_m"]) * 1000
	for start, stop, _radius in tree["segments"]:
		axis.plot(*nodes[[start, stop]].T, color=color, lw=2.2)
	axis.scatter(*nodes[tree["terminals"]].T, color=color, s=16)


def plot_tissue_flow(axis, output):
	reader = vtk.vtkXMLPUnstructuredGridReader()
	reader.SetFileName(str(output / "fields/tissue_fsi/snapshot.pvtu"))
	reader.Update()
	grid = reader.GetOutput()
	flux = vtk_to_numpy(grid.GetCellData().GetArray("darcy_rt0_centroid_flux_m_s"))
	if len(flux) != grid.GetNumberOfCells() or not np.isfinite(flux).all():
		raise ValueError("invalid Darcy RT0 field")
	magnitude = np.linalg.norm(flux, axis=1)
	indices = np.flatnonzero(magnitude > np.percentile(magnitude, 85))
	step = max(1, len(indices) // 28)
	for index in indices[::step][:28]:
		bounds = grid.GetCell(int(index)).GetBounds()
		center = np.array([(bounds[0]+bounds[1])/2,
			(bounds[2]+bounds[3])/2, (bounds[4]+bounds[5])/2]) * 1000
		direction = flux[index] / magnitude[index]
		axis.quiver(*center, *(direction*1.2), color="#32834c",
			arrow_length_ratio=.25, linewidth=.8)


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("output_directory", type=Path)
	parser.add_argument("--case", type=Path,
		default=Path("cases/idealized_cube_dual_tree.json"))
	parser.add_argument("--magnification", type=float, default=1000.0)
	parser.add_argument("--output", type=Path)
	args = parser.parse_args()
	if not np.isfinite(args.magnification) or args.magnification <= 0:
		raise ValueError("magnification must be positive and finite")
	output = args.output_directory.resolve()
	case = json.loads(args.case.read_text(encoding="utf-8"))
	summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
	if summary["case_sha256"] != __import__("hashlib").sha256(args.case.read_bytes()).hexdigest():
		raise ValueError("case hash differs from solved case")
	if not summary["vessel_wall_fsi"]:
		raise ValueError("wall FSI result is absent")
	walls = {name: wall_points(output, name + "_wall") for name in ("arterial", "venous")}
	fig = plt.figure(figsize=(13, 6.5), layout="constrained")
	left = fig.add_subplot(121, projection="3d")
	right = fig.add_subplot(122)
	origin = np.asarray(case["cube_origin_m"]) * 1000
	size = np.asarray(case["cube_size_m"]) * 1000
	cube_edges(left, origin, size)
	left.set(xlim=(origin[0], origin[0] + size[0]),
		ylim=(origin[1], origin[1] + size[1]),
		zlim=(origin[2], origin[2] + size[2]),
		xlabel="x (mm)", ylabel="y (mm)", zlabel="z (mm)")
	left.set_box_aspect(size.copy())
	left.view_init(elev=22, azim=-64)
	plot_tree(left, case["arterial_tree"], "#c43c39")
	plot_tree(left, case["venous_tree"], "#315eaf")
	plot_tissue_flow(left, output)
	left.set_title("Centerline/radius trees and sampled Darcy RT0 directions\n"
		"red: arterial; blue: venous; green arrows: direction only")
	for name, color in (("arterial", "#c43c39"), ("venous", "#315eaf")):
		position, displacement = walls[name]
		if not np.isfinite(displacement).all() or np.max(np.linalg.norm(displacement, axis=1)) <= 0:
			raise ValueError(f"{name} wall did not move")
		step = max(1, len(position) // 1400)
		sample = slice(None, None, step)
		ref = position[sample] * 1000
		moved = (position[sample] + args.magnification * displacement[sample]) * 1000
		right.scatter(ref[:, 0], ref[:, 1], color="0.75", alpha=0.25, s=1)
		right.scatter(moved[:, 0], moved[:, 1], color=color, alpha=0.65, s=2)
	right.set(xlim=(origin[0], origin[0] + size[0]),
		ylim=(origin[1], origin[1] + size[1]),
		xlabel="x (mm)", ylabel="y (mm)")
	right.set_aspect("equal")
	maximum = {name: np.max(np.linalg.norm(field[1], axis=1))
		for name, field in walls.items()}
	right.set_title(f"Wall FEM displacement, x–y projection (grey: reference)\n"
		f"shown ×{args.magnification:g}; actual max: A {maximum['arterial']*1e6:.3g} µm, "
		f"V {maximum['venous']*1e6:.3g} µm")
	flow = summary["flow"]
	fig.suptitle("Idealized artery–Darcy–vein coupling with quasi-steady wall feedback\n"
		f"conservative flow transfer: {flow['fsi_artery_to_tissue_m3_s']:.3e} → "
		f"{flow['fsi_tissue_to_vein_m3_s']:.3e} m³/s", fontsize=11)
	target = args.output.resolve() if args.output else output / "overview.png"
	if target.exists():
		raise ValueError(f"refusing to overwrite {target}")
	fig.savefig(target, dpi=180)
	plt.close(fig)
	print(target)


if __name__ == "__main__":
	main()
