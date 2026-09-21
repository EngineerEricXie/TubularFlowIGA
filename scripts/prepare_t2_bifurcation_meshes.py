#!/usr/bin/env python3
"""Build reproducible coupled surface/volume refinements for the native T2 Y case."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess


LEVELS = ((0.006, 0.004), (0.004, 0.002), (0.002, 0.001), (0.0015, 0.00075))


def run(command):
	completed = subprocess.run(command, text=True, capture_output=True)
	if completed.returncode != 0:
		raise RuntimeError((completed.stderr or completed.stdout).strip())
	return completed.stdout.strip()


def relative_change(value, reference):
	return abs(value-reference)/abs(reference)


def main():
	root = Path(__file__).resolve().parents[1]
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("output_directory")
	parser.add_argument("--maximum-level", type=int, default=len(LEVELS))
	parser.add_argument("--generator",
		default=str(root/"scripts"/"generate_y_pipe_surface.py"))
	parser.add_argument("--converter",
		default=str(root/"scripts"/"surface_to_fem_volume.py"))
	parser.add_argument("--surface-preflight",
		default=str(root/"solvers"/"cpu"/"surface_fem_preflight"))
	args = parser.parse_args()
	if not 1 <= args.maximum_level <= len(LEVELS):
		raise ValueError(f"maximum level must be in [1,{len(LEVELS)}]")
	output = Path(args.output_directory).resolve()
	output.mkdir(parents=True, exist_ok=True)
	levels = []
	for index, (surface_target, volume_target) in enumerate(LEVELS[:args.maximum_level]):
		name = f"level-{index+1}"
		directory = output/name
		directory.mkdir(parents=True, exist_ok=True)
		surface = directory/"surface.vtp"
		mesh = directory/"volume.msh"
		manifest = directory/"volume.json"
		run(["python3", str(Path(args.generator).resolve()), str(surface),
			"--target-size-m", str(surface_target)])
		run(["python3", str(Path(args.converter).resolve()), str(surface), str(mesh),
			"--manifest", str(manifest), "--surface-preflight",
			str(Path(args.surface_preflight).resolve()), "--target-size-m", str(volume_target)])
		volume = json.loads(manifest.read_text(encoding="utf-8"))
		mesh_data = volume["volume_mesh"]
		preflight = volume["surface_preflight"]
		if set(mesh_data["boundary_labels"]) != {"0", "1", "2", "3"}:
			raise RuntimeError(f"{name} boundary labels differ from wall/inlet/outlet contract")
		if not mesh_data["boundary_matches_tetrahedra"] or not volume["gates"]["passed"]:
			raise RuntimeError(f"{name} failed exact-boundary or volume-mesh gates")
		levels.append({
			"name": name,
			"surface_target_size_m": surface_target,
			"volume_target_size_m": volume_target,
			"surface": str(surface.relative_to(output)),
			"volume_mesh": str(mesh.relative_to(output)),
			"volume_manifest": str(manifest.relative_to(output)),
			"surface_sha256": hashlib.sha256(surface.read_bytes()).hexdigest(),
			"mesh_sha256": mesh_data["sha256"],
			"surface_vertices": preflight["canonical_vertices"],
			"surface_triangles": preflight["triangles"],
			"surface_area_m2": preflight["area_m2"],
			"enclosed_volume_m3": preflight["volume_m3"],
			"volume_nodes": mesh_data["nodes"],
			"tetrahedra": mesh_data["elements"],
			"minimum_scaled_jacobian": mesh_data["minimum_scaled_jacobian"],
			"boundary_triangles": mesh_data["boundary_labels"]})
	for previous, current in zip(levels, levels[1:]):
		if not (current["surface_triangles"] > previous["surface_triangles"]
				and current["tetrahedra"] > previous["tetrahedra"]):
			raise RuntimeError("surface and volume refinement must both increase resolution")
	finest = levels[-1]
	for level in levels:
		level["geometry_change_to_finest_discrete_reference"] = {
			"surface_area_relative_change": relative_change(
				level["surface_area_m2"], finest["surface_area_m2"]),
			"enclosed_volume_relative_change": relative_change(
				level["enclosed_volume_m3"], finest["enclosed_volume_m3"])}
	manifest = {
		"schema_version": 1,
		"kind": "t2_bifurcation_coupled_mesh_preparation_not_solver_evidence",
		"case_id": "native_y_bifurcation",
		"geometry_policy": {
			"classification": "coupled_surface_and_volume_refinement",
			"reference_warning": "finest discrete surface is not an exact smooth-geometry error reference",
			"wall_label": 0, "inlet_label": 1, "upper_outlet_label": 2,
			"lower_outlet_label": 3},
		"levels": levels}
	(output/"mesh-series.json").write_text(
		json.dumps(manifest, indent=2, sort_keys=True)+"\n", encoding="utf-8")
	print("prepare_t2_bifurcation_meshes: PASS "
		f"levels={len(levels)} tetrahedra={','.join(str(level['tetrahedra']) for level in levels)}")
	return 0


if __name__ == "__main__":
	try:
		raise SystemExit(main())
	except (KeyError, OSError, RuntimeError, ValueError, ZeroDivisionError) as error:
		print(f"prepare_t2_bifurcation_meshes: ERROR: {error}")
		raise SystemExit(2)
