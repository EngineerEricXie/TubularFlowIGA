#!/usr/bin/env python3
"""Build all frozen T2 circular-pipe FEM mesh levels and preparation evidence."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

from t2_contract import exact_contract_sha256, mesh_contract_sha256


def run(command):
	completed = subprocess.run(command, text=True, capture_output=True)
	if completed.returncode != 0:
		raise RuntimeError((completed.stderr or completed.stdout).strip())
	return completed.stdout.strip()


def main():
	root = Path(__file__).resolve().parents[1]
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("output_directory")
	parser.add_argument("--contract", default=str(root/"benchmarks"/"t2_fixed_flow_contract.json"))
	parser.add_argument("--generator", default=str(root/"scripts"/"generate_circular_pipe_surface.py"))
	parser.add_argument("--converter", default=str(root/"scripts"/"surface_to_fem_volume.py"))
	parser.add_argument("--surface-preflight",
		default=str(root/"solvers"/"cpu"/"surface_fem_preflight"))
	args = parser.parse_args()
	contract_path = Path(args.contract).resolve()
	contract_bytes = contract_path.read_bytes()
	contract = json.loads(contract_bytes)
	length = float(contract["geometry"]["length_m"])
	radius = float(contract["geometry"]["radius_m"])
	output = Path(args.output_directory).resolve()
	output.mkdir(parents=True, exist_ok=True)
	levels = []
	for index, target in enumerate(
			contract["discretization_contract"]["mesh_levels_target_size_m"]):
		name = f"level-{index+1}"
		directory = output/name
		directory.mkdir(parents=True, exist_ok=True)
		# Use half the volume target around the circumference so polygonal
		# geometry error does not dominate the PDE convergence series.
		circumferential = max(8, math.ceil(4.0*math.pi*radius/target))
		axial = max(1, math.ceil(length/target))
		radial = max(1, math.ceil(radius/target))
		surface = directory/"surface.vtp"
		mesh = directory/"volume.msh"
		manifest = directory/"volume.json"
		run(["python3", str(Path(args.generator).resolve()), str(surface),
			"--length-m", str(length), "--radius-m", str(radius),
			"--circumferential-segments", str(circumferential),
			"--axial-segments", str(axial), "--radial-segments", str(radial)])
		run(["python3", str(Path(args.converter).resolve()), str(surface), str(mesh),
			"--manifest", str(manifest), "--surface-preflight",
			str(Path(args.surface_preflight).resolve()), "--target-size-m", str(target)])
		volume = json.loads(manifest.read_text(encoding="utf-8"))
		boundary_counts = volume["volume_mesh"]["boundary_labels"]
		if set(boundary_counts) != {"0", "1", "2"}:
			raise RuntimeError(f"{name} boundary labels differ from wall/inlet/outlet contract")
		polygon_area = 0.5*circumferential*radius**2 \
			*math.sin(2.0*math.pi/circumferential)
		polygon_volume = polygon_area*length
		exact_volume = math.pi*radius**2*length
		if not math.isclose(volume["surface_preflight"]["volume_m3"], polygon_volume,
				rel_tol=5e-13):
			raise RuntimeError(f"{name} preflight volume does not match polygonal cylinder")
		levels.append({
			"name": name, "target_size_m": target,
			"circumferential_segments": circumferential,
			"axial_segments": axial, "radial_segments": radial,
			"surface": str(surface.relative_to(output)),
			"volume_mesh": str(mesh.relative_to(output)),
			"volume_manifest": str(manifest.relative_to(output)),
			"surface_sha256": hashlib.sha256(surface.read_bytes()).hexdigest(),
			"mesh_sha256": volume["volume_mesh"]["sha256"],
			"geometry_error": {
				"maximum_boundary_distance_m": radius*(1.0-math.cos(math.pi/circumferential)),
				"relative_volume_error": (exact_volume-polygon_volume)/exact_volume,
				"boundary_label_counts_match": bool(
					volume["volume_mesh"]["boundary_matches_tetrahedra"])},
			"nodes": volume["volume_mesh"]["nodes"],
			"tetrahedra": volume["volume_mesh"]["elements"],
			"minimum_scaled_jacobian": volume["volume_mesh"]["minimum_scaled_jacobian"]})
		levels[-1]["boundary_triangles"] = boundary_counts
	for previous, current in zip(levels, levels[1:]):
		for metric in ("maximum_boundary_distance_m", "relative_volume_error"):
			if not current["geometry_error"][metric] < previous["geometry_error"][metric]:
				raise RuntimeError(f"geometry metric {metric} did not decrease under refinement")
	index_manifest = {
		"schema_version": 1,
		"kind": "t2_mesh_preparation_not_solver_evidence",
		"case_id": contract["case_id"],
		"contract_sha256": exact_contract_sha256(contract_bytes),
		"mesh_contract_sha256": mesh_contract_sha256(contract),
		"levels": levels}
	(output/"mesh-series.json").write_text(
		json.dumps(index_manifest, indent=2, sort_keys=True)+"\n", encoding="utf-8")
	print(f"prepare_t2_poiseuille_meshes: PASS levels={len(levels)} "
		f"tetrahedra={','.join(str(level['tetrahedra']) for level in levels)}")
	return 0


if __name__ == "__main__":
	try:
		raise SystemExit(main())
	except (KeyError, OSError, RuntimeError, ValueError) as error:
		print(f"prepare_t2_poiseuille_meshes: ERROR: {error}")
		raise SystemExit(2)
