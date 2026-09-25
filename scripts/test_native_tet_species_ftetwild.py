#!/usr/bin/env python3
"""Rebuild a labelled fTetWild pipe and exercise native FEM/0D species flow."""

import argparse
import json
import math
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path


def checked(command, timeout=120):
	result = subprocess.run(command, capture_output=True, text=True,
		timeout=timeout, check=False)
	if result.returncode:
		raise RuntimeError(f"command failed: {' '.join(map(str, command))}\n"
			f"{result.stdout}\n{result.stderr}")
	return result.stdout


def species_snapshot(directory, step, expected, moving, expected_nodes):
	root = directory/f"step_{step}"
	index = ET.parse(root/"snapshot.pvtu").getroot()
	fields = {item.attrib.get("Name") for item in index.findall(".//PPointData/PDataArray")}
	if not {"concentration_mol_m3", "reference_position_m", "displacement_m"}.issubset(fields):
		raise RuntimeError("native species PVTU field schema is incomplete")
	concentrations = {}
	cells = set()
	for item in index.findall(".//Piece"):
		piece = ET.parse(root/item.attrib["Source"]).getroot().find(".//Piece")
		if piece is None:
			raise RuntimeError("native species VTU piece is absent")
		def values(path, kind=float):
			array = piece.find(path)
			if array is None:
				raise RuntimeError("native species VTU array is absent")
			return [kind(value) for value in (array.text or "").split()]
		ids = values("PointData/DataArray[@Name='GlobalPointIds']", int)
		current = values("Points/DataArray")
		reference = values("PointData/DataArray[@Name='reference_position_m']")
		displacement = values("PointData/DataArray[@Name='displacement_m']")
		field = values("PointData/DataArray[@Name='concentration_mol_m3']")
		cell_ids = values("CellData/DataArray[@Name='GlobalCellIds']", int)
		cell_types = values("Cells/DataArray[@Name='types']", int)
		if (len(current) != 3*len(ids) or len(reference) != len(current)
				or len(displacement) != len(current) or len(field) != len(ids)
				or len(cell_ids) != len(cell_types) or any(kind != 10 for kind in cell_types)):
			raise RuntimeError("native species VTU P1 topology or fields are incomplete")
		for cell_id in cell_ids:
			if cell_id in cells:
				raise RuntimeError("native species VTU cell is owned twice")
			cells.add(cell_id)
		for local, global_id in enumerate(ids):
			if global_id >= expected_nodes or global_id < 0:
				raise RuntimeError("native species VTU point id is invalid")
			if not math.isfinite(field[local]) or (expected is not None
					and abs(field[local]-expected[global_id]) > 1e-10):
				raise RuntimeError("native species VTU concentration differs from accepted field")
			for axis in range(3):
				coordinate_index = 3*local+axis
				shift = -1e-5 if moving and step == 2 and axis == 0 else 0.0
				if (not all(math.isfinite(value[coordinate_index]) for value in
						(current, reference, displacement))
						or abs(current[coordinate_index]-reference[coordinate_index]
							-displacement[coordinate_index]) > 1e-10
						or abs(displacement[coordinate_index]-shift) > 1e-10):
					raise RuntimeError("native species VTU ALE coordinate differs")
			if global_id in concentrations and abs(concentrations[global_id]-field[local]) > 1e-12:
				raise RuntimeError("native species shared VTU point differs")
			concentrations[global_id] = field[local]
	if len(concentrations) != expected_nodes or not cells:
		raise RuntimeError("native species VTU lost vertices or cells")
	return len(cells)


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--ftetwild", type=Path, required=True)
	parser.add_argument("--binary", type=Path, required=True)
	parser.add_argument("--ranks", type=int, default=2)
	parser.add_argument("--compare-ranks", action="store_true",
		help="compare full native concentration fields on one mesh at 1/2/4 MPI ranks")
	args = parser.parse_args()
	if not 1 <= args.ranks <= 16:
		raise ValueError("MPI ranks must be between 1 and 16")
	root = Path(__file__).resolve().parents[1]
	ftetwild = args.ftetwild.resolve(strict=True)
	binary = args.binary.resolve(strict=True)
	with tempfile.TemporaryDirectory(prefix="iga-ftetwild-native-species-") as temporary:
		work = Path(temporary)
		surface, mesh = work/"pipe.vtp", work/"pipe.msh"
		manifest = work/"mesh_manifest.json"
		checked([sys.executable, str(root/"scripts/generate_circular_pipe_surface.py"),
			str(surface), "--length-m", "0.03", "--radius-m", "0.005",
			"--circumferential-segments", "8", "--axial-segments", "2"])
		checked([sys.executable, str(root/"scripts/ftetwild_to_fem_volume.py"),
			str(surface), str(mesh), "--manifest", str(manifest),
			"--ftetwild", str(ftetwild), "--target-size-m", "0.007",
			"--envelope-m", "0.0001", "--stop-energy", "12",
			"--max-optimization-passes", "40", "--max-threads", "2"])
		data = json.loads(manifest.read_text())
		volume = data["volume_mesh"]
		if (not data["gates"]["passed"]
				or data["mesher"]["name"] != "fTetWild"
				or not data["mesher"]["executable_sha256"]
				or volume["element_type"] != "tetrahedron_p1"
				or not volume["boundary_matches_tetrahedra"]
				or set(volume["boundary_labels"]) != {"0", "1", "2"}
				or volume["elements"] <= 24
				or volume["minimum_scaled_jacobian"] < 1e-3):
			raise RuntimeError("fTetWild pipe did not pass labelled native tetra gates")
		def solve(variant, ranks, sharp_front=False, moving=False):
			command = ["mpiexec", "-np", str(ranks), str(binary),
				"--ftetwild-pipe", str(mesh)]
			if variant != "baseline":
				command.append(f"--transport-variant={variant}")
			if sharp_front:
				command.extend(("--sharp-front", "--monotone-species"))
			if moving:
				command.append("--moving-two-step")
			output_dir = (work/f"species_{variant}_{ranks}_{'moving' if moving else 'fixed'}") if (
				variant == "baseline" and (not sharp_front or moving)) else None
			if output_dir is not None:
				command.append(f"--species-output-dir={output_dir}")
			if args.compare_ranks or output_dir is not None:
				command.append("--emit-concentration")
			output = checked(command)
			if "fTetWild native source/FEM/RCR species graph passed" not in output:
				raise RuntimeError(f"native FEM species graph did not report acceptance: {output}")
			values = {}
			steps = re.search(r"\bsteps=([0-9]+)", output)
			if not steps or int(steps.group(1)) != (2 if moving else 1):
				raise RuntimeError(f"native FEM species graph accepted-step count differs: {output}")
			for name in ("native_final_mol", "native_minimum_mol_m3",
					"native_rms_delta_mol_m3",
					"source_mol", "global_residual_mol"):
				match = re.search(rf"{name}=([-+0-9.eE]+)", output)
				if not match or not math.isfinite(float(match.group(1))):
					raise RuntimeError(f"native FEM species graph {name} is invalid: {output}")
				values[name] = float(match.group(1))
			if abs(values["global_residual_mol"]) > 1e-12:
				raise RuntimeError(f"native FEM species graph balance is invalid: {output}")
			if sharp_front and (values["native_minimum_mol_m3"] < -1e-12
					or abs(values["source_mol"]) > 1e-12
					or values["native_final_mol"] <= 0):
				raise RuntimeError(f"native FEM sharp-front acceptance is invalid: {output}")
			if args.compare_ranks or output_dir is not None:
				match = re.search(r"^native_concentration_mol_m3=([^\n]+)$",
					output, re.MULTILINE)
				if not match:
					raise RuntimeError("native FEM concentration field is missing")
				values["concentration"] = [float(value) for value in
					match.group(1).split(",")]
				if (len(values["concentration"]) < 24
						or not all(math.isfinite(value) for value in values["concentration"])):
					raise RuntimeError("native FEM concentration field is invalid")
			if output_dir is not None:
				for step in range(1, 3 if moving else 2):
					cell_count = species_snapshot(output_dir, step,
						values["concentration"] if step == (2 if moving else 1) else None,
						moving, len(values["concentration"]))
					if cell_count != volume["elements"]:
						raise RuntimeError("native species VTU lost tetrahedra")
			return values

		def require_unstabilized_front_rejected(ranks):
			command = ["mpiexec", "-np", str(ranks), str(binary),
				"--ftetwild-pipe", str(mesh), "--sharp-front"]
			result = subprocess.run(command, capture_output=True, text=True,
				timeout=120, check=False)
			if result.returncode == 0 or "concentration is invalid" not in (
					result.stdout + result.stderr):
				raise RuntimeError("unstabilized fTetWild front did not reject negative concentration")

		rank_results = {}
		for ranks in ((1, 2, 4) if args.compare_ranks else (args.ranks,)):
			results = {variant: solve(variant, ranks) for variant in
				("baseline", "diffusion", "source", "both")}
			baseline, diffusion, source, both = (results[name] for name in
				("baseline", "diffusion", "source", "both"))
			if (abs(baseline["source_mol"]) > 1e-12
					or abs(diffusion["source_mol"]) > 1e-12
					or source["source_mol"] <= 1e-10
					or both["source_mol"] <= 1e-10
					or source["native_final_mol"]-baseline["native_final_mol"]
						<= 0.25*source["source_mol"]
					or both["native_final_mol"]-diffusion["native_final_mol"]
						<= 0.25*both["source_mol"]):
				raise RuntimeError("fTetWild isolated source variant did not change native inventory")
			diffusion_rms_change = abs(diffusion["native_rms_delta_mol_m3"]
				- baseline["native_rms_delta_mol_m3"])
			if diffusion_rms_change <= 1e-8:
				raise RuntimeError(f"fTetWild isolated diffusion effect is unresolved: {diffusion_rms_change}")
			results["sharp_front"] = solve("baseline", ranks, sharp_front=True)
			results["moving_front"] = solve("baseline", ranks,
				sharp_front=True, moving=True)
			require_unstabilized_front_rejected(ranks)
			rank_results[ranks] = results
		maximum_relative_l2 = 0.0
		if args.compare_ranks:
			for variant in ("baseline", "diffusion", "source", "both",
					"sharp_front", "moving_front"):
				reference = rank_results[1][variant]["concentration"]
				for ranks in (2, 4):
					candidate = rank_results[ranks][variant]["concentration"]
					if len(candidate) != len(reference):
						raise RuntimeError("fTetWild native field sizes differ across ranks")
					relative_l2 = math.sqrt(sum((a-b)**2 for a, b in
						zip(reference, candidate))/sum(a*a for a in reference))
					maximum_relative_l2 = max(maximum_relative_l2, relative_l2)
					gate = 1e-8 if variant.endswith("front") else 1e-10
					if relative_l2 > gate:
						raise RuntimeError(f"fTetWild {variant} rank {ranks} relative L2 {relative_l2} exceeds {gate}")
		results = rank_results[args.ranks if not args.compare_ranks else 1]
		baseline, diffusion, source, both = (results[name] for name in
			("baseline", "diffusion", "source", "both"))
		diffusion_rms_change = abs(diffusion["native_rms_delta_mol_m3"]
			- baseline["native_rms_delta_mol_m3"])
		print(f"native_tet_species_ftetwild: PASS ranks={'1,2,4' if args.compare_ranks else args.ranks} "
			f"tetrahedra={volume['elements']} "
			f"minimum_scaled_jacobian={volume['minimum_scaled_jacobian']:.6g} "
			f"baseline_native_mol={baseline['native_final_mol']:.6g} "
			f"diffusion_rms_change_mol_m3={diffusion_rms_change:.6g} "
			f"source_native_mol={source['native_final_mol']:.6g} "
			f"source_mol={source['source_mol']:.6g} "
			f"both_native_mol={both['native_final_mol']:.6g} "
			f"sharp_front_minimum_mol_m3={results['sharp_front']['native_minimum_mol_m3']:.6g} "
			f"moving_front_minimum_mol_m3={results['moving_front']['native_minimum_mol_m3']:.6g} "
			f"moving_front_native_mol={results['moving_front']['native_final_mol']:.6g} "
			f"unstabilized_front_rejected=1 "
			f"maximum_global_residual_mol={max(abs(value['global_residual_mol']) for group in rank_results.values() for value in group.values()):.6g} "
			f"maximum_rank_field_relative_l2={maximum_relative_l2:.6g} "
			f"mesh_sha256={volume['sha256']}")


if __name__ == "__main__":
	main()
