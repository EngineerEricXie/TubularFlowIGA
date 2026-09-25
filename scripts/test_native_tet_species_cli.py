#!/usr/bin/env python3
"""Exercise the external-mesh native P1 species CLI without a flow solver."""

import argparse
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path


def invoke(binary, ranks, case, *options, success=True):
	result = subprocess.run(["mpiexec", "-np", str(ranks), str(binary),
		str(case), *map(str, options)], capture_output=True, text=True, timeout=120)
	if success and result.returncode:
		raise RuntimeError(f"native species CLI failed: {result.stdout}\n{result.stderr}")
	if not success and result.returncode == 0:
		raise RuntimeError("native species CLI accepted invalid input or overwrite")
	return result.stdout + result.stderr


def snapshot(directory, step, expected_vertices, expected_cells, shift_x):
	root = directory/f"step_{step}"
	index = ET.parse(root/"snapshot.pvtu").getroot()
	fields = {field.attrib.get("Name") for field in index.findall(".//PPointData/PDataArray")}
	if not {"concentration_mol_m3", "reference_position_m", "displacement_m"}.issubset(fields):
		raise RuntimeError("native species CLI PVTU schema is incomplete")
	points, cells = {}, set()
	for item in index.findall(".//Piece"):
		piece = ET.parse(root/item.attrib["Source"]).getroot().find(".//Piece")
		if piece is None:
			raise RuntimeError("native species CLI VTU piece is missing")
		def values(path, kind=float):
			field = piece.find(path)
			if field is None:
				raise RuntimeError("native species CLI VTU array is missing")
			return [kind(value) for value in (field.text or "").split()]
		ids = values("PointData/DataArray[@Name='GlobalPointIds']", int)
		current = values("Points/DataArray")
		reference = values("PointData/DataArray[@Name='reference_position_m']")
		displacement = values("PointData/DataArray[@Name='displacement_m']")
		concentration = values("PointData/DataArray[@Name='concentration_mol_m3']")
		cell_ids = values("CellData/DataArray[@Name='GlobalCellIds']", int)
		cell_types = values("Cells/DataArray[@Name='types']", int)
		if (len(current) != 3*len(ids) or len(reference) != len(current)
				or len(displacement) != len(current) or len(concentration) != len(ids)
				or len(cell_ids) != len(cell_types) or any(kind != 10 for kind in cell_types)):
			raise RuntimeError("native species CLI VTU fields or P1 cells are incomplete")
		for cell_id in cell_ids:
			if cell_id in cells:
				raise RuntimeError("native species CLI duplicate cell ownership")
			cells.add(cell_id)
		for local, identifier in enumerate(ids):
			if not 0 <= identifier < expected_vertices:
				raise RuntimeError("native species CLI global vertex id is invalid")
			if not math.isfinite(concentration[local]):
				raise RuntimeError("native species CLI concentration is nonfinite")
			for axis in range(3):
				coordinate_index = 3*local+axis
				expected_shift = shift_x if axis == 0 else 0.0
				if (not all(math.isfinite(values[coordinate_index]) for values in
						(current, reference, displacement))
						or abs(current[coordinate_index]-reference[coordinate_index]
							-displacement[coordinate_index]) > 1e-12
						or abs(displacement[coordinate_index]-expected_shift) > 1e-12):
					raise RuntimeError("native species CLI reference/current ALE fields differ")
			if identifier in points and abs(points[identifier]-concentration[local]) > 1e-12:
				raise RuntimeError("native species CLI shared concentration differs")
			points[identifier] = concentration[local]
	if len(points) != expected_vertices or len(cells) != expected_cells:
		raise RuntimeError("native species CLI lost tetrahedra or vertices")
	return [points[index] for index in range(expected_vertices)]


def summary(directory, case, mesh, ranks, steps):
	data = json.loads((directory/"run_summary.json").read_text())
	if (data.get("schema_version") != 1
			or data.get("kind") != "native_tet_species_prescribed_velocity"
			or data.get("species_id") != "tracer"
			or data.get("case_sha256") != hashlib.sha256(case.read_bytes()).hexdigest()
			or data.get("mesh_sha256") != hashlib.sha256(mesh.read_bytes()).hexdigest()
			or data.get("mpi_ranks") != ranks or data.get("accepted_steps") != steps
			or not math.isfinite(data.get("last_inventory_mol", math.nan))
			or not math.isfinite(data.get("last_reaction_sink_mol_s", math.nan))
			or data["last_reaction_sink_mol_s"] < 0
			or not math.isfinite(data.get("last_outward_wall_exchange_mol_s", math.nan))
			or abs(data.get("last_balance_defect_mol_s", math.inf)) > 1e-10):
		raise RuntimeError("native species CLI summary or provenance is invalid")


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--binary", type=Path, required=True)
	parser.add_argument("--ftetwild", type=Path,
		help="also run the external-mesh CLI on a newly generated labelled pipe")
	args = parser.parse_args()
	binary = args.binary.resolve(strict=True)
	mesh = (Path(__file__).resolve().parents[1]/
		"solvers/cpu/tests/data/native_tet_hydraulic_star.msh").resolve(strict=True)
	with tempfile.TemporaryDirectory(prefix="native-species-cli-") as temporary:
		work = Path(temporary)
		base = {
			"schema_version": 1, "mesh_file": str(mesh), "species_id": "tracer",
			"initial_concentration_mol_m3": 2.0, "diffusivity_m2_s": 0.0,
			"source_mol_m3_s": 0.0, "fluid_velocity_m_s": [0.0, 0.0, 0.0],
			"mesh_velocity_m_s": [0.0, 0.0, 0.0],
			"inflow_concentration_by_label_mol_m3": {}, "monotone": False,
			"time": {"dt_s": 0.1, "steps": 2},
		}
		def case_file(name, changes):
			case = {**base, **changes}
			path = work/f"{name}.json"
			path.write_text(json.dumps(case), encoding="utf-8")
			return path
		maximum_restart_field_difference = 0.0
		all_restart_checkpoints_bitwise = True
		wall_fields = {}
		reservoir_fields = {}
		reservoir_amounts = {}
		multi_region_fields = {}
		multi_region_amounts = {}
		for ranks in (1, 2, 4):
			fixed = case_file(f"fixed_{ranks}", {})
			if "native_species_input: PASS" not in invoke(binary, ranks, fixed, "--check-input"):
				raise RuntimeError("native species CLI input preflight did not report success")
			output = work/f"fixed_{ranks}"
			log = invoke(binary, ranks, fixed, "--output-dir", output)
			summary(output, fixed, mesh, ranks, 2)
			if len(re.findall(r"native_species_step step=", log)) != 2:
				raise RuntimeError("native species CLI accepted-step count differs")
			for step in (1, 2):
				field = snapshot(output, step, 5, 4, 0.0)
				if max(abs(value-2.0) for value in field) > 1e-10:
					raise RuntimeError("native species CLI stationary constant field changed")
			invoke(binary, ranks, fixed, "--output-dir", output, success=False)
			moving = case_file(f"moving_{ranks}", {
				"fluid_velocity_m_s": [1.0, 0.0, 0.0],
				"mesh_velocity_m_s": [1.0, 0.0, 0.0], "monotone": True,
			})
			moving_output = work/f"moving_{ranks}"
			invoke(binary, ranks, moving, "--output-dir", moving_output)
			summary(moving_output, moving, mesh, ranks, 2)
			for step in (1, 2):
				field = snapshot(moving_output, step, 5, 4, 0.1*step)
				if max(abs(value-2.0) for value in field) > 1e-10:
					raise RuntimeError("native species CLI moving u=w field changed")
			if ranks == 2:
				moving_split = work/"moving_split_2"
				invoke(binary, ranks, moving, "--output-dir", moving_split,
					"--stop-after-step", "1")
				invoke(binary, ranks, moving, "--resume", moving_split)
				summary(moving_split, moving, mesh, ranks, 2)
				all_restart_checkpoints_bitwise &= (
					(moving_split/"checkpoint_step_2.bin").read_bytes() ==
					(moving_output/"checkpoint_step_2.bin").read_bytes())
				for step in (1, 2):
					restarted = snapshot(moving_split, step, 5, 4, 0.1*step)
					uninterrupted = snapshot(moving_output, step, 5, 4, 0.1*step)
					maximum_restart_field_difference = max(maximum_restart_field_difference,
						max(abs(a-b) for a, b in zip(restarted, uninterrupted)))
			decay = case_file(f"decay_{ranks}", {
				"fluid_velocity_m_s": [1.0, 0.0, 0.0],
				"mesh_velocity_m_s": [1.0, 0.0, 0.0],
				"source_mol_m3_s": 0.5, "first_order_decay_rate_s_inv": 2.0})
			decay_output = work/f"decay_{ranks}"
			invoke(binary, ranks, decay, "--check-input")
			invoke(binary, ranks, decay, "--output-dir", decay_output)
			summary(decay_output, decay, mesh, ranks, 2)
			expected_decay = (2.0+0.1*0.5)/1.2
			for step in (1, 2):
				if step == 2:
					expected_decay = (expected_decay+0.1*0.5)/1.2
				field = snapshot(decay_output, step, 5, 4, 0.1*step)
				if max(abs(value-expected_decay) for value in field) > 1e-10:
					raise RuntimeError("native species implicit decay/source field differs")
			decay_summary = json.loads((decay_output/"run_summary.json").read_text())
			if abs(decay_summary["last_reaction_sink_mol_s"]
					-2.0*decay_summary["last_inventory_mol"]) > 1e-12:
				raise RuntimeError("native species reaction sink budget differs")
			if ranks == 2:
				decay_split = work/"decay_split_2"
				invoke(binary, ranks, decay, "--output-dir", decay_split,
					"--stop-after-step", "1")
				invoke(binary, ranks, decay, "--resume", decay_split)
				summary(decay_split, decay, mesh, ranks, 2)
				a = snapshot(decay_output, 2, 5, 4, 0.2)
				b = snapshot(decay_split, 2, 5, 4, 0.2)
				maximum_restart_field_difference = max(maximum_restart_field_difference,
					max(abs(x-y) for x, y in zip(a, b)))
			wall = case_file(f"wall_{ranks}", {"monotone": True,
				"wall_exchange_by_label": {
				"1": {"transfer_coefficient_m_s": 0.5,
					"external_concentration_mol_m3": 0.0}}})
			wall_output = work/f"wall_{ranks}"
			invoke(binary, ranks, wall, "--check-input")
			invoke(binary, ranks, wall, "--output-dir", wall_output)
			summary(wall_output, wall, mesh, ranks, 2)
			wall_summary = json.loads((wall_output/"run_summary.json").read_text())
			if (wall_summary["last_outward_wall_exchange_mol_s"] <= 0
					or wall_summary["last_inventory_mol"] <= 0):
				raise RuntimeError("native species wall sink did not remove inventory")
			wall_fields[ranks] = snapshot(wall_output, 2, 5, 4, 0.0)
			if (min(wall_fields[ranks]) < 0
					or min(wall_fields[ranks]) >= 2.
					or max(wall_fields[ranks]) > 2.+1e-10):
				raise RuntimeError(f"native species wall sink field is invalid: "
					f"ranks={ranks} field={wall_fields[ranks]}")
			if ranks == 2:
				wall_split = work/"wall_split_2"
				invoke(binary, ranks, wall, "--output-dir", wall_split,
					"--stop-after-step", "1")
				invoke(binary, ranks, wall, "--resume", wall_split)
				summary(wall_split, wall, mesh, ranks, 2)
				restarted = snapshot(wall_split, 2, 5, 4, 0.0)
				maximum_restart_field_difference = max(maximum_restart_field_difference,
					max(abs(a-b) for a, b in zip(wall_fields[ranks], restarted)))
				donor = case_file("wall_donor_2", {"monotone": True,
					"wall_exchange_by_label": {
					"1": {"transfer_coefficient_m_s": 0.5,
						"external_concentration_mol_m3": 3.0}}})
				donor_output = work/"wall_donor_2"
				invoke(binary, ranks, donor, "--output-dir", donor_output)
				summary(donor_output, donor, mesh, ranks, 2)
				donor_summary = json.loads((donor_output/"run_summary.json").read_text())
				if donor_summary["last_outward_wall_exchange_mol_s"] >= 0:
					raise RuntimeError("native species wall donor did not supply inventory")
			reservoir = case_file(f"reservoir_{ranks}", {"monotone": True,
				"finite_wall_reservoir": {"boundary_label": 1,
					"transfer_coefficient_m_s": 0.5, "volume_m3": 1.0,
					"initial_concentration_mol_m3": 0.0}})
			reservoir_output = work/f"reservoir_{ranks}"
			invoke(binary, ranks, reservoir, "--check-input")
			invoke(binary, ranks, reservoir, "--output-dir", reservoir_output)
			summary(reservoir_output, reservoir, mesh, ranks, 2)
			reservoir_summary = json.loads((reservoir_output/"run_summary.json").read_text())
			reservoir_amounts[ranks] = reservoir_summary["last_tissue_amount_mol"]
			if (not 0 < reservoir_amounts[ranks] < 1.0
					or abs(reservoir_summary["last_tissue_concentration_mol_m3"]
						-reservoir_amounts[ranks]) > 1e-12
					or abs(reservoir_summary["last_combined_balance_defect_mol_s"]) > 1e-10):
				raise RuntimeError("native species finite reservoir budget is invalid")
			reservoir_fields[ranks] = snapshot(reservoir_output, 2, 5, 4, 0.0)
			for step in (1, 2):
				state = json.loads((reservoir_output/"run_state.json").read_text())
				if not (reservoir_output/f"tissue_step_{step}.bin").is_file():
					raise RuntimeError("native species reservoir checkpoint is missing")
			if step == 2 and state.get("tissue_checkpoint_sha256") != hashlib.sha256(
					(reservoir_output/"tissue_step_2.bin").read_bytes()).hexdigest():
				raise RuntimeError("native species reservoir checkpoint hash differs")
			if ranks == 2:
				reservoir_donor = case_file("reservoir_donor_2", {"monotone": True,
					"finite_wall_reservoir": {"boundary_label": 1,
						"transfer_coefficient_m_s": 0.5, "volume_m3": 1.0,
						"initial_concentration_mol_m3": 3.0}})
				donor_output = work/"reservoir_donor_2"
				invoke(binary, ranks, reservoir_donor, "--output-dir", donor_output)
				summary(donor_output, reservoir_donor, mesh, ranks, 2)
				donor_summary = json.loads((donor_output/"run_summary.json").read_text())
				if (donor_summary["last_tissue_amount_mol"] >= 3.0
						or donor_summary["last_outward_wall_exchange_mol_s"] >= 0
						or abs(donor_summary["last_combined_balance_defect_mol_s"]) > 1e-10):
					raise RuntimeError("native species finite reservoir donor is invalid")
				reservoir_split = work/"reservoir_split_2"
				invoke(binary, ranks, reservoir, "--output-dir", reservoir_split,
					"--stop-after-step", "1")
				for kind in ("missing_tissue", "corrupt_tissue", "orphan_tissue", "missing_hash"):
					bad = work/f"bad_reservoir_{kind}"
					shutil.copytree(reservoir_split, bad)
					if kind == "missing_tissue":
						(bad/"tissue_step_1.bin").unlink()
					elif kind == "corrupt_tissue":
						path = bad/"tissue_step_1.bin"
						contents = bytearray(path.read_bytes())
						contents[len(contents)//2] ^= 1
						path.write_bytes(contents)
					elif kind == "orphan_tissue":
						(bad/"tissue_step_2.bin").write_bytes(b"orphan")
					else:
						state = json.loads((bad/"run_state.json").read_text())
						state.pop("tissue_checkpoint_sha256")
						(bad/"run_state.json").write_text(json.dumps(state))
					invoke(binary, ranks, reservoir, "--resume", bad, success=False)
				invoke(binary, ranks, reservoir, "--resume", reservoir_split)
				summary(reservoir_split, reservoir, mesh, ranks, 2)
				if ((reservoir_split/"tissue_step_2.bin").read_bytes()
						!= (reservoir_output/"tissue_step_2.bin").read_bytes()
						or abs(json.loads((reservoir_split/"run_summary.json").read_text())[
							"last_tissue_amount_mol"]-reservoir_amounts[ranks]) > 1e-12):
					raise RuntimeError("native species reservoir restart amount changed")
				all_restart_checkpoints_bitwise &= (
					(reservoir_split/"checkpoint_step_2.bin").read_bytes() ==
					(reservoir_output/"checkpoint_step_2.bin").read_bytes())
				maximum_restart_field_difference = max(maximum_restart_field_difference,
					max(abs(a-b) for a, b in zip(reservoir_fields[ranks],
						snapshot(reservoir_split, 2, 5, 4, 0.0))))
			multi = case_file(f"multi_reservoir_{ranks}", {"monotone": True,
				"finite_wall_reservoirs_by_label": {
					"1": {"transfer_coefficient_m_s": 0.5, "volume_m3": 1.0,
						"initial_concentration_mol_m3": 0.0},
					"2": {"transfer_coefficient_m_s": 0.5, "volume_m3": 1.0,
						"initial_concentration_mol_m3": 3.0}}})
			multi_output = work/f"multi_reservoir_{ranks}"
			invoke(binary, ranks, multi, "--check-input")
			invoke(binary, ranks, multi, "--output-dir", multi_output)
			summary(multi_output, multi, mesh, ranks, 2)
			multi_summary = json.loads((multi_output/"run_summary.json").read_text())
			multi_region_amounts[ranks] = {
				label: data["amount_mol"] for label, data in
				multi_summary["tissue_regions_by_label"].items()}
			if (set(multi_region_amounts[ranks]) != {"1", "2"}
					or not multi_region_amounts[ranks]["1"] > 0
					or not multi_region_amounts[ranks]["2"] < 3.0
					or abs(multi_summary["last_combined_balance_defect_mol_s"]) > 1e-10):
				raise RuntimeError("native species multi-region budget or direction differs")
			multi_region_fields[ranks] = snapshot(multi_output, 2, 5, 4, 0.0)
			if ranks == 2:
				multi_split = work/"multi_reservoir_split_2"
				invoke(binary, ranks, multi, "--output-dir", multi_split,
					"--stop-after-step", "1")
				for kind in ("missing", "corrupt", "orphan"):
					bad = work/f"bad_multi_{kind}"
					shutil.copytree(multi_split, bad)
					if kind == "missing":
						(bad/"tissue_step_1.bin").unlink()
					elif kind == "corrupt":
						path = bad/"tissue_step_1.bin"
						contents = bytearray(path.read_bytes())
						contents[len(contents)//2] ^= 1
						path.write_bytes(contents)
					else:
						(bad/"tissue_step_2.bin").write_bytes(b"orphan")
					invoke(binary, ranks, multi, "--resume", bad, success=False)
				invoke(binary, ranks, multi, "--resume", multi_split)
				summary(multi_split, multi, mesh, ranks, 2)
				if ((multi_split/"tissue_step_2.bin").read_bytes()
						!= (multi_output/"tissue_step_2.bin").read_bytes()
						or json.loads((multi_split/"run_summary.json").read_text())[
							"tissue_regions_by_label"]
						!= multi_summary["tissue_regions_by_label"]):
					raise RuntimeError("native species multi-region restart amount changed")
				maximum_restart_field_difference = max(maximum_restart_field_difference,
					max(abs(a-b) for a, b in zip(multi_region_fields[ranks],
						snapshot(multi_split, 2, 5, 4, 0.0))))
			sourced = case_file(f"source_{ranks}", {"source_mol_m3_s": 0.5})
			source_output = work/f"source_{ranks}"
			invoke(binary, ranks, sourced, "--output-dir", source_output)
			summary(source_output, sourced, mesh, ranks, 2)
			field = snapshot(source_output, 2, 5, 4, 0.0)
			if max(abs(value-2.1) for value in field) > 1e-10:
				raise RuntimeError("native species CLI body source did not change concentration")
			split_output = work/f"source_split_{ranks}"
			first_log = invoke(binary, ranks, sourced, "--output-dir", split_output,
				"--stop-after-step", "1")
			if (len(re.findall(r"native_species_step step=", first_log)) != 1
					or (split_output/"run_summary.json").exists()
					or not (split_output/"checkpoint_step_1.bin").is_file()
					or json.loads((split_output/"run_state.json").read_text())["accepted_steps"] != 1):
				raise RuntimeError("native species CLI planned stop did not publish one step")
			state_before = (split_output/"run_state.json").read_bytes()
			changed = case_file(f"changed_source_{ranks}", {"source_mol_m3_s": 0.6})
			invoke(binary, ranks, changed, "--resume", split_output, success=False)
			if (split_output/"run_state.json").read_bytes() != state_before:
				raise RuntimeError("native species CLI wrong-case restart mutated the state")
			if ranks == 2:
				invoke(binary, 1, sourced, "--resume", split_output, success=False)
				for kind in ("missing_checkpoint", "corrupt_checkpoint",
						"missing_snapshot", "missing_piece", "orphan_snapshot"):
					bad = work/f"bad_{kind}"
					shutil.copytree(split_output, bad)
					if kind == "missing_checkpoint":
						(bad/"checkpoint_step_1.bin").unlink()
					elif kind == "corrupt_checkpoint":
						path = bad/"checkpoint_step_1.bin"
						contents = bytearray(path.read_bytes())
						contents[len(contents)//2] ^= 1
						path.write_bytes(contents)
					elif kind == "missing_snapshot":
						(bad/"step_1/snapshot.pvtu").unlink()
					elif kind == "missing_piece":
						(bad/"step_1/rank1.vtu").unlink()
					else:
						(bad/"step_2").mkdir()
					invoke(binary, ranks, sourced, "--resume", bad, success=False)
			resume_log = invoke(binary, ranks, sourced, "--resume", split_output)
			if len(re.findall(r"native_species_step step=", resume_log)) != 1:
				raise RuntimeError("native species CLI resume reran an accepted step")
			summary(split_output, sourced, mesh, ranks, 2)
			all_restart_checkpoints_bitwise &= (
				(split_output/"checkpoint_step_2.bin").read_bytes() ==
				(source_output/"checkpoint_step_2.bin").read_bytes())
			for step in (1, 2):
				restarted = snapshot(split_output, step, 5, 4, 0.0)
				uninterrupted = snapshot(source_output, step, 5, 4, 0.0)
				maximum_restart_field_difference = max(maximum_restart_field_difference,
					max(abs(a-b) for a, b in zip(restarted, uninterrupted)))
			invoke(binary, ranks, sourced, "--resume", split_output, success=False)
		for ranks in (1, 4):
			if max(abs(a-b) for a, b in zip(wall_fields[2], wall_fields[ranks])) > 1e-10:
				raise RuntimeError("native species wall exchange MPI field differs")
			if (max(abs(a-b) for a, b in zip(reservoir_fields[2],
					reservoir_fields[ranks])) > 1e-10
					or abs(reservoir_amounts[2]-reservoir_amounts[ranks]) > 1e-10):
				raise RuntimeError("native species reservoir MPI solution differs")
			if (max(abs(a-b) for a, b in zip(multi_region_fields[2],
					multi_region_fields[ranks])) > 1e-10
					or any(abs(multi_region_amounts[2][label]
						-multi_region_amounts[ranks][label]) > 1e-10
						for label in ("1", "2"))):
				raise RuntimeError("native species multi-region MPI solution differs")
		if maximum_restart_field_difference > 1e-10:
			raise RuntimeError(f"native species CLI restart field changed: "
			f"{maximum_restart_field_difference}")
		front = case_file("front", {
			"initial_concentration_mol_m3": 0.0,
			"fluid_velocity_m_s": [1.0, 0.0, 0.0],
			"inflow_concentration_by_label_mol_m3": {"1": 1.0, "2": 1.0},
			"monotone": True, "time": {"dt_s": 0.1, "steps": 1},
		})
		front_output = work/"front"
		invoke(binary, 2, front, "--output-dir", front_output)
		summary(front_output, front, mesh, 2, 1)
		front_field = snapshot(front_output, 1, 5, 4, 0.0)
		if min(front_field) < -1e-12 or max(front_field) <= 0:
			raise RuntimeError("native species CLI monotone front is invalid")
		diffusive_front = case_file("diffusive_front", {
			"initial_concentration_mol_m3": 0.0,
			"fluid_velocity_m_s": [1.0, 0.0, 0.0],
			"inflow_concentration_by_label_mol_m3": {"1": 1.0, "2": 1.0},
			"diffusivity_m2_s": 1.0, "monotone": True,
			"time": {"dt_s": 0.1, "steps": 1},
		})
		diffusive_output = work/"diffusive_front"
		invoke(binary, 2, diffusive_front, "--output-dir", diffusive_output)
		summary(diffusive_output, diffusive_front, mesh, 2, 1)
		diffusive_field = snapshot(diffusive_output, 1, 5, 4, 0.0)
		diffusion_change = max(abs(a-b) for a, b in zip(front_field, diffusive_field))
		if diffusion_change <= 1e-8:
			raise RuntimeError(f"native species CLI diffusivity did not change the field: "
			f"difference={diffusion_change}, front={front_field}, diffusive={diffusive_field}")
		missing_donor = case_file("missing_donor", {
			"initial_concentration_mol_m3": 0.0,
			"fluid_velocity_m_s": [1.0, 0.0, 0.0],
		})
		invoke(binary, 2, missing_donor, "--check-input", success=False)
		wrong_label = case_file("wrong_label", {
			"inflow_concentration_by_label_mol_m3": {"99": 1.0},
		})
		invoke(binary, 2, wrong_label, "--check-input", success=False)
		unknown = case_file("unknown", {"organ": "liver"})
		invoke(binary, 2, unknown, "--check-input", success=False)
		negative_decay = case_file("negative_decay", {
			"first_order_decay_rate_s_inv": -1.0})
		invoke(binary, 2, negative_decay, "--check-input", success=False)
		for name, changes in (
			("empty_multi_reservoir", {"finite_wall_reservoirs_by_label": {}}),
			("mixed_single_multi_reservoir", {"finite_wall_reservoir": {
				"boundary_label": 1, "transfer_coefficient_m_s": 0.5,
				"volume_m3": 1.0, "initial_concentration_mol_m3": 0.0},
				"finite_wall_reservoirs_by_label": {"2": {
					"transfer_coefficient_m_s": 0.5, "volume_m3": 1.0,
					"initial_concentration_mol_m3": 0.0}}}),
			("mixed_wall_models", {"wall_exchange_by_label": {"1": {
				"transfer_coefficient_m_s": 0.5,
				"external_concentration_mol_m3": 0.0}},
				"finite_wall_reservoir": {"boundary_label": 1,
				"transfer_coefficient_m_s": 0.5, "volume_m3": 1.0,
				"initial_concentration_mol_m3": 0.0}}),
			("mixed_empty_wall_models", {"wall_exchange_by_label": {},
				"finite_wall_reservoir": {"boundary_label": 1,
				"transfer_coefficient_m_s": 0.5, "volume_m3": 1.0,
				"initial_concentration_mol_m3": 0.0}}),
			("negative_reservoir_volume", {"finite_wall_reservoir": {
				"boundary_label": 1, "transfer_coefficient_m_s": 0.5,
				"volume_m3": -1.0, "initial_concentration_mol_m3": 0.0}}),
			("missing_wall", {"wall_exchange_by_label": {
				"99": {"transfer_coefficient_m_s": 0.5,
					"external_concentration_mol_m3": 0.0}}}),
			("negative_wall", {"wall_exchange_by_label": {
				"1": {"transfer_coefficient_m_s": -0.5,
					"external_concentration_mol_m3": 0.0}}}),
			("flowing_wall", {"fluid_velocity_m_s": [1.0, 0.0, 0.0],
				"wall_exchange_by_label": {"1": {"transfer_coefficient_m_s": 0.5,
					"external_concentration_mol_m3": 0.0}}})):
			invalid_wall = case_file(name, changes)
			invoke(binary, 2, invalid_wall, "--check-input", success=False)
		if args.ftetwild is not None:
			root = Path(__file__).resolve().parents[1]
			surface, generated = work/"pipe.vtp", work/"pipe.msh"
			manifest = work/"pipe_manifest.json"
			for command in (
				[sys.executable, str(root/"scripts/generate_circular_pipe_surface.py"),
					str(surface), "--length-m", "0.03", "--radius-m", "0.005",
					"--circumferential-segments", "8", "--axial-segments", "2"],
				[sys.executable, str(root/"scripts/ftetwild_to_fem_volume.py"),
					str(surface), str(generated), "--manifest", str(manifest),
					"--ftetwild", str(args.ftetwild.resolve(strict=True)),
					"--target-size-m", "0.007", "--envelope-m", "0.0001",
					"--stop-energy", "12", "--max-optimization-passes", "40",
					"--max-threads", "2"],
			):
				result = subprocess.run(command, capture_output=True, text=True, timeout=120)
				if result.returncode:
					raise RuntimeError(f"fTetWild species CLI mesh failed: {result.stdout}\n{result.stderr}")
			mesh_data = json.loads(manifest.read_text())["volume_mesh"]
			generated_case = case_file("ftetwild", {"mesh_file": str(generated)})
			generated_output = work/"ftetwild_output"
			invoke(binary, 2, generated_case, "--check-input")
			invoke(binary, 2, generated_case, "--output-dir", generated_output)
			summary(generated_output, generated_case, generated, 2, 2)
			for step in (1, 2):
				field = snapshot(generated_output, step, mesh_data["nodes"],
					mesh_data["elements"], 0.0)
				if max(abs(value-2.0) for value in field) > 1e-10:
					raise RuntimeError("fTetWild native species CLI constant field changed")
		print("native_tet_species_cli: PASS external mesh, 1/2/4 MPI, "
			"moving ALE, source, implicit decay, wall exchange, one/multi-region finite reservoirs, diffusion, monotone front, split/restart, preflight, "
			f"no-overwrite; diffusion_field_change={diffusion_change:.6g}"
			f" restart_max_field_difference={maximum_restart_field_difference:.6g}"
			f" restart_checkpoints_bitwise={int(all_restart_checkpoints_bitwise)}"
			+ (f" ftetwild_tetrahedra={mesh_data['elements']}" if args.ftetwild is not None else ""))


if __name__ == "__main__":
	main()
