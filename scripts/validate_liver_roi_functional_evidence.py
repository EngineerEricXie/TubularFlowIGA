#!/usr/bin/env python3
"""Check scoped ROI functional evidence; never upgrade it to organ validation."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import xml.etree.ElementTree as ET


SOURCE_FILES = (
	"scripts/audit_dicom_seg_region_overlap.py",
	"scripts/audit_seg_exposed_vessel_faces.py",
	"scripts/dicom_seg_to_multiregion_tet.py",
	"scripts/split_multiregion_tet_for_native.py",
	"scripts/run_native_matching_roi_functional.py",
	"scripts/time_native_mpi_rank.py",
	"scripts/run_liver_roi_functional_case.py",
	"scripts/run_liver_roi_from_raw.py",
	"cases/liver_roi_functional.json",
	"solvers/cpu/src/native_tet_matching_solved_roi_smoke.cpp",
	"solvers/cpu/include/NativeTetMatchingInterfaceSource.hpp",
	"solvers/cpu/include/NativeTetMeshComponents.hpp",
	"solvers/cpu/include/NativeTetAlePetscRuntime.hpp",
	"solvers/cpu/include/NativeTetDarcyPetsc.hpp",
	"solvers/cpu/include/NativeTetDarcyVisualization.hpp",
	"solvers/cpu/include/NativeTetHydraulicVisualization.hpp",
	"solvers/cpu/include/ParallelVtkOutput.hpp",
	"include/PartitionedVtkOutput.hpp",
)


def reject(message):
	raise ValueError(message)


def source_hash(root):
	digest = hashlib.sha256()
	for name in SOURCE_FILES:
		digest.update(name.encode("utf-8")+b"\0")
		digest.update((root/name).read_bytes())
	return digest.hexdigest()


def check_metrics(metrics):
	keys = {"facets", "fluid_inlet_m3_s", "interface_source_m3_s",
		"darcy_outward_m3_s", "fluid_newton_iterations", "fluid_linear_reason",
		"darcy_linear_reason", "darcy_max_cell_defect_m3_s"}
	if set(metrics) != keys or any(not isinstance(value, (int, float))
			or not math.isfinite(value) for value in metrics.values()):
		reject("ROI metrics are incomplete or nonfinite")
	source = metrics["interface_source_m3_s"]
	scale = max(abs(source), 1e-18)
	if metrics["facets"] != 158 or source <= 0 \
			or metrics["fluid_inlet_m3_s"] >= 0 \
			or abs(metrics["fluid_inlet_m3_s"]+source) > 1e-6*scale \
			or abs(source-metrics["darcy_outward_m3_s"]) > 1e-8*scale \
			or metrics["darcy_max_cell_defect_m3_s"] > 1e-8*scale \
			or metrics["fluid_newton_iterations"] <= 0 \
			or metrics["fluid_linear_reason"] <= 0 \
			or metrics["darcy_linear_reason"] <= 0:
		reject("ROI convergence, interface, or conservation gate failed")


def field_bundle_hash(hashes):
	digest = hashlib.sha256()
	for name, value in sorted(hashes.items()):
		digest.update(name.encode("utf-8")+b"\0"+value.encode("ascii")+b"\0")
	return digest.hexdigest()


def read_rank_owned_fields(directory, ranks, role, cell_count, fields):
	"""Reassemble fields by explicit global IDs, rejecting overlap and holes."""
	index = ET.parse(directory/role/"snapshot.pvtu")
	if [piece.get("Source") for piece in index.findall(".//Piece")] != [
			f"rank{rank}.vtu" for rank in range(ranks)]:
		reject(f"{role} PVTU pieces are incomplete")
	assembled = {name: {} for name in fields}
	cells = set()
	for rank in range(ranks):
		piece = ET.parse(directory/role/f"rank{rank}.vtu").find(".//Piece")
		if piece is None:
			reject(f"{role} rank-owned VTU piece is absent")
		for kind, id_name in (("PointData", "GlobalPointIds"),
				("CellData", "GlobalCellIds")):
			parent = piece.find(kind)
			if parent is None:
				reject(f"{role} rank {rank} lacks {kind}")
			ids = parent.find(f"DataArray[@Name='{id_name}']")
			if ids is None:
				reject(f"{role} rank {rank} lacks {id_name}")
			keys = [int(value) for value in (ids.text or "").split()]
			if len(keys) != len(set(keys)) or len(keys) != int(piece.get(
					"NumberOfPoints" if kind == "PointData" else "NumberOfCells", "-1")):
				reject(f"{role} rank {rank} has malformed global IDs")
			if kind == "CellData":
				if cells.intersection(keys):
					reject(f"{role} rank-owned cells overlap")
				cells.update(keys)
			for name, location in fields.items():
				if location != kind:
					continue
				array = parent.find(f"DataArray[@Name='{name}']")
				if array is None:
					reject(f"{role} rank {rank} lacks {name}")
				components = int(array.get("NumberOfComponents", "1"))
				values = [float(value) for value in (array.text or "").split()]
				if len(values) != components*len(keys) or not all(map(math.isfinite, values)):
					reject(f"{role} rank {rank} has malformed {name}")
				for offset, key in enumerate(keys):
					value = tuple(values[components*offset:components*(offset+1)])
					old = assembled[name].setdefault(key, value)
					if old != value:
						reject(f"{role} shared point {name} differs across pieces")
	if len(cells) != cell_count:
		reject(f"{role} rank-owned cells do not cover the mesh")
	return assembled


def compare_fields(reference, candidate, label):
	for name, values in reference.items():
		other = candidate[name]
		if values.keys() != other.keys():
			reject(f"{label} {name} global IDs differ across MPI ranks")
		scale = max((abs(number) for value in values.values() for number in value),
			default=0.)
		maximum = max((abs(a-b) for key, value in values.items()
			for a, b in zip(value, other[key])), default=0.)
		if maximum > 1e-6*scale+1e-14:
			reject(f"{label} {name} differs across MPI ranks")


def read_facet_ledger(directory, expected_facets, expected_source):
	ledger = json.loads((directory/"interface_facets.json").read_text(encoding="utf-8"))
	if ledger.get("schema_version") != 1 \
			or ledger.get("kind") != "native_matching_facet_flow_functional_only" \
			or ledger.get("flow_sign") != "vessel_outward_positive_tissue_source_positive" \
			or len(ledger.get("facets", [])) != expected_facets:
		reject("matching-facet ledger schema or size is invalid")
	result = {}
	for facet in ledger["facets"]:
		if set(facet) != {"port_name", "vessel_cell_id", "tissue_cell_id",
				"vessel_outward_m3_s", "triangle_m"} \
				or facet["port_name"] != "solved_roi_functional" \
				or not math.isfinite(facet["vessel_outward_m3_s"]):
			reject("matching-facet ledger entry is invalid")
		triangle = tuple(tuple(point) for point in facet["triangle_m"])
		if len(triangle) != 3 or any(len(point) != 3 for point in triangle) \
				or any(not math.isfinite(value) for point in triangle for value in point) \
				or triangle in result:
			reject("matching-facet ledger geometry is invalid or duplicate")
		result[triangle] = (facet["vessel_cell_id"], facet["tissue_cell_id"],
			facet["vessel_outward_m3_s"])
	flow = sum(value[2] for value in result.values())
	if abs(flow-expected_source) > 1e-8*max(abs(expected_source), 1e-18):
		reject("matching-facet ledger does not close to source")
	return result


def compare_facet_ledgers(reference, candidate):
	if reference.keys() != candidate.keys():
		reject("matching-facet geometry differs across MPI ranks")
	for triangle, first in reference.items():
		other = candidate[triangle]
		if first[:2] != other[:2] or abs(first[2]-other[2]) > 1e-6*max(
				abs(first[2]), abs(other[2]), 1e-18):
			reject("matching-facet ownership or signed flow differs across MPI ranks")


def require_v2_summary_bindings(summary):
	if summary.get("schema_version") == 2 and (
			"source_face_audit_sha256" not in summary
			or "stage_metrics" not in summary):
		reject("fresh series version 2 requires source-face audit and stage metrics")


def validate_fresh_series(case_path, series_dir, solver_path,
		require_summary=True):
	"""Validate a newly rebuilt raw-SEG-to-fields run without archived hashes."""
	try:
		from scripts.run_liver_roi_functional_case import validate_case
	except ModuleNotFoundError:
		from run_liver_roi_functional_case import validate_case
	case, paths = validate_case(case_path, series_dir)
	case_hash = hashlib.sha256(case_path.read_bytes()).hexdigest()
	manifest = json.loads(paths["geometry_manifest"].read_text(encoding="utf-8"))
	split = json.loads(paths["split_manifest"].read_text(encoding="utf-8"))
	solver_hash = hashlib.sha256(solver_path.read_bytes()).hexdigest()
	if manifest.get("mesh_sha256") != hashlib.sha256(paths[
			"multiregion_mesh"].read_bytes()).hexdigest() \
			or split.get("source_mesh_sha256") != manifest["mesh_sha256"]:
			reject("fresh series geometry manifest does not match generated mesh")
	reference = None
	rank_resource_presence = {}
	for rank in (1, 2, 4):
		result = json.loads((series_dir/f"rank{rank}"/"result.json").read_text(
			encoding="utf-8"))
		if result.get("schema_version") != 1 \
				or result.get("kind") != "patient_derived_roi_solved_fluid_darcy_functional_only" \
				or result.get("full_liver_case") is not False \
				or result.get("physiological_validation") is not False \
				or result.get("case_file_sha256") != case_hash \
				or result.get("source_seg_sha256") != case[
					"files_relative_to_data_dir"]["source_seg"]["sha256"] \
				or result.get("multiregion_mesh_sha256") != manifest["mesh_sha256"] \
				or result.get("multiregion_contract_sha256") != hashlib.sha256(
					paths["multiregion_contract"].read_bytes()).hexdigest() \
				or result.get("solver_binary_sha256") != solver_hash \
				or result.get("explicit_functional_inputs") != case[
					"explicit_functional_inputs"] \
				or result.get("mpi_ranks") != rank:
			reject(f"fresh series rank {rank} provenance or classification differs")
		if result.get("regions") != {"vessel": case["regions"]["vessel"],
				"tissue": case["regions"]["tissue"]} \
				or result.get("labels") != {
					"interface": split["interface_boundary_labels"][
						case["regions"]["matching_interface"]],
					"inlet": case["boundary_labels"]["artificial_roi_cut_fluid_inlet"],
					"tissue_exits": case["boundary_labels"][
						"artificial_tissue_pressure_exits"]}:
			reject(f"fresh series rank {rank} region or boundary labels differ")
		for role, region in result["regions"].items():
			key = "vessel_mesh_sha256" if role == "vessel" else "tissue_mesh_sha256"
			if result.get(key) != split["submeshes"][region]["sha256"]:
				reject(f"fresh series rank {rank} native submesh hash differs")
		metrics = result["metrics"]
		check_metrics(metrics)
		directory = (series_dir/f"rank{rank}"/"fields").resolve()
		if result.get("rank_owned_field_output_directory") != str(directory):
			reject(f"fresh series rank {rank} field directory differs")
		hashes = result.get("rank_owned_field_files_sha256")
		expected_names = {f"{role}/rank{owner}.vtu"
			for role in ("fluid", "tissue") for owner in range(rank)} \
			| {"fluid/snapshot.pvtu", "tissue/snapshot.pvtu",
				"interface_facets.json"}
		if not isinstance(hashes, dict) or set(hashes) != expected_names:
			reject(f"fresh series rank {rank} field file inventory differs")
		for name, expected in hashes.items():
			if hashlib.sha256((directory/name).read_bytes()).hexdigest() != expected:
				reject(f"fresh series rank {rank} field file hash differs: {name}")
		resources = result.get("rank_solver_resources")
		rank_resource_presence[rank] = resources is not None
		if resources is not None:
			metrics_dir = directory.parent/(directory.name+".rank-metrics")
			if resources.get("rank_metrics_directory") != str(metrics_dir) \
					or "per rank, not summed" not in resources.get("scope", "") \
					or set(resources.get("ranks", {})) != {
						str(owner) for owner in range(rank)} \
					or {path.name for path in metrics_dir.iterdir()} != {
						f"rank{owner}.time.txt" for owner in range(rank)}:
				reject(f"fresh series rank {rank} solver resource inventory differs")
			for owner in range(rank):
				path = metrics_dir/f"rank{owner}.time.txt"
				values = path.read_text(encoding="utf-8").split()
				entry = resources["ranks"][str(owner)]
				if len(values) != 2 or len(entry) != 3 \
						or not math.isfinite(entry["wall_seconds"]) \
						or entry["wall_seconds"] < 0 \
						or entry["gnu_time_maximum_rss_kb"] <= 0 \
						or float(values[0]) != entry["wall_seconds"] \
						or int(values[1]) != entry["gnu_time_maximum_rss_kb"] \
						or hashlib.sha256(path.read_bytes()).hexdigest() != entry[
							"timing_log_sha256"]:
					reject(f"fresh series rank {rank} solver resource file differs")
		fluid = read_rank_owned_fields(directory, rank, "fluid",
			split["submeshes"][case["regions"]["vessel"]]["audit"]["tetrahedra"],
			{"velocity_m_s": "PointData", "pressure_pa": "PointData"})
		tissue = read_rank_owned_fields(directory, rank, "tissue",
			split["submeshes"][case["regions"]["tissue"]]["audit"]["tetrahedra"],
			{"pressure_pa": "PointData", "darcy_rt0_centroid_flux_m_s": "CellData",
				"conservative_outward_face_flow_m3_s": "CellData", "source_s_inv": "CellData"})
		ledger = read_facet_ledger(directory,int(metrics["facets"]),
			metrics["interface_source_m3_s"])
		if reference is None:
			reference = {"metrics": metrics, "fluid": fluid, "tissue": tissue,
				"ledger": ledger}
		else:
			for name in ("fluid_inlet_m3_s", "interface_source_m3_s",
					"darcy_outward_m3_s"):
				first = reference["metrics"][name]
				if abs(first-metrics[name]) > 1e-8*max(abs(first), 1e-18):
					reject(f"fresh series {name} differs across MPI ranks")
			compare_fields(reference["fluid"], fluid, "fluid")
			compare_fields(reference["tissue"], tissue, "tissue")
			compare_facet_ledgers(reference["ledger"], ledger)
	summary_path = series_dir/"summary.json"
	if require_summary and not summary_path.is_file():
		reject("fresh series completed summary is missing")
	if summary_path.exists():
		summary = json.loads(summary_path.read_text(encoding="utf-8"))
		if summary.get("schema_version") not in (1, 2) \
				or summary.get("kind") != "patient_derived_liver_roi_raw_to_native_fem_functional_only" \
				or summary.get("full_liver_case") is not False \
				or summary.get("physiological_validation") is not False \
				or summary.get("case_file_sha256") != case_hash \
				or summary.get("source_seg_sha256") != case[
					"files_relative_to_data_dir"]["source_seg"]["sha256"] \
				or summary.get("multiregion_mesh_sha256") != manifest["mesh_sha256"] \
				or summary.get("solver_binary_sha256") != solver_hash \
				or summary.get("mpi_ranks_tested") != [1, 2, 4]:
				reject("fresh series summary provenance or classification differs")
		require_v2_summary_bindings(summary)
		if summary["schema_version"] == 2 and not all(
				rank_resource_presence.values()):
			reject("fresh series version 2 requires every native solver rank resource record")
		if "source_face_audit_sha256" in summary:
			face_path = series_dir/"source-face-audit.json"
			if hashlib.sha256(face_path.read_bytes()).hexdigest() != summary[
					"source_face_audit_sha256"]:
				reject("fresh series source-face audit hash differs")
			faces = json.loads(face_path.read_text(encoding="utf-8"))
			if faces.get("kind") != "seg_vessel_face_inventory_not_ports_or_boundary_conditions" \
					or faces.get("source_seg_sha256") != summary["source_seg_sha256"] \
					or faces.get("tissue_segment_number") != 1 \
					or faces.get("vessel_segment_numbers") != [3, 4] \
					or [len(vessel["components"]) for vessel in faces["vessels"]] != [1, 116]:
				reject("fresh series source-face audit provenance differs")
			for vessel, expected in zip(faces["vessels"],
					((1876, 3048, 0, 0), (33862, 6376, 0, 0))):
				if tuple(sum(component[role]["total_faces"] for component in
						vessel["components"]) for role in
						("tissue", "background", "other_vessel", "grid_extent")) != expected:
					reject("fresh series source-face category totals differ")
		if "stage_metrics" in summary:
			expected_stages = {"mesh", "split", "solve_rank1", "solve_rank2",
				"solve_rank4"}
			if "source_face_audit_sha256" in summary:
				expected_stages.add("source_face_audit")
			if set(summary["stage_metrics"]) != expected_stages \
					or "not the sum of MPI rank memory" not in summary.get(
						"resource_measurement_scope", ""):
				reject("fresh series stage-resource scope or inventory differs")
			for stage in summary["stage_metrics"].values():
				if set(stage) != {"wall_seconds", "gnu_time_maximum_rss_kb"} \
						or not math.isfinite(stage["wall_seconds"]) \
						or stage["wall_seconds"] < 0 \
						or stage["gnu_time_maximum_rss_kb"] <= 0:
					reject("fresh series stage-resource value is invalid")
	return {"case_file_sha256": case_hash,
		"source_seg_sha256": case["files_relative_to_data_dir"]["source_seg"]["sha256"],
		"multiregion_mesh_sha256": manifest["mesh_sha256"],
		"solver_binary_sha256": solver_hash,
		"reference": reference}


def compare_fresh_series(first, second):
	for key in ("case_file_sha256", "source_seg_sha256",
			"multiregion_mesh_sha256", "solver_binary_sha256"):
		if first[key] != second[key]:
			reject(f"fresh reruns differ in {key}")
	a = first["reference"]
	b = second["reference"]
	for name in ("fluid_inlet_m3_s", "interface_source_m3_s",
			"darcy_outward_m3_s"):
		left, right = a["metrics"][name], b["metrics"][name]
		if abs(left-right) > 1e-8*max(abs(left), abs(right), 1e-18):
			reject(f"fresh reruns differ in {name}")
	compare_fields(a["fluid"], b["fluid"], "fluid rerun")
	compare_fields(a["tissue"], b["tissue"], "tissue rerun")
	compare_facet_ledgers(a["ledger"], b["ledger"])
	return True


def validate(root, card_path, results_dir=None, field_results_dir=None,
		case_result=None, component_roi_audit=None):
	data = json.loads(card_path.read_text(encoding="utf-8"))
	if data.get("schema_version") != 1 \
			or data.get("kind") != "patient_derived_liver_roi_functional_not_physiological" \
			or data.get("full_liver_case") is not False \
			or data.get("physiological_validation") is not False \
			or data.get("source_files_sha256") != source_hash(root):
		reject("ROI card classification or current source hash is invalid")
	if data.get("case_file_sha256") != hashlib.sha256((root/
			"cases/liver_roi_functional.json").read_bytes()).hexdigest():
		reject("ROI versioned case file differs from evidence card")
	if data.get("mpi_ranks_tested") != [1, 2, 4] \
			or set(data.get("metrics_by_rank", {})) != {"1", "2", "4"}:
		reject("ROI card lacks 1/2/4-rank measurements")
	inputs = data.get("explicit_functional_inputs", {})
	if inputs != {"inlet_velocity_m_s": [0, 0, -0.01],
			"density_kg_m3": 1000, "viscosity_pa_s": 0.004,
			"mobility_m2_pa_s": 1, "interface_pressure_pa": 0,
			"tissue_exit_pressure_pa": 0}:
		reject("ROI artificial functional inputs changed without new evidence")
	for rank in (1, 2, 4):
		metrics = data["metrics_by_rank"][str(rank)]
		check_metrics(metrics)
		if results_dir is not None:
			path = results_dir/f"solved-functional-r{rank}.json"
			current = json.loads(path.read_text(encoding="utf-8"))
			if current.get("mpi_ranks") != rank \
					or current.get("kind") != "patient_derived_roi_solved_fluid_darcy_functional_only" \
					or current.get("full_liver_case") is not False \
					or current.get("physiological_validation") is not False \
					or current.get("metrics") != metrics:
				reject(f"rank {rank} current result differs from evidence card")
			for key in ("source_seg_sha256", "multiregion_mesh_sha256",
					"multiregion_contract_sha256", "vessel_mesh_sha256",
					"tissue_mesh_sha256"):
				if current.get(key) != data.get(key):
					reject(f"rank {rank} source identity differs from evidence card")
	first = data["metrics_by_rank"]["1"]["interface_source_m3_s"]
	for rank in (2, 4):
		other = data["metrics_by_rank"][str(rank)]["interface_source_m3_s"]
		if abs(first-other) > 1e-8*abs(first):
			reject("ROI source flow is inconsistent across ranks")
	if field_results_dir is not None:
		reference = {}
		for rank in (1, 2, 4):
			result = json.loads((field_results_dir/
				f"solved-functional-ledger-final-r{rank}.json").read_text(encoding="utf-8"))
			if result.get("mpi_ranks") != rank \
					or result.get("kind") != "patient_derived_roi_solved_fluid_darcy_functional_only" \
					or result.get("full_liver_case") is not False \
					or result.get("physiological_validation") is not False \
					or result.get("solver_binary_sha256") != data.get(
						"solver_binary_sha256_for_field_evidence") \
					or result.get("source_seg_sha256") != data["source_seg_sha256"] \
					or result.get("vessel_mesh_sha256") != data["vessel_mesh_sha256"] \
					or result.get("tissue_mesh_sha256") != data["tissue_mesh_sha256"]:
				reject(f"rank {rank} field provenance differs from evidence card")
			check_metrics(result["metrics"])
			for name in ("facets", "fluid_newton_iterations", "fluid_linear_reason",
					"darcy_linear_reason", "fluid_inlet_m3_s",
					"interface_source_m3_s", "darcy_outward_m3_s"):
				value = data["metrics_by_rank"][str(rank)][name]
				if abs(value-result["metrics"][name]) > 1e-8*max(abs(value), 1e-18):
					reject(f"rank {rank} field-run metric {name} differs from evidence card")
			hashes = result.get("rank_owned_field_files_sha256")
			if not isinstance(hashes, dict) or field_bundle_hash(hashes) != data[
					"rank_owned_field_bundle_sha256_by_rank"][str(rank)]:
				reject(f"rank {rank} field bundle differs from evidence card")
			directory = Path(result["rank_owned_field_output_directory"])
			for name, expected in hashes.items():
				if hashlib.sha256((directory/name).read_bytes()).hexdigest() != expected:
					reject(f"rank {rank} field file hash differs: {name}")
			fluid = read_rank_owned_fields(directory,rank,"fluid",data["portal_tetrahedra"],
				{"velocity_m_s": "PointData", "pressure_pa": "PointData"})
			tissue = read_rank_owned_fields(directory,rank,"tissue",
				data["tissue_candidate_tetrahedra"],
				{"pressure_pa": "PointData", "darcy_rt0_centroid_flux_m_s": "CellData",
					"conservative_outward_face_flow_m3_s": "CellData", "source_s_inv": "CellData"})
			ledger = read_facet_ledger(directory,data["interface_facets"],
				result["metrics"]["interface_source_m3_s"])
			if rank == 1:
				reference = {"fluid": fluid, "tissue": tissue, "ledger": ledger}
			else:
				compare_fields(reference["fluid"], fluid, "fluid")
				compare_fields(reference["tissue"], tissue, "tissue")
				compare_facet_ledgers(reference["ledger"], ledger)
	if case_result is not None:
		result = json.loads(case_result.read_text(encoding="utf-8"))
		if result.get("case_file_sha256") != data["case_file_sha256"] \
				or result.get("source_seg_sha256") != data["source_seg_sha256"] \
				or result.get("multiregion_mesh_sha256") != data["multiregion_mesh_sha256"] \
				or result.get("full_liver_case") is not False \
				or result.get("physiological_validation") is not False \
				or result.get("explicit_functional_inputs") != data["explicit_functional_inputs"]:
			reject("versioned case result provenance or classification differs")
		check_metrics(result["metrics"])
		directory = Path(result["rank_owned_field_output_directory"])
		hashes = result["rank_owned_field_files_sha256"]
		if not isinstance(hashes, dict):
			reject("versioned case has no rank-owned fields")
		for name, expected in hashes.items():
			if hashlib.sha256((directory/name).read_bytes()).hexdigest() != expected:
				reject(f"versioned case field hash differs: {name}")
		for role, count, fields in (
			("fluid", data["portal_tetrahedra"],
				{"velocity_m_s": "PointData", "pressure_pa": "PointData"}),
			("tissue", data["tissue_candidate_tetrahedra"],
				{"pressure_pa": "PointData", "source_s_inv": "CellData"})):
			read_rank_owned_fields(directory,result["mpi_ranks"],role,count,fields)
		read_facet_ledger(directory,data["interface_facets"],
			result["metrics"]["interface_source_m3_s"])
	if component_roi_audit is not None:
		if hashlib.sha256(component_roi_audit.read_bytes()).hexdigest() != data[
				"component_roi_audit_sha256"]:
			reject("component ROI audit hash differs from evidence card")
		inventory = json.loads(component_roi_audit.read_text(encoding="utf-8"))
		if inventory.get("source_sha256") != data["source_seg_sha256"] \
				or inventory.get("kind") != "same_seg_voxel_overlap_and_face_contact_not_mesh_or_flow" \
				or inventory["roi_component_membership"] != {
					"roi_zyx_half_open": data["roi_zyx_half_open"],
					"tissue_candidate": [{"component_id": 1, "roi_voxels": 2188}],
					"vessels": {"3": [], "4": [{"component_id": 1,
						"roi_voxels": 32}]},
					"role": "full-grid component IDs intersecting ROI only; no port, BC, or component-selection policy"}:
			reject("component ROI audit source or membership differs")
	return True


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--results-dir", type=Path)
	parser.add_argument("--field-results-dir", type=Path)
	parser.add_argument("--case-result", type=Path)
	parser.add_argument("--component-roi-audit", type=Path)
	parser.add_argument("--fresh-series-dir", type=Path)
	parser.add_argument("--compare-fresh-series-dir", type=Path)
	parser.add_argument("--fresh-case", type=Path,
		default=Path("cases/liver_roi_functional.json"))
	parser.add_argument("--fresh-solver", type=Path,
		default=Path("solvers/cpu/native_tet_matching_solved_roi_smoke"))
	args = parser.parse_args()
	root = Path(__file__).resolve().parents[1]
	card = root/"benchmarks"/"liver_roi_functional_evidence.json"
	try:
		if args.compare_fresh_series_dir is not None and args.fresh_series_dir is None:
			reject("--compare-fresh-series-dir requires --fresh-series-dir")
		if args.fresh_series_dir is not None:
			first = validate_fresh_series(args.fresh_case,
				args.fresh_series_dir.resolve(),args.fresh_solver)
			if args.compare_fresh_series_dir is not None:
				second = validate_fresh_series(args.fresh_case,
					args.compare_fresh_series_dir.resolve(),args.fresh_solver)
				compare_fresh_series(first, second)
		else:
			validate(root, card, args.results_dir, args.field_results_dir,
				args.case_result, args.component_roi_audit)
	except (OSError, ValueError, KeyError, json.JSONDecodeError) as error:
		print(f"liver ROI functional evidence: ERROR: {error}", file=sys.stderr)
		return 2
	print("liver ROI functional evidence: PASS (local ROI function only; no full-organ or physiological validation)")
	return 0


if __name__ == "__main__":
	sys.exit(main())
