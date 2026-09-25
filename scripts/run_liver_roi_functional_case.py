#!/usr/bin/env python3
"""Run a versioned, expressly nonphysiological liver ROI FEM/Darcy case."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys


INPUT_NAMES = ("source_seg", "geometry_manifest", "multiregion_mesh",
	"multiregion_contract", "split_manifest")
PHYSICAL_NAMES = ("inlet_velocity_m_s", "density_kg_m3", "viscosity_pa_s",
	"mobility_m2_pa_s", "interface_pressure_pa", "tissue_exit_pressure_pa")


def reject(message):
	raise ValueError(message)


def sha256(path):
	return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_case(path, data_dir):
	data = json.loads(path.read_text(encoding="utf-8"))
	if set(data) != {"schema_version", "kind", "full_liver_case",
			"physiological_validation", "source", "files_relative_to_data_dir",
			"regions", "boundary_labels", "explicit_functional_inputs",
			"input_origin", "limitations"} \
			or data["schema_version"] != 1 \
			or data["kind"] != "patient_derived_liver_roi_functional_only" \
			or data["full_liver_case"] is not False \
			or data["physiological_validation"] is not False:
		reject("case schema or nonphysiological classification is invalid")
	source = data["source"]
	if source.get("collection") != "TCIA Colorectal-Liver-Metastases" \
			or source.get("case_id") != "CRLM-CT-1072" \
			or source.get("license") != "CC BY 4.0" \
			or source.get("provenance_document") != "docs/LIVER_GEOMETRY_CANDIDATES.md" \
			or source.get("roi_zyx_half_open") != [[10, 20], [210, 226], [101, 117]]:
		reject("case provenance or ROI differs from the audited specimen")
	if set(data["files_relative_to_data_dir"]) != set(INPUT_NAMES):
		reject("case input file set is incomplete")
	root = data_dir.resolve()
	paths = {}
	for name in INPUT_NAMES:
		entry = data["files_relative_to_data_dir"][name]
		if set(entry) != {"path", "sha256"} or not isinstance(entry["path"], str):
			reject(f"{name} file entry is malformed")
		relative = Path(entry["path"])
		if relative.is_absolute() or not relative.parts or ".." in relative.parts:
			reject(f"{name} path escapes the case data directory")
		candidate = (root/relative).resolve()
		if not candidate.is_relative_to(root) or not candidate.is_file() \
				or sha256(candidate) != entry["sha256"]:
			reject(f"{name} file is missing, outside the case directory, or has wrong hash")
		paths[name] = candidate
	regions = data["regions"]
	if set(regions) != {"vessel", "tissue", "matching_interface"} \
			or regions["vessel"] != "seg_4_portal" \
			or regions["tissue"] != "tissue_candidate" \
			or regions["matching_interface"] != "interface_seg_4_portal_tissue_candidate":
		reject("case region/interface semantics are invalid")
	labels = data["boundary_labels"]
	if set(labels) != {"artificial_roi_cut_fluid_inlet",
			"artificial_tissue_pressure_exits", "matching_interface"} \
			or labels["artificial_roi_cut_fluid_inlet"] != 1 \
			or labels["artificial_tissue_pressure_exits"] != [2, 3] \
			or labels["matching_interface"] != 4:
		reject("case boundary labels differ from the audited ROI")
	inputs = data["explicit_functional_inputs"]
	if set(inputs) != set(PHYSICAL_NAMES) \
			or not isinstance(inputs["inlet_velocity_m_s"], list) \
			or len(inputs["inlet_velocity_m_s"]) != 3:
		reject("case physical inputs are incomplete")
	for name in PHYSICAL_NAMES:
		values = inputs[name] if name == "inlet_velocity_m_s" else [inputs[name]]
		if any(isinstance(value, bool) or not isinstance(value, (int, float))
				or not math.isfinite(value) for value in values):
			reject(f"{name} must contain finite SI values")
	if any(inputs[name] <= 0 for name in ("density_kg_m3", "viscosity_pa_s",
			"mobility_m2_pa_s")) or not any(inputs["inlet_velocity_m_s"]):
		reject("case material values and inlet speed must be nonzero/positive")
	if "not measured or sourced patient physiology" not in data["input_origin"]:
		reject("case must explicitly classify artificial functional inputs")
	return data, paths


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--case", type=Path,
		default=Path("cases/liver_roi_functional.json"))
	parser.add_argument("--data-dir", type=Path, required=True)
	parser.add_argument("--ranks", type=int, default=1)
	parser.add_argument("--solver", type=Path,
		default=Path("solvers/cpu/native_tet_matching_solved_roi_smoke"))
	parser.add_argument("--output", type=Path)
	parser.add_argument("--field-output-dir", type=Path)
	parser.add_argument("--check-only", action="store_true")
	args = parser.parse_args()
	try:
		if args.ranks < 1 or args.ranks > 64:
			reject("MPI ranks must be between 1 and 64")
		data, paths = validate_case(args.case, args.data_dir)
		case_hash = sha256(args.case)
		if args.check_only:
			print(f"liver ROI functional case input: PASS case_sha256={case_hash}")
			return 0
		if args.output is None or args.field_output_dir is None:
			reject("solving requires --output and --field-output-dir")
		labels = data["boundary_labels"]
		inputs = data["explicit_functional_inputs"]
		command = [sys.executable,
			str(Path(__file__).with_name("run_native_matching_roi_functional.py")),
			"--source-seg", str(paths["source_seg"]),
			"--geometry-manifest", str(paths["geometry_manifest"]),
			"--multiregion-mesh", str(paths["multiregion_mesh"]),
			"--multiregion-contract", str(paths["multiregion_contract"]),
			"--split-manifest", str(paths["split_manifest"]),
			"--vessel-region", data["regions"]["vessel"],
			"--tissue-region", data["regions"]["tissue"],
			"--interface-name", data["regions"]["matching_interface"],
			"--inlet-label", str(labels["artificial_roi_cut_fluid_inlet"]),
			"--tissue-exit-labels", *map(str, labels["artificial_tissue_pressure_exits"]),
			"--inlet-velocity-m-s", *map(str, inputs["inlet_velocity_m_s"]),
			"--density-kg-m3", str(inputs["density_kg_m3"]),
			"--viscosity-pa-s", str(inputs["viscosity_pa_s"]),
			"--mobility-m2-pa-s", str(inputs["mobility_m2_pa_s"]),
			"--interface-pressure-pa", str(inputs["interface_pressure_pa"]),
			"--tissue-exit-pressure-pa", str(inputs["tissue_exit_pressure_pa"]),
			"--ranks", str(args.ranks), "--solver", str(args.solver),
			"--output", str(args.output), "--field-output-dir",
			str(args.field_output_dir), "--case-file-sha256", case_hash]
		return subprocess.run(command, check=False).returncode
	except (KeyError, TypeError, OSError, ValueError, json.JSONDecodeError) as error:
		print(f"liver ROI functional case rejected: {error}", file=sys.stderr)
		return 2


if __name__ == "__main__":
	sys.exit(main())
