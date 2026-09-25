#!/usr/bin/env python3
"""Rebuild the audited local liver ROI and run native FEM/Darcy on 1/2/4 ranks.

This is an artificial-BC functional workflow, never a full-organ or
physiological perfusion claim. Keep its patient-derived output outside git.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

from run_liver_roi_functional_case import validate_case
from validate_liver_roi_functional_evidence import validate_fresh_series


def reject(message):
	raise ValueError(message)


def sha256(path):
	return hashlib.sha256(path.read_bytes()).hexdigest()


def run(stage, command, output):
	"""Record GNU time's stage/launcher RSS, not summed MPI worker memory."""
	log = output/f"{stage}.time.txt"
	if log.exists():
		reject(f"stage timing log already exists: {log}")
	completed = subprocess.run(["/usr/bin/time", "-f", "%e %M", "-o",
		str(log), *command], check=False)
	if completed.returncode:
		reject(f"stage exited {completed.returncode}: {' '.join(map(str, command[:3]))}")
	values = log.read_text(encoding="utf-8").split()
	if len(values) != 2:
		reject(f"stage timing log is malformed: {stage}")
	wall, rss = float(values[0]), int(values[1])
	if not math.isfinite(wall) or wall < 0 or rss <= 0:
		reject(f"stage timing or memory is invalid: {stage}")
	return {"wall_seconds": wall, "gnu_time_maximum_rss_kb": rss}


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--source-seg", type=Path, required=True)
	parser.add_argument("--output-dir", type=Path, required=True)
	parser.add_argument("--case", type=Path,
		default=Path("cases/liver_roi_functional.json"))
	parser.add_argument("--solver", type=Path,
		default=Path("solvers/cpu/native_tet_matching_solved_roi_smoke"))
	args = parser.parse_args()
	try:
		root = Path(__file__).resolve().parents[1]
		output = args.output_dir.resolve()
		case = json.loads(args.case.read_text(encoding="utf-8"))
		if case.get("schema_version") != 1 \
				or case.get("kind") != "patient_derived_liver_roi_functional_only" \
				or case.get("full_liver_case") is not False \
				or case.get("physiological_validation") is not False:
			reject("raw-to-result case classification is invalid")
		seg = args.source_seg.resolve()
		if not seg.is_file() or sha256(seg) != case[
				"files_relative_to_data_dir"]["source_seg"]["sha256"]:
			reject("original DICOM SEG is missing or has a different SHA-256")
		solver = args.solver.resolve()
		if not solver.is_file():
			reject("native matching FEM solver binary is absent")
		if output.exists() or output.is_relative_to(root):
			reject("output directory already exists or is inside the source repository")
		roi = case["source"]["roi_zyx_half_open"]
		if roi != [[10, 20], [210, 226], [101, 117]]:
			reject("raw-to-result ROI differs from the audited case")
		output.mkdir(parents=True, exist_ok=False)
		shutil.copyfile(seg, output/"seg.dcm")
		if sha256(output/"seg.dcm") != sha256(seg):
			reject("copied original SEG hash differs")
		mesh = output/"roi-portal-connected.msh"
		contract = output/"roi-portal-connected.contract.json"
		split = output/"roi-native-connected-split"
		mesh_script = Path(__file__).with_name("dicom_seg_to_multiregion_tet.py")
		face_script = Path(__file__).with_name("audit_seg_exposed_vessel_faces.py")
		face_audit = output/"source-face-audit.json"
		split_script = Path(__file__).with_name("split_multiregion_tet_for_native.py")
		case_script = Path(__file__).with_name("run_liver_roi_functional_case.py")
		stages = {}
		stages["source_face_audit"] = run("source_face_audit",
			[sys.executable, str(face_script), str(output/"seg.dcm"),
			"1", "3", "4", "--manifest", str(face_audit)], output)
		stages["mesh"] = run("mesh", [sys.executable, str(mesh_script),
			str(output/"seg.dcm"), "1", "3", "4",
			str(mesh), "--roi", "10:20,210:226,101:117", "--max-voxels", "50000"], output)
		stages["split"] = run("split", [sys.executable, str(split_script),
			str(mesh), str(contract), str(split),
			"--max-tetrahedra", "500000"], output)
		validate_case(args.case, output)
		for rank in (1, 2, 4):
			folder = output/f"rank{rank}"
			stages[f"solve_rank{rank}"] = run(f"solve_rank{rank}",
				[sys.executable, str(case_script), "--case", str(args.case),
				"--data-dir", str(output), "--ranks", str(rank),
				"--solver", str(solver), "--output", str(folder/"result.json"),
				"--field-output-dir", str(folder/"fields")], output)
		validate_fresh_series(args.case, output, solver, require_summary=False)
		summary = {"schema_version": 2,
			"kind": "patient_derived_liver_roi_raw_to_native_fem_functional_only",
			"full_liver_case": False, "physiological_validation": False,
			"case_file_sha256": sha256(args.case),
			"source_seg_sha256": sha256(seg),
			"source_face_audit_sha256": sha256(face_audit),
			"multiregion_mesh_sha256": sha256(mesh),
			"split_manifest_sha256": sha256(split/"manifest.json"),
			"solver_binary_sha256": sha256(solver),
			"mpi_ranks_tested": [1, 2, 4],
			"stage_metrics": stages,
			"resource_measurement_scope": "GNU time peak RSS for each stage command/launcher; not the sum of MPI rank memory or a GPU allocation measurement",
			"limitations": "Local original-SEG ROI with deliberately artificial BC and material; no full liver or physiological validation"}
		with (output/"summary.json").open("x", encoding="utf-8") as stream:
			stream.write(json.dumps(summary, indent=2, sort_keys=True)+"\n")
		validate_fresh_series(args.case, output, solver)
	except (AttributeError, IndexError, KeyError, OSError, ValueError,
			json.JSONDecodeError) as error:
		print(f"liver ROI raw-to-result rejected: {error}", file=sys.stderr)
		return 2
	print(f"liver ROI raw-to-result: PASS output={output} (functional only)")
	return 0


if __name__ == "__main__":
	sys.exit(main())
