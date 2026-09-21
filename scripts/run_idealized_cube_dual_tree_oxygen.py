#!/usr/bin/env python3
"""Rerun passive-tracer transport on a checked artificial frozen-FSI flow run."""

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def digest(path):
	result = hashlib.sha256()
	with path.open("rb") as stream:
		for block in iter(lambda: stream.read(1024*1024), b""):
			result.update(block)
	return result.hexdigest()


def require(condition, message):
	if not condition:
		raise ValueError(message)


def validate(case, path):
	required = {"schema_version", "kind", "physiological_validation",
		"flow_case_file", "frozen_post_fsi_velocity", "tissue_porosity_assumed",
		"tissue_consumption_mol_m3_s", "inlet_concentration_mol_m3",
		"diffusivity_m2_s", "time_step_s", "minimum_steps", "maximum_steps",
		"venous_breakthrough_fraction", "limitations"}
	require(set(case) == required and case["schema_version"] == 1 and
		case["kind"] == "idealized_cube_passive_oxygen_functional_only" and
		case["physiological_validation"] is False and
		case["frozen_post_fsi_velocity"] is True and
		case["tissue_porosity_assumed"] == 1.0 and
		case["tissue_consumption_mol_m3_s"] == 0.0,
		"oxygen case schema or functional-only assumptions differ")
	for key in ("inlet_concentration_mol_m3", "diffusivity_m2_s",
		"time_step_s", "venous_breakthrough_fraction"):
		require(type(case[key]) in (int, float) and math.isfinite(case[key]),
			f"oxygen {key} is invalid")
	require(case["inlet_concentration_mol_m3"] > 0 and
		case["diffusivity_m2_s"] >= 0 and case["time_step_s"] > 0 and
		0 < case["venous_breakthrough_fraction"] < 1 and
		type(case["minimum_steps"]) is int and
		type(case["maximum_steps"]) is int and
		2 <= case["minimum_steps"] <= case["maximum_steps"] <= 1000,
		"oxygen artificial concentration/diffusion/time contract is invalid")
	require(type(case["flow_case_file"]) is str and
		Path(case["flow_case_file"]).name == case["flow_case_file"] and
		type(case["limitations"]) is str and case["limitations"],
		"oxygen source case or limitations are invalid")
	flow_case = path.parent/case["flow_case_file"]
	require(flow_case.is_file(), "oxygen referenced flow case is absent")
	return flow_case


def run(command, log):
	completed = subprocess.run(command, stdout=subprocess.PIPE,
		stderr=subprocess.STDOUT, text=True, check=False)
	log.write_text(completed.stdout, encoding="utf-8")
	if completed.returncode:
		raise ValueError(f"stage {command[0]} failed with exit {completed.returncode}: "
			f"{completed.stdout[-1600:]}")
	return completed.stdout


def check_ledger(path, summary, case):
	with path.open(newline="", encoding="utf-8") as stream:
		rows = list(csv.DictReader(stream))
	require(len(rows) == summary["accepted_steps"] and
		case["minimum_steps"] <= len(rows) <= case["maximum_steps"],
		"oxygen step ledger count differs")
	for index, row in enumerate(rows, 1):
		require(int(row["step"]) == index and
			abs(float(row["time_s"])-index*case["time_step_s"]) <= 1e-10,
			"oxygen step sequence or clock differs")
		for key in ("artery_inventory_mol", "tissue_inventory_mol",
			"vein_inventory_mol", "venous_outlet_concentration_mol_m3"):
			require(math.isfinite(float(row[key])) and float(row[key]) >= 0,
				f"oxygen {key} is invalid")
		for key in ("artery_reason", "tissue_reason", "vein_reason"):
			require(int(row[key]) > 0, f"oxygen {key} did not converge")
		incoming = max(abs(float(row["artery_to_tissue_mol_s"])), 1e-18)
		require(abs(float(row["global_balance_defect_mol_s"])) <=
			1e-14+1e-7*incoming, "oxygen global amount does not balance")
	require(float(rows[-1]["venous_outlet_concentration_mol_m3"]) >=
		case["venous_breakthrough_fraction"]*case["inlet_concentration_mol_m3"],
		"oxygen venous outlet breakthrough was not reached")
	return rows


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("flow_run", type=Path)
	parser.add_argument("output_directory", type=Path)
	parser.add_argument("--case", type=Path,
		default=Path("cases/idealized_cube_dual_tree_oxygen.json"))
	parser.add_argument("--split-directory", type=Path,
		help="reuse matching generated submeshes after hash verification")
	parser.add_argument("--solver", type=Path,
		default=Path("solvers/cpu/native_tet_cube_oxygen"))
	parser.add_argument("--ranks", type=int, default=8)
	parser.add_argument("--check-only", action="store_true")
	args = parser.parse_args()
	try:
		case_path = args.case.resolve()
		case = json.loads(case_path.read_text(encoding="utf-8"))
		flow_case = validate(case, case_path)
		run_directory = args.flow_run.resolve()
		flow_summary_path = run_directory/"summary.json"
		flow = json.loads(flow_summary_path.read_text(encoding="utf-8"))
		require(flow.get("kind") ==
			"idealized_cube_dual_tree_quasi_steady_fsi_functional_only" and
			flow.get("physiological_validation") is False and
			flow.get("vessel_wall_fsi") is True and
			flow.get("case_sha256") == digest(flow_case),
			"oxygen flow run or referenced case identity differs")
		require(1 <= args.ranks <= 16, "oxygen local MPI ranks must be 1..16")
		if args.split_directory:
			split = args.split_directory.resolve()
			require(digest(split/"manifest.json") == flow["split_manifest_sha256"],
				"oxygen split manifest differs from frozen flow run")
		if args.check_only:
			print("artificial oxygen transport input: PASS")
			return 0
		solver = args.solver.resolve()
		require(solver.is_file(), "native oxygen solver binary is absent")
		output = args.output_directory.resolve()
		require(not output.exists() and output.parent.is_dir(),
			"oxygen output exists or parent directory is missing")
		with tempfile.TemporaryDirectory(prefix="cube-oxygen-", dir=output.parent) as temp:
			stage = Path(temp)
			submeshes = stage/"submeshes"
			if args.split_directory:
				shutil.copytree(split, submeshes)
			else:
				generator = Path(__file__).with_name("generate_idealized_cube_dual_tree.py")
				geometry = stage/"geometry"
				run([sys.executable, str(generator), str(geometry), "--case",
					str(flow_case)], stage/"geometry.log")
				splitter = Path(__file__).with_name("split_multiregion_tet_for_native.py")
				run([sys.executable, str(splitter), str(geometry)+".volume.msh",
					str(geometry)+".contract.json", str(submeshes)], stage/"split.log")
			require(digest(submeshes/"manifest.json") == flow["split_manifest_sha256"],
				"oxygen rebuilt mesh does not match frozen flow run")
			prepare = Path(__file__).with_name("prepare_idealized_cube_oxygen_fields.py")
			fields = stage/"frozen_fields"
			run([sys.executable, str(prepare), str(run_directory),
				str(submeshes), str(fields)], stage/"field_preparation.log")
			manifest = json.loads((submeshes/"manifest.json").read_text(encoding="utf-8"))
			labels = manifest["interface_boundary_labels"]
			command = [str(solver), str(submeshes/"arterial_lumen.msh"),
				str(submeshes/"fixed_darcy_tissue.msh"),
				str(submeshes/"venous_lumen.msh"),
				str(fields/"arterial_fsi_state.txt"),
				str(fields/"tissue_rt0_face_flow.txt"),
				str(fields/"venous_fsi_state.txt"), "207",
				*[str(labels[f"arterial_terminal_{i}"]) for i in range(4)],
				*[str(labels[f"venous_terminal_{i}"]) for i in range(4)],
				"208", str(case["diffusivity_m2_s"]),
				str(case["inlet_concentration_mol_m3"]), str(case["time_step_s"]),
				str(case["minimum_steps"]), str(case["maximum_steps"]),
				str(case["venous_breakthrough_fraction"]),
				str(stage/"transport")]
			if args.ranks > 1:
				command = ["mpiexec", "-n", str(args.ranks), *command]
			stdout = run(command, stage/"transport.log")
			require(stdout.count("native_cube_oxygen: PASS") == 1,
				"oxygen solver did not report one PASS result")
			summary = json.loads((stage/"transport"/"summary.json").read_text(
				encoding="utf-8"))
			require(summary.get("kind") == case["kind"] and
				summary.get("physiological_validation") is False and
				summary.get("ranks") == args.ranks and
				summary.get("inlet_concentration_mol_m3") ==
					case["inlet_concentration_mol_m3"] and
				summary.get("diffusivity_m2_s") == case["diffusivity_m2_s"],
				"oxygen native summary differs from input")
			rows = check_ledger(stage/"transport"/"ledger.csv", summary, case)
			(stage/"provenance.json").write_text(json.dumps({
				"schema_version": 1,
				"kind": case["kind"],
				"physiological_validation": False,
				"case_sha256": digest(case_path),
				"flow_summary_sha256": digest(flow_summary_path),
				"split_manifest_sha256": digest(submeshes/"manifest.json"),
				"frozen_field_manifest_sha256": digest(fields/"source_manifest.json"),
				"solver_sha256": digest(solver),
				"ranks": args.ranks,
				"accepted_steps": len(rows),
				"maximum_abs_global_balance_defect_mol_s": max(abs(float(
					row["global_balance_defect_mol_s"])) for row in rows),
				"limitations": case["limitations"],
			}, indent=2, sort_keys=True)+"\n", encoding="utf-8")
			os.replace(stage, output)
		print("artificial oxygen transport: PASS "
			f"steps={summary['accepted_steps']} "
			f"outlet={summary['venous_outlet_concentration_mol_m3']:.12g} "
			f"output={output}")
		return 0
	except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as error:
		print(f"artificial oxygen transport rejected: {error}", file=sys.stderr)
		return 2


if __name__ == "__main__":
	sys.exit(main())
