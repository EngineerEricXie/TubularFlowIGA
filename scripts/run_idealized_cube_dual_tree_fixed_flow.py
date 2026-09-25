#!/usr/bin/env python3
"""Rerun the conservative artificial dual-tree small-strain FSI case."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from generate_idealized_cube_dual_tree import validate


def reject(message):
	raise ValueError(message)


def sha256(path):
	digest = hashlib.sha256()
	with path.open("rb") as stream:
		for block in iter(lambda: stream.read(1024*1024), b""):
			digest.update(block)
	return digest.hexdigest()


def run(command, log):
	completed = subprocess.run(command, text=True, stdout=subprocess.PIPE,
		stderr=subprocess.STDOUT, check=False)
	log.write_text(completed.stdout, encoding="utf-8")
	if completed.returncode:
		reject(f"stage {command[0]} failed with exit {completed.returncode}: "
			f"{completed.stdout[-1200:]}")
	return completed.stdout


def parse_fields(line):
	result = {}
	for token in line.split()[1:]:
		if "=" not in token:
			continue
		name, value = token.split("=", 1)
		result[name] = float(value)
	return result


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("output_directory", type=Path)
	parser.add_argument("--case", type=Path,
		default=Path("cases/idealized_cube_dual_tree.json"))
	parser.add_argument("--solver", type=Path,
		default=Path("solvers/cpu/native_tet_cube_dual_tree_fixed_flow"))
	parser.add_argument("--ranks", type=int, default=1)
	parser.add_argument("--check-only", action="store_true")
	args = parser.parse_args()
	try:
		case = validate(json.loads(args.case.read_text(encoding="utf-8")))
		if not 1 <= args.ranks <= 16:
			reject("local MPI ranks must be 1..16")
		if args.check_only:
			print("idealized dual-tree fixed-flow input: PASS")
			return 0
		solver = args.solver.resolve()
		if not solver.is_file():
			reject("native solver binary missing; build native_tet_cube_dual_tree_fixed_flow")
		output = args.output_directory.resolve()
		if output.exists() or not output.parent.is_dir():
			reject("output exists or its parent directory does not exist")
		with tempfile.TemporaryDirectory(prefix="cube-dual-tree-", dir=output.parent) as temporary:
			work = Path(temporary)
			generator = Path(__file__).with_name("generate_idealized_cube_dual_tree.py")
			mesh_prefix = work/"geometry"
			run([sys.executable, str(generator), str(mesh_prefix),
				"--case", str(args.case.resolve())], work/"geometry.log")
			mesh = Path(str(mesh_prefix)+".volume.msh")
			contract = Path(str(mesh_prefix)+".contract.json")
			splitter = Path(__file__).with_name("split_multiregion_tet_for_native.py")
			split_dir = work/"submeshes"
			run([sys.executable, str(splitter), str(mesh), str(contract),
				str(split_dir)], work/"split.log")
			split = json.loads((split_dir/"manifest.json").read_text(encoding="utf-8"))
			labels = split["interface_boundary_labels"]
			params = case["functional_parameters"]
			command = [str(solver), str(split_dir/"arterial_lumen.msh"),
				str(split_dir/"fixed_darcy_tissue.msh"),
				str(split_dir/"venous_lumen.msh"),
				"207", str(labels["arterial_fluid_wall"]),
				*[str(labels[f"arterial_terminal_{i}"]) for i in range(4)],
				str(labels["venous_fluid_wall"]),
				*[str(labels[f"venous_terminal_{i}"]) for i in range(4)],
				"208", str(params["inlet_speed_m_s"]),
				str(params["fluid_density_kg_m3"]),
				str(params["fluid_viscosity_pa_s"]),
				str(params["darcy_mobility_m2_pa_s"]),
				str(split_dir/"arterial_wall.msh"),
				str(split_dir/"venous_wall.msh"),
				str(labels["arterial_wall_tissue"]),
				str(labels["venous_wall_tissue"]), "210", "218",
				str(params["wall_young_modulus_pa"]),
				str(params["wall_poisson_ratio"]), str(work/"fields")]
			if args.ranks > 1:
				command = ["mpiexec", "-n", str(args.ranks), *command]
			stdout = run(command, work/"fixed_flow.log")
			lines = stdout.splitlines()
			pass_lines = [line for line in lines
				if line.startswith("native_cube_dual_tree_fixed_flow: PASS ")]
			if len(pass_lines) != 1:
				reject("fixed-flow solver did not report exactly one PASS line")
			fields = parse_fields(pass_lines[0])
			terminals = {}
			fsi_terminals = {}
			for line in lines:
				for prefix, destination in (("terminal_fsi_", fsi_terminals),
					("terminal_", terminals)):
					if line.startswith(prefix):
						index = line.split()[0][len(prefix):]
						if not index.isdigit() or index in destination:
							reject("duplicate or invalid terminal ledger index")
						destination[index] = parse_fields(line)
						break
			if set(terminals) != {str(i) for i in range(4)} or \
				set(fsi_terminals) != {str(i) for i in range(4)}:
				reject("initial or FSI per-terminal ledger is incomplete")
			for ledger in (terminals, fsi_terminals):
				for entry in ledger.values():
					for field in ("artery_to_tissue_m3_s", "tissue_to_vein_m3_s",
						"vein_inward_m3_s"):
						if entry.get(field, 0) <= 0:
							reject("terminal flow reversed or absent")
					if abs(entry["tissue_to_vein_m3_s"] - entry["vein_inward_m3_s"]) > \
						1e-6 * max(entry["tissue_to_vein_m3_s"], 1e-18):
						reject("terminal tissue-to-vein ledger does not close")
			if fields.get("coupled_pressure_continuity") != 0 or \
				fields.get("vessel_wall_fsi") != 1 or \
				any(fields.get(reason, 0) <= 0 for reason in
					("artery_linear_reason", "darcy_linear_reason", "vein_linear_reason",
					"fsi_artery_linear_reason", "fsi_darcy_linear_reason",
					"fsi_vein_linear_reason")) or \
				not 0 <= fields.get("fsi_displacement_relative_residual", -1) <= 1e-3 or \
				fields.get("arterial_wall_max_displacement_m", 0) <= 0 or \
				fields.get("venous_wall_max_displacement_m", 0) <= 0:
				reject("fixed-wall stage capability classification changed")
			summary = {"schema_version": 1,
				"kind": "idealized_cube_dual_tree_quasi_steady_fsi_functional_only",
				"physiological_validation": False, "vessel_wall_fsi": True,
				"fsi_model": "small_strain_quasi_steady_partitioned_reference",
				"coupled_pressure_continuity": False,
				"ranks": args.ranks,
				"case_sha256": sha256(args.case.resolve()),
				"solver_sha256": sha256(solver),
				"geometry_manifest_sha256": sha256(Path(str(mesh_prefix)+".json")),
				"split_manifest_sha256": sha256(split_dir/"manifest.json"),
				"flow": fields, "terminals": terminals,
				"fsi_terminals": fsi_terminals,
				"limitations": "Two-way small-strain quasi-steady vessel-wall displacement/traction and moving fluid meshes, with fixed Darcy tissue and conservative one-way hydraulic terminal flow transfer; no cross-domain pressure-continuity iteration, transient FSI or physiological calibration."}
			(work/"summary.json").write_text(json.dumps(summary, indent=2,
				sort_keys=True)+"\n", encoding="utf-8")
			os.replace(work, output)
		print("idealized dual-tree fixed flow: PASS "
			f"ranks={args.ranks} inlet={fields['artery_inlet_m3_s']:.12g} "
			f"outlet={fields['vein_outlet_m3_s']:.12g}")
		return 0
	except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
		print(f"idealized dual-tree fixed flow rejected: {error}", file=sys.stderr)
		return 2


if __name__ == "__main__":
	sys.exit(main())
