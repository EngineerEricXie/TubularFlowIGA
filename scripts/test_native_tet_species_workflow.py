#!/usr/bin/env python3
"""Exercise the explicit native prescribed-velocity species workflow route."""

import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from test_native_tet_species_cli import snapshot, summary


def invoke(command, success=True):
	result = subprocess.run(command, capture_output=True, text=True, timeout=180)
	if (result.returncode == 0) != success:
		raise RuntimeError(f"unexpected workflow result: {command}\n"
			f"{result.stdout}\n{result.stderr}")
	return result


def status(path, state, accepted):
	value = json.loads((path/"workflow_status.json").read_text())
	if (value.get("state") != state or value.get("accepted_steps") != accepted
			or value.get("backend") != "native_tet_p1_prescribed_species"
			or value.get("evidence_classification") != {
				"numerical_acceptance": "functional_checks_only" if state == "passed"
					else "not_established",
				"physical_acceptance": "not_established",
				"demonstration_visualization": "available"}):
		raise RuntimeError("species workflow status or classification differs")
	return value


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--solver", type=Path, required=True)
	parser.add_argument("--ftetwild", type=Path)
	parser.add_argument("--gmsh", action="store_true")
	args = parser.parse_args()
	if args.gmsh and not args.ftetwild:
		parser.error("--gmsh requires --ftetwild for the shared surface fixture")
	root = Path(__file__).resolve().parents[1]
	solver = args.solver.resolve(strict=True)
	mesh = (root/"solvers/cpu/tests/data/native_tet_hydraulic_star.msh").resolve()
	workflow = root/"scripts/run_native_tet_workflow.py"
	with tempfile.TemporaryDirectory(prefix="native-species-workflow-") as temporary:
		work = Path(temporary)
		case = {"schema_version": 1, "mesh_file": str(mesh), "species_id": "tracer",
			"initial_concentration_mol_m3": 2.0, "diffusivity_m2_s": 0.0,
			"source_mol_m3_s": 0.5, "first_order_decay_rate_s_inv": 1.0,
			"fluid_velocity_m_s": [1.0, 0.0, 0.0],
			"mesh_velocity_m_s": [1.0, 0.0, 0.0],
			"inflow_concentration_by_label_mol_m3": {}, "monotone": False,
			"time": {"dt_s": 0.1, "steps": 2}}
		(work/"species.json").write_text(json.dumps(case))
		config = {"schema_version": 1, "backend": "native_tet_p1_prescribed_species",
			"input_route": "volume", "input_file": str(mesh),
			"case_file": "species.json", "output_directory": "full", "mpi_ranks": 2}
		def command(name, options=()):
			path = work/f"{name}.json"
			path.write_text(json.dumps({**config, "output_directory": name}))
			return [sys.executable, str(workflow), str(path), "--solver", str(solver),
				*options]
		full = command("full")
		invoke(full)
		status(work/"full", "passed", 2)
		summary(work/"full/fields", work/"full/case.json", work/"full/volume.msh", 2, 2)
		expected = ((2.0+0.1*0.5)/1.1+0.1*0.5)/1.1
		if max(abs(value-expected) for value in snapshot(work/"full/fields", 2, 5, 4,
				0.2)) > 1e-10:
			raise RuntimeError("species workflow source/ALE field differs")
		invoke(full, success=False)
		paused = command("paused", ("--stop-after-step", "1"))
		invoke(paused)
		status(work/"paused", "partial", 1)
		if (work/"paused/fields/run_summary.json").exists():
			raise RuntimeError("partial species workflow has a final summary")
		changed = dict(config, output_directory="paused", mpi_ranks=1)
		(work/"changed.json").write_text(json.dumps(changed))
		invoke([sys.executable, str(workflow), str(work/"changed.json"),
			"--solver", str(solver), "--resume"], success=False)
		state = work/"paused/fields/run_state.json"
		state.rename(work/"paused/fields/run_state.json.hold")
		invoke(command("paused", ("--resume",)), success=False)
		(work/"paused/fields/run_state.json.hold").rename(state)
		invoke(command("paused", ("--resume",)))
		status(work/"paused", "passed", 2)
		summary(work/"paused/fields", work/"paused/case.json",
			work/"paused/volume.msh", 2, 2)
		left = snapshot(work/"full/fields", 2, 5, 4, 0.2)
		right = snapshot(work/"paused/fields", 2, 5, 4, 0.2)
		if max(abs(a-b) for a, b in zip(left, right)) > 1e-10:
			raise RuntimeError("species workflow cross-process field differs")
		invoke(command("paused", ("--resume",)), success=False)
		reservoir_case = dict(case, source_mol_m3_s=0.0,
			first_order_decay_rate_s_inv=0.0, monotone=True,
			finite_wall_reservoir={"boundary_label": 1,
				"transfer_coefficient_m_s": 0.5, "volume_m3": 1.0,
				"initial_concentration_mol_m3": 0.0})
		(work/"reservoir.json").write_text(json.dumps(reservoir_case))
		reservoir_config = dict(config, case_file="reservoir.json")
		def reservoir_command(name, options=()):
			path = work/f"{name}.json"
			path.write_text(json.dumps(dict(reservoir_config, output_directory=name)))
			return [sys.executable, str(workflow), str(path), "--solver", str(solver),
				*options]
		invoke(reservoir_command("reservoir_full"))
		invoke(reservoir_command("reservoir_pause", ("--stop-after-step", "1")))
		status(work/"reservoir_pause", "partial", 1)
		invoke(reservoir_command("reservoir_pause", ("--resume",)))
		for name in ("reservoir_full", "reservoir_pause"):
			output = work/name
			status(output, "passed", 2)
			summary(output/"fields", output/"case.json", output/"volume.msh", 2, 2)
			if not (output/"fields/tissue_step_2.bin").is_file():
				raise RuntimeError("species workflow did not retain tissue checkpoint")
		full_tissue = json.loads((work/"reservoir_full/fields/run_summary.json").read_text())
		pause_tissue = json.loads((work/"reservoir_pause/fields/run_summary.json").read_text())
		if (not 0 < full_tissue["last_tissue_amount_mol"] < 1.0
				or abs(full_tissue["last_tissue_amount_mol"]
					-pause_tissue["last_tissue_amount_mol"]) > 1e-12
				or abs(full_tissue["last_combined_balance_defect_mol_s"]) > 1e-10
				or max(abs(a-b) for a, b in zip(
					snapshot(work/"reservoir_full/fields", 2, 5, 4, 0.2),
					snapshot(work/"reservoir_pause/fields", 2, 5, 4, 0.2))) > 1e-10):
			raise RuntimeError("species workflow paired reservoir restart differs")
		multi_case = dict(case, source_mol_m3_s=0.0,
			first_order_decay_rate_s_inv=0.0, monotone=True,
			finite_wall_reservoirs_by_label={
				"1": {"transfer_coefficient_m_s": 0.5, "volume_m3": 1.0,
					"initial_concentration_mol_m3": 0.0},
				"2": {"transfer_coefficient_m_s": 0.5, "volume_m3": 1.0,
					"initial_concentration_mol_m3": 3.0}})
		(work/"multi.json").write_text(json.dumps(multi_case))
		multi_config = dict(config, case_file="multi.json")
		def multi_command(name, options=()):
			path = work/f"{name}.json"
			path.write_text(json.dumps(dict(multi_config, output_directory=name)))
			return [sys.executable, str(workflow), str(path), "--solver", str(solver),
				*options]
		invoke(multi_command("multi_full"))
		invoke(multi_command("multi_pause", ("--stop-after-step", "1")))
		invoke(multi_command("multi_pause", ("--resume",)))
		for name in ("multi_full", "multi_pause"):
			output = work/name
			status(output, "passed", 2)
			summary(output/"fields", output/"case.json", output/"volume.msh", 2, 2)
		left_multi = json.loads((work/"multi_full/fields/run_summary.json").read_text())
		right_multi = json.loads((work/"multi_pause/fields/run_summary.json").read_text())
		if (left_multi["tissue_regions_by_label"] != right_multi["tissue_regions_by_label"]
				or set(left_multi["tissue_regions_by_label"]) != {"1", "2"}
				or abs(left_multi["last_combined_balance_defect_mol_s"]) > 1e-10
				or max(abs(a-b) for a, b in zip(
					snapshot(work/"multi_full/fields", 2, 5, 4, 0.2),
					snapshot(work/"multi_pause/fields", 2, 5, 4, 0.2))) > 1e-10):
			raise RuntimeError("species workflow multi-region paired restart differs")
		wrong = dict(config, backend="dolfin", output_directory="external")
		(work/"external.json").write_text(json.dumps(wrong))
		invoke([sys.executable, str(workflow), str(work/"external.json"),
			"--solver", str(solver)], success=False)
		if (work/"external").exists():
			raise RuntimeError("unsupported FEM backend created output")
		wrong_solver = dict(config, output_directory="wrong_solver")
		(work/"wrong_solver.json").write_text(json.dumps(wrong_solver))
		invoke([sys.executable, str(workflow), str(work/"wrong_solver.json"),
			"--solver", "/bin/true"], success=False)
		wrong_status = json.loads((work/"wrong_solver/workflow_status.json").read_text())
		if (wrong_status["state"] != "failed" or wrong_status["phase"] != "input_preflight"
				or (work/"wrong_solver/simulation.log").exists()):
			raise RuntimeError("wrong solver was not rejected before simulation")
		if args.ftetwild:
			surface = work/"pipe.vtp"
			invoke([sys.executable, str(root/"scripts/generate_circular_pipe_surface.py"),
				str(surface), "--length-m", "0.03", "--radius-m", "0.005",
				"--circumferential-segments", "8", "--axial-segments", "2"])
			pipe_case = dict(case, fluid_velocity_m_s=[0.0, 0.0, 0.0],
				mesh_velocity_m_s=[0.0, 0.0, 0.0], source_mol_m3_s=0.0,
				monotone=True,
				wall_exchange_by_label={"0": {"transfer_coefficient_m_s": 0.5,
					"external_concentration_mol_m3": 0.0}},
				time={"dt_s": 0.1, "steps": 1})
			(work/"pipe.json").write_text(json.dumps(pipe_case))
			pipe = dict(config, input_route="surface", input_file="pipe.vtp",
				case_file="pipe.json", output_directory="surface",
				mesher={"kind": "ftetwild", "target_size_m": 0.007,
					"envelope_m": 0.0001,
					"executable": str(args.ftetwild.resolve(strict=True))})
			(work/"surface.json").write_text(json.dumps(pipe))
			invoke([sys.executable, str(workflow), str(work/"surface.json"),
				"--solver", str(solver)])
			status(work/"surface", "passed", 1)
			summary(work/"surface/fields", work/"surface/case.json",
				work/"surface/volume.msh", 2, 1)
			if json.loads((work/"surface/fields/run_summary.json").read_text())[
					"last_outward_wall_exchange_mol_s"] <= 0:
				raise RuntimeError("fTetWild surface wall sink was not active")
			pipe_multi_case = dict(pipe_case)
			pipe_multi_case.pop("wall_exchange_by_label")
			pipe_multi_case["finite_wall_reservoirs_by_label"] = {
				"0": {"transfer_coefficient_m_s": 1e-4,
					"volume_m3": 1e-6, "initial_concentration_mol_m3": 0.0},
				"1": {"transfer_coefficient_m_s": 1e-4,
					"volume_m3": 1e-6, "initial_concentration_mol_m3": 3.0}}
			(work/"pipe_multi.json").write_text(json.dumps(pipe_multi_case))
			pipe_multi = dict(pipe, case_file="pipe_multi.json",
				output_directory="surface_multi")
			(work/"surface_multi.json").write_text(json.dumps(pipe_multi))
			invoke([sys.executable, str(workflow), str(work/"surface_multi.json"),
				"--solver", str(solver)])
			status(work/"surface_multi", "passed", 1)
			summary(work/"surface_multi/fields", work/"surface_multi/case.json",
				work/"surface_multi/volume.msh", 2, 1)
			pipe_multi_summary = json.loads((work/"surface_multi/fields/run_summary.json").read_text())
			if (set(pipe_multi_summary["tissue_regions_by_label"]) != {"0", "1"}
					or pipe_multi_summary["tissue_regions_by_label"]["0"]["amount_mol"] <= 0
					or pipe_multi_summary["tissue_regions_by_label"]["1"]["amount_mol"] >= 3e-6
					or abs(pipe_multi_summary["last_combined_balance_defect_mol_s"]) > 1e-10):
				raise RuntimeError("fTetWild two-region tissue exchange differs")
			if args.gmsh:
				gmsh = dict(pipe, output_directory="surface_gmsh",
					mesher={"kind": "gmsh", "target_size_m": 0.007})
				(work/"surface_gmsh.json").write_text(json.dumps(gmsh))
				invoke([sys.executable, str(workflow), str(work/"surface_gmsh.json"),
					"--solver", str(solver)])
				status(work/"surface_gmsh", "passed", 1)
				summary(work/"surface_gmsh/fields", work/"surface_gmsh/case.json",
					work/"surface_gmsh/volume.msh", 2, 1)
				if json.loads((work/"surface_gmsh/fields/run_summary.json").read_text())[
						"last_outward_wall_exchange_mol_s"] <= 0:
					raise RuntimeError("Gmsh surface wall sink was not active")
	print("native_tet_species_workflow: PASS volume, ALE/source/decay, one/multi-region reservoir paired restart, pause/resume, "
		"backend rejection"+("; fTetWild surface" if args.ftetwild else "")
		+("; Gmsh surface" if args.gmsh else "")
		+"; labelled wall and fTetWild multi-region exchange")


if __name__ == "__main__":
	main()
