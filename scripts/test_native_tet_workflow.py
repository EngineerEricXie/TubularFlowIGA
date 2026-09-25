#!/usr/bin/env python3
"""Smoke the explicit volume and optional surface routes of the native workflow."""

import argparse
import fcntl
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from test_native_tet_hydraulic_cli import compare_snapshots, vtk_snapshot


def invoke(command, expected_success=True):
	result = subprocess.run(command, capture_output=True, text=True,
		timeout=120, check=False)
	if (result.returncode == 0) != expected_success:
		raise RuntimeError(f"unexpected workflow result: {' '.join(map(str, command))}\n"
			f"{result.stdout}\n{result.stderr}")
	return result


def check_evidence_classification(status, numerical, visualization):
	if status.get("evidence_classification") != {
		"numerical_acceptance": numerical,
		"physical_acceptance": "not_established",
		"demonstration_visualization": visualization}:
		raise RuntimeError("native workflow evidence classification is misleading")


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--solver", type=Path, required=True)
	parser.add_argument("--ftetwild", type=Path)
	parser.add_argument("--gmsh", action="store_true",
		help="also exercise the Gmsh surface fallback (requires --ftetwild)")
	args = parser.parse_args()
	if args.gmsh and not args.ftetwild:
		parser.error("--gmsh currently requires --ftetwild for the shared surface fixture")
	root = Path(__file__).resolve().parents[1]
	solver = args.solver.resolve(strict=True)
	workflow_script = root/"scripts/run_native_tet_workflow.py"
	fixture = root/"solvers/cpu/tests/data/native_tet_hydraulic_star.json"
	with tempfile.TemporaryDirectory(prefix="iga-native-workflow-") as temporary:
		work = Path(temporary)
		volume_config = {"schema_version": 1,
			"backend": "native_tet_p2p1_ale_hydraulic", "input_route": "volume",
			"input_file": str(fixture.with_suffix(".msh")),
			"case_file": str(fixture), "output_directory": "volume_result",
			"mpi_ranks": 1}
		config_path = work/"volume_workflow.json"
		config_path.write_text(json.dumps(volume_config))
		command = [sys.executable, str(workflow_script), str(config_path),
			"--solver", str(solver)]
		invoke(command)
		status = json.loads((work/"volume_result/workflow_status.json").read_text())
		if status["state"] != "passed" or "mesh_sha256" not in status:
			raise RuntimeError("volume workflow did not record a passed result")
		check_evidence_classification(status, "functional_checks_only", "available")
		vtk_snapshot(work/"volume_result/fields", 2,
			expected_displacement_x=0.0001)
		paused_config = dict(volume_config, output_directory="paused_result")
		paused_path = work/"paused_workflow.json"
		paused_path.write_text(json.dumps(paused_config))
		paused_command = [sys.executable, str(workflow_script), str(paused_path),
			"--solver", str(solver)]
		invoke([*paused_command, "--stop-after-step", "1"])
		paused = work/"paused_result"
		paused_status = json.loads((paused/"workflow_status.json").read_text())
		if (paused_status["state"] != "partial" or paused_status["accepted_steps"] != 1
				or (paused/"fields/step_2").exists()):
			raise RuntimeError("planned partial workflow did not stop at step 1")
		check_evidence_classification(paused_status, "not_established", "available")
		changed_ranks = dict(paused_config, mpi_ranks=2)
		(work/"changed_ranks.json").write_text(json.dumps(changed_ranks))
		invoke([sys.executable, str(workflow_script), str(work/"changed_ranks.json"),
			"--solver", str(solver), "--resume"], expected_success=False)
		if json.loads((paused/"workflow_status.json").read_text())["state"] != "partial":
			raise RuntimeError("rejected resume modified the accepted partial result")
		paused_mesh = paused/"volume.msh"
		original_mesh = paused_mesh.read_bytes()
		paused_mesh.write_bytes(original_mesh+b"\n")
		invoke([*paused_command, "--resume"], expected_success=False)
		paused_mesh.write_bytes(original_mesh)
		if json.loads((paused/"workflow_status.json").read_text())["state"] != "partial":
			raise RuntimeError("tampered-mesh resume modified the partial result")
		with (paused/"workflow.lock").open("r+") as lock_file:
			fcntl.flock(lock_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
			invoke([*paused_command, "--resume"], expected_success=False)
		orphan = paused/"fields/step_2"
		orphan.mkdir()
		(orphan/"snapshot.pvtu").write_text("orphan")
		invoke([*paused_command, "--resume"], expected_success=False)
		(orphan/"snapshot.pvtu").unlink()
		orphan.rmdir()
		if json.loads((paused/"workflow_status.json").read_text())["state"] != "partial":
			raise RuntimeError("concurrent/orphan resume modified the partial result")
		interrupted = json.loads((paused/"workflow_status.json").read_text())
		interrupted["state"] = "running"
		interrupted["phase"] = "simulation"
		(paused/"workflow_status.json").write_text(json.dumps(interrupted))
		invoke([*paused_command, "--resume"])
		resumed_status = json.loads((paused/"workflow_status.json").read_text())
		if (resumed_status["state"] != "passed"
				or resumed_status["accepted_steps"] != 2
				or resumed_status["resume_runs"] != 1):
			raise RuntimeError("planned partial workflow did not finish on resume")
		check_evidence_classification(resumed_status, "functional_checks_only", "available")
		compare_snapshots(work/"volume_result/fields", paused/"fields", 2,
			expected_displacement_x=0.0001)
		invoke([*paused_command, "--resume"], expected_success=False)
		invoke(command, expected_success=False)
		if json.loads((work/"volume_result/workflow_status.json").read_text())["state"] != "passed":
			raise RuntimeError("workflow overwrote an existing result")
		unsupported = dict(volume_config, input_route="centerline",
			output_directory="unsupported_result")
		(work/"unsupported.json").write_text(json.dumps(unsupported))
		invoke([sys.executable, str(workflow_script), str(work/"unsupported.json"),
			"--solver", str(solver)], expected_success=False)
		if (work/"unsupported_result").exists():
			raise RuntimeError("unsupported route created output")
		wrong_backend = dict(volume_config, backend="dolfin",
			output_directory="external_backend_result")
		(work/"external_backend.json").write_text(json.dumps(wrong_backend))
		invoke([sys.executable, str(workflow_script), str(work/"external_backend.json"),
			"--solver", str(solver)], expected_success=False)
		if (work/"external_backend_result").exists():
			raise RuntimeError("external FEM backend was not rejected")
		bad_case = json.loads(fixture.read_text())
		bad_case["boundaries"]["outlet_label"] = 3
		(work/"bad_case.json").write_text(json.dumps(bad_case))
		bad_case_config = dict(volume_config, case_file="bad_case.json",
			output_directory="failed_result")
		(work/"bad_case_workflow.json").write_text(json.dumps(bad_case_config))
		invoke([sys.executable, str(workflow_script),
			str(work/"bad_case_workflow.json"), "--solver", str(solver)],
			expected_success=False)
		failed = work/"failed_result"
		failed_status = json.loads((failed/"workflow_status.json").read_text())
		if (failed_status["state"] != "failed"
				or failed_status["phase"] != "input_preflight"
				or not (failed/"input_preflight.log").is_file()
				or (failed/"simulation.log").exists()):
			raise RuntimeError("failed native workflow lost its diagnostic evidence")
		check_evidence_classification(failed_status, "not_established", "not_available")
		invoke([sys.executable, str(workflow_script),
			str(work/"bad_case_workflow.json"), "--solver", str(solver), "--resume"],
			expected_success=False)
		if args.ftetwild:
			ftetwild = args.ftetwild.resolve(strict=True)
			surface = work/"pipe.vtp"
			invoke([sys.executable,
				str(root/"scripts/generate_circular_pipe_surface.py"),
				str(surface), "--length-m", "0.03", "--radius-m", "0.005",
				"--circumferential-segments", "8", "--axial-segments", "2"])
			case = json.loads(fixture.read_text())
			case["boundaries"]["wall_labels"] = [0]
			case["fluid"].update({"density_kg_m3": 1.0,
				"dynamic_viscosity_pa_s": 0.001,
				"initial_velocity_x_m_s": -0.01,
				"nonlinear_tolerance": 1e-8, "maximum_iterations": 20})
			case["source"] = {"capacitance_m3_pa": 0.01,
				"resistance_pa_s_m3": 10.0,
				"prescribed_flow_m3_s": 1e-6, "initial_pressure_pa": 0.0}
			case["terminal_rcr"] = {"proximal_resistance_pa_s_m3": 1.0,
				"distal_resistance_pa_s_m3": 10.0,
				"capacitance_m3_pa": 0.01, "distal_pressure_pa": 0.0,
				"initial_pressure_pa": 0.0}
			case["motion"] = {"kind": "fixed", "speed_x_m_s": 0.0}
			case["coupling"]["method"] = "explicit"
			case["coupling"]["maximum_iterations"] = 1
			case["time"] = {"dt_s": 0.01, "steps": 1}
			(work/"pipe_case.json").write_text(json.dumps(case))
			surface_config = dict(volume_config, input_route="surface",
				input_file=surface.name, case_file="pipe_case.json",
				output_directory="surface_result", mesher={"kind": "ftetwild",
					"target_size_m": 0.007, "envelope_m": 0.0001,
					"executable": str(ftetwild)})
			(work/"surface_workflow.json").write_text(json.dumps(surface_config))
			invoke([sys.executable, str(workflow_script),
				str(work/"surface_workflow.json"), "--solver", str(solver)])
			result = work/"surface_result"
			status = json.loads((result/"workflow_status.json").read_text())
			manifest = json.loads((result/"mesh_manifest.json").read_text())
			if status["state"] != "passed" or not manifest["gates"]["passed"]:
				raise RuntimeError("surface workflow did not pass")
			vtk_snapshot(result/"fields", 1, expected_points=None,
				expected_cells=manifest["volume_mesh"]["elements"],
				expected_first_cell=None, expected_displacement_x=0.0)
			if args.gmsh:
				gmsh_config = dict(surface_config, output_directory="gmsh_result",
					mesher={"kind": "gmsh", "target_size_m": 0.007})
				(work/"gmsh_workflow.json").write_text(json.dumps(gmsh_config))
				invoke([sys.executable, str(workflow_script),
					str(work/"gmsh_workflow.json"), "--solver", str(solver)])
				gmsh_result = work/"gmsh_result"
				gmsh_status = json.loads((gmsh_result/"workflow_status.json").read_text())
				gmsh_manifest = json.loads((gmsh_result/"mesh_manifest.json").read_text())
				if gmsh_status["state"] != "passed" or not gmsh_manifest["gates"]["passed"]:
					raise RuntimeError("Gmsh surface workflow did not pass")
				vtk_snapshot(gmsh_result/"fields", 1, expected_points=None,
					expected_cells=gmsh_manifest["volume_mesh"]["elements"],
					expected_first_cell=None, expected_displacement_x=0.0)
		print("native_tet_workflow: PASS volume, interrupted restart, fail-closed routing, no overwrite"
			+(", fTetWild surface" if args.ftetwild else "")
			+(", Gmsh fallback" if args.gmsh else ""))


if __name__ == "__main__":
	main()
