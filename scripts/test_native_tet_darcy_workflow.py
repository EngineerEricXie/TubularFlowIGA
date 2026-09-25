#!/usr/bin/env python3
"""Exercise the steady native Darcy surface/volume workflow route."""

import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from test_native_tet_darcy_cli import fields, summary


def invoke(command, success=True):
	result = subprocess.run(command, capture_output=True, text=True, timeout=180)
	if (result.returncode == 0) != success:
		raise RuntimeError(f"unexpected Darcy workflow result: {command}\n"
			f"{result.stdout}\n{result.stderr}")
	return result


def status(directory, state, accepted):
	data = json.loads((directory/"workflow_status.json").read_text())
	if (data.get("state") != state or data.get("accepted_steps") != accepted
			or data.get("backend") != "native_tet_p1_darcy_steady"
			or data.get("evidence_classification") != {
				"numerical_acceptance": "functional_checks_only" if state == "passed"
					else "not_established",
				"physical_acceptance": "not_established",
				"demonstration_visualization": "available" if accepted else "not_available"}):
		raise RuntimeError("Darcy workflow status classification differs")
	return data


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
	with tempfile.TemporaryDirectory(prefix="native-darcy-workflow-") as temporary:
		work = Path(temporary)
		case = {"schema_version": 1, "mesh_file": str(mesh),
			"mobility_m2_pa_s": 2.0, "source_s_inv": 1.0,
			"pressure_by_boundary_label_pa": {"1": 1.0},
			"outward_flux_by_boundary_label_m_s": {"2": 0.0}}
		(work/"darcy.json").write_text(json.dumps(case))
		config = {"schema_version": 1, "backend": "native_tet_p1_darcy_steady",
			"input_route": "volume", "input_file": str(mesh),
			"case_file": "darcy.json", "output_directory": "volume",
			"mpi_ranks": 2}
		def command(name, options=()):
			path = work/f"{name}.json"
			path.write_text(json.dumps(dict(config, output_directory=name)))
			return [sys.executable, str(workflow), str(path), "--solver", str(solver),
				*options]
		invoke(command("volume"))
		status(work/"volume", "passed", 1)
		summary(work/"volume/fields", work/"volume/case.json",
			work/"volume/volume.msh", 2, 4)
		fields(work/"volume/fields", 5, 4)
		invoke(command("volume"), success=False)
		invoke(command("volume", ("--resume",)), success=False)
		if status(work/"volume", "passed", 1).get("resume_runs"):
			raise RuntimeError("steady Darcy workflow attempted restart")
		invoke(command("partial", ("--stop-after-step", "1")), success=False)
		if (work/"partial").exists():
			raise RuntimeError("Darcy planned stop created output")
		wrong = dict(config, backend="dolfin", output_directory="external")
		(work/"external.json").write_text(json.dumps(wrong))
		invoke([sys.executable, str(workflow), str(work/"external.json"),
			"--solver", str(solver)], success=False)
		if (work/"external").exists():
			raise RuntimeError("unsupported Darcy FEM backend created output")
		wrong_solver = dict(config, output_directory="wrong_solver")
		(work/"wrong_solver.json").write_text(json.dumps(wrong_solver))
		invoke([sys.executable, str(workflow), str(work/"wrong_solver.json"),
			"--solver", "/bin/true"], success=False)
		bad = json.loads((work/"wrong_solver/workflow_status.json").read_text())
		if (bad["state"] != "failed" or bad["phase"] != "input_preflight"
				or (work/"wrong_solver/simulation.log").exists()):
			raise RuntimeError("wrong Darcy solver was not rejected before simulation")
		if args.ftetwild:
			surface = work/"pipe.vtp"
			invoke([sys.executable, str(root/"scripts/generate_circular_pipe_surface.py"),
				str(surface), "--length-m", "0.03", "--radius-m", "0.005",
				"--circumferential-segments", "8", "--axial-segments", "2"])
			pipe_case = dict(case, mobility_m2_pa_s=1e-4, source_s_inv=0.0,
				pressure_by_boundary_label_pa={"1": 1.0, "2": 0.0},
				outward_flux_by_boundary_label_m_s={})
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
			manifest = json.loads((work/"surface/mesh_manifest.json").read_text())
			volume = manifest["volume_mesh"]
			data = summary(work/"surface/fields", work/"surface/case.json",
				work/"surface/volume.msh", 2, volume["elements"])
			fields(work/"surface/fields", volume["nodes"], volume["elements"])
			if (data["outward_boundary_flow_m3_s"]["1"] >= 0
					or data["outward_boundary_flow_m3_s"]["2"] <= 0):
				raise RuntimeError("fTetWild Darcy workflow flow direction differs")
			if args.gmsh:
				gmsh = dict(pipe, output_directory="surface_gmsh",
					mesher={"kind": "gmsh", "target_size_m": 0.007})
				(work/"surface_gmsh.json").write_text(json.dumps(gmsh))
				invoke([sys.executable, str(workflow), str(work/"surface_gmsh.json"),
					"--solver", str(solver)])
				status(work/"surface_gmsh", "passed", 1)
				data = json.loads((work/"surface_gmsh/mesh_manifest.json").read_text())[
					"volume_mesh"]
				summary(work/"surface_gmsh/fields", work/"surface_gmsh/case.json",
					work/"surface_gmsh/volume.msh", 2, data["elements"])
				fields(work/"surface_gmsh/fields", data["nodes"], data["elements"])
		print("native_tet_darcy_workflow: PASS volume, status/provenance, "
			"steady restart rejection, backend rejection"
			+ ("; fTetWild surface" if args.ftetwild else "")
			+ ("; Gmsh surface" if args.gmsh else ""))


if __name__ == "__main__":
	main()
