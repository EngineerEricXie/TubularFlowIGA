#!/usr/bin/env python3
"""Explicit surface/volume to project-owned tetra hydraulic, species or Darcy FEM."""

import argparse
import fcntl
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time


def require_keys(value, allowed, required, context):
	if not isinstance(value, dict):
		raise ValueError(f"{context} must be an object")
	unknown = set(value)-set(allowed)
	missing = set(required)-set(value)
	if unknown or missing:
		raise ValueError(f"{context} unknown={sorted(unknown)} missing={sorted(missing)}")


def positive_number(value, name):
	if isinstance(value, bool) or not isinstance(value, (int, float)) \
			or not math.isfinite(value) or value <= 0:
		raise ValueError(f"{name} must be finite and positive")
	return value


def existing_file(base, value, name):
	if not isinstance(value, str) or not value:
		raise ValueError(f"{name} must be a nonempty path")
	path = (base/value).resolve(strict=True)
	if not path.is_file():
		raise ValueError(f"{name} is not a file: {path}")
	return path


def sha256(path):
	return hashlib.sha256(path.read_bytes()).hexdigest()


def invoke(command, logfile):
	start = time.monotonic()
	with logfile.open("w", encoding="utf-8") as output:
		result = subprocess.run(command, stdout=output, stderr=subprocess.STDOUT,
			text=True, check=False)
	seconds = time.monotonic()-start
	if result.returncode:
		raise RuntimeError(f"command exited {result.returncode}; see {logfile}")
	return seconds


def accepted_steps(logfile, backend):
	if backend == "native_tet_p1_darcy_steady":
		if re.search(r"^native_darcy: PASS\b", logfile.read_text(encoding="utf-8"),
				flags=re.MULTILINE) is None:
			raise RuntimeError(f"native Darcy solver did not report completion: {logfile}")
		return 1
	pattern = (r"^native_hydraulic_complete accepted_steps=(\d+)\b"
		if backend == "native_tet_p2p1_ale_hydraulic"
		else r"^native_species_step step=(\d+)\b")
	steps = re.findall(pattern,
		logfile.read_text(encoding="utf-8"), flags=re.MULTILINE)
	if not steps or (backend == "native_tet_p2p1_ale_hydraulic" and len(steps) != 1):
		raise RuntimeError(f"native solver did not report accepted steps: {logfile}")
	return int(steps[-1])


def published_checkpoint_step(checkpoints):
	if not checkpoints.is_dir():
		raise RuntimeError("native workflow has no checkpoint directory")
	steps = []
	for manifest in checkpoints.glob("*/manifest"):
		matches = re.findall(r"^accepted_steps (\d+)$",
			manifest.read_text(encoding="utf-8"), flags=re.MULTILINE)
		if len(matches) != 1:
			raise RuntimeError(f"checkpoint manifest step is invalid: {manifest}")
		steps.append(int(matches[0]))
	if not steps:
		raise RuntimeError("native workflow has no published checkpoint")
	return max(steps)


def published_species_step(fields):
	state_path = fields/"run_state.json"
	if not state_path.is_file():
		raise RuntimeError("native species workflow has no published state")
	state = json.loads(state_path.read_text(encoding="utf-8"))
	step = state.get("accepted_steps")
	if type(step) is not int or step < 1:
		raise RuntimeError("native species workflow state has no accepted step")
	return step


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("configuration", type=Path)
	parser.add_argument("--solver", type=Path,
		help="native executable matching the explicitly selected backend")
	parser.add_argument("--stop-after-step", type=int,
		help="publish a resumable planned partial run after this accepted step")
	parser.add_argument("--resume", action="store_true",
		help="continue a planned partial run in its existing output directory")
	args = parser.parse_args()
	config_path = args.configuration.resolve(strict=True)
	base = config_path.parent
	config = json.loads(config_path.read_text(encoding="utf-8"))
	require_keys(config,
		{"schema_version", "backend", "input_route", "input_file", "case_file",
		"output_directory", "mpi_ranks", "mesher"},
		{"schema_version", "backend", "input_route", "input_file", "case_file",
		"output_directory", "mpi_ranks"}, "workflow")
	if type(config["schema_version"]) is not int or config["schema_version"] != 1:
		raise ValueError("workflow schema_version must be 1")
	backends = {
		"native_tet_p2p1_ale_hydraulic": "native_tet_hydraulic_graph",
		"native_tet_p1_prescribed_species": "native_tet_species_transport",
		"native_tet_p1_darcy_steady": "native_tet_darcy"}
	if config["backend"] not in backends:
		raise ValueError("only project-owned native tetra hydraulic/species/Darcy backends are supported")
	species = config["backend"] == "native_tet_p1_prescribed_species"
	darcy = config["backend"] == "native_tet_p1_darcy_steady"
	if darcy and (args.stop_after_step is not None or args.resume):
		raise ValueError("steady Darcy workflow has no step horizon or restart")
	if config["input_route"] not in ("surface", "volume"):
		raise ValueError("only surface and volume tetra input routes are implemented")
	if type(config["mpi_ranks"]) is not int or not 1 <= config["mpi_ranks"] <= 16:
		raise ValueError("mpi_ranks must be an integer between 1 and 16")
	input_path = existing_file(base, config["input_file"], "input_file")
	case_template = existing_file(base, config["case_file"], "case_file")
	solver = (args.solver or Path(__file__).resolve().parents[1]/
		"solvers/cpu"/backends[config["backend"]]).resolve(strict=True)
	if not solver.is_file():
		raise ValueError("native solver is not a file")
	case = json.loads(case_template.read_text(encoding="utf-8"))
	if not isinstance(case, dict) or "mesh_file" not in case:
		raise ValueError("native case template must contain mesh_file")
	coupled_species = (config["backend"] == "native_tet_p2p1_ale_hydraulic"
		and bool(case.get("species")))
	if coupled_species and args.resume:
		raise ValueError("coupled flow/species workflow has no paired checkpoint")
	if darcy:
		if "time" in case:
			raise ValueError("steady Darcy case must not define a time horizon")
		final_steps = 1
	else:
		time_config = case.get("time")
		if not isinstance(time_config, dict) or type(time_config.get("steps")) is not int \
				or time_config["steps"] < 1:
			raise ValueError("native case template must contain a positive step count")
		final_steps = time_config["steps"]
	if args.stop_after_step is not None and not 1 <= args.stop_after_step <= final_steps:
		raise ValueError("stop-after-step must be within the case step horizon")
	if coupled_species and args.stop_after_step not in (None, final_steps):
		raise ValueError("coupled flow/species workflow cannot publish a resumable partial run")
	mesher = config.get("mesher")
	if config["input_route"] == "surface":
		require_keys(mesher, {"kind", "target_size_m", "envelope_m",
			"executable"}, {"kind", "target_size_m"}, "mesher")
		if mesher["kind"] not in ("ftetwild", "gmsh"):
			raise ValueError("surface mesher must be ftetwild or gmsh")
		positive_number(mesher["target_size_m"], "target_size_m")
		if mesher["kind"] == "ftetwild":
			if "envelope_m" not in mesher or "executable" not in mesher:
				raise ValueError("ftetwild requires envelope_m and executable")
			positive_number(mesher["envelope_m"], "envelope_m")
			mesher_binary = existing_file(base, mesher["executable"], "mesher.executable")
		elif "envelope_m" in mesher or "executable" in mesher:
			raise ValueError("gmsh does not accept fTetWild envelope/executable")
	elif mesher is not None:
		raise ValueError("volume route must not specify a mesher")
	if not isinstance(config["output_directory"], str) or not config["output_directory"]:
		raise ValueError("output_directory must be a nonempty path")
	output = (base/config["output_directory"]).resolve()
	if output.exists() and not args.resume:
		raise FileExistsError(f"workflow output already exists: {output}")
	status_path = output/"workflow_status.json"
	identity = {"schema_version": 1,
		"acceptance": "functional_only_not_physical_validation",
		"backend": config["backend"], "input_route": config["input_route"],
		"input_file": str(input_path), "case_template_file": str(case_template),
		"solver_file": str(solver),
		"input_sha256": sha256(input_path), "case_template_sha256": sha256(case_template),
		"configuration_sha256": sha256(config_path), "solver_sha256": sha256(solver),
		"mpi_ranks": config["mpi_ranks"]}
	if config["input_route"] == "surface" and mesher["kind"] == "ftetwild":
		identity["mesher_binary_sha256"] = sha256(mesher_binary)
	def publish():
		status["evidence_classification"] = {
			"numerical_acceptance": "functional_checks_only"
				if status["state"] == "passed" else "not_established",
			"physical_acceptance": "not_established",
			"demonstration_visualization": "available"
				if status.get("accepted_steps", 0) > 0 else "not_available"}
		temporary = status_path.with_name("workflow_status.json.tmp")
		temporary.write_text(json.dumps(status, indent=2, sort_keys=True)+"\n",
			encoding="utf-8")
		temporary.replace(status_path)
	if args.resume:
		if not output.is_dir() or not status_path.is_file():
			raise RuntimeError("resume requires an existing native workflow output")
		lock_path = output/"workflow.lock"
		if not lock_path.is_file():
			raise RuntimeError("resume requires the original workflow lock file")
		lock_file = lock_path.open("r+")
		fcntl.flock(lock_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
		status = json.loads(status_path.read_text(encoding="utf-8"))
		if status.get("state") not in ("partial", "failed", "running") \
				or status.get("phase") not in ("simulation", "simulation_resume") \
				or any(status.get(key) != value for key, value in identity.items()):
			raise RuntimeError("resume state or workflow identity differs")
		mesh = output/"volume.msh"
		case_path = output/"case.json"
		if (not mesh.is_file() or not case_path.is_file()
				or sha256(mesh) != status.get("mesh_sha256")
				or sha256(case_path) != status.get("effective_case_sha256")):
			raise RuntimeError("resume mesh or effective case differs")
		last = (published_species_step(output/"fields") if species
			else published_checkpoint_step(output/"checkpoints"))
		if last < 1 or last > final_steps \
				or (status.get("accepted_steps") is not None
				and last < status["accepted_steps"]):
			raise RuntimeError("resume checkpoint is inconsistent with workflow status")
		if args.stop_after_step is not None and args.stop_after_step <= last:
			raise ValueError("resume stop step must advance past the checkpoint")
		for step in range(1, last+1):
			if not (output/"fields"/f"step_{step}"/"snapshot.pvtu").is_file():
				raise RuntimeError("resume is missing an accepted visualization step")
		for entry in (output/"fields").glob("step_*"):
			if entry.name[5:].isdigit() and int(entry.name[5:]) > last:
				raise RuntimeError("resume found an uncheckpointed visualization step")
		status["state"] = "running"
		status["phase"] = "simulation_resume"
		status["accepted_steps"] = last
		status.pop("error", None)
		resume_run = status.get("resume_runs", 0)+1
		status["resume_runs"] = resume_run
		publish()
		try:
			logfile = output/f"simulation_resume_{resume_run}.log"
			command = ["mpiexec", "-np", str(config["mpi_ranks"]), str(solver),
				str(case_path)]
			if species:
				command += ["--resume", str(output/"fields")]
			else:
				command += ["--checkpoint-dir", str(output/"checkpoints"),
					"--restart-dir", str(output/"checkpoints"),
					"--output-dir", str(output/"fields")]
			if args.stop_after_step is not None:
				command += ["--stop-after-step", str(args.stop_after_step)]
			status["timings_s"][f"simulation_resume_{resume_run}"] = invoke(command,
				logfile)
			status["accepted_steps"] = accepted_steps(logfile, config["backend"])
			target = args.stop_after_step or final_steps
			if status["accepted_steps"] != target:
				raise RuntimeError("resumed solver stopped at the wrong step")
			status["state"] = "passed" if status["accepted_steps"] == final_steps else "partial"
			publish()
			print(f"native_tet_workflow: {status['state'].upper()} "
				f"accepted_steps={status['accepted_steps']} output={output}")
		except Exception as error:
			status["state"] = "failed"
			status["error"] = str(error)
			publish()
			raise
		return
	output.mkdir(parents=True)
	lock_file = (output/"workflow.lock").open("x+")
	fcntl.flock(lock_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
	status = dict(identity, state="running", phase="mesh_preparation", timings_s={})
	publish()
	try:
		mesh = output/"volume.msh"
		if config["input_route"] == "volume":
			start = time.monotonic()
			shutil.copyfile(input_path, mesh)
			status["timings_s"]["volume_copy"] = time.monotonic()-start
		else:
			script = (Path(__file__).resolve().parent/
				("ftetwild_to_fem_volume.py" if mesher["kind"] == "ftetwild"
				else "surface_to_fem_volume.py"))
			command = [sys.executable, str(script), str(input_path), str(mesh),
				"--manifest", str(output/"mesh_manifest.json"),
				"--target-size-m", str(mesher["target_size_m"])]
			if mesher["kind"] == "ftetwild":
				command += ["--ftetwild", str(mesher_binary),
					"--envelope-m", str(mesher["envelope_m"])]
			status["mesher"] = mesher["kind"]
			status["timings_s"]["volume_meshing"] = invoke(command,
				output/"meshing.log")
		case["mesh_file"] = mesh.name
		case_path = output/"case.json"
		case_path.write_text(json.dumps(case, indent=2)+"\n",
			encoding="utf-8")
		status["mesh_sha256"] = sha256(mesh)
		status["effective_case_sha256"] = sha256(case_path)
		status["phase"] = "input_preflight"
		publish()
		status["timings_s"]["input_preflight"] = invoke(
			["mpiexec", "-np", str(config["mpi_ranks"]), str(solver),
				str(case_path), "--check-input"], output/"input_preflight.log")
		expected_preflight = ("native_darcy_input: PASS" if darcy else
			"native_species_input: PASS" if species else "native_hydraulic_input: PASS")
		if expected_preflight not in (output/"input_preflight.log").read_text(
				encoding="utf-8"):
			raise RuntimeError("native solver preflight did not match selected backend")
		status["phase"] = "simulation"
		publish()
		command = ["mpiexec", "-np", str(config["mpi_ranks"]), str(solver),
			str(case_path)]
		if not species and not darcy and not coupled_species:
			command += ["--checkpoint-dir", str(output/"checkpoints")]
		command += ["--output-dir", str(output/"fields")]
		if args.stop_after_step is not None:
			command += ["--stop-after-step", str(args.stop_after_step)]
		logfile = output/"simulation.log"
		status["timings_s"]["simulation"] = invoke(command, logfile)
		status["accepted_steps"] = accepted_steps(logfile, config["backend"])
		target = args.stop_after_step or final_steps
		if status["accepted_steps"] != target:
			raise RuntimeError("native solver stopped at the wrong step")
		status["state"] = "passed" if status["accepted_steps"] == final_steps else "partial"
		publish()
		print(f"native_tet_workflow: {status['state'].upper()} route={config['input_route']} "
			f"mesher={status.get('mesher', 'none')} ranks={config['mpi_ranks']} "
			f"output={output}")
	except Exception as error:
		status["state"] = "failed"
		status["error"] = str(error)
		publish()
		raise


if __name__ == "__main__":
	main()
