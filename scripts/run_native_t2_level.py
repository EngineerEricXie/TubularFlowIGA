#!/usr/bin/env python3
"""Run one frozen T2 mesh level through the repository-owned native FEM solver."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

from t2_contract import exact_contract_sha256, mesh_contract_sha256


SOURCE_FILES = (
	"solvers/cpu/include/NativeTetFem.hpp",
	"solvers/cpu/src/native_tet_flow.cpp",
)


def sha256(path):
	return hashlib.sha256(path.read_bytes()).hexdigest()


def combined_sha256(paths):
	digest = hashlib.sha256()
	for label, path in paths:
		contents = path.read_bytes()
		digest.update(label.encode("utf-8")+b"\0")
		digest.update(str(len(contents)).encode("ascii")+b"\0")
		digest.update(contents)
	return digest.hexdigest()


def git(root, *arguments):
	completed = subprocess.run(["git", *arguments], cwd=root, text=True,
		capture_output=True, check=True)
	return completed.stdout.strip()


def solver_options(mode, profile):
	strict = mode == "validation"
	options = [
		"-native_tet_validation", "true" if strict else "false",
		"-native_tet_ksp_rtol", "1e-10" if strict else "1e-8",
		"-native_tet_ksp_atol", "1e-12" if strict else "1e-10",
		"-native_tet_ksp_max_it", "600",
	]
	if profile == "direct":
		return options+["-native_tet_ksp_type", "preonly",
			"-native_tet_pc_type", "lu",
			"-native_tet_pc_factor_mat_solver_type", "mumps"]
	return options+[
		"-native_tet_ksp_type", "fgmres",
		"-native_tet_ksp_gmres_restart", "300",
		"-native_tet_pc_type", "fieldsplit",
		"-native_tet_pc_fieldsplit_type", "schur",
		"-native_tet_pc_fieldsplit_schur_fact_type", "full",
		"-native_tet_pc_fieldsplit_schur_precondition", "user",
		"-native_tet_fieldsplit_velocity_ksp_type", "preonly",
		"-native_tet_fieldsplit_velocity_pc_type", "lu",
		"-native_tet_fieldsplit_velocity_pc_factor_mat_solver_type", "mumps",
		"-native_tet_fieldsplit_pressure_ksp_type", "preonly",
		"-native_tet_fieldsplit_pressure_pc_type", "lu",
	]


def main():
	root = Path(__file__).resolve().parents[1]
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("mesh_series")
	parser.add_argument("level", type=int, choices=(1, 2, 3))
	parser.add_argument("output")
	parser.add_argument("--contract", default=str(root/"benchmarks"/"t2_fixed_flow_contract.json"))
	parser.add_argument("--executable", default=str(root/"solvers"/"cpu"/"native_tet_flow"))
	parser.add_argument("--mode", choices=("functional", "validation"), default="functional")
	parser.add_argument("--solver-profile", choices=("fieldsplit", "direct"), default="fieldsplit")
	parser.add_argument("--mpi-ranks", type=int, default=1,
		help="launch this many MPI ranks (default: 1)")
	parser.add_argument("--mpiexec", default="mpiexec",
		help="MPI launcher used when --mpi-ranks is greater than one")
	parser.add_argument("--mpi-launcher-option", action="append", default=[],
		help="option placed before -np for the MPI launcher; repeat as needed")
	parser.add_argument("--timeout-s", type=float)
	parser.add_argument("--petsc-option", action="append", default=[],
		help="append one PETSc option token; repeat for option values")
	args = parser.parse_args()

	series_path = Path(args.mesh_series).resolve()
	series = json.loads(series_path.read_text(encoding="utf-8"))
	contract_path = Path(args.contract).resolve()
	contract_bytes = contract_path.read_bytes()
	contract = json.loads(contract_bytes)
	if series.get("kind") != "t2_mesh_preparation_not_solver_evidence":
		raise RuntimeError("mesh-series kind is not the frozen T2 preparation artifact")
	if series.get("case_id") != contract["case_id"]:
		raise RuntimeError("mesh-series case id differs from the contract")
	if series.get("mesh_contract_sha256") != mesh_contract_sha256(contract):
		raise RuntimeError("mesh-series geometry contract hash is stale")
	levels = series.get("levels")
	if not isinstance(levels, list) or len(levels) != 3:
		raise RuntimeError("mesh-series must contain exactly three levels")
	level = levels[args.level-1]
	targets = contract["discretization_contract"]["mesh_levels_target_size_m"]
	if level.get("target_size_m") != targets[args.level-1]:
		raise RuntimeError("selected mesh target differs from the frozen contract")
	mesh = (series_path.parent/level["volume_mesh"]).resolve()
	manifest = (series_path.parent/level["volume_manifest"]).resolve()
	if sha256(mesh) != level.get("mesh_sha256"):
		raise RuntimeError("selected mesh SHA-256 differs from mesh-series")
	volume_manifest = json.loads(manifest.read_text(encoding="utf-8"))
	if volume_manifest["volume_mesh"]["sha256"] != level["mesh_sha256"]:
		raise RuntimeError("volume manifest mesh SHA-256 differs from mesh-series")

	executable = Path(args.executable).resolve()
	if not executable.is_file():
		raise RuntimeError(f"native FEM executable does not exist: {executable}")
	if args.mpi_ranks < 1:
		raise RuntimeError("--mpi-ranks must be positive")
	timeout = args.timeout_s
	if timeout is None:
		timeout = float(contract["local_resource_budget"]["maximum_wall_s"])
	with tempfile.TemporaryDirectory(prefix="native-t2-level-") as directory:
		raw = Path(directory)/"raw-result.json"
		command = [str(executable), str(mesh), str(contract_path),
			*solver_options(args.mode, args.solver_profile),
			"-native_tet_level_result", str(raw), *args.petsc_option]
		if args.mpi_ranks > 1:
			launcher = shutil.which(args.mpiexec)
			if launcher is None:
				raise RuntimeError(f"MPI launcher does not exist: {args.mpiexec}")
			command = [launcher, *args.mpi_launcher_option, "-np", str(args.mpi_ranks), *command]
		try:
			completed = subprocess.run(command, text=True, capture_output=True, timeout=timeout)
		except subprocess.TimeoutExpired as error:
			raise RuntimeError(
				f"native FEM level solve exceeded {timeout:g} s; no result was promoted") from error
		if completed.returncode != 0:
			detail = (completed.stderr or completed.stdout).strip()
			raise RuntimeError(f"native FEM level solve failed ({completed.returncode}): {detail}")
		if not raw.is_file():
			raise RuntimeError("native FEM solver did not write its raw result")
		result = json.loads(raw.read_text(encoding="utf-8"))

	expected_class = "physical_validation_candidate" if args.mode == "validation" \
		else "functional_smoke"
	if result.get("result_classification") != expected_class:
		raise RuntimeError("native FEM result classification differs from requested mode")
	if result.get("validation_gates_enforced") != (args.mode == "validation"):
		raise RuntimeError("native FEM validation-gate flag differs from requested mode")
	if result.get("mpi_ranks") != args.mpi_ranks:
		raise RuntimeError("native FEM reported an MPI rank count different from the launcher request")
	source_paths = [(name, root/name) for name in SOURCE_FILES]
	for name, path in source_paths:
		if not path.is_file():
			raise RuntimeError(f"missing native FEM source file: {name}")
	commit = git(root, "rev-parse", "HEAD")
	dirty = bool(git(root, "status", "--porcelain", "--", *SOURCE_FILES))
	result.update({
		"target_size_m": level["target_size_m"],
		"source_commit": commit,
		"source_tree_state": "dirty_with_source_manifest" if dirty else "clean",
		"source_files_sha256": combined_sha256(source_paths),
		"executable_sha256": sha256(executable),
		"input_sha256": combined_sha256((("contract", contract_path),
			("mesh-series", series_path), ("mesh", mesh), ("volume-manifest", manifest))),
		"mesh_sha256": level["mesh_sha256"],
		"geometry_error": level["geometry_error"],
		"contract_sha256": exact_contract_sha256(contract_bytes),
		"mesh_contract_sha256": series["mesh_contract_sha256"],
		"mesh_level": args.level,
		"solver_profile": args.solver_profile,
		"requested_mpi_ranks": args.mpi_ranks,
		"mpi_launcher_options": args.mpi_launcher_option,
	})
	output = Path(args.output).resolve()
	output.parent.mkdir(parents=True, exist_ok=True)
	output.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n", encoding="utf-8")
	print(f"run_native_t2_level: PASS mode={args.mode} level={args.level} "
		f"classification={expected_class} output={output}")
	return 0


if __name__ == "__main__":
	try:
		raise SystemExit(main())
	except (KeyError, OSError, RuntimeError, subprocess.SubprocessError, ValueError) as error:
		print(f"run_native_t2_level: ERROR: {error}")
		raise SystemExit(2)
