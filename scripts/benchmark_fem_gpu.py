#!/usr/bin/env python3
"""Compare native Tet FEM solvers on CPU and GPU for growing box meshes.

For each box mesh (scripts/generate_tet_box_mesh.py), P1 Darcy and P1 species
run with the CPU PETSc solver and the CUDA solver; an optional second CUDA
build (for example the former dense-LU solvers) runs too while its dense
matrix fits. The table reports the solver's own solve time (PETSc KSPSolve
from -log_view on the CPU), Krylov iterations, host peak RSS, CUDA peak
allocation, and the relative L2 difference of the final field from the CPU.
With --flow-divisions, P2/P1 channel flow runs on the GPU only (it has no CPU
CLI with the same case); its velocity is compared with the reference build.
Needs the vtk Python module to read VTU output.
"""

import argparse
import array
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import vtk

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"scripts"))
from generate_tet_box_mesh import box_mesh  # noqa: E402
from generate_fsi_channel_fixture import write_msh  # noqa: E402


def read_point_field(path, name):
	reader = vtk.vtkXMLUnstructuredGridReader()
	reader.SetFileName(str(path))
	reader.Update()
	data = reader.GetOutput().GetPointData()
	values = data.GetArray(name)
	if values is None:
		raise RuntimeError(f"{path}: no point array {name}")
	ids = data.GetArray("GlobalPointIds")
	order = [int(ids.GetTuple1(i)) for i in range(values.GetNumberOfTuples())] if ids \
		else list(range(values.GetNumberOfTuples()))
	field = [0.0]*len(order)
	for index, point in enumerate(order):
		field[point] = values.GetTuple1(index)
	return field


def relative_l2(reference, value):
	numerator = sum((a-b)**2 for a, b in zip(reference, value))
	denominator = sum(a*a for a in reference)
	return (numerator/denominator)**0.5 if denominator else numerator**0.5


def run(command, log, environment=None):
	started = time.monotonic()
	with log.open("w") as stream:
		result = subprocess.run([str(part) for part in command], stdout=stream, stderr=subprocess.STDOUT,
			env=environment)
	if result.returncode:
		raise RuntimeError(f"{' '.join(map(str, command))} failed; see {log}")
	return time.monotonic()-started, log.read_text()


def petsc_solve_seconds(text):
	"""Total KSPSolve time from a PETSc -log_view report."""
	line = next(line for line in text.splitlines() if line.startswith("KSPSolve "))
	return float(line.split()[3])


def completion(text, prefix):
	line = next(line for line in text.splitlines() if line.startswith(prefix))
	return {key: value for key, value in re.findall(r"(\w+)=(\S+)", line)}


def darcy_case(directory, mesh):
	case = directory/"darcy.json"
	case.write_text(json.dumps({"schema_version": 1, "mesh_file": mesh.name, "mobility_m2_pa_s": 1e-8,
		"source_s_inv": 0.0, "pressure_by_boundary_label_pa": {"1": 1.0, "2": 0.0},
		"outward_flux_by_boundary_label_m_s": {}}))
	return case


def species_case(directory, mesh, steps):
	case = directory/"species.json"
	case.write_text(json.dumps({"schema_version": 1, "mesh_file": mesh.name, "species_id": "tracer",
		"initial_concentration_mol_m3": 1.0, "diffusivity_m2_s": 1e-5, "source_mol_m3_s": 0.0,
		"fluid_velocity_m_s": [1e-3, 0.0, 0.0], "mesh_velocity_m_s": [0.0, 0.0, 0.0],
		"inflow_concentration_by_label_mol_m3": {"1": 2.0}, "monotone": False,
		"time": {"dt_s": 0.5, "steps": steps}}))
	return case


def flow_case(directory, mesh):
	case = directory/"flow.json"
	case.write_text(json.dumps({"schema_version": 1, "mesh_file": mesh.name,
		"fluid": {"density_kg_m3": 1000.0, "dynamic_viscosity_pa_s": 0.004,
			"initial_velocity_m_s": [0.0, 0.0, 0.0], "nonlinear_tolerance": 1e-8, "maximum_iterations": 12},
		"boundaries": {"velocity_by_label_m_s": {"1": [0.01, 0.0, 0.0], "3": [0.0, 0.0, 0.0]},
			"pressure_by_label_pa": {"2": 0.0}},
		"time": {"dt_s": 0.05, "steps": 2}}))
	return case


def binary_field(path):
	values = array.array("d")
	values.frombytes(path.read_bytes())
	return list(values)


def gpu_row(name, binary, case, output, field, prefix, reference):
	wall, text = run([binary, case, output], output.with_suffix(".log"))
	values = completion(text, prefix)
	flow = prefix.startswith("native_tet_flow")
	row = {"solver": name, "wall_s": wall, "assembly_s": float(values["assembly_s"]),
		"solve_s": float(values["solve_s"]), "host_peak_rss_mib": int(values["host_peak_rss_kib"])/1024,
		"cuda_peak_mib": int(values["cuda_peak_bytes"])/2**20}
	if flow:
		row["iterations"] = text.count("_iteration step=")
	for key in ("cg_iterations", "krylov_iterations"):
		if key in values:
			row["iterations"] = int(values[key])
	result = binary_field(field(output)) if flow else read_point_field(*field(output))
	if reference is not None:
		row["relative_l2_vs_reference" if flow else "relative_l2_vs_cpu"] = relative_l2(reference, result)
	row["field"] = result
	return row


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--divisions", type=int, nargs="+", default=[8, 16, 32])
	parser.add_argument("--output-dir", type=Path, required=True)
	parser.add_argument("--species-steps", type=int, default=4)
	parser.add_argument("--reference-gpu-dir", type=Path,
		help="directory with another native_tet_{darcy,species}_cuda build to compare")
	parser.add_argument("--reference-max-dofs", type=int, default=20000,
		help="skip the reference GPU build above this many vertices (dense LU memory)")
	parser.add_argument("--flow-divisions", type=int, nargs="*", default=[],
		help="also run P2/P1 GPU channel flow on these box sizes")
	parser.add_argument("--flow-reference-max-tetrahedra", type=int, default=1500,
		help="skip the reference GPU flow build above this many tetrahedra")
	args = parser.parse_args()
	cpu, gpu = ROOT/"solvers/cpu", ROOT/"solvers/cuda"
	rows = []
	for divisions in args.divisions:
		directory = args.output_dir.resolve()/f"box{divisions}"
		directory.mkdir(parents=True, exist_ok=False)
		points, cells, triangles = box_mesh(divisions, 0.01)
		mesh = directory/"box.msh"
		write_msh(mesh, points, cells, triangles)
		size = {"divisions": divisions, "vertices": len(points), "tetrahedra": len(cells)}
		darcy = darcy_case(directory, mesh)
		logged = dict(os.environ, PETSC_OPTIONS="-log_view")
		wall, text = run([cpu/"native_tet_darcy", darcy, "--output-dir", directory/"darcy-cpu"],
			directory/"darcy-cpu.log", logged)
		summary = json.loads((directory/"darcy-cpu/run_summary.json").read_text())
		reference = read_point_field(directory/"darcy-cpu/fields/rank0.vtu", "pressure_pa")
		rows.append(dict(size, case="darcy", solver="cpu", wall_s=wall, solve_s=petsc_solve_seconds(text),
			iterations=summary["linear_iterations"]))
		builds = [("gpu", gpu)]
		if args.reference_gpu_dir and len(points) <= args.reference_max_dofs:
			builds.append(("reference-gpu", args.reference_gpu_dir))
		for name, build in builds:
			rows.append(dict(size, case="darcy", **gpu_row(name, build/"native_tet_darcy_cuda", darcy,
				directory/f"darcy-{name}", lambda out: (out/"tissue.vtu", "pressure_pa"),
				"native_tet_darcy_cuda_complete", reference)))
		species = species_case(directory, mesh, args.species_steps)
		wall, text = run([cpu/"native_tet_species_transport", species, "--output-dir", directory/"species-cpu",
			"--visualization-format", "pvtu"], directory/"species-cpu.log", logged)
		final = f"step_{args.species_steps}/rank0.vtu"
		reference = read_point_field(directory/"species-cpu"/final, "concentration_mol_m3")
		iterations = sum(int(value) for value in re.findall(r"linear_iterations=(\d+)", text))
		rows.append(dict(size, case="species", solver="cpu", wall_s=wall, solve_s=petsc_solve_seconds(text),
			iterations=iterations))
		for name, build in builds:
			rows.append(dict(size, case="species", **gpu_row(name, build/"native_tet_species_cuda", species,
				directory/f"species-{name}",
				lambda out: (out/f"species_step_{args.species_steps}.vtu", "concentration_mol_m3"),
				"native_tet_species_cuda_complete", reference)))
	for divisions in args.flow_divisions:
		directory = args.output_dir.resolve()/f"flow{divisions}"
		directory.mkdir(parents=True, exist_ok=False)
		points, cells, triangles = box_mesh(divisions, 0.01)
		mesh = directory/"box.msh"
		write_msh(mesh, points, cells, triangles)
		size = {"divisions": divisions, "vertices": len(points), "tetrahedra": len(cells)}
		flow = flow_case(directory, mesh)
		field = lambda out: out/"flow_velocity_step_2.bin"
		reference = None
		if args.reference_gpu_dir and len(cells) <= args.flow_reference_max_tetrahedra:
			row = gpu_row("reference-gpu", args.reference_gpu_dir/"native_tet_flow_cuda", flow,
				directory/"flow-reference", field, "native_tet_flow_cuda_complete", None)
			reference = row.pop("field")
			rows.append(dict(size, case="flow", **row))
		row = gpu_row("gpu", gpu/"native_tet_flow_cuda", flow, directory/"flow-gpu", field,
			"native_tet_flow_cuda_complete", reference)
		row.pop("field")
		rows.append(dict(size, case="flow", **row))
	for row in rows:
		row.pop("field", None)
	(args.output_dir/"benchmark.json").write_text(json.dumps(rows, indent=2)+"\n")
	columns = ("case", "vertices", "tetrahedra", "solver", "wall_s", "assembly_s", "solve_s", "iterations",
		"host_peak_rss_mib", "cuda_peak_mib", "relative_l2_vs_cpu", "relative_l2_vs_reference")
	print(" | ".join(columns))
	for row in rows:
		print(" | ".join(f"{row[c]:.3g}" if isinstance(row.get(c), float) else str(row.get(c, "")) for c in columns))


if __name__ == "__main__":
	main()
