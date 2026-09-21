#!/usr/bin/env python3
"""Exercise the project-owned labelled-tetra Darcy CLI and rank-owned VTU."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET


def invoke(command, success=True):
	result = subprocess.run(command, capture_output=True, text=True, timeout=180)
	if (result.returncode == 0) != success:
		raise RuntimeError(f"unexpected Darcy result: {command}\n"
			f"{result.stdout}\n{result.stderr}")
	return result.stdout + result.stderr


def fields(directory, vertices, tetrahedra):
	index = ET.parse(directory/"fields/snapshot.pvtu").getroot()
	point_fields = {field.attrib.get("Name") for field in
		index.findall(".//PPointData/PDataArray")}
	cell_fields = {field.attrib.get("Name") for field in
		index.findall(".//PCellData/PDataArray")}
	if "pressure_pa" not in point_fields or not {
			"darcy_flux_m_s", "darcy_rt0_centroid_flux_m_s",
			"conservative_outward_face_flow_m3_s",
			"mobility_m2_pa_s", "source_s_inv"}.issubset(cell_fields):
		raise RuntimeError("Darcy PVTU schema is incomplete")
	pressures, fluxes, mobilities, sources, face_flows, rt0_fluxes = (
		{}, {}, {}, {}, {}, {})
	for entry in index.findall(".//Piece"):
		piece = ET.parse(directory/"fields"/entry.attrib["Source"]).getroot().find(
			".//Piece")
		if piece is None:
			raise RuntimeError("Darcy VTU piece is missing")
		def values(path, kind=float):
			array = piece.find(path)
			if array is None:
				raise RuntimeError(f"Darcy VTU field is missing: {path}")
			return [kind(value) for value in (array.text or "").split()]
		ids = values("PointData/DataArray[@Name='GlobalPointIds']", int)
		concentration = values("PointData/DataArray[@Name='pressure_pa']")
		cell_ids = values("CellData/DataArray[@Name='GlobalCellIds']", int)
		flux = values("CellData/DataArray[@Name='darcy_flux_m_s']")
		rt0_flux = values("CellData/DataArray[@Name='darcy_rt0_centroid_flux_m_s']")
		face_flow = values("CellData/DataArray[@Name='conservative_outward_face_flow_m3_s']")
		mobility = values("CellData/DataArray[@Name='mobility_m2_pa_s']")
		source = values("CellData/DataArray[@Name='source_s_inv']")
		cell_types = values("Cells/DataArray[@Name='types']", int)
		if (len(ids) != len(concentration) or len(flux) != 3*len(cell_ids)
				or len(rt0_flux) != 3*len(cell_ids)
				or len(face_flow) != 4*len(cell_ids)
				or len(mobility) != len(cell_ids) or len(source) != len(cell_ids)
				or any(kind != 10 for kind in cell_types)):
			raise RuntimeError("Darcy VTU array sizes or cell types differ")
		for identifier, pressure in zip(ids, concentration):
			if not math.isfinite(pressure):
				raise RuntimeError("Darcy pressure is nonfinite")
			if identifier in pressures and abs(pressures[identifier]-pressure) > 1e-10:
				raise RuntimeError("Darcy shared pressure differs across pieces")
			pressures[identifier] = pressure
		for local, identifier in enumerate(cell_ids):
			if identifier in fluxes:
				raise RuntimeError("Darcy cell has duplicate ownership")
			value = tuple(flux[3*local:3*local+3])
			if not all(math.isfinite(number) for number in value):
				raise RuntimeError("Darcy cell flux is nonfinite")
			fluxes[identifier] = value
			rt0_fluxes[identifier] = tuple(rt0_flux[3*local:3*local+3])
			if not all(math.isfinite(number) for number in rt0_fluxes[identifier]):
				raise RuntimeError("Darcy RT0 centroid flux is nonfinite")
			face_flows[identifier] = tuple(face_flow[4*local:4*local+4])
			if not all(math.isfinite(number) for number in face_flows[identifier]):
				raise RuntimeError("Darcy conservative face flux is nonfinite")
			mobilities[identifier] = mobility[local]
			sources[identifier] = source[local]
	if len(pressures) != vertices or len(fluxes) != tetrahedra:
		raise RuntimeError("Darcy VTU lost vertices or tetrahedra")
	return pressures, fluxes, mobilities, sources, face_flows, rt0_fluxes


def summary(directory, case, mesh, ranks, tetrahedra):
	result = json.loads((directory/"run_summary.json").read_text())
	if (result.get("schema_version") != 1
			or result.get("kind") != "native_tet_p1_darcy_steady"
			or result.get("case_sha256") != hashlib.sha256(case.read_bytes()).hexdigest()
			or result.get("mesh_sha256") != hashlib.sha256(mesh.read_bytes()).hexdigest()
			or result.get("mpi_ranks") != ranks
			or result.get("tetrahedra") != tetrahedra
			or result.get("boundary_flow_kind") != "P1_cell_gradient_diagnostic"
			or result.get("conservative_face_flow_kind")
				!= "least_squares_cell_balance_recovery"
			or result.get("recovered_vector_field_kind")
				!= "rt0_hdiv_from_conservative_faces"
			or result.get("maximum_cell_balance_defect_m3_s", math.inf) > 1e-10
			or not isinstance(result.get("conservative_outward_boundary_flow_m3_s"), dict)
			or not math.isfinite(result.get("linear_residual_norm", math.nan))):
		raise RuntimeError("Darcy summary or provenance differs")
	return result


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--binary", type=Path, required=True)
	parser.add_argument("--ftetwild", type=Path)
	args = parser.parse_args()
	binary = args.binary.resolve(strict=True)
	root = Path(__file__).resolve().parents[1]
	mesh = (root/"solvers/cpu/tests/data/native_tet_hydraulic_star.msh").resolve()
	with tempfile.TemporaryDirectory(prefix="native-darcy-cli-") as temporary:
		work = Path(temporary)
		base = {"schema_version": 1, "mesh_file": str(mesh),
			"mobility_m2_pa_s": 2.0, "source_s_inv": 1.0,
			"mobility_by_cell_id_m2_pa_s": {"5": 4.0},
			"source_by_cell_id_s_inv": {"6": 2.0},
			"pressure_by_boundary_label_pa": {"1": 1.0},
			"outward_flux_by_boundary_label_m_s": {"2": 0.0}}
		def case_file(name, changes):
			path = work/f"{name}.json"
			path.write_text(json.dumps({**base, **changes}))
			return path
		fields_by_rank = {}
		for ranks in (1, 2, 4):
			case = case_file(f"case_{ranks}", {})
			if "native_darcy_input: PASS" not in invoke(
					["mpiexec", "-np", str(ranks), str(binary), str(case), "--check-input"]):
				raise RuntimeError("Darcy CLI input preflight did not report success")
			output = work/f"run_{ranks}"
			invoke(["mpiexec", "-np", str(ranks), str(binary), str(case),
				"--output-dir", str(output)])
			data = summary(output, case, mesh, ranks, 4)
			if not 0 < data["volume_source_m3_s"] < 1:
				raise RuntimeError("Darcy CLI volume source budget is invalid")
			fields_by_rank[ranks] = fields(output, 5, 4)
			if (fields_by_rank[ranks][2] != {5: 4.0, 6: 2.0, 7: 2.0, 8: 2.0}
					or fields_by_rank[ranks][3] != {5: 1.0, 6: 2.0, 7: 1.0, 8: 1.0}):
				raise RuntimeError("Darcy cell material/source override was not published")
			invoke(["mpiexec", "-np", str(ranks), str(binary), str(case),
				"--output-dir", str(output)], success=False)
		for ranks in (1, 4):
			for expected, actual in zip(fields_by_rank[2], fields_by_rank[ranks]):
				if expected.keys() != actual.keys():
					raise RuntimeError("Darcy MPI ids differ")
				for identifier in expected:
					left, right = expected[identifier], actual[identifier]
					if isinstance(left, tuple):
						if max(abs(a-b) for a, b in zip(left, right)) > 1e-10:
							raise RuntimeError("Darcy MPI flux differs")
					elif abs(left-right) > 1e-10:
						raise RuntimeError("Darcy MPI pressure differs")
		for name, changes in (
			("unknown", {"organ": "liver"}),
			("bad_label", {"pressure_by_boundary_label_pa": {"99": 0.0}}),
			("conflicting", {"outward_flux_by_boundary_label_m_s": {"1": 0.0}}),
			("negative_mobility", {"mobility_m2_pa_s": -1.0}),
			("missing_cell", {"mobility_by_cell_id_m2_pa_s": {"99": 1.0}})):
			bad = case_file(name, changes)
			invoke(["mpiexec", "-np", "2", str(binary), str(bad), "--check-input"],
				success=False)
		if args.ftetwild:
			surface, generated, manifest = work/"pipe.vtp", work/"pipe.msh", work/"mesh.json"
			invoke([sys.executable, str(root/"scripts/generate_circular_pipe_surface.py"),
				str(surface), "--length-m", "0.03", "--radius-m", "0.005",
				"--circumferential-segments", "8", "--axial-segments", "2"])
			invoke([sys.executable, str(root/"scripts/ftetwild_to_fem_volume.py"),
				str(surface), str(generated), "--manifest", str(manifest),
				"--ftetwild", str(args.ftetwild.resolve(strict=True)),
				"--target-size-m", "0.007", "--envelope-m", "0.0001",
				"--stop-energy", "12", "--max-optimization-passes", "40",
				"--max-threads", "2"])
			volume = json.loads(manifest.read_text())["volume_mesh"]
			pipe = case_file("pipe", {"mesh_file": str(generated),
				"mobility_m2_pa_s": 1e-4, "source_s_inv": 0.0,
				"mobility_by_cell_id_m2_pa_s": {},
				"source_by_cell_id_s_inv": {},
				"pressure_by_boundary_label_pa": {"1": 1.0, "2": 0.0},
				"outward_flux_by_boundary_label_m_s": {}})
			pipe_output = work/"pipe_output"
			invoke(["mpiexec", "-np", "2", str(binary), str(pipe),
				"--output-dir", str(pipe_output)])
			data = summary(pipe_output, pipe, generated, 2, volume["elements"])
			if (data["outward_boundary_flow_m3_s"]["1"] >= 0
					or data["outward_boundary_flow_m3_s"]["2"] <= 0):
				raise RuntimeError("Darcy fTetWild pipe flow direction differs")
			fields(pipe_output, volume["nodes"], volume["elements"])
		print("native_tet_darcy_cli: PASS 1/2/4 MPI, preflight, rank-owned pressure/flux VTU, "
			"case/mesh provenance, no-overwrite, fail-closed inputs"
			+ ("; fTetWild pipe" if args.ftetwild else ""))


if __name__ == "__main__":
	main()
