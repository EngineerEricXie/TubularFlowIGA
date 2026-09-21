#!/usr/bin/env python3
"""Re-run a labelled ROI solved-flow-to-Darcy functional case with hash gates.

All physical values are explicit command-line inputs. This is not a patient
boundary-condition or physiological-validation workflow.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import resource
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

from audit_multiregion_tet_mesh import audit as audit_multiregion
from audit_multiregion_tet_mesh import parse_msh, read_contract


def reject(message):
	raise ValueError(message)


def sha256(path):
	return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
	try:
		return json.loads(path.read_text(encoding="utf-8"))
	except (OSError, json.JSONDecodeError) as error:
		reject(f"cannot read JSON {path}: {error}")


def positive(value, name):
	if not math.isfinite(value) or value <= 0:
		reject(f"{name} must be positive and finite")


def parse_result(stdout):
	lines = [line for line in stdout.splitlines()
		if line.startswith("native_tet_matching_solved_roi_smoke: PASS ")]
	if len(lines) != 1:
		reject("solver did not emit exactly one PASS metrics line")
	values = {}
	for field in lines[0].split()[2:]:
		if "=" not in field:
			reject("solver metrics line has an invalid field")
		key, raw = field.split("=", 1)
		if key in values:
			reject("solver metrics line has duplicate fields")
		value = float(raw)
		if not math.isfinite(value):
			reject("solver metric is nonfinite")
		values[key] = value
	required = {"facets", "fluid_inlet_m3_s", "interface_source_m3_s",
		"darcy_outward_m3_s", "fluid_newton_iterations", "fluid_linear_reason",
		"darcy_linear_reason", "darcy_max_cell_defect_m3_s"}
	if set(values) != required or values["facets"] < 1 \
			or values["fluid_linear_reason"] <= 0 or values["darcy_linear_reason"] <= 0:
		reject("solver metrics are incomplete or unconverged")
	return values


def tetra_volume(nodes, cell):
	a, b, c, d = (nodes[tag] for tag in cell)
	ab = [b[axis]-a[axis] for axis in range(3)]
	ac = [c[axis]-a[axis] for axis in range(3)]
	ad = [d[axis]-a[axis] for axis in range(3)]
	determinant = (ab[0]*(ac[1]*ad[2]-ac[2]*ad[1])
		-ab[1]*(ac[0]*ad[2]-ac[2]*ad[0])
		+ab[2]*(ac[0]*ad[1]-ac[1]*ad[0]))
	if not math.isfinite(determinant) or determinant <= 0:
		reject("rank-owned source audit found nonpositive tetra Jacobian")
	return determinant/6


def validate_field_output(root, ranks, region_cells, expected_facets,
		expected_source_m3_s, mesh_by_role, interface_label):
	"""Check complete, non-overlapping rank-owned cell pieces before publication."""
	files = {}
	owned_cells = {}
	tissue_source_by_cell = {}
	for role, expected_cells, arrays in (
		("fluid", region_cells["vessel"], ("velocity_m_s", "pressure_pa")),
		("tissue", region_cells["tissue"], ("pressure_pa", "darcy_rt0_centroid_flux_m_s",
			"conservative_outward_face_flow_m3_s", "source_s_inv"))):
		directory = root/role
		index = directory/"snapshot.pvtu"
		tree = ET.parse(index)
		pieces = tree.findall(".//Piece")
		if [piece.get("Source") for piece in pieces] != [
				f"rank{rank}.vtu" for rank in range(ranks)]:
			reject(f"{role} PVTU rank pieces are incomplete")
		files[f"{role}/snapshot.pvtu"] = sha256(index)
		seen = set()
		for rank in range(ranks):
			name = f"rank{rank}.vtu"
			path = directory/name
			piece = ET.parse(path).find(".//Piece")
			if piece is None:
				reject(f"{role} rank {rank} VTU piece is absent")
			ids = piece.find("./CellData/DataArray[@Name='GlobalCellIds']")
			if ids is None or ids.text is None:
				reject(f"{role} rank {rank} has no global cell IDs")
			values = [int(value) for value in ids.text.split()]
			if len(values) != int(piece.get("NumberOfCells", "-1")) \
					or len(values) != len(set(values)) or seen.intersection(values):
				reject(f"{role} rank-owned cell IDs overlap or are malformed")
			seen.update(values)
			available = {field.get("Name") for field in piece.findall("./PointData/DataArray")
				+ piece.findall("./CellData/DataArray")}
			if not set(arrays).issubset(available):
				reject(f"{role} rank {rank} solution arrays are incomplete")
			if role == "tissue":
				array = piece.find("./CellData/DataArray[@Name='source_s_inv']")
				if array is None:
					reject("tissue source field is absent")
				sources = [float(value) for value in (array.text or "").split()]
				if len(sources) != len(values) or not all(map(math.isfinite, sources)):
					reject("tissue source field is malformed")
				tissue_source_by_cell.update(zip(values, sources))
			files[f"{role}/{name}"] = sha256(path)
		if len(seen) != expected_cells:
			reject(f"{role} rank-owned cells do not cover the submesh")
		owned_cells[role] = seen
	mesh_cells = {}
	interface_triangles = {}
	for role in ("fluid", "tissue"):
		nodes, elements = parse_msh(mesh_by_role[role])
		mesh_cells[role] = {cell_id: tetra_volume(nodes, connectivity)
			for cell_id, dimension, _, _, connectivity in elements if dimension == 3}
		if set(mesh_cells[role]) != owned_cells[role]:
			reject(f"{role} rank-owned global cell IDs differ from native mesh")
		interface_triangles[role] = {
			tuple(sorted(nodes[tag] for tag in connectivity))
			for _, dimension, _, physical, connectivity in elements
			if dimension == 2 and physical == f"boundary_label_{interface_label}"}
		if len(interface_triangles[role]) != expected_facets:
			reject(f"{role} native mesh interface triangle count differs")
	if interface_triangles["fluid"] != interface_triangles["tissue"]:
		reject("native mesh interface geometry is not matching")
	ledger_path = root/"interface_facets.json"
	ledger = read_json(ledger_path)
	if ledger.get("schema_version") != 1 \
			or ledger.get("kind") != "native_matching_facet_flow_functional_only" \
			or ledger.get("flow_sign") != "vessel_outward_positive_tissue_source_positive" \
			or not isinstance(ledger.get("facets"), list) \
			or len(ledger["facets"]) != expected_facets:
		reject("interface facet ledger schema or facet count differs")
	seen_triangles = set()
	flow_sum = 0.0
	flow_by_tissue_cell = {}
	for facet in ledger["facets"]:
		if set(facet) != {"port_name", "vessel_cell_id", "tissue_cell_id",
				"vessel_outward_m3_s", "triangle_m"} \
				or facet["port_name"] != "solved_roi_functional" \
				or facet["vessel_cell_id"] not in owned_cells["fluid"] \
				or facet["tissue_cell_id"] not in owned_cells["tissue"] \
				or not isinstance(facet["vessel_outward_m3_s"], (int, float)) \
				or not math.isfinite(facet["vessel_outward_m3_s"]):
			reject("interface facet ledger owner or flow is invalid")
		triangle = facet["triangle_m"]
		if not isinstance(triangle, list) or len(triangle) != 3 \
				or any(not isinstance(point, list) or len(point) != 3
					or any(not isinstance(value, (int, float)) or not math.isfinite(value)
						for value in point) for point in triangle):
			reject("interface facet ledger triangle is invalid")
		key = tuple(tuple(point) for point in triangle)
		if key in seen_triangles:
			reject("interface facet ledger has duplicate geometry")
		seen_triangles.add(key)
		flow_sum += facet["vessel_outward_m3_s"]
		flow_by_tissue_cell[facet["tissue_cell_id"]] = flow_by_tissue_cell.get(
			facet["tissue_cell_id"], 0.)+facet["vessel_outward_m3_s"]
	if seen_triangles != interface_triangles["fluid"]:
		reject("interface facet ledger differs from native mesh geometry")
	if abs(flow_sum-expected_source_m3_s) > 1e-8*max(
			abs(expected_source_m3_s), 1e-18):
		reject("interface per-facet ledger does not sum to tissue source")
	for cell_id, volume in mesh_cells["tissue"].items():
		mapped_flow = flow_by_tissue_cell.get(cell_id, 0.)
		field_flow = tissue_source_by_cell[cell_id]*volume
		if abs(mapped_flow-field_flow) > 1e-8*max(abs(mapped_flow),
				abs(field_flow), 1e-18):
			reject("interface facet ledger differs from tissue cell source field")
	files["interface_facets.json"] = sha256(ledger_path)
	return files


def validate_rank_solver_resources(directory, ranks):
	if {path.name for path in directory.iterdir()} != {
			f"rank{rank}.time.txt" for rank in range(ranks)}:
		reject("native solver per-rank timing files are incomplete")
	measurements = {}
	for rank in range(ranks):
		path = directory/f"rank{rank}.time.txt"
		values = path.read_text(encoding="utf-8").split()
		if len(values) != 2:
			reject("native solver per-rank timing record is malformed")
		wall, rss = float(values[0]), int(values[1])
		if not math.isfinite(wall) or wall < 0 or rss <= 0:
			reject("native solver per-rank timing or RSS is invalid")
		measurements[str(rank)] = {"wall_seconds": wall,
			"gnu_time_maximum_rss_kb": rss, "timing_log_sha256": sha256(path)}
	return {"scope": "GNU time for each native MPI solver rank process; peak RSS is per rank, not summed or simultaneous total",
		"rank_metrics_directory": str(directory), "ranks": measurements}


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--source-seg", type=Path, required=True)
	parser.add_argument("--geometry-manifest", type=Path, required=True)
	parser.add_argument("--multiregion-mesh", type=Path, required=True)
	parser.add_argument("--multiregion-contract", type=Path, required=True)
	parser.add_argument("--split-manifest", type=Path, required=True)
	parser.add_argument("--vessel-region", required=True)
	parser.add_argument("--tissue-region", required=True)
	parser.add_argument("--interface-name", required=True)
	parser.add_argument("--inlet-label", type=int, required=True)
	parser.add_argument("--tissue-exit-labels", type=int, nargs=2, required=True)
	parser.add_argument("--inlet-velocity-m-s", type=float, nargs=3, required=True)
	parser.add_argument("--density-kg-m3", type=float, required=True)
	parser.add_argument("--viscosity-pa-s", type=float, required=True)
	parser.add_argument("--mobility-m2-pa-s", type=float, required=True)
	parser.add_argument("--interface-pressure-pa", type=float, required=True)
	parser.add_argument("--tissue-exit-pressure-pa", type=float, required=True)
	parser.add_argument("--ranks", type=int, required=True)
	parser.add_argument("--solver", type=Path,
		default=Path("solvers/cpu/native_tet_matching_solved_roi_smoke"))
	parser.add_argument("--output", type=Path, required=True)
	parser.add_argument("--field-output-dir", type=Path,
		help="optional immutable root for rank-owned fluid and tissue VTU snapshots")
	parser.add_argument("--case-file-sha256", help="validated case file SHA-256")
	args = parser.parse_args()
	try:
		for name in ("density_kg_m3", "viscosity_pa_s", "mobility_m2_pa_s"):
			positive(getattr(args, name), name)
		if args.case_file_sha256 is not None and (len(args.case_file_sha256) != 64
				or any(char not in "0123456789abcdef" for char in args.case_file_sha256)):
			reject("case file SHA-256 is malformed")
		if args.ranks < 1 or args.ranks > 64 or args.inlet_label < 0 \
				or any(label < 0 for label in args.tissue_exit_labels) \
				or args.tissue_exit_labels[0] == args.tissue_exit_labels[1] \
				or all(value == 0 for value in args.inlet_velocity_m_s) \
				or not all(math.isfinite(value) for value in [
					*args.inlet_velocity_m_s, args.interface_pressure_pa,
					args.tissue_exit_pressure_pa]):
			reject("rank, label, velocity, or pressure input is invalid")
		geometry = read_json(args.geometry_manifest)
		split = read_json(args.split_manifest)
		if geometry.get("kind") != "voxel_faithful_multiregion_roi_geometry_only" \
				or split.get("kind") != "native_single_region_split_geometry_only":
			reject("geometry or split manifest kind is invalid")
		if sha256(args.source_seg) != geometry.get("source_sha256") \
				or sha256(args.multiregion_mesh) != geometry.get("mesh_sha256") \
				or geometry["mesh_sha256"] != split.get("source_mesh_sha256") \
				or sha256(args.multiregion_contract) != split.get("source_contract_sha256"):
			reject("source SEG, multi-region mesh, or contract hash mismatch")
		_, region_map, boundary_map, interface_map = read_contract(args.multiregion_contract)
		current_geometry_audit = audit_multiregion(*parse_msh(args.multiregion_mesh),
			region_map, boundary_map, interface_map, 1.0e-3)
		if current_geometry_audit["tetrahedra"] != geometry.get("tetrahedra") \
				or not current_geometry_audit["exterior_and_interface_partition_exact"]:
			reject("current multi-region geometry audit differs from manifest")
		if args.vessel_region == args.tissue_region \
				or args.tissue_region != "tissue_candidate":
			reject("vessel and tissue region roles are invalid")
		submeshes = split.get("submeshes", {})
		if args.vessel_region not in submeshes or args.tissue_region not in submeshes:
			reject("requested region is absent from split manifest")
		paths = {}
		for region in (args.vessel_region, args.tissue_region):
			path = args.split_manifest.resolve().parent/submeshes[region]["file"]
			if path.resolve().parent != args.split_manifest.resolve().parent \
					or sha256(path) != submeshes[region]["sha256"]:
				reject("native submesh path or hash mismatch")
			paths[region] = path
		labels = split.get("interface_boundary_labels", {})
		if args.interface_name not in labels:
			reject("requested interface is absent from split manifest")
		interface = labels[args.interface_name]
		if interface in (args.inlet_label, *args.tissue_exit_labels) \
				or args.inlet_label in args.tissue_exit_labels:
			reject("interface, inlet, and exit labels must differ")
		solver = args.solver.resolve()
		if not solver.is_file():
			reject("native solver binary is absent")
		output = args.output.resolve()
		fields = args.field_output_dir.resolve() if args.field_output_dir else None
		if output.exists() or output in {args.source_seg.resolve(),
				args.geometry_manifest.resolve(), args.multiregion_mesh.resolve(),
				args.multiregion_contract.resolve(), args.split_manifest.resolve(),
				*map(Path.resolve, paths.values())}:
			reject("result output exists or overlaps an input")
		if fields is not None and (fields.exists() or fields == output
				or output.is_relative_to(fields) or fields.is_relative_to(output)):
			reject("field output directory exists or overlaps result output")
		native_command = [str(solver),
			str(paths[args.vessel_region]), str(paths[args.tissue_region]),
			str(interface), str(args.inlet_label), *map(str, args.tissue_exit_labels),
			*map(str, args.inlet_velocity_m_s), str(args.density_kg_m3),
			str(args.viscosity_pa_s), str(args.mobility_m2_pa_s),
			str(args.interface_pressure_pa), str(args.tissue_exit_pressure_pa)]
		if fields is not None:
			native_command.append(str(fields))
			metrics_dir = fields.parent/(fields.name+".rank-metrics")
			if metrics_dir.exists():
				reject("native solver rank metrics directory already exists")
			metrics_dir.mkdir(parents=True, exist_ok=False)
			command = ["mpiexec", "-np", str(args.ranks), sys.executable,
				str(Path(__file__).with_name("time_native_mpi_rank.py")),
				"--metrics-dir", str(metrics_dir), "--", *native_command]
		else:
			metrics_dir = None
			command = ["mpiexec", "-np", str(args.ranks), *native_command]
		start = time.perf_counter()
		completed = subprocess.run(command, text=True, stdout=subprocess.PIPE,
			stderr=subprocess.PIPE)
		wall = time.perf_counter()-start
		if completed.returncode:
			reject(f"native solver rejected case (exit {completed.returncode}): "
				+completed.stderr.strip()[-1000:])
		metrics = parse_result(completed.stdout)
		scale = max(abs(metrics["interface_source_m3_s"]), 1e-18)
		interface_pair_key = "_".join(sorted((args.vessel_region, args.tissue_region)))
		if metrics["facets"] != geometry["interfaces"][interface_pair_key] \
				or metrics["facets"] != current_geometry_audit["interfaces"][args.interface_name] \
				or abs(metrics["fluid_inlet_m3_s"]+metrics["interface_source_m3_s"]) > 1e-6*scale \
				or abs(metrics["interface_source_m3_s"]-metrics["darcy_outward_m3_s"]) > 1e-8*scale \
				or metrics["darcy_max_cell_defect_m3_s"] > 1e-8*scale:
			reject("reproducibility wrapper conservation or interface facet gate failed")
		field_hashes = None
		rank_resources = None
		if fields is not None:
			field_hashes = validate_field_output(fields, args.ranks, {
				"vessel": submeshes[args.vessel_region]["audit"]["tetrahedra"],
				"tissue": submeshes[args.tissue_region]["audit"]["tetrahedra"]},
				int(metrics["facets"]), metrics["interface_source_m3_s"],
				{"fluid": paths[args.vessel_region], "tissue": paths[args.tissue_region]},
				interface)
			rank_resources = validate_rank_solver_resources(metrics_dir, args.ranks)
		result = {"schema_version": 1,
			"kind": "patient_derived_roi_solved_fluid_darcy_functional_only",
			"case_file_sha256": args.case_file_sha256,
			"physiological_validation": False, "full_liver_case": False,
			"source_seg_sha256": geometry["source_sha256"],
			"multiregion_mesh_sha256": geometry["mesh_sha256"],
			"multiregion_contract_sha256": split["source_contract_sha256"],
			"vessel_mesh_sha256": submeshes[args.vessel_region]["sha256"],
			"tissue_mesh_sha256": submeshes[args.tissue_region]["sha256"],
			"solver_binary_sha256": sha256(solver),
			"regions": {"vessel": args.vessel_region, "tissue": args.tissue_region},
			"labels": {"interface": interface, "inlet": args.inlet_label,
				"tissue_exits": args.tissue_exit_labels},
			"explicit_functional_inputs": {"inlet_velocity_m_s": args.inlet_velocity_m_s,
				"density_kg_m3": args.density_kg_m3,
				"viscosity_pa_s": args.viscosity_pa_s,
				"mobility_m2_pa_s": args.mobility_m2_pa_s,
				"interface_pressure_pa": args.interface_pressure_pa,
				"tissue_exit_pressure_pa": args.tissue_exit_pressure_pa},
			"mpi_ranks": args.ranks, "metrics": metrics,
			"rank_owned_field_output_directory": str(fields) if fields is not None else None,
			"rank_owned_field_files_sha256": field_hashes,
			"rank_solver_resources": rank_resources,
			"wall_seconds": wall,
			"launcher_child_peak_rss_kb": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
			"gates": {"fluid_linear_converged": True,
				"darcy_linear_converged": True, "boundary_flow_closed": True,
				"interface_to_darcy_closed": True, "cell_balance_closed": True},
			"limitations": "ROI crop and user-supplied artificial functional BC/material; no physiological or full-liver validation"}
		output.parent.mkdir(parents=True, exist_ok=True)
		with output.open("x", encoding="utf-8") as stream:
			stream.write(json.dumps(result, indent=2, sort_keys=True)+"\n")
	except (AttributeError, IndexError, KeyError, OSError, ValueError) as error:
		print(f"native matching ROI functional run rejected: {error}", file=sys.stderr)
		return 2
	print(f"native matching ROI functional run: PASS ranks={args.ranks} output={output}")
	return 0


if __name__ == "__main__":
	sys.exit(main())
