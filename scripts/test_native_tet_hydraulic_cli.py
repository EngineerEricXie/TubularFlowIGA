#!/usr/bin/env python3
"""Integration test for the file-backed native tetra/0D command-line runner."""

import argparse
import json
import math
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path


def run(binary, ranks, case, *options, expect_success=True):
	command = ["mpiexec", "-np", str(ranks), str(binary), str(case), *map(str, options)]
	result = subprocess.run(command, capture_output=True, text=True,
		timeout=120, check=False)
	if expect_success and result.returncode:
		raise RuntimeError(f"native CLI failed: {' '.join(command)}\n"
			f"{result.stdout}\n{result.stderr}")
	if not expect_success and result.returncode == 0:
		raise RuntimeError(f"native CLI unexpectedly accepted: {' '.join(command)}")
	return result.stdout


def step_values(output, step):
	for line in output.splitlines():
		if not line.startswith(f"native_hydraulic_step step={step} "):
			continue
		return {key: float(value) for key, value in
			(token.split("=", 1) for token in line.split()[1:])}
	raise RuntimeError(f"native CLI did not publish step {step}")


def vtk_snapshot(directory, step, expected_points=15, expected_cells=4,
		expected_first_cell=(4, 1, 2, 3, 5, 8, 6, 7, 9, 10),
		expected_displacement_x=None, return_stress=False):
	snapshot = directory/f"step_{step}"
	index = ET.parse(snapshot/"snapshot.pvtu").getroot()
	pieces = [item.attrib["Source"] for item in index.findall(".//Piece")]
	if not pieces:
		raise RuntimeError("native CLI PVTU index has no pieces")
	parallel_fields = {item.attrib.get("Name") for item in index.findall(".//PPointData/PDataArray")}
	if not {"velocity_m_s", "pressure_pa", "reference_position_m",
			"displacement_m"}.issubset(parallel_fields):
		raise RuntimeError("native CLI PVTU field schema is incomplete")
	parallel_cell_fields = {item.attrib.get("Name"): item.attrib for item in index.findall(".//PCellData/PDataArray")}
	if parallel_cell_fields.get("cauchy_stress_pa", {}).get("NumberOfComponents") != "9":
		raise RuntimeError("native CLI PVTU stress field is missing")
	points = {}
	cells = {}
	stresses = {}
	for name in pieces:
		piece = ET.parse(snapshot/name).getroot().find(".//Piece")
		if piece is None:
			raise RuntimeError("native CLI VTU piece is missing")
		data = lambda path: piece.find(path).text or ""
		ids = [int(value) for value in data("PointData/DataArray[@Name='GlobalPointIds']").split()]
		coordinates = [float(value) for value in data("Points/DataArray").split()]
		velocity = [float(value) for value in data("PointData/DataArray[@Name='velocity_m_s']").split()]
		pressure = [float(value) for value in data("PointData/DataArray[@Name='pressure_pa']").split()]
		reference = [float(value) for value in data("PointData/DataArray[@Name='reference_position_m']").split()]
		displacement = [float(value) for value in data("PointData/DataArray[@Name='displacement_m']").split()]
		if (len(coordinates) != 3*len(ids) or len(velocity) != 3*len(ids)
				or len(pressure) != len(ids) or len(reference) != 3*len(ids)
				or len(displacement) != 3*len(ids)):
			raise RuntimeError("native CLI VTU point arrays are incomplete")
		for i, identifier in enumerate(ids):
			value = (tuple(coordinates[3*i:3*i+3]),
				tuple(velocity[3*i:3*i+3]), pressure[i],
				tuple(reference[3*i:3*i+3]),
				tuple(displacement[3*i:3*i+3]))
			if not all(math.isfinite(number) for field in (value[0], value[1],
					(value[2],), value[3], value[4]) for number in field):
				raise RuntimeError("native CLI VTU field contains a nonfinite value")
			for axis in range(3):
				if abs(value[0][axis]-value[3][axis]-value[4][axis]) > 1e-10:
					raise RuntimeError("native CLI reference/current displacement differs")
				if expected_displacement_x is not None:
					expected = expected_displacement_x if axis == 0 else 0.0
					if abs(value[4][axis]-expected) > 1e-10:
						raise RuntimeError("native CLI prescribed ALE displacement differs")
			if identifier in points and points[identifier] != value:
				raise RuntimeError("native CLI shared VTU point differs between ranks")
			points[identifier] = value
		connectivity = [int(value) for value in data("Cells/DataArray[@Name='connectivity']").split()]
		offsets = [int(value) for value in data("Cells/DataArray[@Name='offsets']").split()]
		types = [int(value) for value in data("Cells/DataArray[@Name='types']").split()]
		cell_ids = [int(value) for value in data("CellData/DataArray[@Name='GlobalCellIds']").split()]
		stress = [float(value) for value in data("CellData/DataArray[@Name='cauchy_stress_pa']").split()]
		if len(offsets) != len(types) or len(types) != len(cell_ids):
			raise RuntimeError("native CLI VTU cell arrays are incomplete")
		if len(stress) != 9*len(cell_ids) or not all(math.isfinite(value) for value in stress):
			raise RuntimeError("native CLI VTU Cauchy stress is incomplete or nonfinite")
		begin = 0
		for cell_index, (cell_id, end, cell_type) in enumerate(zip(cell_ids, offsets, types)):
			if cell_type != 24 or end-begin != 10:
				raise RuntimeError("native CLI VTU cell is not a quadratic tetrahedron")
			cells[cell_id] = tuple(ids[local] for local in connectivity[begin:end])
			stresses[cell_id] = tuple(stress[9*cell_index:9*(cell_index+1)])
			begin = end
		if begin != len(connectivity):
			raise RuntimeError("native CLI VTU connectivity has trailing entries")
	if ((expected_points is not None and len(points) != expected_points)
			or len(cells) != expected_cells):
		raise RuntimeError("native CLI VTU lost native P2 nodes or tetrahedra")
	if expected_first_cell is not None and cells.get(5) != expected_first_cell:
		raise RuntimeError("native P2 edge nodes were not reordered for VTK")
	for cell in cells.values():
		for mid, (first, second) in zip(cell[4:],
			((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))):
			expected_pressure = 0.5*(points[cell[first]][2]+points[cell[second]][2])
			if abs(points[mid][2]-expected_pressure) > 1e-10*max(1.0, abs(expected_pressure)):
				raise RuntimeError("native CLI P1 pressure is wrong at a P2 midside node")
	return (points, cells, stresses) if return_stress else (points, cells)


def compare_snapshots(reference_directory, resumed_directory, step,
		expected_points=15, expected_cells=4, expected_first_cell=(4, 1, 2, 3, 5, 8, 6, 7, 9, 10),
		expected_displacement_x=None):
	reference_points, reference_cells, reference_stress = vtk_snapshot(reference_directory, step,
		expected_points, expected_cells, expected_first_cell, expected_displacement_x, True)
	resumed_points, resumed_cells, resumed_stress = vtk_snapshot(resumed_directory, step,
		expected_points, expected_cells, expected_first_cell, expected_displacement_x, True)
	if reference_cells != resumed_cells or reference_points.keys() != resumed_points.keys():
		raise RuntimeError("native CLI restart changed quadratic tetra topology")
	if reference_stress.keys() != resumed_stress.keys():
		raise RuntimeError("native CLI restart changed stress cell ids")
	for cell_id, values in reference_stress.items():
		for actual, reference in zip(resumed_stress[cell_id], values):
			if abs(actual-reference) > 1e-10*max(1.0, abs(reference)):
				raise RuntimeError("native CLI restart changed Cauchy stress")
	for identifier, values in reference_points.items():
		actual = resumed_points[identifier]
		actual_values = (*actual[0], *actual[1], actual[2], *actual[3], *actual[4])
		reference_values = (*values[0], *values[1], values[2], *values[3], *values[4])
		for actual, reference in zip(actual_values, reference_values):
			if abs(actual-reference) > 1e-10*max(1.0, abs(reference)):
				raise RuntimeError("native CLI restart changed the full VTU field")
	return reference_points


def channel_case(directory, reference_case, split_outlet=False):
	"""Write the 24-tetra, three-boundary channel used by the native graph tests."""
	def node(i, j, k):
		return (i*3+j)*3+k
	points = [(float(i), 0.5*j, 0.5*k)
		for i in range(2) for j in range(3) for k in range(3)]
	cells = []
	for j in range(2):
		for k in range(2):
			a, b = node(0, j, k), node(1, j, k)
			c, d = node(0, j+1, k), node(1, j+1, k)
			e, f = node(0, j, k+1), node(1, j, k+1)
			g, h = node(0, j+1, k+1), node(1, j+1, k+1)
			for tetra in ((a, b, d, h), (a, d, c, h), (a, c, g, h),
				(a, g, e, h), (a, e, f, h), (a, f, b, h)):
				x0, x1, x2, x3 = (points[index] for index in tetra)
				u = [x1[axis]-x0[axis] for axis in range(3)]
				v = [x2[axis]-x0[axis] for axis in range(3)]
				w = [x3[axis]-x0[axis] for axis in range(3)]
				determinant = (u[0]*(v[1]*w[2]-v[2]*w[1])
					-u[1]*(v[0]*w[2]-v[2]*w[0])
					+u[2]*(v[0]*w[1]-v[1]*w[0]))
				if determinant < 0:
					tetra = (tetra[0], tetra[2], tetra[1], tetra[3])
				elif determinant == 0:
					raise RuntimeError("native channel has a degenerate tetrahedron")
				cells.append(tetra)
	faces = {}
	for tetra in cells:
		for opposite in range(4):
			face = tuple(tetra[local] for local in range(4) if local != opposite)
			key = tuple(sorted(face))
			faces.setdefault(key, []).append(face)
	by_label = {label: [] for label in ((1, 2, 3, 4) if split_outlet
		else (1, 2, 3))}
	for owners in faces.values():
		if len(owners) != 1:
			continue
		face = owners[0]
		x = [points[index][0] for index in face]
		label = 1 if all(value == 1.0 for value in x) else (
			(4 if split_outlet and sum(points[index][1] for index in face)/3 > 0.5
				else 2) if all(value == 0.0 for value in x) else 3)
		by_label[label].append(face)
	if len(points) != 18 or len(cells) != 24 or any(not faces for faces in by_label.values()):
		raise RuntimeError("native channel fixture topology differs")
	count = sum(map(len, by_label.values()))+len(cells)
	lines = ["$MeshFormat", "4.1 0 8", "$EndMeshFormat",
		"$PhysicalNames", str(len(by_label)+1)]
	lines += [f'2 {label} "boundary_label_{label}"' for label in by_label]
	lines += ['3 1 "fluid"', "$EndPhysicalNames",
		"$Entities", f"0 0 {len(by_label)} 1"]
	lines += [f"{label} 0 0 0 1 1 1 1 {label} 0" for label in by_label]
	lines += ["1 0 0 0 1 1 1 1 1 0", "$EndEntities",
		"$Nodes", "1 18 1 18", "3 1 0 18"]
	lines += [str(index+1) for index in range(len(points))]
	lines += [" ".join(map(str, point)) for point in points]
	lines += ["$EndNodes", "$Elements", f"{len(by_label)+1} {count} 1 {count}"]
	element_id = 1
	for label, surface in by_label.items():
		lines.append(f"2 {label} 2 {len(surface)}")
		for face in surface:
			lines.append(" ".join(map(str, (element_id, *(index+1 for index in face)))))
			element_id += 1
	lines.append(f"3 1 4 {len(cells)}")
	for tetra in cells:
		lines.append(" ".join(map(str, (element_id, *(index+1 for index in tetra)))))
		element_id += 1
	lines.append("$EndElements")
	mesh = directory/("native_channel_multi.msh" if split_outlet else "native_channel.msh")
	mesh.write_text("\n".join(lines)+"\n")
	case = json.loads(reference_case.read_text())
	case["mesh_file"] = str(mesh)
	case["boundaries"]["wall_labels"] = [3]
	case["fluid"]["density_kg_m3"] = 1.0
	case["fluid"]["dynamic_viscosity_pa_s"] = 0.1
	case["fluid"]["initial_velocity_x_m_s"] = -0.01
	case["source"] = {"capacitance_m3_pa": 0.01,
		"resistance_pa_s_m3": 10.0, "prescribed_flow_m3_s": 0.01,
		"initial_pressure_pa": 0.0}
	case["terminal_rcr"] = {"proximal_resistance_pa_s_m3": 1.0,
		"distal_resistance_pa_s_m3": 10.0, "capacitance_m3_pa": 0.01,
		"distal_pressure_pa": 0.0, "initial_pressure_pa": 0.0}
	case["coupling"]["method"] = "explicit"
	case["coupling"]["maximum_iterations"] = 1
	if split_outlet:
		case["schema_version"] = 2
		del case["boundaries"]["outlet_label"]
		terminal = case.pop("terminal_rcr")
		case["terminal_rcrs"] = [dict(terminal, boundary_label=label)
			for label in (2, 4)]
		case["terminal_rcrs"][1]["distal_resistance_pa_s_m3"] = 20.0
	case_path = directory/("native_channel_multi.json" if split_outlet
		else "native_channel.json")
	case_path.write_text(json.dumps(case))
	return case_path


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--ranks", type=int, required=True)
	parser.add_argument("--binary", type=Path, required=True)
	parser.add_argument("--case", type=Path, required=True)
	args = parser.parse_args()
	if args.ranks < 1 or args.ranks > 16:
		raise ValueError("native CLI test ranks must be between 1 and 16")
	binary = args.binary.resolve(strict=True)
	case = args.case.resolve(strict=True)
	with tempfile.TemporaryDirectory(prefix="iga-native-hydraulic-cli-") as directory:
		root = Path(directory)
		checked = run(binary, args.ranks, case, "--check-input")
		if "native_hydraulic_input: PASS" not in checked:
			raise RuntimeError("native CLI input preflight did not inspect its case")
		run(binary, args.ranks, case, "--check-input",
			"--output-dir", root/"forbidden_preflight_output", expect_success=False)
		if (root/"forbidden_preflight_output").exists():
			raise RuntimeError("native CLI input preflight created solver output")
		baseline = run(binary, args.ranks, case,
			"--output-dir", root/"baseline")
		published_piece = root/"baseline"/"step_1"/"rank0.vtu"
		published_bytes = published_piece.read_bytes()
		run(binary, args.ranks, case,
			"--output-dir", root/"baseline", expect_success=False)
		if published_piece.read_bytes() != published_bytes:
			raise RuntimeError("native CLI overwrote a published VTU snapshot")
		first = run(binary, args.ranks, case,
			"--checkpoint-dir", root/"checkpoints", "--output-dir", root/"split",
			"--stop-after-step", 1)
		if "native_hydraulic_step step=1 " not in first:
			raise RuntimeError("native CLI did not save the first accepted step")
		resumed = run(binary, args.ranks, case,
			"--restart-dir", root/"checkpoints", "--output-dir", root/"split")
		if "native_hydraulic_step step=1 " in resumed:
			raise RuntimeError("native CLI replayed an accepted step after restart")
		for key, reference in step_values(baseline, 2).items():
			actual = step_values(resumed, 2)[key]
			if not math.isfinite(actual) or abs(actual-reference) > 1e-10*max(1.0, abs(reference)):
				raise RuntimeError(f"native CLI restart {key} differs: {actual} vs {reference}")
		compare_snapshots(root/"baseline", root/"split", 2,
			expected_displacement_x=0.0001)
		vtk_snapshot(root/"split", 1, expected_displacement_x=0.00005)
		changed = json.loads(case.read_text())
		changed["mesh_file"] = str((case.parent/changed["mesh_file"]).resolve())
		changed["motion"]["speed_x_m_s"] *= 2
		wrong_case = root/"changed_case.json"
		wrong_case.write_text(json.dumps(changed))
		rejected = run(binary, args.ranks, wrong_case,
			"--restart-dir", root/"checkpoints", expect_success=False)
		if "native_hydraulic_step" in rejected:
			raise RuntimeError("wrong native CLI case published a step")
		wrong_labels = json.loads(case.read_text())
		wrong_labels["mesh_file"] = str((case.parent/wrong_labels["mesh_file"]).resolve())
		wrong_labels["boundaries"]["outlet_label"] = 3
		wrong_label_case = root/"wrong_labels.json"
		wrong_label_case.write_text(json.dumps(wrong_labels))
		run(binary, args.ranks, wrong_label_case, expect_success=False)
		run(binary, args.ranks, wrong_label_case, "--check-input",
			expect_success=False)
		mesh_text = (case.parent/json.loads(case.read_text())["mesh_file"]).read_text()
		if "8 1 2 3 5" not in mesh_text:
			raise RuntimeError("native tetra test fixture orientation is unknown")
		bad_mesh = root/"inverted.msh"
		bad_mesh.write_text(mesh_text.replace("8 1 2 3 5", "8 1 3 2 5"))
		inverted = json.loads(case.read_text())
		inverted["mesh_file"] = str(bad_mesh)
		inverted_case = root/"inverted_case.json"
		inverted_case.write_text(json.dumps(inverted))
		run(binary, args.ranks, inverted_case, expect_success=False)
		run(binary, args.ranks, inverted_case, "--check-input",
			expect_success=False)
		channel = channel_case(root, case)
		channel_result = run(binary, args.ranks, channel,
			"--output-dir", root/"channel_fields")
		channel_step = step_values(channel_result, 2)
		if not (channel_step["outlet_m3_s"] > 1e-5
				and channel_step["inlet_m3_s"] < -1e-5):
			raise RuntimeError("native no-slip channel has no throughflow")
		channel_points, _ = vtk_snapshot(root/"channel_fields", 2,
			expected_points=None, expected_cells=24, expected_first_cell=None,
			expected_displacement_x=0.0001)
		if len(channel_points) <= 18:
			raise RuntimeError("native no-slip channel has no P2 midside nodes")
		for coordinate, velocity, *_ in channel_points.values():
			if any(abs(axis-edge) < 1e-12 for axis in coordinate[1:]
				for edge in (0.0, 1.0)):
				if (abs(velocity[0]-0.001) > 1e-8
						or abs(velocity[1]) > 1e-8 or abs(velocity[2]) > 1e-8):
					raise RuntimeError("native no-slip channel wall velocity differs")
		run(binary, args.ranks, channel,
			"--checkpoint-dir", root/"channel_checkpoint",
			"--output-dir", root/"channel_split", "--stop-after-step", 1)
		channel_resumed = run(binary, args.ranks, channel,
			"--restart-dir", root/"channel_checkpoint",
			"--output-dir", root/"channel_split")
		for key, reference in channel_step.items():
			actual = step_values(channel_resumed, 2)[key]
			if abs(actual-reference) > 1e-10*max(1.0, abs(reference)):
				raise RuntimeError(f"native channel restart {key} differs")
		compare_snapshots(root/"channel_fields", root/"channel_split", 2,
			expected_points=None, expected_cells=24, expected_first_cell=None,
			expected_displacement_x=0.0001)
		missing_wall = json.loads(channel.read_text())
		missing_wall["boundaries"]["wall_labels"] = []
		missing_wall_case = root/"missing_wall.json"
		missing_wall_case.write_text(json.dumps(missing_wall))
		run(binary, args.ranks, missing_wall_case, expect_success=False)
		multi = channel_case(root, case, split_outlet=True)
		multi_output = run(binary, args.ranks, multi,
			"--output-dir", root/"multi_fields")
		multi_step = step_values(multi_output, 2)
		if not (multi_step["inlet_m3_s"] < -1e-5
				and multi_step["outlet_2_m3_s"] > 1e-5
				and multi_step["outlet_4_m3_s"] > 1e-5
				and abs(multi_step["inlet_m3_s"]+multi_step["outlet_2_m3_s"]
					+multi_step["outlet_4_m3_s"]) < 1e-6):
			raise RuntimeError("native multi-outlet channel flow does not balance")
		vtk_snapshot(root/"multi_fields", 2, expected_points=None,
			expected_cells=24, expected_first_cell=None,
			expected_displacement_x=0.0001)
		run(binary, args.ranks, multi,
			"--checkpoint-dir", root/"multi_checkpoint",
			"--output-dir", root/"multi_split", "--stop-after-step", 1)
		multi_resumed = run(binary, args.ranks, multi,
			"--restart-dir", root/"multi_checkpoint",
			"--output-dir", root/"multi_split")
		for key, reference in multi_step.items():
			actual = step_values(multi_resumed, 2)[key]
			if abs(actual-reference) > 1e-10*max(1.0, abs(reference)):
				raise RuntimeError(f"native multi-outlet restart {key} differs")
		compare_snapshots(root/"multi_fields", root/"multi_split", 2,
			expected_points=None, expected_cells=24, expected_first_cell=None,
			expected_displacement_x=0.0001)
		changed_terminal = json.loads(multi.read_text())
		changed_terminal["terminal_rcrs"][1]["distal_resistance_pa_s_m3"] *= 2
		changed_terminal_case = root/"changed_terminal.json"
		changed_terminal_case.write_text(json.dumps(changed_terminal))
		run(binary, args.ranks, changed_terminal_case,
			"--restart-dir", root/"multi_checkpoint", expect_success=False)
		duplicate_outlet = json.loads(multi.read_text())
		duplicate_outlet["terminal_rcrs"][1]["boundary_label"] = 2
		duplicate_case = root/"duplicate_outlet.json"
		duplicate_case.write_text(json.dumps(duplicate_outlet))
		run(binary, args.ranks, duplicate_case, expect_success=False)
		multi_strong = json.loads(multi.read_text())
		multi_strong["coupling"]["method"] = "fixed"
		multi_strong["coupling"]["maximum_iterations"] = 20
		strong_multi_case = root/"strong_multi.json"
		strong_multi_case.write_text(json.dumps(multi_strong))
		strong_multi_step = step_values(run(binary, args.ranks, strong_multi_case), 2)
		if not (strong_multi_step["iterations"] >= 2
				and strong_multi_step["outlet_2_m3_s"] > 1e-5
				and strong_multi_step["outlet_4_m3_s"] > 1e-5):
			raise RuntimeError("native strongly coupled two-outlet channel failed")
		multi_nonconverged = json.loads(strong_multi_case.read_text())
		multi_nonconverged["coupling"]["maximum_iterations"] = 1
		multi_nonconverged_case = root/"multi_nonconverged.json"
		multi_nonconverged_case.write_text(json.dumps(multi_nonconverged))
		nonconverged = subprocess.run(["mpiexec", "-np", str(args.ranks),
			str(binary), str(multi_nonconverged_case)], capture_output=True,
			text=True, timeout=120, check=False)
		if nonconverged.returncode == 0 or "edge=source_tet" not in nonconverged.stderr:
			raise RuntimeError("native fixed-point failure omitted edge diagnostics")
		multi_reverse = json.loads(multi.read_text())
		multi_reverse["source"]["prescribed_flow_m3_s"] = -0.01
		multi_reverse["fluid"]["initial_velocity_x_m_s"] = 0.01
		multi_reverse_case = root/"multi_reverse.json"
		multi_reverse_case.write_text(json.dumps(multi_reverse))
		multi_reverse_step = step_values(run(binary, args.ranks, multi_reverse_case), 2)
		if not (multi_reverse_step["inlet_m3_s"] > 1e-5
				and multi_reverse_step["outlet_2_m3_s"] < -1e-5
				and multi_reverse_step["outlet_4_m3_s"] < -1e-5):
			raise RuntimeError("native two-outlet reverse flow is absent")
		multi_reverse["coupling"]["method"] = "fixed"
		multi_reverse["coupling"]["maximum_iterations"] = 20
		multi_reverse_fixed_case = root/"multi_reverse_fixed.json"
		multi_reverse_fixed_case.write_text(json.dumps(multi_reverse))
		multi_reverse_fixed = step_values(run(binary, args.ranks,
			multi_reverse_fixed_case), 2)
		if not (multi_reverse_fixed["iterations"] >= 2
				and multi_reverse_fixed["inlet_m3_s"] > 1e-5
				and multi_reverse_fixed["outlet_2_m3_s"] < -1e-5
				and multi_reverse_fixed["outlet_4_m3_s"] < -1e-5):
			raise RuntimeError("native strongly coupled two-outlet reverse flow failed")
		reverse = json.loads(channel.read_text())
		reverse["source"]["prescribed_flow_m3_s"] = -0.01
		reverse["fluid"]["initial_velocity_x_m_s"] = 0.01
		reverse_case = root/"reverse_channel.json"
		reverse_case.write_text(json.dumps(reverse))
		reverse_output = run(binary, args.ranks, reverse_case)
		reverse_step = step_values(reverse_output, 2)
		if not (reverse_step["inlet_m3_s"] > 1e-5
				and reverse_step["outlet_m3_s"] < -1e-5):
			raise RuntimeError("native no-slip channel reverse flow is absent")
		reverse["coupling"]["method"] = "fixed"
		reverse["coupling"]["maximum_iterations"] = 20
		strong_reverse_case = root/"strong_reverse_channel.json"
		strong_reverse_case.write_text(json.dumps(reverse))
		strong_reverse = step_values(run(binary, args.ranks, strong_reverse_case), 2)
		if not (strong_reverse["iterations"] >= 2
				and strong_reverse["inlet_m3_s"] > 1e-5
				and strong_reverse["outlet_m3_s"] < -1e-5):
			raise RuntimeError("native strongly coupled reverse channel failed")
	print(f"native tetra/0D CLI quadratic VTU, restart, multi-outlet and reverse-flow channel: PASS ranks={args.ranks}")


if __name__ == "__main__":
	main()
