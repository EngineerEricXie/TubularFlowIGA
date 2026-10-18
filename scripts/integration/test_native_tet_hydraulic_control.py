#!/usr/bin/env python3
"""Integration test for the opt-in control features of native_tet_hydraulic_graph.

It covers the implicit 0D predictor (one coupling iteration, same answer as a
tightly converged fixed point), steady momentum, the --summary-json record,
controlled restarts that change declared boundary values (and every rejection
around them), and more than 16 terminal RCRs. The fixtures are the small
channels of test_native_tet_hydraulic_cli.py plus an n-by-n outlet grid.
"""

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_native_tet_hydraulic_cli import channel_case, step_values  # noqa: E402


def run(binary, ranks, case, *options, expect_success=True, environment=None):
	command = ["mpiexec", "-np", str(ranks), str(binary), str(case), *map(str, options)]
	result = subprocess.run(command, capture_output=True, text=True, timeout=300,
		check=False, env=environment)
	if expect_success and result.returncode:
		raise RuntimeError(f"native CLI failed: {' '.join(command)}\n{result.stdout}\n{result.stderr}")
	if not expect_success and result.returncode == 0:
		raise RuntimeError(f"native CLI unexpectedly accepted: {' '.join(command)}")
	return result.stdout+result.stderr


def solver_values(output, step):
	for line in output.splitlines():
		if line.startswith(f"native_hydraulic_solver step={step} "):
			return {key: float(value) for key, value in
				(token.split("=", 1) for token in line.split()[1:])}
	raise RuntimeError(f"native CLI did not report solver work for step {step}")


def write_case(path, case):
	path.write_text(json.dumps(case, indent=1))
	return path


def close(actual, reference, tolerance):
	return abs(actual-reference) <= tolerance*max(abs(reference), 1e-30)


def write_msh(path, points, tetrahedra, labels):
	"""ASCII Gmsh 4.1 volume with one physical surface per boundary label."""
	count = sum(map(len, labels.values()))+len(tetrahedra)
	lines = ["$MeshFormat", "4.1 0 8", "$EndMeshFormat", "$PhysicalNames", str(len(labels)+1)]
	lines += [f'2 {label} "boundary_label_{label}"' for label in sorted(labels)]
	lines += ['3 1 "fluid"', "$EndPhysicalNames", "$Entities", f"0 0 {len(labels)} 1"]
	lines += [f"{label} 0 0 0 1 1 1 1 {label} 0" for label in sorted(labels)]
	lines += ["1 0 0 0 1 1 1 1 1 0", "$EndEntities", "$Nodes",
		f"1 {len(points)} 1 {len(points)}", f"3 1 0 {len(points)}"]
	lines += [str(index+1) for index in range(len(points))]
	lines += [" ".join(map(repr, point)) for point in points]
	lines += ["$EndNodes", "$Elements", f"{len(labels)+1} {count} 1 {count}"]
	element = 1
	for label in sorted(labels):
		lines.append(f"2 {label} 2 {len(labels[label])}")
		for face in labels[label]:
			lines.append(" ".join(map(str, (element, *(index+1 for index in face)))))
			element += 1
	lines.append(f"3 1 4 {len(tetrahedra)}")
	for tetra in tetrahedra:
		lines.append(" ".join(map(str, (element, *(index+1 for index in tetra)))))
		element += 1
	lines.append("$EndElements")
	path.write_text("\n".join(lines)+"\n")


def wall_locked_case(directory, reference_case):
	"""Two tetrahedra; the lower one has three no-slip wall faces, so its apex
	pressure has no free P2 velocity node and P2/P1 is singular."""
	points = [(0., 0., 0.), (1., 0., 0.), (0., 1., 0.), (.2, .2, -.5), (.3, .3, .6)]
	tetrahedra = [(0, 2, 1, 3), (0, 1, 2, 4)]
	labels = {0: [(0, 3, 1), (1, 3, 2), (0, 2, 3), (0, 1, 4)], 1: [(1, 2, 4)], 2: [(0, 4, 2)]}
	mesh = directory/"wall_locked.msh"
	write_msh(mesh, points, tetrahedra, labels)
	case = json.loads(reference_case.read_text())
	case["mesh_file"] = str(mesh)
	case["boundaries"]["wall_labels"] = [0]
	case["motion"] = {"kind": "fixed", "speed_x_m_s": 0.0}
	return write_case(directory/"wall_locked.json", case)


def outlet_grid_case(directory, reference_case, cells):
	"""Unit cube, one cell deep in x: inlet label 1 at x=1, wall label 3, and
	one outlet label per (y,z) square at x=0, numbered 10, 11, ..."""
	def node(i, j, k):
		return (i*(cells+1)+j)*(cells+1)+k
	h = 1.0/cells
	points = [(float(i), h*j, h*k) for i in range(2) for j in range(cells+1) for k in range(cells+1)]
	tetrahedra = []
	for j in range(cells):
		for k in range(cells):
			a, b = node(0, j, k), node(1, j, k)
			c, d = node(0, j+1, k), node(1, j+1, k)
			e, f = node(0, j, k+1), node(1, j, k+1)
			g, hh = node(0, j+1, k+1), node(1, j+1, k+1)
			for tetra in ((a, b, d, hh), (a, d, c, hh), (a, c, g, hh),
					(a, g, e, hh), (a, e, f, hh), (a, f, b, hh)):
				x0, x1, x2, x3 = (points[index] for index in tetra)
				u = [x1[axis]-x0[axis] for axis in range(3)]
				v = [x2[axis]-x0[axis] for axis in range(3)]
				w = [x3[axis]-x0[axis] for axis in range(3)]
				determinant = (u[0]*(v[1]*w[2]-v[2]*w[1])-u[1]*(v[0]*w[2]-v[2]*w[0])
					+u[2]*(v[0]*w[1]-v[1]*w[0]))
				if determinant < 0:
					tetra = (tetra[0], tetra[2], tetra[1], tetra[3])
				tetrahedra.append(tetra)
	faces = {}
	for tetra in tetrahedra:
		for opposite in range(4):
			face = tuple(tetra[local] for local in range(4) if local != opposite)
			faces.setdefault(tuple(sorted(face)), []).append(face)
	labels = {}
	for owners in faces.values():
		if len(owners) != 1:
			continue
		face = owners[0]
		x = [points[index][0] for index in face]
		if all(value == 1.0 for value in x):
			label = 1
		elif all(value == 0.0 for value in x):
			y = sum(points[index][1] for index in face)/3
			z = sum(points[index][2] for index in face)/3
			label = 10+int(y/h)*cells+int(z/h)
		else:
			label = 3
		labels.setdefault(label, []).append(face)
	mesh = directory/f"outlet_grid_{cells}.msh"
	write_msh(mesh, points, tetrahedra, labels)
	case = json.loads(reference_case.read_text())
	case["schema_version"] = 2
	case["mesh_file"] = str(mesh)
	case["boundaries"] = {"inlet_label": 1, "wall_labels": [3]}
	case["fluid"].update({"density_kg_m3": 1.0, "dynamic_viscosity_pa_s": 0.1,
		"initial_velocity_x_m_s": 0.0})
	case["source"] = {"capacitance_m3_pa": 0.01, "resistance_pa_s_m3": 10.0,
		"prescribed_flow_m3_s": 0.01, "initial_pressure_pa": 0.0}
	case.pop("terminal_rcr")
	outlets = sorted(label for label in labels if label >= 10)
	case["terminal_rcrs"] = [{"boundary_label": label, "proximal_resistance_pa_s_m3": 1.0,
		"distal_resistance_pa_s_m3": 10.0*(1+index % 5), "capacitance_m3_pa": 0.01,
		"distal_pressure_pa": 0.0, "initial_pressure_pa": 0.0}
		for index, label in enumerate(outlets)]
	case["motion"] = {"kind": "fixed", "speed_x_m_s": 0.0}
	case["coupling"] = {"method": "fixed", "maximum_iterations": 4,
		"pressure_relative_tolerance": 1e-6, "flow_relative_tolerance": 1e-8,
		"relaxation_factor": 0.5, "predictor": "implicit_0d"}
	return write_case(directory/f"outlet_grid_{cells}.json", case), outlets


def summary_checks(summary, output, outlets, steps):
	if summary.get("kind") != "native_tet_hydraulic_summary" or len(summary["steps"]) != steps:
		raise RuntimeError("summary JSON has the wrong kind or step count")
	for record in summary["steps"]:
		step = record["step"]
		printed = step_values(output, step)
		if [item["label"] for item in record["outlets"]] != outlets:
			raise RuntimeError("summary JSON outlet labels differ")
		for item in record["outlets"]:
			key = f"outlet_{item['label']}_m3_s"
			if key in printed and not close(item["outward_flow_m3_s"], printed[key], 1e-12):
				raise RuntimeError("summary JSON outlet flow differs from the step line")
		inflow = record["inlet"]["total_inflow_m3_s"]
		if not close(record["total_outlet_flow_m3_s"], inflow, 1e-6):
			raise RuntimeError("summary JSON flows do not balance")
		solver = record["solver"]
		if not (record["peak_rss_bytes_max_rank"] > 0 and record["wall_time_s"] > 0
				and solver["trial_solves"] >= 1 and math.isfinite(record["inlet"]["mean_pressure_pa"])):
			raise RuntimeError("summary JSON solver or resource fields are incomplete")


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--ranks", type=int, required=True)
	parser.add_argument("--binary", type=Path, required=True)
	parser.add_argument("--case", type=Path, required=True)
	args = parser.parse_args()
	binary = args.binary.resolve(strict=True)
	reference = args.case.resolve(strict=True)
	ranks = args.ranks
	with tempfile.TemporaryDirectory(prefix="iga-native-hydraulic-control-") as directory:
		root = Path(directory)
		multi = json.loads(channel_case(root, reference, split_outlet=True).read_text())
		multi["motion"] = {"kind": "fixed", "speed_x_m_s": 0.0}
		multi["time"]["steps"] = 3
		multi["fluid"]["nonlinear_tolerance"] = 1e-12

		# 1. Predictor: one coupling iteration, same state as the converged fixed point.
		fixed = dict(multi, coupling={"method": "fixed", "maximum_iterations": 80,
			"pressure_relative_tolerance": 1e-11, "flow_relative_tolerance": 1e-12,
			"relaxation_factor": 0.5})
		fixed_output = run(binary, ranks, write_case(root/"fixed.json", fixed))
		predicted = json.loads(json.dumps(fixed))
		predicted["coupling"]["predictor"] = "implicit_0d"
		predicted_case = write_case(root/"predicted.json", predicted)
		predicted_output = run(binary, ranks, predicted_case, "--summary-json", root/"predicted_summary.json")
		for step in (1, 2, 3):
			reference_values = step_values(fixed_output, step)
			values = step_values(predicted_output, step)
			if reference_values["iterations"] < 3 or values["iterations"] != 1:
				raise RuntimeError("predictor did not replace the outer fixed-point iterations")
			for key, value in reference_values.items():
				if key != "iterations" and not close(values[key], value, 1e-8):
					raise RuntimeError(f"predictor {key} differs: {values[key]} vs {value}")
			work = solver_values(predicted_output, step)
			if not (work["predictor"] == 1 and work["coupling_iterations"] == 1
					and work["newton_iterations"] == 0):
				raise RuntimeError("predictor seed was not accepted without extra Newton work")
		summary_checks(json.loads((root/"predicted_summary.json").read_text()), predicted_output,
			[2, 4], 3)
		run(binary, ranks, predicted_case, "--summary-json", root/"x.json", "--check-input",
			expect_success=False)

		# 2. Steady momentum with the predictor.
		steady = json.loads(json.dumps(predicted))
		steady["fluid"]["steady"] = True
		steady_output = run(binary, ranks, write_case(root/"steady.json", steady))
		if step_values(steady_output, 3)["iterations"] != 1:
			raise RuntimeError("steady predictor run needed outer iterations")
		moving = json.loads(json.dumps(steady))
		moving["motion"] = {"kind": "rigid_translation_x", "speed_x_m_s": 0.001}
		run(binary, ranks, write_case(root/"steady_moving.json", moving), expect_success=False)

		# 3. Controlled restart with declared mutable boundary values.
		controlled = json.loads(json.dumps(predicted))
		controlled["boundary_updates"] = {"mutable": ["source.prescribed_flow_m3_s",
			"terminal_rcrs[].proximal_resistance_pa_s_m3",
			"terminal_rcrs[].distal_resistance_pa_s_m3"]}
		base_case = write_case(root/"controlled.json", controlled)
		uninterrupted = run(binary, ranks, base_case)
		run(binary, ranks, base_case, "--checkpoint-dir", root/"ck1", "--stop-after-step", 1)
		same = run(binary, ranks, base_case, "--restart-dir", root/"ck1")
		for step in (2, 3):
			for key, value in step_values(uninterrupted, step).items():
				if not close(step_values(same, step)[key], value, 1e-10):
					raise RuntimeError("controlled restart with unchanged values differs")
		changed = json.loads(json.dumps(controlled))
		changed["source"]["prescribed_flow_m3_s"] *= 1.5
		changed["terminal_rcrs"][1]["distal_resistance_pa_s_m3"] *= 2.0
		changed["terminal_rcrs"][0]["proximal_resistance_pa_s_m3"] *= 3.0
		changed_case = write_case(root/"changed.json", changed)
		rejected = run(binary, ranks, changed_case, "--restart-dir", root/"ck1", expect_success=False)
		if "--accept-boundary-update" not in rejected or "native_hydraulic_step" in rejected:
			raise RuntimeError("changed boundary values were not rejected without the flag")
		updated_output = run(binary, ranks, changed_case, "--restart-dir", root/"ck1",
			"--accept-boundary-update", "--checkpoint-dir", root/"ck2",
			"--summary-json", root/"updated_summary.json")
		summary = json.loads((root/"updated_summary.json").read_text())
		summary_checks(summary, updated_output, [2, 4], 2)
		if (summary["restart"]["restored_accepted_steps"] != 1
				or sorted(summary["restart"]["boundary_updated_domains"]) != ["source", "terminal_2", "terminal_4"]
				or not summary["steps"][0]["boundary_values"]["updated_at_this_step"]
				or summary["steps"][1]["boundary_values"]["updated_at_this_step"]):
			raise RuntimeError("summary JSON does not record the boundary update")
		history = summary["boundary_history"]
		if ([(item["first_step"], item["last_step"]) for item in history] != [(1, 1), (2, 3)]
				or history[0]["values_sha256"] == history[1]["values_sha256"]
				or history[1]["values_sha256"] != summary["boundary_values_sha256"]):
			raise RuntimeError(f"checkpoint boundary-value history is wrong: {history}")
		for record in summary["steps"]:
			# The terminal port pressure uses the new proximal resistance.
			for item in record["outlets"]:
				resistance = next(value["proximal_resistance_pa_s_m3"]
					for value in record["boundary_values"]["terminals"] if value["label"] == item["label"])
				expected = item["terminal_stored_pressure_pa"]+resistance*item["outward_flow_m3_s"]
				if not close(item["terminal_port_pressure_pa"], expected, 1e-10):
					raise RuntimeError("updated terminal does not use the new proximal resistance")
			if not close(record["boundary_values"]["source_prescribed_flow_m3_s"], 0.015, 1e-15):
				raise RuntimeError("updated source flow is not recorded")
		# The source capacitance stores part of the step change, so only require
		# a clearly larger inflow (outward inlet flow is negative).
		if not abs(step_values(updated_output, 2)["inlet_m3_s"]) > 1.05*abs(step_values(uninterrupted, 2)["inlet_m3_s"]):
			raise RuntimeError("updated source flow did not change the inflow")
		# The updated checkpoint matches the updated case, so no flag is needed;
		# the step horizon stays part of the identity (use --stop-after-step).
		resumed = run(binary, ranks, changed_case, "--restart-dir", root/"ck2", "--stop-after-step", 3)
		if "native_hydraulic_step" in resumed:
			raise RuntimeError("restart at the stop step published a step")
		continued = json.loads(json.dumps(changed))
		continued["time"]["steps"] = 4
		run(binary, ranks, write_case(root/"continued.json", continued), "--restart-dir", root/"ck2",
			expect_success=False)
		for name, mutate in (
				("capacitance", lambda case: case["terminal_rcrs"][0].update(capacitance_m3_pa=0.02)),
				("viscosity", lambda case: case["fluid"].update(dynamic_viscosity_pa_s=0.2)),
				("undeclared", lambda case: case["terminal_rcrs"][0].update(distal_pressure_pa=1.0)),
				("policy", lambda case: case["boundary_updates"]["mutable"].pop()),
				("labels", lambda case: case["terminal_rcrs"].reverse())):
			bad = json.loads(json.dumps(changed))
			mutate(bad)
			output = run(binary, ranks, write_case(root/f"bad_{name}.json", bad), "--restart-dir",
				root/"ck1", "--accept-boundary-update", expect_success=False)
			if "native_hydraulic_step" in output:
				raise RuntimeError(f"{name} change published a step")
		plain = json.loads(json.dumps(predicted))
		plain["source"]["prescribed_flow_m3_s"] *= 1.5
		run(binary, ranks, write_case(root/"plain.json", plain), "--restart-dir", root/"ck1",
			"--accept-boundary-update", expect_success=False)
		bad_path = json.loads(json.dumps(controlled))
		bad_path["boundary_updates"]["mutable"].append("fluid.dynamic_viscosity_pa_s")
		run(binary, ranks, write_case(root/"bad_path.json", bad_path), "--check-input",
			expect_success=False)
		run(binary, ranks, base_case, "--accept-boundary-update", expect_success=False)

		# 4. More than 16 terminal RCRs, including checkpoint/restart.
		grid_case, outlets = outlet_grid_case(root, reference, 5)
		grid_output = run(binary, ranks, grid_case, "--summary-json", root/"grid_summary.json")
		summary = json.loads((root/"grid_summary.json").read_text())
		summary_checks(summary, grid_output, outlets, 2)
		last = summary["steps"][-1]
		if len(outlets) != 25 or last["coupling_iterations"] != 1 or not all(
				item["outward_flow_m3_s"] > 0 for item in last["outlets"]):
			raise RuntimeError("25-outlet grid did not solve with the predictor")
		run(binary, ranks, grid_case, "--checkpoint-dir", root/"grid_ck", "--stop-after-step", 1)
		grid_resumed = run(binary, ranks, grid_case, "--restart-dir", root/"grid_ck")
		for key, value in step_values(grid_output, 2).items():
			if not close(step_values(grid_resumed, 2)[key], value, 1e-10):
				raise RuntimeError("25-outlet restart differs")
		# 5. A wall-locked P2/P1 pressure vertex is rejected before any solve.
		locked = run(binary, ranks, wall_locked_case(root, reference), "--check-input",
			expect_success=False)
		if "pressure vertices whose tetrahedra have every P2 velocity node" not in locked:
			raise RuntimeError("wall-locked pressure vertex was not diagnosed")
	print(f"native hydraulic predictor, steady, summary, controlled restart, 25-outlet "
		f"and wall-locked check: PASS ranks={ranks}")


if __name__ == "__main__":
	main()
