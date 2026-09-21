#!/usr/bin/env python3
"""Generate a real Y surface and solve its two native RCR outlets."""

import argparse
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

from test_native_tet_hydraulic_cli import step_values, vtk_snapshot


def run(command):
	result = subprocess.run(command, capture_output=True, text=True,
		timeout=180, check=False)
	if result.returncode:
		raise RuntimeError(f"Y case failed: {' '.join(map(str, command))}\n"
			f"{result.stdout}\n{result.stderr}")
	return result.stdout


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--solver", type=Path, required=True)
	parser.add_argument("--ranks", type=int, default=2)
	parser.add_argument("--coupling", choices=("explicit", "fixed"),
		default="explicit")
	parser.add_argument("--ftetwild", type=Path,
		help="also test fTetWild in addition to the Gmsh fallback")
	args = parser.parse_args()
	if not 1 <= args.ranks <= 16:
		raise ValueError("Y case ranks must be between 1 and 16")
	root = Path(__file__).resolve().parents[1]
	solver = args.solver.resolve(strict=True)
	with tempfile.TemporaryDirectory(prefix="iga-native-y-rcr-") as temporary:
		work = Path(temporary)
		surface = work/"y.vtp"
		run([sys.executable, str(root/"scripts/generate_y_pipe_surface.py"),
			str(surface), "--radius-m", "0.005", "--trunk-length-m", "0.025",
			"--branch-x-m", "0.02", "--branch-y-m", "0.015",
			"--target-size-m", "0.006"])
		case = json.loads((root/"solvers/cpu/tests/data/native_tet_hydraulic_star.json").read_text())
		case["schema_version"] = 2
		case["boundaries"] = {"inlet_label": 1, "wall_labels": [0]}
		case["fluid"].update({"density_kg_m3": 1.0,
			"dynamic_viscosity_pa_s": 0.001, "initial_velocity_x_m_s": 0.01,
			"nonlinear_tolerance": 1e-8, "maximum_iterations": 20})
		case["source"] = {"capacitance_m3_pa": 1e-11,
			"resistance_pa_s_m3": 1e8, "prescribed_flow_m3_s": 1e-6,
			"initial_pressure_pa": 0.0}
		terminal = {"proximal_resistance_pa_s_m3": 1000.0,
			"distal_resistance_pa_s_m3": 10000.0, "capacitance_m3_pa": 1e-6,
			"distal_pressure_pa": 0.0, "initial_pressure_pa": 0.0}
		case.pop("terminal_rcr")
		case["terminal_rcrs"] = [dict(terminal, boundary_label=2),
			dict(terminal, boundary_label=3,
				proximal_resistance_pa_s_m3=2000.0,
				distal_resistance_pa_s_m3=20000.0)]
		case["motion"] = {"kind": "fixed", "speed_x_m_s": 0.0}
		case["coupling"]["method"] = args.coupling
		case["coupling"]["maximum_iterations"] = 1 if args.coupling == "explicit" else 20
		if args.coupling == "fixed":
			case["coupling"]["pressure_relative_tolerance"] = 1e-4
			case["coupling"]["flow_relative_tolerance"] = 1e-6
		case["time"] = {"dt_s": 0.01, "steps": 2}
		(work/"y_case.json").write_text(json.dumps(case, indent=2)+"\n")
		meshers = [("gmsh", {"kind": "gmsh", "target_size_m": 0.007})]
		if args.ftetwild:
			meshers.append(("ftetwild", {"kind": "ftetwild",
				"target_size_m": 0.007, "envelope_m": 0.0001,
				"executable": str(args.ftetwild.resolve(strict=True))}))
		for name, mesher in meshers:
			configuration = {"schema_version": 1,
				"backend": "native_tet_p2p1_ale_hydraulic",
				"input_route": "surface", "input_file": surface.name,
				"case_file": "y_case.json", "output_directory": name+"_result",
				"mpi_ranks": args.ranks, "mesher": mesher}
			config_path = work/(name+"_workflow.json")
			config_path.write_text(json.dumps(configuration, indent=2)+"\n")
			result = work/(name+"_result")
			try:
				run([sys.executable, str(root/"scripts/run_native_tet_workflow.py"),
					str(config_path), "--solver", str(solver)])
			except RuntimeError as error:
				log = result/"simulation.log"
				detail = log.read_text() if log.is_file() else "no simulation log"
				raise RuntimeError(f"{name} {args.coupling} Y run failed: {detail}") from error
			status = json.loads((result/"workflow_status.json").read_text())
			mesh = json.loads((result/"mesh_manifest.json").read_text())
			volume = mesh["volume_mesh"]
			if (status["state"] != "passed" or status["accepted_steps"] != 2
					or not mesh["gates"]["passed"]
					or set(volume["boundary_labels"]) != {"0", "1", "2", "3"}):
				raise RuntimeError(f"{name} Y geometry or run did not pass")
			output = (result/"simulation.log").read_text()
			for step_number in (1, 2):
				step = step_values(output, step_number)
				if args.coupling == "fixed" and step["iterations"] < 2:
					raise RuntimeError(f"{name} Y fixed-point did not iterate: {step}")
				inlet, upper, lower = (step["inlet_m3_s"],
					step["outlet_2_m3_s"], step["outlet_3_m3_s"])
				scale = abs(inlet)+abs(upper)+abs(lower)
				if (not all(math.isfinite(value) for value in step.values())
						or inlet >= -1e-9 or upper <= 1e-9 or lower <= 1e-9
						or abs(step["terminal_2_pa"]-step["terminal_3_pa"])
							<= 0.05*max(abs(step["terminal_2_pa"]),
								abs(step["terminal_3_pa"]))
						or abs(inlet+upper+lower) > 1e-4*scale):
					raise RuntimeError(f"{name} Y outlet flow or mass balance differs: {step}")
			vtk_snapshot(result/"fields", 2, expected_points=None,
				expected_cells=volume["elements"], expected_first_cell=None,
				expected_displacement_x=0.0)
			final = step_values(output, 2)
			print(f"native_tet_y_rcr: PASS mesher={name} coupling={args.coupling} ranks={args.ranks} "
				f"tetrahedra={volume['elements']} "
				f"iterations={final['iterations']:.0f} "
				f"inlet_m3_s={final['inlet_m3_s']:.6g} "
				f"outlet_2_m3_s={final['outlet_2_m3_s']:.6g} "
				f"outlet_3_m3_s={final['outlet_3_m3_s']:.6g}")


if __name__ == "__main__":
	main()
