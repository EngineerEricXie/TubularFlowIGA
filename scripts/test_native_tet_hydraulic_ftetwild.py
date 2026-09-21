#!/usr/bin/env python3
"""Real fTetWild pipe to project-owned P2/P1 FEM and RCR smoke test."""

import argparse
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

from test_native_tet_hydraulic_cli import run, step_values, vtk_snapshot


def checked(command):
	result = subprocess.run(command, capture_output=True, text=True,
		timeout=120, check=False)
	if result.returncode:
		raise RuntimeError(f"command failed: {' '.join(map(str, command))}\n"
			f"{result.stdout}\n{result.stderr}")
	return result.stdout


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--ftetwild", type=Path, required=True)
	parser.add_argument("--binary", type=Path, required=True)
	parser.add_argument("--ranks", type=int, default=2)
	args = parser.parse_args()
	if not 1 <= args.ranks <= 16:
		raise ValueError("MPI ranks must be between 1 and 16")
	root = Path(__file__).resolve().parents[1]
	ftetwild = args.ftetwild.resolve(strict=True)
	binary = args.binary.resolve(strict=True)
	with tempfile.TemporaryDirectory(prefix="iga-ftetwild-native-flow-") as temporary:
		work = Path(temporary)
		surface, mesh = work/"pipe.vtp", work/"pipe.msh"
		manifest, case_path = work/"mesh_manifest.json", work/"case.json"
		print(checked([sys.executable, str(root/"scripts/generate_circular_pipe_surface.py"),
			str(surface), "--length-m", "0.03", "--radius-m", "0.005",
			"--circumferential-segments", "8", "--axial-segments", "2"]).strip())
		print(checked([sys.executable, str(root/"scripts/ftetwild_to_fem_volume.py"),
			str(surface), str(mesh), "--manifest", str(manifest),
			"--ftetwild", str(ftetwild), "--target-size-m", "0.007",
			"--envelope-m", "0.0001", "--stop-energy", "12",
			"--max-optimization-passes", "40", "--max-threads", "2"]).strip())
		mesh_data = json.loads(manifest.read_text())
		volume = mesh_data["volume_mesh"]
		if (not mesh_data["gates"]["passed"]
				or volume["element_type"] != "tetrahedron_p1"
				or not volume["boundary_matches_tetrahedra"]
				or set(volume["boundary_labels"]) != {"0", "1", "2"}
				or volume["elements"] < 1):
			raise RuntimeError("fTetWild mesh did not pass labelled native tetra gates")
		case = json.loads((root/"solvers/cpu/tests/data/native_tet_hydraulic_star.json").read_text())
		case["mesh_file"] = mesh.name
		case["boundaries"]["wall_labels"] = [0]
		case["fluid"].update({"density_kg_m3": 1.0,
			"dynamic_viscosity_pa_s": 0.001,
			"initial_velocity_x_m_s": -0.01,
			"nonlinear_tolerance": 1e-8, "maximum_iterations": 20})
		case["source"] = {"capacitance_m3_pa": 0.01,
			"resistance_pa_s_m3": 10.0, "prescribed_flow_m3_s": 1e-6,
			"initial_pressure_pa": 0.0}
		case["terminal_rcr"] = {"proximal_resistance_pa_s_m3": 1.0,
			"distal_resistance_pa_s_m3": 10.0, "capacitance_m3_pa": 0.01,
			"distal_pressure_pa": 0.0, "initial_pressure_pa": 0.0}
		case["motion"] = {"kind": "fixed", "speed_x_m_s": 0.0}
		case["coupling"].update({"method": "explicit", "maximum_iterations": 1,
			"pressure_relative_tolerance": 0.001,
			"flow_relative_tolerance": 0.001})
		case["time"] = {"dt_s": 0.01, "steps": 1}
		case_path.write_text(json.dumps(case, indent=2)+"\n")
		output = run(binary, args.ranks, case_path,
			"--output-dir", work/"fields")
		step = step_values(output, 1)
		if not (step["inlet_m3_s"] < -1e-9
				and step["outlet_m3_s"] > 1e-9
				and abs(step["inlet_m3_s"]+step["outlet_m3_s"])
					< 1e-4*(abs(step["inlet_m3_s"])+abs(step["outlet_m3_s"]))
				and all(math.isfinite(value) for value in step.values())):
			raise RuntimeError(f"fTetWild pipe has invalid native hydraulic flow: {step}")
		points, cells = vtk_snapshot(work/"fields", 1,
			expected_points=None, expected_cells=volume["elements"],
			expected_first_cell=None, expected_displacement_x=0.0)
		if len(points) <= volume["nodes"]:
			raise RuntimeError("native P2 field has no midside nodes")
		if not all(math.isfinite(value) for coordinate, velocity, pressure, *_ in points.values()
				for value in (*coordinate, *velocity, pressure)):
			raise RuntimeError("native P2 field contains nonfinite values")
		print(f"native_tet_hydraulic_ftetwild: PASS ranks={args.ranks} "
			f"tetrahedra={len(cells)} p2_nodes={len(points)} "
			f"minimum_scaled_jacobian={volume['minimum_scaled_jacobian']:.6g} "
			f"inlet_m3_s={step['inlet_m3_s']:.6g} "
			f"outlet_m3_s={step['outlet_m3_s']:.6g}")


if __name__ == "__main__":
	main()
