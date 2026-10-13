#!/usr/bin/env python3
"""Compare the cuDSS and iterative linear solvers of native_tet_flow_cuda.

Channel flow on a small box mesh runs with fluid.linear_solver "direct" and
"iterative"; the final velocities must agree to 1e-10 relative L2.
"""
import argparse
import array
import json
import math
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/"scripts"))
from generate_tet_box_mesh import box_mesh  # noqa: E402
from generate_fsi_channel_fixture import write_msh  # noqa: E402


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("root", type=Path, help="new output directory")
	args = parser.parse_args()
	root = args.root.resolve()
	root.mkdir(parents=True)
	points, cells, triangles = box_mesh(4, 0.01)
	write_msh(root/"box.msh", points, cells, triangles)
	velocities = {}
	for solver in ("direct", "iterative"):
		case = {"schema_version": 1, "mesh_file": "box.msh",
			"fluid": {"density_kg_m3": 1000.0, "dynamic_viscosity_pa_s": 0.004,
				"initial_velocity_m_s": [0.0, 0.0, 0.0], "nonlinear_tolerance": 1e-8,
				"maximum_iterations": 12, "linear_solver": solver},
			"boundaries": {"velocity_by_label_m_s": {"1": [0.01, 0.0, 0.0], "3": [0.0, 0.0, 0.0]},
				"pressure_by_label_pa": {"2": 0.0}},
			"time": {"dt_s": 0.05, "steps": 2}}
		(root/f"{solver}.json").write_text(json.dumps(case))
		subprocess.run([str(HERE.parent/"native_tet_flow_cuda"), str(root/f"{solver}.json"), str(root/solver)],
			check=True, stdout=subprocess.DEVNULL)
		values = array.array("d")
		values.frombytes((root/solver/"flow_velocity_step_2.bin").read_bytes())
		velocities[solver] = values
	direct, iterative = velocities["direct"], velocities["iterative"]
	difference = math.sqrt(sum((a-b)**2 for a, b in zip(direct, iterative))/sum(a*a for a in direct))
	if not difference < 1e-10:
		raise SystemExit(f"iterative and direct velocities differ: relative L2 {difference:.3e}")
	print(f"flow linear solvers agree: relative L2 {difference:.3e}")


if __name__ == "__main__":
	main()
