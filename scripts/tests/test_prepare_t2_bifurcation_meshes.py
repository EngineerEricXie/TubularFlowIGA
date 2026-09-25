#!/usr/bin/env python3
"""Exercise the reproducible coupled Y-bifurcation mesh preparation path."""

import json
from pathlib import Path
import subprocess
import tempfile


def main():
	root = Path(__file__).resolve().parents[2]
	with tempfile.TemporaryDirectory(prefix="tubularflow-y-series-") as directory:
		completed = subprocess.run(["python3",
			str(root/"scripts"/"prepare_t2_bifurcation_meshes.py"), directory,
			"--maximum-level", "2"], text=True, capture_output=True)
		assert completed.returncode == 0, completed.stdout+completed.stderr
		data = json.loads((Path(directory)/"mesh-series.json").read_text(encoding="utf-8"))
		assert data["kind"] == "t2_bifurcation_coupled_mesh_preparation_not_solver_evidence"
		assert len(data["levels"]) == 2
		coarse, fine = data["levels"]
		assert coarse["surface_triangles"] < fine["surface_triangles"]
		assert coarse["tetrahedra"] < fine["tetrahedra"]
		assert set(coarse["boundary_triangles"]) == {"0", "1", "2", "3"}
		assert coarse["minimum_scaled_jacobian"] >= 1e-3
		assert fine["minimum_scaled_jacobian"] >= 1e-3
	print("prepare_t2_bifurcation_meshes_test: PASS")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
