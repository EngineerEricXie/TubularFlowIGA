#!/usr/bin/env python3
"""Exercise circular-pipe generation through the shared C++ surface preflight."""

import json
import math
from pathlib import Path
import subprocess
import tempfile


def main():
	root = Path(__file__).resolve().parents[2]
	generator = root/"scripts"/"generate_circular_pipe_surface.py"
	preflight = root/"solvers"/"cpu"/"surface_fem_preflight"
	with tempfile.TemporaryDirectory(prefix="tubularflow-circular-pipe-") as directory:
		temporary = Path(directory)
		surface = temporary/"pipe.vtp"
		completed = subprocess.run(["python3", str(generator), str(surface), "--length-m",
			"0.1", "--radius-m", "0.005", "--circumferential-segments", "16",
			"--axial-segments", "4"], text=True, capture_output=True)
		assert completed.returncode == 0, completed.stderr
		manifest = temporary/"preflight.json"
		completed = subprocess.run([str(preflight), str(surface), str(temporary/"surface.msh"),
			"--manifest", str(manifest)], text=True, capture_output=True)
		assert completed.returncode == 0, completed.stderr
		data = json.loads(manifest.read_text(encoding="utf-8"))
		assert data["checks"] == {"closed": True, "oriented": True, "manifold": True,
			"connected": True, "self_intersection_free": True, "positive_volume": True}
		assert {label: item["triangles"] for label, item in data["boundary_labels"].items()} \
			== {"0": 128, "1": 16, "2": 16}
		expected_volume = 0.5*16*0.005**2*math.sin(2*math.pi/16)*0.1
		assert math.isclose(data["volume_m3"], expected_volume, rel_tol=2e-14)
		bad = subprocess.run(["python3", str(generator), str(temporary/"bad.vtp"),
			"--length-m", "0.1", "--radius-m", "0.005", "--circumferential-segments",
			"3", "--axial-segments", "1"], text=True, capture_output=True)
		assert bad.returncode != 0
	print("circular_pipe_surface_test: PASS")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
