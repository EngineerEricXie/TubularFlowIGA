#!/usr/bin/env python3
"""Validate the reproducible labelled bent-pipe comparison surface."""

import json
from pathlib import Path
import subprocess
import tempfile


def main():
	root = Path(__file__).resolve().parents[2]
	generator = root/"scripts"/"generate_bent_pipe_surface.py"
	preflight = root/"solvers"/"cpu"/"surface_fem_preflight"
	with tempfile.TemporaryDirectory(prefix="tubularflow-bent-pipe-") as directory:
		directory = Path(directory); surface = directory/"bent.vtp"
		generated = subprocess.run(["python3", str(generator), str(surface),
			"--bend-radius-m", "0.04", "--tube-radius-m", "0.005",
			"--bend-angle-deg", "90", "--circumferential-segments", "16",
			"--axial-segments", "12", "--radial-segments", "2"],
			text=True, capture_output=True)
		assert generated.returncode == 0, generated.stdout+generated.stderr
		manifest = directory/"surface.json"
		checked = subprocess.run([str(preflight), str(surface), str(directory/"surface.msh"),
			"--manifest", str(manifest)], text=True, capture_output=True)
		assert checked.returncode == 0, checked.stdout+checked.stderr
		data = json.loads(manifest.read_text(encoding="utf-8"))
		assert data["checks"]["closed"] and data["checks"]["self_intersection_free"]
		assert {label: value["triangles"] for label, value in data["boundary_labels"].items()} \
			== {"0": 384, "1": 48, "2": 48}
	print("bent_pipe_surface_test: PASS")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
