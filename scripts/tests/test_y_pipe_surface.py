#!/usr/bin/env python3
"""Validate the reproducible labelled Y-pipe comparison surface."""

import json
from pathlib import Path
import subprocess
import tempfile


def main():
	root = Path(__file__).resolve().parents[2]
	generator = root/"scripts"/"generate_y_pipe_surface.py"
	preflight = root/"solvers"/"cpu"/"surface_fem_preflight"
	with tempfile.TemporaryDirectory(prefix="tubularflow-y-pipe-") as directory:
		directory = Path(directory); surface = directory/"y.vtp"
		generated = subprocess.run(["python3", str(generator), str(surface)],
			text=True, capture_output=True)
		assert generated.returncode == 0, generated.stdout+generated.stderr
		manifest = directory/"surface.json"
		checked = subprocess.run([str(preflight), str(surface), str(directory/"surface.msh"),
			"--manifest", str(manifest)], text=True, capture_output=True)
		assert checked.returncode == 0, checked.stdout+checked.stderr
		data = json.loads(manifest.read_text(encoding="utf-8"))
		assert data["checks"]["closed"] and data["checks"]["self_intersection_free"]
		counts = {label: value["triangles"] for label, value in data["boundary_labels"].items()}
		assert set(counts) == {"0","1","2","3"} and all(value > 0 for value in counts.values())
	print("y_pipe_surface_test: PASS")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
