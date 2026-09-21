#!/usr/bin/env python3
"""Generate a labelled fTetWild pipe and solve native P1 Darcy flow on it."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile


def invoke(command):
	result = subprocess.run(command, capture_output=True, text=True, timeout=180)
	if result.returncode:
		raise RuntimeError(f"command failed: {command}\n{result.stdout}\n{result.stderr}")
	return result.stdout + result.stderr


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--binary", type=Path, required=True)
	parser.add_argument("--ftetwild", type=Path, required=True)
	parser.add_argument("--ranks", type=int, default=2)
	args = parser.parse_args()
	if not 1 <= args.ranks <= 16:
		raise ValueError("MPI ranks must be between 1 and 16")
	root = Path(__file__).resolve().parents[1]
	with tempfile.TemporaryDirectory(prefix="native-darcy-ftetwild-") as temporary:
		work = Path(temporary)
		surface, mesh, manifest = work/"pipe.vtp", work/"pipe.msh", work/"mesh.json"
		invoke([sys.executable, str(root/"scripts/generate_circular_pipe_surface.py"),
			str(surface), "--length-m", "0.03", "--radius-m", "0.005",
			"--circumferential-segments", "8", "--axial-segments", "2"])
		invoke([sys.executable, str(root/"scripts/ftetwild_to_fem_volume.py"),
			str(surface), str(mesh), "--manifest", str(manifest),
			"--ftetwild", str(args.ftetwild.resolve(strict=True)),
			"--target-size-m", "0.007", "--envelope-m", "0.0001",
			"--stop-energy", "12", "--max-optimization-passes", "40",
			"--max-threads", "2"])
		volume = json.loads(manifest.read_text())["volume_mesh"]
		output = invoke(["mpiexec", "-np", str(args.ranks),
			str(args.binary.resolve(strict=True)), str(mesh)])
		match = re.search(r"native_tet_darcy: PASS ranks=(\d+)", output)
		if match is None or int(match.group(1)) != args.ranks:
			raise RuntimeError("native Darcy did not solve the generated mesh")
		diagnostic = re.search(r"pipe diagnostic pressure_error=([^ ]+) "
			r"flux_error=([^ ]+) inlet=([^ ]+) outlet=([^ ]+) wall=([^\n]+)", output)
		if diagnostic is None:
			raise RuntimeError("native Darcy omitted generated-pipe diagnostics")
		print("native Darcy fTetWild diagnostic:", diagnostic.group(0))
		print(f"native_tet_darcy_ftetwild: PASS ranks={args.ranks} "
			f"tetrahedra={volume['elements']} mesh_sha256="
			f"{hashlib.sha256(mesh.read_bytes()).hexdigest()}")


if __name__ == "__main__":
	main()
