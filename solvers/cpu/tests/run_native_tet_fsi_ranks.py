#!/usr/bin/env python3
"""Check native_tet_fsi_steps across solvers and MPI rank counts.

On the generated FSI channel fixture, the PETSc solvers on one and on two
ranks must give the same step ledger to round-off (every rank holds the whole
problem), and with tight Newton tolerances they must agree with the dense
one-rank solvers to 1e-6 relative.
"""
import argparse
import re
import shlex
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TIGHT = ["-fsi_fluid_tolerance", "1e-11", "-fsi_solid_rtol", "1e-11"]


def ledger(text):
	steps = [line for line in text.splitlines() if line.startswith("native_fsi_step ")]
	if not steps:
		raise SystemExit("no native_fsi_step ledger lines")
	values = []
	for line in steps:
		for key, value in re.findall(r"(\w+)=(\S+)", line):
			if key in ("step", "iterations"):
				values.append((key, int(value)))
			else:
				values.append((key, float(value)))
	return values


def compare(name, reference, other, tolerance):
	if [key for key, _ in reference] != [key for key, _ in other]:
		raise SystemExit(f"{name}: ledgers have different fields")
	for (key, a), (_, b) in zip(reference, other):
		if isinstance(a, int):
			if a != b:
				raise SystemExit(f"{name}: {key} differs ({a} vs {b})")
		elif abs(a-b) > tolerance*max(abs(a), abs(b), 1e-20) and abs(a-b) > 1e-18:
			raise SystemExit(f"{name}: {key} differs ({a!r} vs {b!r})")


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("root", type=Path, help="new output directory")
	parser.add_argument("--launcher", default="mpiexec --oversubscribe")
	args = parser.parse_args()
	root = args.root.resolve()
	root.mkdir(parents=True)
	subprocess.run([sys.executable, str(ROOT/"scripts/generate_fsi_channel_fixture.py"), str(root)],
		check=True, stdout=subprocess.DEVNULL)
	binary = str(HERE.parent/"native_tet_fsi_steps")
	launcher = shlex.split(args.launcher)

	def run(name, ranks, options):
		command = [*launcher, "-np", str(ranks), binary, str(root/"case.json"), str(root/name), *options]
		result = subprocess.run(command, capture_output=True, text=True, cwd=root)
		if result.returncode:
			raise SystemExit(f"{name} failed:\n{result.stdout[-2000:]}{result.stderr[-2000:]}")
		return ledger(result.stdout)

	one = run("petsc-1", 1, TIGHT)
	two = run("petsc-2", 2, TIGHT)
	dense = run("dense", 1, TIGHT+["-fsi_dense_solvers"])
	compare("PETSc 1 vs 2 ranks", one, two, 1e-10)
	compare("dense vs PETSc", dense, one, 1e-6)
	print(f"native FSI ranks: PETSc 1/2 ranks and dense agree over {len(one)} ledger values")


if __name__ == "__main__":
	main()
