#!/usr/bin/env python3
"""Run the FSI channel fixture as a multidomain graph case and as a plain case.

native_tet_fsi_graph (graph case: two native domains joined by one FSI
coupling) and native_tet_fsi_steps (the combined case) solve the same problem,
so their step ledgers must match to round-off on one and on two ranks.
"""
import argparse
import re
import shlex
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def ledger(text):
	steps = [line for line in text.splitlines() if line.startswith("native_fsi_step ")]
	if not steps:
		raise SystemExit("no native_fsi_step ledger lines")
	return [(key, float(value)) for line in steps for key, value in re.findall(r"(\w+)=(\S+)", line)]


def compare(name, reference, other, tolerance):
	if [key for key, _ in reference] != [key for key, _ in other]:
		raise SystemExit(f"{name}: ledgers have different fields")
	for (key, a), (_, b) in zip(reference, other):
		if abs(a-b) > tolerance*max(abs(a), abs(b), 1e-20) and abs(a-b) > 1e-18:
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
	launcher = shlex.split(args.launcher)

	def run(name, ranks, command):
		result = subprocess.run([*launcher, "-np", str(ranks), *command], capture_output=True, text=True, cwd=root)
		if result.returncode:
			raise SystemExit(f"{name} failed:\n{result.stdout[-2000:]}{result.stderr[-2000:]}")
		return ledger(result.stdout)

	steps = str(ROOT/"solvers/cpu/native_tet_fsi_steps")
	graph = str(HERE.parent/"native_tet_fsi_graph")
	plain = run("plain", 1, [steps, str(root/"case.json"), str(root/"plain")])
	one = run("graph-1", 1, [graph, "--graph-case", str(root/"graph"), "--output-dir", str(root/"graph-1")])
	two = run("graph-2", 2, [graph, "--graph-case", str(root/"graph"), "--output-dir", str(root/"graph-2")])
	compare("graph vs plain case", plain, one, 1e-10)
	compare("graph 1 vs 2 ranks", one, two, 1e-10)
	for name in ("flow.pvd", "wall.pvd", "graph_binding_manifest.json"):
		if not (root/"graph-1"/name).is_file():
			raise SystemExit(f"graph run did not write {name}")
	print(f"native FSI graph: graph and plain ledgers agree over {len(one)} values on 1 and 2 ranks")


if __name__ == "__main__":
	main()
