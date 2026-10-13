#!/usr/bin/env python3
"""Run native_graph_checkpoint_groups_test on two MPI groups and compare outputs.

Three ranks split into a one-rank and a two-rank group, each running its own
native graph. Each group runs uninterrupted, saves after step 3, and resumes;
the resumed outputs must match the uninterrupted ones byte for byte.
"""
import argparse
import filecmp
import shlex
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]/"scripts/hpc"))
from hpc_native_graph_checkpoint import fixture  # noqa: E402


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("root", type=Path, help="new output directory")
	parser.add_argument("--launcher", default="mpiexec --oversubscribe")
	args = parser.parse_args()
	binaries = HERE.parent
	root = args.root.resolve()
	root.mkdir(parents=True)
	for group, ranks in ((0, 1), (1, 2)):
		fixture(binaries/"native_graph_checkpoint_fixture_test", root/f"group-{group}"/"fixture", "flow", ranks)
	for action in ("full", "save", "resume"):
		subprocess.run([*shlex.split(args.launcher), "-np", "3",
			str(binaries/"native_graph_checkpoint_groups_test"), str(root), action], check=True)
	for group in (0, 1):
		full, resumed = root/f"group-{group}/full-output", root/f"group-{group}/resume-output"
		names = sorted(path.name for path in full.iterdir())
		_, mismatch, errors = filecmp.cmpfiles(full, resumed, names, shallow=False)
		if not names or mismatch or errors:
			raise SystemExit(f"group {group}: resumed output differs from the uninterrupted run: {mismatch+errors}")
	print(f"native graph checkpoint groups: both groups resumed identically ({len(names)} files each)")


if __name__ == "__main__":
	main()
