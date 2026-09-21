#!/usr/bin/env python3
"""Exercise native tetra/0D file restart across two MPI process lifetimes."""

import argparse
import subprocess
import tempfile
from pathlib import Path


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--ranks", type=int, required=True)
	parser.add_argument("--binary", type=Path, required=True)
	args = parser.parse_args()
	if args.ranks < 1 or args.ranks > 16:
		raise ValueError("cross-process test ranks must be between 1 and 16")
	binary = args.binary.resolve(strict=True)
	with tempfile.TemporaryDirectory(prefix="iga-native-hydraulic-cross-") as directory:
		for phase in ("--cross-process-save", "--cross-process-resume"):
			command = ("mpiexec", "-np", str(args.ranks), str(binary), phase, directory)
			result = subprocess.run(command, capture_output=True, text=True,
				timeout=120, check=False)
			if result.returncode:
				raise RuntimeError(
					f"{phase} failed with exit {result.returncode}:\n"
					f"{result.stdout}\n{result.stderr}")
			print(result.stdout.strip())
	print(f"native tetra/0D cross-process file restart: PASS ranks={args.ranks}")


if __name__ == "__main__":
	main()
