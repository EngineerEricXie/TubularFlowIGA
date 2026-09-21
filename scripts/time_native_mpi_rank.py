#!/usr/bin/env python3
"""Measure one native MPI solver rank with GNU time; never aggregate RSS."""

import argparse
import os
from pathlib import Path
import subprocess
import sys


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--metrics-dir", type=Path, required=True)
	parser.add_argument("command", nargs=argparse.REMAINDER)
	args = parser.parse_args()
	command = args.command[1:] if args.command[:1] == ["--"] else args.command
	try:
		if not command:
			raise ValueError("native MPI rank command is empty")
		rank_text = next((os.environ[name] for name in (
			"OMPI_COMM_WORLD_RANK", "PMI_RANK", "PMIX_RANK")
			if name in os.environ), None)
		if rank_text is None or not rank_text.isdecimal():
			raise ValueError("MPI rank environment is absent or malformed")
		rank = int(rank_text)
		path = args.metrics_dir/f"rank{rank}.time.txt"
		if not args.metrics_dir.is_dir() or path.exists():
			raise ValueError("rank metrics directory is absent or output already exists")
		return subprocess.run(["/usr/bin/time", "-f", "%e %M", "-o",
			str(path), *command], check=False).returncode
	except (OSError, ValueError) as error:
		print(f"native MPI rank timing rejected: {error}", file=sys.stderr)
		return 2


if __name__ == "__main__":
	sys.exit(main())
