#!/usr/bin/env python3
"""Hash every regular file below a directory for immutable-input checks."""

import argparse
import json
from pathlib import Path

from hpc_inventory import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    if not root.is_dir() or args.output.exists():
        parser.error("root must be a directory and output must not exist")
    values = {str(path.resolve()): digest(path)
              for path in sorted(root.rglob("*")) if path.is_file()}
    if not values:
        parser.error("tree contains no regular files")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(values, indent=2) + "\n")


if __name__ == "__main__":
    main()
