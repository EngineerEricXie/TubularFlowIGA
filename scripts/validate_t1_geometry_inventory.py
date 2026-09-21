#!/usr/bin/env python3
"""Validate T1 route/class coverage and repository evidence paths."""

import argparse
import json
from pathlib import Path
import sys


REQUIRED_CLASSES = {"tube", "bifurcation", "chamber"}
REQUIRED_ROUTES = {"surface_to_fem_volume", "centerline_to_iga_volume",
	"surface_to_immersed_background"}
CASE_KEYS = {"id", "geometry_class", "route", "inputs", "implementation",
	"evidence", "local_command", "gates", "limitations"}


def error(message):
	print(f"T1 geometry inventory error: {message}", file=sys.stderr)
	return 1


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--inventory", default="benchmarks/t1_geometry_inventory.json")
	parser.add_argument("--root")
	args = parser.parse_args()
	path = Path(args.inventory).resolve()
	root = Path(args.root).resolve() if args.root else path.parent.parent
	try:
		data = json.loads(path.read_text(encoding="utf-8"))
	except (OSError, json.JSONDecodeError) as failure:
		return error(f"cannot read inventory: {failure}")
	if data.get("schema_version") != 1 or not isinstance(data.get("cases"), list):
		return error("schema_version must be 1 and cases must be an array")
	classes = set()
	routes = set()
	identifiers = set()
	for index, case in enumerate(data["cases"]):
		if not isinstance(case, dict) or set(case) != CASE_KEYS:
			return error(f"cases[{index}] has missing or unknown keys")
		identifier = case["id"]
		if not isinstance(identifier, str) or not identifier or identifier in identifiers:
			return error(f"cases[{index}].id is empty or duplicated")
		identifiers.add(identifier)
		classes.add(case["geometry_class"])
		routes.add(case["route"])
		for field in ("inputs", "implementation", "evidence", "gates"):
			if not isinstance(case[field], list) or not case[field]:
				return error(f"{identifier}: {field} must be a nonempty array")
		for field in ("inputs", "implementation", "evidence"):
			for relative in case[field]:
				if not isinstance(relative, str) or Path(relative).is_absolute() \
						or not (root/relative).exists():
					return error(f"{identifier}: missing or invalid path {relative!r}")
		for field in ("local_command", "limitations"):
			if not isinstance(case[field], str) or not case[field].strip():
				return error(f"{identifier}: {field} must be a nonempty string")
	if classes != REQUIRED_CLASSES:
		return error(f"geometry classes are {sorted(classes)}, expected {sorted(REQUIRED_CLASSES)}")
	if routes != REQUIRED_ROUTES:
		return error(f"routes are {sorted(routes)}, expected {sorted(REQUIRED_ROUTES)}")
	print(f"T1 geometry inventory: PASS ({len(identifiers)} cases, three routes)")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
