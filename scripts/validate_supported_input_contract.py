#!/usr/bin/env python3
"""Fail-closed validation for the supported input/unit/label contract."""

import argparse
import json
from pathlib import Path
import sys


ROUTES = {
	"centerline_to_iga_volume": {".swc", ".obj"},
	"surface_to_fem_volume": {".vtp", ".stl"},
	"surface_to_immersed_background": {".vtp"},
	"network_1d": {".swc", ".obj"},
}
RUNTIMES = {"native_tetrahedral_fem", "packed_body_fitted_iga"}


def fail(message):
	print(f"supported input contract error: {message}", file=sys.stderr)
	return 1


def nonempty_strings(value):
	return isinstance(value, list) and value and all(isinstance(x, str) and x.strip() for x in value)


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--contract", default="benchmarks/supported_input_contract.json")
	parser.add_argument("--root")
	args = parser.parse_args()
	path = Path(args.contract).resolve()
	root = Path(args.root).resolve() if args.root else path.parent.parent
	try:
		data = json.loads(path.read_text(encoding="utf-8"))
	except (OSError, json.JSONDecodeError) as error:
		return fail(f"cannot read contract: {error}")
	if data.get("schema_version") != 1:
		return fail("schema_version must be 1")
	global_units = data.get("global_conventions", {})
	if global_units.get("solver_length_unit") != "m" \
			or global_units.get("label_semantics_must_be_explicit") is not True \
			or global_units.get("organ_name_inference_for_route_or_labels") is not False:
		return fail("global SI and explicit-label policy changed")
	routes = data.get("routes")
	if not isinstance(routes, list) or {x.get("id") for x in routes if isinstance(x, dict)} != set(ROUTES):
		return fail("route catalog is incomplete or contains an unknown route")
	for route in routes:
		identifier = route["id"]
		if route.get("status") != "implemented" or not nonempty_strings(route.get("minimum_case_inputs")):
			return fail(f"{identifier}: status or minimum input set is incomplete")
		geometry = route.get("accepted_geometry")
		if not isinstance(geometry, list) or {x.get("extension") for x in geometry if isinstance(x, dict)} != ROUTES[identifier]:
			return fail(f"{identifier}: accepted extensions changed")
		for item in geometry:
			if set(item) != {"extension", "kind", "contract"} \
					or not all(isinstance(item[k], str) and item[k].strip() for k in item):
				return fail(f"{identifier}: malformed geometry entry")
		for field in ("source_units", "labels"):
			if not isinstance(route.get(field), str) or not route[field].strip():
				return fail(f"{identifier}: missing {field}")
		if not nonempty_strings(route.get("limitations")) or not nonempty_strings(route.get("implementation_evidence")):
			return fail(f"{identifier}: missing limitations or evidence")
		for relative in route["implementation_evidence"]:
			if Path(relative).is_absolute() or not (root/relative).is_file():
				return fail(f"{identifier}: evidence path does not exist: {relative}")
	fem = next(x for x in routes if x["id"] == "surface_to_fem_volume")
	if set(fem.get("external_ftetwild_cli_formats_not_project_inputs", [])) != {".off", ".obj", ".stl", ".ply"}:
		return fail("fTetWild upstream formats are not explicitly separated from project inputs")
	runtimes = data.get("runtime_interfaces")
	if not isinstance(runtimes, list) or {x.get("id") for x in runtimes if isinstance(x, dict)} != RUNTIMES:
		return fail("runtime interface catalog is incomplete")
	for runtime in runtimes:
		identifier = runtime["id"]
		if not nonempty_strings(runtime.get("minimum_inputs")) \
				or not nonempty_strings(runtime.get("implementation_evidence")):
			return fail(f"{identifier}: minimum input set or evidence is incomplete")
		for field in ("format", "labels", "discretization"):
			if not isinstance(runtime.get(field), str) or not runtime[field].strip():
				return fail(f"{identifier}: missing {field}")
		for relative in runtime["implementation_evidence"]:
			if Path(relative).is_absolute() or not (root/relative).is_file():
				return fail(f"{identifier}: evidence path does not exist: {relative}")
	if not nonempty_strings(data.get("explicit_rejections")):
		return fail("explicit rejection policy is empty")

	# Bind high-risk claims to the actual readers so documentation cannot drift silently.
	surface_preflight = (root/"solvers/cpu/src/surface_fem_preflight.cpp").read_text(encoding="utf-8")
	for token in ('extension==".vtp"', 'extension==".stl"',
			'surface input extension must be .vtp or .stl'):
		if token not in surface_preflight:
			return fail(f"surface preflight no longer proves {token!r}")
	immersed = (root/"include/ImmersedFlowCase.hpp").read_text(encoding="utf-8")
	if "SurfaceReaders::ReadVtpPath" not in immersed or "ValidateLabelPartition" not in immersed:
		return fail("production immersed VTP/label-partition binding changed")
	native = (root/"solvers/cpu/include/NativeTetFem.hpp").read_text(encoding="utf-8")
	for token in ("ASCII Gmsh 4.1", '"boundary_label_"', 'name->second != "fluid"',
			"boundary != labelled"):
		if token not in native:
			return fail(f"native FEM input binding no longer proves {token!r}")
	mesh = (root/"preprocessing/mesh/src/MeshGenerator.cpp").read_text(encoding="utf-8")
	if "mesh.labels.at(offset+index) = 0" not in mesh or "int tip_label = 1" not in mesh:
		return fail("centerline-to-IGA label convention changed")
	print("Supported input contract: PASS (4 routes, 2 runtime interfaces)")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
