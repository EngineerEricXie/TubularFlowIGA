#!/usr/bin/env python3
"""Validate the native-FEM T2 Poiseuille card and analytic values."""

import argparse
import json
import math
from pathlib import Path
import sys


def fail(message):
	print(f"T2 fixed-flow contract error: {message}", file=sys.stderr)
	return 1


def close(actual, expected):
	return math.isfinite(actual) and math.isclose(actual, expected, rel_tol=2e-14, abs_tol=1e-15)


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--contract", default="benchmarks/t2_fixed_flow_contract.json")
	args = parser.parse_args()
	try:
		data = json.loads(Path(args.contract).read_text(encoding="utf-8"))
	except (OSError, json.JSONDecodeError) as error:
		return fail(f"cannot read contract: {error}")
	if data.get("schema_version") != 1 or data.get("case_id") != "steady-circular-pipe-poiseuille":
		return fail("unexpected schema or case id")
	physics = data["physics"]
	geometry = data["geometry"]
	boundary = data["boundary_conditions"]
	reference = data["analytic_reference"]
	gates = data["qoi_gates"]
	radius = float(geometry["radius_m"])
	length = float(geometry["length_m"])
	viscosity = float(physics["dynamic_viscosity_pa_s"])
	density = float(physics["density_kg_m3"])
	flow = float(boundary["prescribed_flow_m3_s"])
	if not all(math.isfinite(value) and value > 0.0
			for value in (radius, length, viscosity, density, flow)):
		return fail("physical parameters must be finite and positive")
	if geometry.get("boundary_labels") != {"wall": 0, "inlet": 1, "outlet": 2}:
		return fail("canonical boundary labels changed")
	if boundary.get("flow_sign") != "outward_positive" \
			or boundary.get("expected_inlet_outward_flow_m3_s") != -flow \
			or boundary.get("expected_outlet_outward_flow_m3_s") != flow:
		return fail("outward-positive flow convention is inconsistent")
	area = math.pi*radius**2
	expected = {
		"cross_section_area_m2": area,
		"volume_m3": area*length,
		"mean_velocity_m_s": flow/area,
		"centerline_velocity_m_s": 2.0*flow/area,
		"pressure_drop_pa": 8.0*viscosity*length*flow/(math.pi*radius**4),
		"wall_shear_stress_pa": 4.0*viscosity*flow/(math.pi*radius**3),
		"diameter_reynolds_number": 2.0*density*flow/(math.pi*viscosity*radius)}
	for name, value in expected.items():
		if not close(float(reference[name]), value):
			return fail(f"analytic {name} is {reference[name]!r}, expected {value!r}")
	levels = data["discretization_contract"]["mesh_levels_target_size_m"]
	if len(levels) < 3 or any(not (levels[index] > levels[index+1] > 0.0)
			for index in range(len(levels)-1)):
		return fail("mesh convergence requires at least three strictly refined levels")
	if gates.get("minimum_mesh_levels") != len(levels) \
			or gates.get("absolute_qoi_gate_level") != "finest":
		return fail("three levels and a finest-level absolute QoI gate are required")
	if gates.get("maximum_relative_mass_imbalance", 1.0) > 1e-8:
		return fail("mass gate was loosened beyond the frozen limit")
	required = data.get("required_result_fields", [])
	for name in ("source_tree_state", "source_files_sha256", "linear_solver", "velocity_space",
			"pressure_space", "stabilization", "geometry_error",
			"assembly_wall_s", "solve_wall_s", "peak_rss_bytes"):
		if name not in required:
			return fail(f"required result field {name!r} is missing")
	decision = data.get("backend_decision", {})
	if decision != {"architecture": "native",
			"selected_backend": "native_cpp_petsc_tetrahedral_fem",
			"allowed_next_choices": ["native_cpp_petsc_tetrahedral_fem"],
			"selection_required_before_solver_adapter": False,
			"external_discretization_frameworks_allowed": False}:
		return fail("native backend decision changed or permits an external FEM framework")
	print("T2 fixed-flow contract: PASS (native C++/PETSc FEM card)")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
