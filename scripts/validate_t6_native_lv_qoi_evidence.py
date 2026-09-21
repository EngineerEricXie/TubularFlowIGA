#!/usr/bin/env python3
"""Audit gauge-relative prescribed-LV numerical QoI without claiming physiology."""

import json
import math
from pathlib import Path

from validate_t4_native_ale_evidence import source_hash


def require(condition, message):
	if not condition:
		raise RuntimeError(message)


def close(actual, expected, absolute=1e-14, relative=1e-10):
	return math.isfinite(actual) and math.isfinite(expected) and math.isclose(
		actual, expected, abs_tol=absolute, rel_tol=relative)


def pressure_volume_integral(pressures, volumes, initial_volume):
	previous_pressure = 0.0
	previous_volume = initial_volume
	result = 0.0
	for pressure, volume in zip(pressures, volumes):
		result += 0.5*(previous_pressure+pressure)*(volume-previous_volume)
		previous_pressure, previous_volume = pressure, volume
	return result


def main():
	root = Path(__file__).resolve().parents[1]
	card = json.loads((root/"benchmarks/t6_native_lv_qoi_evidence.json").read_text())
	native = json.loads((root/"benchmarks/t4_native_ale_reference_evidence.json").read_text())
	require(card.get("schema_version") == 1
		and card.get("kind") ==
		"native_prescribed_left_ventricle_qoi_diagnostic_not_physical_validation"
		and card.get("external_fem_framework") is False,
		"native LV QoI card misclassifies its solver or evidence")
	require(card.get("native_ale_source_files_sha256") == source_hash(root)
		== native.get("source_files_sha256"),
		"native LV QoI source identity is stale")
	case = card["case"]
	require(case.get("fixture") == "IdealizedLeftVentricleFixture"
		and case.get("steps") == 16 and close(case.get("period_s"), 0.8)
		and close(case.get("density_kg_m3"), 1050.0)
		and close(case.get("dynamic_viscosity_pa_s"), 0.012)
		and close(case.get("inlet_target_outward_flow_m3_s"), -1e-6)
		and case.get("outlet_pressure_pa") == 0.0
		and case.get("outlet_backflow_beta") == 0.5
		and close(case.get("initial_closed_surface_volume_m3"),
			5.4822741234586634e-5),
		"native LV QoI case is not the matched prescribed-motion case")
	definition = card["definitions"]
	require("P1 pressure" in definition.get("mean_pressure_pa", "")
		and "P2 velocity" in definition.get("kinetic_energy_j", "")
		and "p_previous+p_current" in definition.get("pressure_volume_integral_j", "")
		and "not myocardial" in definition.get("interpretation", ""),
		"native LV QoI definitions or limitations are missing")
	levels = card["levels"]
	require(len(levels) == 3
		and [(level.get("name"), level.get("tetrahedra")) for level in levels]
		== [("coarse_center_fan", 96), ("radial_internal_refinement", 384),
			("two_shell_internal_refinement", 672)],
		"native LV QoI mesh levels are missing")
	for level in levels:
		volumes = level["volume_m3"]
		pressures = level["mean_pressure_pa"]
		kinetic = level["kinetic_energy_j"]
		require(level.get("mpi_ranks") == 1 and level.get("accepted_steps") == 16
			and level.get("status") == "passed_numerical_qoi_diagnostic"
			and len(volumes) == len(pressures) == len(kinetic) == 16
			and all(math.isfinite(value) and value > 0 for value in volumes+kinetic)
			and all(math.isfinite(value) for value in pressures),
			"native LV QoI series is incomplete or nonfinite")
		require(close(volumes[-1], case["initial_closed_surface_volume_m3"])
			and close(level.get("periodic_volume_relative_error"),
				abs(volumes[-1]-case["initial_closed_surface_volume_m3"])
				/case["initial_closed_surface_volume_m3"], absolute=1e-13)
			and level["periodic_volume_relative_error"] < 1e-12
			and close(level.get("final_kinetic_energy_j"), kinetic[-1])
			and level["final_kinetic_energy_j"] > level["initial_kinetic_energy_j"]
			and close(level.get("maximum_absolute_mean_pressure_pa"),
				max(abs(value) for value in pressures))
			and pressures[8] < 0.0 and pressures[8] == min(pressures)
			and close(level.get("pressure_volume_integral_j"),
				pressure_volume_integral(pressures, volumes,
					case["initial_closed_surface_volume_m3"]), absolute=1e-12)
			and level["pressure_volume_integral_j"] < 0.0,
			"native LV QoI values do not reproduce the cycle diagnostics")
		require(level.get("maximum_relative_inlet_flow_error", 1) < 1e-12
			and level.get("maximum_relative_total_boundary_flow_defect", 1) < 1e-10
			and level.get("maximum_relative_gcl_error", 1) < 2e-12,
			"native LV QoI flow or geometric conservation failed")
	for index in range(16):
		require(all(close(levels[0]["volume_m3"][index],
			level["volume_m3"][index], absolute=1e-14) for level in levels[1:]),
			"native LV QoI meshes do not retain the same moving volume")
	return_diagnostic = card["two_cycle_return_diagnostic"]
	require("P2 velocity L2" in return_diagnostic.get("definition", "")
		and return_diagnostic.get("same_prescribed_geometry_at_cycle_ends") is True
		and return_diagnostic.get("periodic_flow_state_verified") is False
		and len(return_diagnostic.get("levels", [])) == 2,
		"native LV two-cycle return definition is incomplete")
	for one_cycle, two_cycle in zip(levels[:2], return_diagnostic["levels"]):
		first = two_cycle["first_cycle_end"]
		second = two_cycle["second_cycle_end"]
		require(two_cycle.get("name") == one_cycle["name"]
			and two_cycle.get("tetrahedra") == one_cycle["tetrahedra"]
			and two_cycle.get("accepted_steps") == 32
			and two_cycle.get("status") == "passed_two_cycles_not_periodic_flow"
			and close(first.get("time_s"), 0.8)
			and close(second.get("time_s"), 1.6)
			and close(first.get("volume_m3"), one_cycle["volume_m3"][-1])
			and close(first.get("mean_pressure_pa"), one_cycle["mean_pressure_pa"][-1])
			and close(first.get("kinetic_energy_j"), one_cycle["kinetic_energy_j"][-1])
			and close(first["volume_m3"], second["volume_m3"])
			and second["kinetic_energy_j"] > first["kinetic_energy_j"]
			and 0.1 < two_cycle.get("cycle_end_velocity_relative_l2", 0) < 1.0
			and two_cycle.get("maximum_relative_inlet_flow_error", 1) < 1e-12
			and two_cycle.get("maximum_relative_total_boundary_flow_defect", 1) < 1e-10
			and two_cycle.get("maximum_relative_gcl_error", 1) < 2e-12,
			"native LV two-cycle field-return evidence is inconsistent or overclaimed")
	comparison = card["comparison"]
	differences = [abs(a-b) for a, b in zip(
		levels[0]["mean_pressure_pa"], levels[1]["mean_pressure_pa"])]
	require(close(comparison.get("largest_absolute_pressure_difference_pa"),
		max(differences), absolute=1e-10)
		and comparison.get("largest_pressure_difference_step") ==
			differences.index(max(differences))+1
		and close(comparison.get("pressure_volume_integral_relative_difference_refined_denominator"),
			abs(levels[0]["pressure_volume_integral_j"]
				-levels[1]["pressure_volume_integral_j"])
				/abs(levels[1]["pressure_volume_integral_j"]))
		and close(comparison.get("two_shell_vs_one_shell_pressure_volume_integral_relative_difference_two_shell_denominator"),
			abs(levels[1]["pressure_volume_integral_j"]
				-levels[2]["pressure_volume_integral_j"])
				/abs(levels[2]["pressure_volume_integral_j"]))
		and close(comparison.get("two_shell_vs_one_shell_largest_absolute_pressure_difference_pa"),
			max(abs(a-b) for a,b in zip(levels[1]["mean_pressure_pa"],
				levels[2]["mean_pressure_pa"])), absolute=1e-10)
		and comparison.get("geometric_volume_periodic_on_all_levels") is True
		and comparison.get("fluid_state_periodicity_verified") is False
		and comparison.get("spatial_convergence_verified") is False
		and comparison.get("radial_grading_same_across_refinements") is False
		and comparison.get("physiological_validation_passed") is False,
		"native LV QoI comparison is inconsistent or overclaimed")
	print("T6 native prescribed LV numerical QoI evidence: PASS (physiology and convergence open)")


if __name__ == "__main__":
	main()
