#!/usr/bin/env python3
"""Bind the local matched-port LV QoI comparison to its implementations."""

import hashlib
import json
import math
from pathlib import Path


IMMERSED_PATHS = (
	"solvers/cpu/include/ImmersedFlowPort.hpp",
	"solvers/cpu/include/ImmersedTransientFlowRuntime.hpp",
	"solvers/cpu/include/ImmersedTransientDistributedOperator.hpp",
	"solvers/cpu/include/MovingImmersedTransientFlowRuntime.hpp",
	"solvers/cpu/include/IdealizedLeftVentricleFixture.hpp",
	"solvers/cpu/include/NativeTetVelocityField.hpp",
	"solvers/cpu/tests/test_phase7_lv_closure.cpp",
	"solvers/cpu/tests/test_native_tet_velocity_field.cpp",
	"solvers/cpu/tests/test_immersed_flow_port.cpp",
	"solvers/cpu/tests/test_immersed_transient_flow.cpp",
)


def require(condition, message):
	if not condition:
		raise RuntimeError(message)


def source_hash(root):
	digest = hashlib.sha256()
	for label in IMMERSED_PATHS:
		contents = (root/label).read_bytes()
		digest.update(label.encode()+b"\0"+str(len(contents)).encode()+b"\0"+contents)
	return digest.hexdigest()


def close(actual, expected, absolute=1e-12, relative=1e-10):
	return math.isfinite(actual) and math.isfinite(expected) and math.isclose(
		actual, expected, abs_tol=absolute, rel_tol=relative)


def main():
	root = Path(__file__).resolve().parents[1]
	card = json.loads((root/"benchmarks/t4_matched_lv_qoi_evidence.json").read_text())
	ale = json.loads((root/"benchmarks/t4_native_ale_reference_evidence.json").read_text())
	require(card.get("schema_version") == 1, "unsupported matched LV QoI evidence schema")
	require(card.get("kind") == "matched_open_port_lv_functional_qoi_comparison_not_spatial_or_physical_validation"
		and card.get("external_fem_framework") is False, "matched LV classification is inaccurate")
	require(card.get("native_ale_source_files_sha256") == ale.get("source_files_sha256"),
		"matched LV card does not bind the native ALE evidence source")
	require(card.get("immersed_source_files_sha256") == source_hash(root),
		"immersed LV source changed; rerun the matched cycle and update evidence")
	case = card["shared_case"]
	require(case.get("steps") == 16 and close(case.get("period_s"), 0.8)
		and close(case.get("density_kg_m3"), 1050.0)
		and close(case.get("dynamic_viscosity_pa_s"), 0.012)
		and close(case.get("inlet_outward_flow_target_m3_s"), -1e-6)
		and close(case.get("outlet_pressure_pa"), 0.0)
		and close(case.get("outlet_backflow_beta"), 0.5)
		and close(case.get("immersed_wall_inertial_gamma0"), 0.0)
		and case.get("same_source_surface_labels_and_motion") is True
		and case.get("same_continuous_open_boundary_term") is True
		and case.get("same_spatial_discretization_or_wall_enforcement") is False,
		"matched LV physics or discrete-model scope changed")
	runtime = card["runtime"]
	require(runtime.get("native_ale_mpi_ranks") == 2
		and runtime.get("native_ale_tetrahedra") == 96
		and runtime.get("immersed_mpi_ranks") == 1
		and runtime.get("immersed_cartesian_grid") == [6, 6, 7]
		and runtime.get("immersed_cut_cell_max_depth") == 2
		and runtime.get("accepted_steps_each") == 16
		and runtime.get("status") == "passed_functional_cycle",
		"matched LV functional cycle evidence is incomplete")
	rms = card["rms_speed_m_s"]
	ale_values = rms["native_ale_values_from_six_digit_console_output"]
	immersed_values = rms["immersed_values"]
	require(len(ale_values) == len(immersed_values) == 16
		and all(math.isfinite(value) and value > 0 for value in ale_values+immersed_values)
		and rms.get("relative_difference_denominator") == "immersed",
		"matched LV RMS samples are invalid")
	errors = [abs(a-i)/i for a, i in zip(ale_values, immersed_values)]
	for key, step in (("end_systolic_step_8_relative_difference", 8),
		("first_expansion_step_9_relative_difference", 9),
		("final_end_diastolic_step_16_relative_difference", 16)):
		require(close(rms.get(key), errors[step-1], absolute=1e-14),
			f"matched LV {key} does not match recorded samples")
	require(close(rms.get("maximum_relative_difference"), max(errors), absolute=1e-14)
		and rms.get("maximum_difference_step") == errors.index(max(errors))+1
		and rms.get("relative_l2_field_error_measured") is True
		and rms.get("relative_l2_field_error_scope") == "full_16_step_cycle"
		and rms.get("spatial_convergence_verified") is False,
		"matched LV field-error conclusion is inaccurate")
	flow = card["selected_outlet_outward_flow_m3_s"]
	require(abs(flow["ale_step_8_from_console"]-flow["immersed_step_8"]) < 5e-11
		and abs(flow["ale_step_16_from_console"]-flow["immersed_step_16"]) < 5e-11,
		"matched LV selected outlet flows do not agree at ALE console precision")
	geometry = card["geometry_integration"]
	for prefix in ("end_diastolic", "end_systolic"):
		closed = geometry[f"shared_closed_surface_{prefix}_volume_m3"]
		quadrature = geometry[f"immersed_quadrature_{prefix}_volume_m3"]
		require(close(geometry[f"{prefix}_relative_quadrature_excess"],
			(quadrature-closed)/closed, absolute=1e-14),
			f"matched LV {prefix} integration discrepancy is inconsistent")
		depth_3 = geometry[f"depth_3_{prefix}_quadrature_volume_m3"]
		require(close(geometry[f"depth_3_{prefix}_relative_quadrature_excess"],
			(depth_3-closed)/closed, absolute=1e-14)
			and abs(depth_3-closed) < abs(quadrature-closed)
			and geometry[f"depth_3_{prefix}_logical_points"] > 100000,
			f"matched LV {prefix} depth-three integration did not improve")
	require(geometry.get("immersed_quadrature_geometry_error_refined_away") is False
		and geometry.get("depth_3_reduces_both_absolute_volume_errors") is True
		and 0 < card.get("maximum_immersed_normalized_moving_continuity_defect", 1) < 1e-10
		and 0 < card.get("maximum_immersed_normalized_wall_relative_leakage", 1) < 0.01
		and card.get("physical_validation_passed") is False,
		"matched LV error limitations are missing")
	first_step = card["depth_3_first_step_flow_probe"]
	field = card["first_step_common_point_field_probe"]
	cycle_field = card["full_cycle_common_point_field_probe"]
	values = cycle_field.get("relative_l2_velocity_error_by_step", [])
	points = cycle_field.get("overlap_points_by_step", [])
	require(cycle_field.get("command") == "phase7_lv_closure_test --matched-field-probe 16"
		and cycle_field.get("mpi_ranks") == 1
		and cycle_field.get("accepted_steps_each") == 16
		and cycle_field.get("status") == "passed_functional_field_comparison"
		and len(values) == len(points) == 16
		and all(math.isfinite(value) and 0 < value < 1 for value in values)
		and all(isinstance(count, int) and count > 0 for count in points)
		and cycle_field.get("excluded_points_each_step") == 0
		and cycle_field.get("excluded_volume_m3_each_step") == 0.0
		and close(cycle_field.get("maximum_relative_l2_error"), max(values))
		and cycle_field.get("maximum_error_step") == values.index(max(values))+1
		and close(cycle_field.get("minimum_relative_l2_error"), min(values))
		and cycle_field.get("minimum_error_step") == values.index(min(values))+1
		and cycle_field.get("spatial_or_temporal_convergence_verified") is False
		and cycle_field.get("wall_enforcement_difference_attributed") is False
		and cycle_field.get("physical_validation_passed") is False,
		"matched LV full-cycle common-point field evidence is missing or overclaimed")
	require(field.get("command") == "phase7_lv_closure_test --matched-field-probe 1"
		and field.get("step") == 1 and close(field.get("time_s"), 0.05)
		and field.get("mpi_ranks") == 1
		and field.get("overlap_points") == 97824
		and field.get("excluded_points") == 0
		and close(field.get("overlap_volume_m3"), 5.3857711532501893e-5)
		and field.get("excluded_volume_m3") == 0.0
		and close(field.get("relative_l2_velocity_error"), 0.43893307288816696)
		and field.get("ale_residual", 1.0) < 1e-9
		and field.get("field_convergence_or_physical_validation_proven") is False,
		"matched LV common-point first-step field evidence is missing or overclaimed")
	require(close(values[0], field["relative_l2_velocity_error"])
		and points[0] == field["overlap_points"],
		"first-step and full-cycle common-point field results disagree")
	depth_three_field = card["depth_3_first_step_common_point_field_probe"]
	require(depth_three_field.get("command") == "phase7_lv_closure_test --matched-field-probe 1 3"
		and depth_three_field.get("accepted_steps") == 1
		and depth_three_field.get("immersed_cut_cell_depth") == 3
		and close(depth_three_field.get("depth_2_relative_l2_velocity_error"), values[0])
		and close(depth_three_field.get("absolute_relative_l2_reduction"),
			values[0]-depth_three_field["relative_l2_velocity_error"])
		and 0.0 < depth_three_field["absolute_relative_l2_reduction"] < 0.01
		and depth_three_field.get("overlap_points", 0) > points[0]
		and depth_three_field.get("excluded_points") == 0
		and depth_three_field.get("excluded_volume_m3") == 0.0
		and depth_three_field.get("ale_residual", 1.0) < 1e-9
		and depth_three_field.get("elapsed_wall_s", 0) > 0
		and depth_three_field.get("peak_rss_kib", 0) > 0
		and depth_three_field.get("status") == "passed_single_step_diagnostic_not_spatial_convergence",
		"depth-three LV common-point first-step field evidence is missing or overclaimed")
	region = card["first_step_cell_region_error_diagnostic"]
	require(region.get("command") == "phase7_lv_closure_test --matched-field-probe 1 2"
		and region.get("interior_points", 0) > 0
		and region.get("cut_points", 0) > 0
		and region["interior_points"]+region["cut_points"] == points[0]
		and close(region.get("interior_cell_relative_l2"), 0.47653989391853041)
		and close(region.get("cut_cell_relative_l2"), 0.42898021692569827)
		and 0 < region.get("cut_volume_fraction", 0) < 1
		and 0 < region.get("cut_difference_squared_integral_fraction", 0) < 1
		and region["interior_cell_relative_l2"] > region["cut_cell_relative_l2"]
		and region.get("wall_enforcement_difference_attributed") is False,
		"first-step cell-region field diagnostic is missing or overclaimed")
	radial = card["radial_ale_first_step_common_point_field_probe"]
	require(radial.get("command") == "phase7_lv_closure_test --matched-field-probe 1 2 radial-ale"
		and radial.get("accepted_steps") == 1
		and radial.get("native_ale_tetrahedra") == 384
		and radial.get("unrefined_native_ale_tetrahedra") == 96
		and radial.get("shared_external_surface_and_wall_motion") is True
		and close(radial.get("unrefined_relative_l2_velocity_error"), values[0])
		and close(radial.get("absolute_relative_l2_reduction"),
			values[0]-radial["relative_l2_velocity_error"])
		and radial["absolute_relative_l2_reduction"] >
			depth_three_field["absolute_relative_l2_reduction"]
		and radial.get("overlap_points") == points[0]
		and radial.get("excluded_points") == 0
		and 0 < radial.get("interior_cell_relative_l2", 0) < region["interior_cell_relative_l2"]
		and 0 < radial.get("cut_cell_relative_l2", 0) < region["cut_cell_relative_l2"]
		and 0 < radial.get("cut_difference_squared_integral_fraction", 0) < 1
		and radial.get("ale_residual", 1.0) < 1e-9
		and radial.get("elapsed_wall_s", 0) > 0
		and radial.get("peak_rss_kib", 0) > 0
		and radial.get("full_cycle_convergence_proven") is False
		and radial.get("known_unstabilized_refined_cycle_failure_step") == 9
		and radial.get("backflow_stabilized_refined_cycle_tested") is True
		and radial.get("status") == "passed_single_step_sensitivity_not_spatial_convergence",
		"radial ALE LV first-step field sensitivity is missing or overclaimed")
	two_shell = card["two_shell_ale_first_step_common_point_field_probe"]
	require(two_shell.get("command") ==
		"phase7_lv_closure_test --matched-field-probe 1 2 two-shell-ale"
		and two_shell.get("accepted_steps") == 1
		and two_shell.get("native_ale_tetrahedra") == 672
		and two_shell.get("shared_external_surface_and_wall_motion") is True
		and two_shell.get("immersed_cut_cell_depth") == 2
		and close(two_shell.get("one_shell_relative_l2_velocity_error"),
			radial["relative_l2_velocity_error"])
		and close(two_shell.get("absolute_relative_l2_reduction_vs_one_shell"),
			radial["relative_l2_velocity_error"]-two_shell["relative_l2_velocity_error"])
		and values[0] > radial["relative_l2_velocity_error"]
			> two_shell["relative_l2_velocity_error"] > 0
		and two_shell.get("overlap_points") == points[0]
		and two_shell.get("excluded_points") == 0
		and two_shell.get("excluded_volume_m3") == 0
		and close(two_shell.get("immersed_overlap_rms_m_s"),
			radial["immersed_overlap_rms_m_s"])
		and 0 < two_shell.get("interior_cell_relative_l2", 0)
			< radial["interior_cell_relative_l2"]
		and 0 < two_shell.get("cut_cell_relative_l2", 0)
			< radial["cut_cell_relative_l2"]
		and two_shell.get("ale_residual", 1.0) < 1e-9
		and two_shell.get("full_cycle_or_spatial_convergence_proven") is False
		and two_shell.get("status") ==
			"passed_single_step_three_mesh_sensitivity_not_convergence",
		"two-shell ALE LV field diagnostic is missing or overclaimed")
	radial_cycle = card["radial_ale_full_cycle_common_point_field_probe"]
	radial_values = radial_cycle.get("relative_l2_velocity_error_by_step", [])
	require(radial_cycle.get("command") == "phase7_lv_closure_test --matched-field-probe 16 2 radial-ale"
		and radial_cycle.get("accepted_steps_each") == 16
		and radial_cycle.get("native_ale_tetrahedra") == 384
		and radial_cycle.get("immersed_cut_cell_depth") == 2
		and len(radial_values) == len(values) == 16
		and all(math.isfinite(value) and 0 < value < baseline for value,baseline in zip(radial_values,values))
		and close(radial_cycle.get("maximum_relative_l2_error"), max(radial_values))
		and radial_cycle.get("maximum_error_step") == radial_values.index(max(radial_values))+1
		and close(radial_cycle.get("minimum_relative_l2_error"), min(radial_values))
		and radial_cycle.get("minimum_error_step") == radial_values.index(min(radial_values))+1
		and close(radial_cycle.get("mean_relative_l2_error"), sum(radial_values)/16)
		and close(radial_cycle.get("mean_relative_l2_reduction_vs_96_tetra"),
			sum(a-b for a,b in zip(values,radial_values))/16)
		and close(radial_cycle.get("first_expansion_step_reduction_vs_96_tetra"),
			values[8]-radial_values[8])
		and close(radial_cycle.get("final_end_diastolic_reduction_vs_96_tetra"),
			values[15]-radial_values[15])
		and close(radial_values[0],radial["relative_l2_velocity_error"])
		and radial_cycle.get("all_steps_lower_error_than_96_tetra") is True
		and radial_cycle.get("overlap_point_counts_match_96_tetra") is True
		and radial_cycle.get("excluded_points_each_step") == 0
		and radial_cycle.get("spatial_convergence_verified") is False
		and radial_cycle.get("immersed_grid_refined_in_this_comparison") is False
		and radial_cycle.get("status") == "passed_full_cycle_mesh_sensitivity_not_convergence",
		"radial ALE full-cycle field sensitivity is missing or overclaimed")
	two_shell_cycle = card["two_shell_ale_full_cycle_common_point_field_probe"]
	two_shell_values = two_shell_cycle.get("relative_l2_velocity_error_by_step", [])
	two_shell_points = two_shell_cycle.get("overlap_points_by_step", [])
	require(two_shell_cycle.get("command") ==
		"phase7_lv_closure_test --matched-field-probe 16 2 two-shell-ale"
		and two_shell_cycle.get("accepted_steps_each") == 16
		and two_shell_cycle.get("native_ale_tetrahedra") == 672
		and two_shell_cycle.get("immersed_cut_cell_depth") == 2
		and len(two_shell_values) == len(two_shell_points) == 16
		and all(0 < value < previous for value,previous in
			zip(two_shell_values,radial_values))
		and all(isinstance(count,int) and count > 0 for count in two_shell_points)
		and two_shell_points == points
		and two_shell_points[0] == two_shell["overlap_points"]
		and close(two_shell_values[0],two_shell["relative_l2_velocity_error"])
		and two_shell_cycle.get("excluded_points_each_step") == 0
		and two_shell_cycle.get("excluded_volume_m3_each_step") == 0
		and two_shell_cycle.get("maximum_ale_residual", 1.0) < 1e-9
		and close(two_shell_cycle.get("maximum_relative_l2_error"),max(two_shell_values))
		and two_shell_cycle.get("maximum_error_step") ==
			two_shell_values.index(max(two_shell_values))+1
		and close(two_shell_cycle.get("minimum_relative_l2_error"),min(two_shell_values))
		and two_shell_cycle.get("minimum_error_step") ==
			two_shell_values.index(min(two_shell_values))+1
		and close(two_shell_cycle.get("mean_relative_l2_error"),
			sum(two_shell_values)/16)
		and close(two_shell_cycle.get("mean_relative_l2_reduction_vs_384_tetra"),
			sum(a-b for a,b in zip(radial_values,two_shell_values))/16)
		and close(two_shell_cycle.get("first_expansion_step_reduction_vs_384_tetra"),
			radial_values[8]-two_shell_values[8])
		and close(two_shell_cycle.get("final_end_diastolic_reduction_vs_384_tetra"),
			radial_values[15]-two_shell_values[15])
		and two_shell_cycle.get("all_steps_lower_error_than_384_tetra") is True
		and two_shell_cycle.get("overlap_point_counts_match_384_tetra") is True
		and two_shell_cycle.get("spatial_convergence_verified") is False
		and two_shell_cycle.get("radial_grading_same_as_384_tetra") is False
		and two_shell_cycle.get("immersed_grid_refined_in_this_comparison") is False
		and two_shell_cycle.get("status") ==
			"passed_full_cycle_three_mesh_sensitivity_not_convergence",
		"two-shell ALE full-cycle field sensitivity is missing or overclaimed")
	require(ale["radial_backflow_stabilized_lv_cycle"].get("status") == "passed"
		and ale["radial_backflow_stabilized_lv_cycle"].get("accepted_steps") == 16,
		"matched radial field card lacks independently accepted native ALE cycle")
	require(first_step.get("accepted_steps") == 1
		and first_step.get("nonlinear_iterations", 0) > 0
		and first_step.get("residual_norm", 1) < 1e-10
		and first_step.get("status") == "passed_single_step_diagnostic"
		and first_step.get("flow_field_convergence_proven") is False,
		"depth-three LV single-step result is missing or overclaimed")
	ale_first = first_step["native_ale_first_step_rms_speed_from_console_m_s"]
	for depth in (2, 3):
		immersed_first = first_step[f"depth_{depth}_first_step_rms_speed_m_s"] if depth == 2 else first_step["rms_speed_m_s"]
		require(close(first_step[f"depth_{depth}_relative_rms_difference_vs_ale"],
			abs(ale_first-immersed_first)/immersed_first, absolute=1e-14),
			f"depth-{depth} first-step RMS difference is inconsistent")
	require(first_step["depth_3_relative_rms_difference_vs_ale"]
		> first_step["depth_2_relative_rms_difference_vs_ale"],
		"depth-three flow field unexpectedly claimed improvement")
	legacy = card["legacy_full_immersed_transient_regression"]
	require(legacy.get("exit_code") == 124 and legacy.get("counted_as_pass") is False,
		"timed-out full immersed regression must not count as passed")
	print("T4 matched LV functional QoI evidence: PASS (spatial/physical validation still open)")


if __name__ == "__main__":
	main()
