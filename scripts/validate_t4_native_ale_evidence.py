#!/usr/bin/env python3
"""Validate native ALE small-domain evidence and bind it to its source files."""

import hashlib
import json
import math
from pathlib import Path


def require(condition, message):
	if not condition:
		raise RuntimeError(message)


def source_hash(root):
	# T5 adapters build on ALE but have their own evidence card. Keep this T4
	# digest bound to the standalone ALE operators/runtimes so an FSI-only edit
	# does not falsely stale otherwise unchanged T4 numerical measurements.
	paths = sorted([path for path in (root/"solvers/cpu/include").glob("NativeTetAle*.hpp")
		if "Fsi" not in path.name]
		+ [root/"solvers/cpu/include/NativeTetFem.hpp",
			root/"solvers/cpu/include/NativeTetStarRadialRefinement.hpp",
			root/"solvers/cpu/include/IdealizedLeftVentricleFixture.hpp",
			root/"solvers/cpu/include/PrescribedSurfaceMotion.hpp",
			root/"solvers/cpu/include/MaterialSurfaceKinematics.hpp"]
		+ [path for path in (root/"solvers/cpu/tests").glob("test_native_tet_ale_*.cpp")
			if "fsi" not in path.name]
		+ [root/"solvers/cpu/tests/test_native_tet_star_radial_refinement.cpp"])
	digest = hashlib.sha256()
	for path in paths:
		label = str(path.relative_to(root))
		contents = path.read_bytes()
		digest.update(label.encode()+b"\0"+str(len(contents)).encode()+b"\0"+contents)
	return digest.hexdigest()


def main():
	root = Path(__file__).resolve().parents[1]
	data = json.loads((root/"benchmarks/t4_native_ale_reference_evidence.json").read_text())
	require(data.get("schema_version") == 1, "unsupported T4 evidence schema")
	require(data.get("kind") ==
		"native_ALE_small_domain_numerical_reference_not_production_validation",
		"T4 evidence classification is inaccurate")
	require(data.get("external_fem_framework") is False,
		"external FEM output cannot count as native ALE evidence")
	require(data.get("source_files_sha256") == source_hash(root),
		"native ALE source hash is stale; rerun tests and update evidence")
	runtime = data["global_runtime"]
	require(runtime.get("tetrahedra") == 4 and runtime.get("free_velocity_components") == 15,
		"global runtime must exercise free cross-element velocity DOFs")
	require(runtime.get("accepted_steps", 0) >= 5
		and runtime.get("paired_geometry_and_fluid_commit") is True,
		"multi-step transactional evidence is missing")
	require(runtime.get("GCL_and_moving_domain_balance_checked_each_step") is True,
		"per-step conservation diagnostics are missing")
	require(runtime.get("fallible_prepare_then_noexcept_geometry_finalize") is True,
		"ALE geometry commit is not safe for paired FSI finalization")
	distributed = data["distributed_assembly"]
	require(distributed.get("mpi_ranks", 0) >= 2 and distributed.get("status") == "passed",
		"distributed PETSc ALE assembly evidence is missing")
	require(distributed.get("maximum_relative_entry_error", 1.0) < 1e-12,
		"distributed PETSc ALE assembly differs from serial native assembly")
	distributed_solve = data["distributed_solve"]
	require(distributed_solve.get("mpi_ranks", 0) >= 2
		and distributed_solve.get("status") == "passed",
		"distributed PETSc ALE solve evidence is missing")
	require(distributed_solve.get("newton_iterations", 0) > 0
		and distributed_solve.get("final_residual_l2", 1.0) < 1e-10,
		"distributed PETSc ALE nonlinear solve failed its gate")
	steady = data["distributed_steady_solve"]
	require(steady.get("mpi_ranks", 0) >= 2 and steady.get("status") == "passed"
		and steady.get("newton_iterations", 0) > 0
		and steady.get("final_residual_l2", 1.0) < 1e-10,
		"distributed steady PETSc solve evidence is missing")
	iterative = data["distributed_iterative_solve"]
	require(iterative.get("classification") ==
		"small_domain_solver_equivalence_not_scalability_validation"
		and iterative.get("mpi_ranks", 0) >= 2 and iterative.get("status") == "passed"
		and iterative.get("physical_row_column_scaling") is True
		and iterative.get("direct_residual_l2", 1.0) < 1e-10
		and iterative.get("iterative_residual_l2", 1.0) < 1e-10,
		"distributed iterative ALE equivalence evidence is missing")
	time_loop = data["distributed_time_loop"]
	require(time_loop.get("mpi_ranks", 0) >= 2 and time_loop.get("accepted_steps", 0) >= 5
		and time_loop.get("status") == "passed",
		"distributed ALE multi-step evidence is missing")
	require(time_loop.get("quality_and_conservation_gates_each_step") is True
		and time_loop.get("paired_geometry_and_fluid_commit") is True,
		"distributed ALE time loop did not preserve the transactional contract")
	coarse = data["coarse_mesh_translation_smoke"]
	require(coarse.get("classification") ==
		"functional_smoke_not_deforming_domain_validation" and coarse.get("status") == "passed",
		"coarse ALE mesh evidence is misclassified or failed")
	require(coarse.get("tetrahedra", 0) >= 5000 and coarse.get("mixed_dofs", 0) >= 29000,
		"coarse ALE mesh evidence is not representative of the prepared mesh")
	require(coarse.get("final_residual_l2", 1.0) < 1e-8
		and coarse.get("maximum_velocity_error_m_s", 1.0)
		< coarse.get("functional_velocity_tolerance_m_s", 0.0),
		"coarse ALE mesh functional gate failed")
	require(coarse.get("minimum_determinant_ratio", 0.0) >= 0.05
		and coarse.get("minimum_scaled_jacobian", 0.0) >= 1e-3,
		"coarse ALE mesh quality gate failed")
	require(abs(coarse.get("gcl_residual_m3_s", 1.0)) < 1e-12
		and abs(coarse.get("moving_domain_balance_residual_m3_s", 1.0)) < 1e-12,
		"coarse ALE mesh conservation gate failed")
	shear = data["coarse_mesh_affine_shear"]
	require(shear.get("classification") ==
		"deforming_domain_functional_verification_not_chamber_validation"
		and shear.get("status") == "passed",
		"coarse deforming ALE evidence is missing or misclassified")
	require(shear.get("tetrahedra", 0) >= 5000 and shear.get("mixed_dofs", 0) >= 29000,
		"deforming ALE evidence did not use the prepared coarse mesh")
	require(shear.get("final_residual_l2", 1.0) < 1e-8
		and shear.get("maximum_velocity_error_m_s", 1.0)
		< shear.get("functional_velocity_tolerance_m_s", 0.0),
		"deforming ALE functional gate failed")
	require(shear.get("minimum_determinant_ratio", 0.0) >= 0.05
		and shear.get("minimum_scaled_jacobian", 0.0) >= 1e-3,
		"deforming ALE mesh quality gate failed")
	require(abs(shear.get("gcl_residual_m3_s", 1.0)) < 1e-12
		and abs(shear.get("moving_domain_balance_residual_m3_s", 1.0)) < 1e-12,
		"deforming ALE conservation gate failed")
	shared=data["shared_immersed_chamber_motion_contract"]
	require(shared.get("classification")=="shared_geometry_and_wall_motion_parity_not_flow_field_comparison"
		and shared.get("time_steps")==16 and shared.get("same_material_vertex_ids") is True
		and shared.get("functional_rms_and_port_qoi_comparison_available") is True
		and shared.get("flow_field_qoi_comparison_completed") is False,
		"ALE/immersed shared chamber contract is missing or overclaimed")
	require(shared.get("maximum_relative_volume_difference",1)>=0
		and shared.get("maximum_relative_volume_difference",1)<2e-14
		and shared.get("maximum_wall_position_difference_m",1)<1e-15
		and shared.get("maximum_wall_velocity_difference_m_s",1)<1e-14
		and shared.get("maximum_relative_gcl_residual",1)<2e-12,
		"ALE/immersed shared chamber motion parity failed")
	cycle=data["optional_backflow_stabilized_lv_cycle"]
	require(cycle.get("classification")==
		"prescribed_motion_functional_cycle_not_physical_or_spatial_validation"
		and cycle.get("external_fem_framework") is False
		and cycle.get("mpi_ranks_tested")==[1,2]
		and cycle.get("accepted_steps")==16
		and cycle.get("pressure_port_backflow_beta")==0.5
		and cycle.get("same_open_boundary_model_as_original_immersed_lv_reference") is False
		and cycle.get("matched_open_boundary_immersed_probe_available") is True
		and cycle.get("periodic_geometry_returned") is True
		and cycle.get("periodic_flow_state_verified") is False
		and cycle.get("status")=="passed",
		"optional-backflow native ALE LV evidence is missing or overclaimed")
	require(cycle.get("maximum_relative_inlet_flow_error",1)<1e-12
		and cycle.get("maximum_relative_total_boundary_flow_defect",1)<1e-10
		and cycle.get("maximum_relative_gcl_residual",1)<2e-12
		and cycle.get("minimum_determinant_ratio",0)>0.05
		and cycle.get("minimum_scaled_jacobian",0)>1e-3
		and cycle.get("end_systolic_outlet_flow_m3_s",0)>0
		and cycle.get("final_end_diastolic_outlet_flow_m3_s",0)<0
		and cycle.get("total_newton_updates",0)>0,
		"optional-backflow native ALE LV numerical gate failed")
	radial=data["radial_backflow_stabilized_lv_cycle"]
	require(radial.get("classification") ==
		"one_rank_prescribed_motion_functional_cycle_not_spatial_convergence"
		and radial.get("command") == "native_tet_ale_idealized_lv_radial_backflow_flow_test"
		and radial.get("tetrahedra") == 384
		and radial.get("same_exterior_triangles_labels_and_wall_motion_as_96_tetra_case") is True
		and radial.get("mpi_ranks") == 1
		and radial.get("accepted_steps") == 16
		and radial.get("pressure_port_backflow_beta") == cycle["pressure_port_backflow_beta"]
		and radial.get("maximum_relative_inlet_flow_error",1) < 1e-12
		and radial.get("maximum_relative_total_boundary_flow_defect",1) < 1e-10
		and radial.get("maximum_relative_gcl_residual",1) < 2e-12
		and radial.get("minimum_determinant_ratio",0) > 0.05
		and radial.get("minimum_scaled_jacobian",0) > 1e-3
		and radial.get("end_systolic_outlet_flow_m3_s",0) > 0
		and radial.get("final_end_diastolic_outlet_flow_m3_s",0) < 0
		and radial.get("total_newton_updates",0) == radial.get("total_linear_iterations")
		and radial.get("total_newton_updates",0) > 0
		and radial.get("unstabilized_radial_cycle_fails_at_step") == 9
		and radial.get("spatial_convergence_verified") is False
		and radial.get("status") == "passed",
		"backflow-stabilized refined native ALE LV cycle is missing or overclaimed")
	controls=data["native_ale_inlet_flow_and_outlet_pressure"]
	require(controls.get("classification")==
		"two_rank_functional_boundary_control_not_chamber_comparison"
		and controls.get("external_fem_framework") is False
		and controls.get("mpi_ranks")>=2 and controls.get("status")=="passed"
		and controls.get("maximum_absolute_flow_error_m3_s",1)<1e-12
		and controls.get("nonzero_outlet_pressure_shift_pa")==5.0
		and controls.get("nonzero_pressure_shift_preserves_velocity") is True
		and controls.get("nonzero_pressure_shift_changes_controller_multiplier_by_pa")==5.0,
		"native ALE pressure/flow boundary control evidence is missing or misclassified")
	attempt=data["idealized_lv_flow_attempt"]
	require(attempt.get("classification")==
		"incomplete_negative_result_not_ale_immersed_flow_comparison"
		and attempt.get("same_source_surface_motion_and_boundary_labels_as_immersed") is True
		and attempt.get("same_density_viscosity_initial_rotation_and_port_targets_as_immersed") is True
		and attempt.get("completed_steps_before_failure")==13
		and attempt.get("failure_step")==14 and attempt.get("requested_steps")==16
		and attempt.get("residual_backtracking_line_search_resolved") is False
		and [(item.get("internal_steps"),item.get("failure_step"))
			for item in attempt.get("time_refinement_trials",[])]==[(32,29),(64,62)]
		and attempt.get("radial_refinement_trial",{}).get("tetrahedra")==384
		and attempt.get("radial_refinement_trial",{}).get("failure_step")==9
		and attempt.get("promoted_to_validation") is False,
		"incomplete ALE LV attempt must remain an explicit negative result")
	failed = data.get("failed_attempts", [])
	require(any(item.get("physical_validation_promoted") is False for item in failed),
		"coarse smoke tolerance change must retain its rejected attempt")
	temporal = data["temporal_convergence"]
	dt = temporal["time_steps_s"]
	errors = temporal["absolute_errors"]
	orders = temporal["observed_orders"]
	require(len(dt) == len(errors) == 4 and len(orders) == 3,
		"temporal refinement evidence has the wrong number of levels")
	for level in range(1, len(dt)):
		require(math.isclose(dt[level-1]/dt[level], 2.0),
			"temporal refinements must halve dt")
		require(errors[level] < errors[level-1], "temporal error does not decrease")
		computed = math.log(errors[level-1]/errors[level], 2.0)
		require(math.isclose(computed, orders[level-1], rel_tol=2e-6),
			"recorded temporal order differs from errors")
		require(0.85 < computed < 1.15, "backward-Euler temporal order is outside its gate")
	print("T4 native ALE reference evidence: PASS")


if __name__ == "__main__":
	try:
		main()
	except (KeyError, OSError, RuntimeError, TypeError, ValueError) as error:
		print(f"T4 native ALE reference evidence: ERROR: {error}")
		raise SystemExit(2)
