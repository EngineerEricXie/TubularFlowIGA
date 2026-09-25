#!/usr/bin/env python3
"""Validate the native FEM Y-bifurcation functional evidence."""

import hashlib
import json
from pathlib import Path

SOURCE_FILES = (
	"solvers/cpu/include/NativeTetFem.hpp",
	"solvers/cpu/include/NativeTetAleTransient.hpp",
	"solvers/cpu/include/NativeTetAlePetscAssembly.hpp",
	"solvers/cpu/include/NativeTetAlePetscRuntime.hpp",
	"solvers/cpu/include/NativeTetBoundaryFlow.hpp",
	"solvers/cpu/include/NativeTetDarcyPetsc.hpp",
	"solvers/cpu/include/NativeTetMatchingInterfaceSource.hpp",
	"solvers/cpu/tests/test_native_tet_bifurcation_runtime.cpp",
	"scripts/generate_y_pipe_surface.py",
	"scripts/prepare_t2_bifurcation_meshes.py",
	"scripts/surface_to_fem_volume.py",
)

def require(condition, message):
	if not condition:
		raise RuntimeError(message)

def source_hash(root):
	digest = hashlib.sha256()
	for name in SOURCE_FILES:
		contents = (root/name).read_bytes()
		digest.update(name.encode()+b"\0"+str(len(contents)).encode()+b"\0"+contents)
	return digest.hexdigest()

def main():
	root = Path(__file__).resolve().parents[1]
	data = json.loads((root/"benchmarks/t2_native_bifurcation_evidence.json").read_text())
	require(data.get("kind") ==
		"native_FEM_bifurcation_numerical_convergence_not_physiological_validation",
		"bifurcation evidence is misclassified")
	require(data.get("external_fem_framework") is False,
		"external FEM output cannot count as native bifurcation evidence")
	require(data.get("current_source_files_sha256") == source_hash(root),
		"native bifurcation current source hash is stale")
	require(data.get("current_source_status") ==
		"historical_natural_outlet_qoi_not_current_zero_pressure_outlet_regression",
		"historical bifurcation QoI must not be claimed as current-source validation")
	mesh = data["mesh"]
	require(mesh.get("tetrahedra", 0) > 1000 and mesh.get("mixed_dofs", 0) > 8000,
		"bifurcation mesh evidence is incomplete")
	require(mesh.get("minimum_scaled_jacobian", 0.0) >= 1e-3
		and mesh.get("relative_volume_error", 1.0) < 5e-3,
		"bifurcation mesh failed geometry gates")
	runtime = data["runtime"]
	require(runtime.get("mpi_ranks", 0) >= 2 and runtime.get("time_steps", 0) >= 5
		and runtime.get("final_residual_l2", 1.0) < 1e-8,
		"bifurcation runtime failed its functional gate")
	qoi = data["qoi"]
	require(qoi.get("inlet_outward_flow_m3_s", 0.0) < 0.0
		and qoi.get("upper_outlet_flow_m3_s", 0.0) > 0.0
		and qoi.get("lower_outlet_flow_m3_s", 0.0) > 0.0,
		"bifurcation flow signs are incorrect")
	require(0.45 < qoi.get("upper_flow_fraction", 0.0) < 0.55,
		"symmetric bifurcation flow split is outside its functional gate")
	require(qoi.get("relative_mass_imbalance", 1.0) < 1e-10,
		"bifurcation mass balance failed")
	refinement = data["fixed_surface_volume_refinement"]
	require(refinement.get("surface_geometry_identical_across_levels") is True,
		"bifurcation refinement changed the surface geometry")
	levels = refinement.get("levels", [])
	require(len(levels) == 3 and all(levels[index]["tetrahedra"]
		< levels[index+1]["tetrahedra"] for index in range(2)),
		"bifurcation refinement levels are invalid")
	require(all(level.get("relative_mass_imbalance", 1.0) < 1e-10 for level in levels),
		"a bifurcation refinement level failed mass balance")
	changes = refinement["medium_to_fine_relative_changes"]
	require(changes.get("upper_flow_fraction", 1.0)
		< refinement.get("flow_split_gate_maximum_relative_change", 0.0)
		and refinement.get("flow_split_status") == "passed",
		"bifurcation flow split did not pass its refinement gate")
	stress_gate = refinement.get("stress_qoi_gate_maximum_relative_change", 0.0)
	require(any(changes.get(key, 0.0) > stress_gate for key in
		("inlet_average_pressure", "wall_mean_wss", "wall_rms_wss",
		 "wall_integrated_traction_norm"))
		and refinement.get("pressure_traction_wss_status") == "not_converged",
		"nonconverged bifurcation stress QoIs must remain explicitly open")
	steady = data["steady_fixed_surface_spatial_convergence"]
	require(steady.get("classification") ==
		"idealized_numerical_convergence_not_physiological_validation"
		and steady.get("time_integration") == "steady"
		and steady.get("surface_geometry_identical_across_levels") is True
		and steady.get("status") == "passed",
		"steady bifurcation convergence evidence is misclassified")
	steady_levels = steady.get("levels", [])
	require(len(steady_levels) == 3 and all(steady_levels[index]["tetrahedra"]
		< steady_levels[index+1]["tetrahedra"] for index in range(2)),
		"steady bifurcation refinement levels are invalid")
	require(all(level.get("relative_mass_imbalance", 1.0)
		< steady.get("mass_gate_maximum", 0.0)
		and level.get("final_residual_l2", 1.0) < 1e-8 for level in steady_levels),
		"steady bifurcation level failed residual or mass gates")
	steady_changes = steady["medium_to_fine_relative_changes"]
	require(steady_changes.get("upper_flow_fraction", 1.0)
		< steady.get("flow_split_gate_maximum_relative_change", 0.0),
		"steady bifurcation flow split did not converge")
	require(all(steady_changes.get(key, 1.0)
		< steady.get("other_qoi_gate_maximum_relative_change", 0.0) for key in
		("inlet_average_pressure", "outlet_referenced_pressure_drop",
		 "volume_mean_speed", "volume_rms_speed", "kinetic_energy",
		 "wall_mean_total_traction", "wall_rms_total_traction",
		 "wall_mean_wss", "wall_rms_wss", "regular_wall_mean_total_traction",
		 "regular_wall_rms_total_traction", "regular_wall_mean_wss",
		 "regular_wall_rms_wss")),
		"a steady bifurcation pressure/velocity/traction/WSS QoI did not converge")
	study = data["one_step_fixed_surface_spatial_study"]
	require(study.get("time_steps") == 1 and study.get("dt_s") == 0.05,
		"bifurcation spatial-study time policy changed")
	study_levels = study.get("levels", [])
	require(len(study_levels) == 4 and all(study_levels[index]["tetrahedra"]
		< study_levels[index+1]["tetrahedra"] for index in range(3)),
		"four-level bifurcation spatial study is incomplete")
	require(all(level.get("relative_mass_imbalance", 1.0) < 1e-10
		and level.get("final_residual_l2", 1.0) < 1e-8 for level in study_levels),
		"a bifurcation spatial-study level failed its functional gates")
	study_changes = study["fine_to_finest_relative_changes"]
	require(study_changes.get("volume_rms_speed", 1.0)
		< study.get("velocity_qoi_gate_maximum_relative_change", 0.0)
		and study_changes.get("upper_flow_fraction", 1.0)
		< study.get("flow_split_gate_maximum_relative_change", 0.0)
		and study.get("velocity_flow_status") == "passed_integral_qoi_only",
		"bifurcation integrated velocity/flow QoIs did not pass their study gates")
	require(any(study_changes.get(key, 0.0)
		> study.get("stress_qoi_gate_maximum_relative_change", 0.0) for key in
		("inlet_average_pressure", "wall_mean_wss", "wall_rms_wss",
		 "wall_integrated_traction_norm"))
		and study.get("pressure_traction_wss_status") == "not_converged",
		"four-level nonconverged stress QoIs must remain explicitly open")
	coupled = data["coupled_surface_volume_refinement"]
	require(coupled.get("mesh_series_generator") ==
		"scripts/prepare_t2_bifurcation_meshes.py",
		"coupled bifurcation mesh series is not reproducible")
	coupled_levels = coupled.get("levels", [])
	require(len(coupled_levels) == 3 and all(
		coupled_levels[index]["surface_triangles"]
		< coupled_levels[index+1]["surface_triangles"]
		and coupled_levels[index]["tetrahedra"]
		< coupled_levels[index+1]["tetrahedra"] for index in range(2)),
		"coupled surface/volume refinement levels are invalid")
	require(all(level.get("minimum_scaled_jacobian", 0.0) >= 1e-3
		and level.get("relative_mass_imbalance", 1.0) < 1e-10
		and level.get("final_residual_l2", 1.0) < 1e-8 for level in coupled_levels),
		"a coupled bifurcation level failed geometry or solver gates")
	coupled_changes = coupled["medium_to_fine_qoi_relative_changes"]
	integral_gate = coupled.get("other_integral_qoi_gate_maximum_relative_change", 0.0)
	require(coupled_changes.get("upper_flow_fraction", 1.0)
		< coupled.get("flow_split_gate_maximum_relative_change", 0.0)
		and all(coupled_changes.get(key, 1.0) < integral_gate for key in
		("inlet_average_pressure", "volume_mean_speed", "volume_rms_speed",
		 "kinetic_energy", "wall_mean_total_traction", "wall_rms_total_traction",
		 "regular_wall_mean_total_traction", "regular_wall_rms_total_traction"))
		and coupled.get("velocity_flow_pressure_total_traction_status") ==
		"passed_integral_qoi_only",
		"coupled integral QoIs did not pass their functional gates")
	wss_gate = coupled.get("wss_gate_maximum_relative_change", 0.0)
	require(all(coupled_changes.get(key, 0.0) > wss_gate for key in
		("wall_mean_wss", "wall_rms_wss", "regular_wall_mean_wss",
		 "regular_wall_rms_wss"))
		and coupled.get("raw_and_regular_wall_wss_status") == "not_converged",
		"coupled nonconverged WSS must remain explicitly open")
	require(all(coupled_changes.get(key, 0.0) > wss_gate for key in
		("recovered_wall_mean_wss", "recovered_wall_rms_wss",
		 "recovered_regular_wall_mean_wss", "recovered_regular_wall_rms_wss"))
		and coupled.get("patch_recovered_wss_status") ==
		"not_converged_not_promoted",
		"unsuccessful patch-recovered WSS must not be promoted")
	timeouts = [item for item in data.get("failed_attempts", [])
		if item.get("case") == "coupled_surface_volume_refinement_level_4"]
	require(len(timeouts) == 1 and timeouts[0].get("tetrahedra") == 92179
		and timeouts[0].get("outcome") == "timeout_exit_124_no_result_promoted"
		and timeouts[0].get("accepted") is False,
		"fourth coupled-level timeout is missing or was promoted")
	iterative = data["iterative_field_split"]
	medium = iterative["medium_level"]
	require(medium.get("status") == "passed" and medium.get("mpi_ranks", 0) >= 4
		and medium.get("total_linear_iterations", 0) > 0
		and medium.get("final_linear_converged_reason", 0) > 0
		and medium.get("final_residual_l2", 1.0) < 1e-8
		and medium.get("relative_mass_imbalance", 1.0) < 1e-10
		and medium.get("maximum_relative_qoi_difference_from_direct", 1.0) < 1e-8,
		"medium scaled field-split evidence failed correctness gates")
	require(iterative.get("large_level_status") ==
		"not_scalable_with_current_selfp_configuration"
		and sum(item.get("case", "").startswith("field_split_selfp_large_level")
			and item.get("accepted") is False for item in data.get("failed_attempts", [])) == 2,
		"large field-split timeout limitations are missing")
	require(any(item.get("accepted") is False for item in data.get("failed_attempts", [])),
		"rejected pressure-gauge formulation was not retained")
	print("T2 native bifurcation evidence: PASS (historical QoI only; current outlet BC differs)")

if __name__ == "__main__":
	try:
		main()
	except (KeyError, OSError, RuntimeError, TypeError, ValueError) as error:
		print(f"T2 native bifurcation evidence: ERROR: {error}")
		raise SystemExit(2)
