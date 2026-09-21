#!/usr/bin/env python3
"""Audit the native moving-tetra species reference path and its scope."""

import hashlib
import json
import math
from pathlib import Path


SOURCE_FILES = (
	"solvers/cpu/include/NativeTetMovingSpeciesTransport.hpp",
	"solvers/cpu/tests/test_native_tet_moving_species_transport.cpp",
	"solvers/cpu/include/NativeTetMovingSpeciesPetscRuntime.hpp",
	"solvers/cpu/include/NativeTetSpeciesVisualization.hpp",
	"solvers/cpu/tests/test_native_tet_moving_species_petsc_runtime.cpp",
	"solvers/cpu/tests/test_native_tet_moving_species_spatial_convergence.cpp",
	"solvers/cpu/include/NativeTetMovingSpeciesCheckpoint.hpp",
	"solvers/cpu/tests/test_native_tet_moving_species_checkpoint.cpp",
	"solvers/cpu/include/NativeTetMovingSpeciesPorts.hpp",
	"solvers/cpu/tests/test_native_tet_moving_species_ports.cpp",
	"solvers/cpu/include/NativeTetAleGraphPorts.hpp",
	"solvers/cpu/tests/test_native_tet_ale_graph_ports.cpp",
	"solvers/cpu/tests/test_native_tet_ale_graph_petsc_runtime.cpp",
	"solvers/cpu/include/NativeTetAleFlowTransportDomainAdapter.hpp",
	"solvers/cpu/tests/test_native_tet_ale_flow_transport_domain_adapter.cpp",
	"solvers/cpu/tests/test_native_tet_ale_species_graph_1d.cpp",
	"scripts/test_native_tet_species_ftetwild.py",
	"include/SpeciesCoupling.hpp",
	"include/ZeroDSpeciesReservoir.hpp",
	"include/ZeroDFlowSpeciesCheckpoint.hpp",
	"include/SpeciesGraphCheckpointIdentity.hpp",
	"include/ZeroDSourceReservoirSpecies.hpp",
	"include/ZeroDSourceReservoirSpeciesDomainRuntime.hpp",
	"include/ZeroDTerminalRcrSpecies.hpp",
	"include/ZeroDTerminalRcrSpeciesDomainRuntime.hpp",
	"solvers/coupling/tests/test_zero_d_species_reservoir.cpp",
	"solvers/coupling/tests/test_zero_d_source_reservoir_species.cpp",
	"solvers/coupling/tests/test_zero_d_source_reservoir_species_runtime.cpp",
	"solvers/coupling/tests/test_zero_d_flow_species_checkpoint.cpp",
	"solvers/coupling/tests/test_zero_d_terminal_rcr_species.cpp",
	"solvers/coupling/tests/test_zero_d_terminal_rcr_species_runtime.cpp",
	"solvers/coupling/tests/test_zero_d_terminal_rcr_species_graph.cpp",
	"include/OneDFlowDomainAdapter.hpp",
	"solvers/one_d/tests/test_one_d_runtime.cpp",
)


def source_hash(root):
	digest = hashlib.sha256()
	for name in SOURCE_FILES:
		contents = (root/name).read_bytes()
		digest.update(name.encode()+b"\0"+str(len(contents)).encode()+b"\0"+contents)
	return digest.hexdigest()


def require(condition, message):
	if not condition:
		raise RuntimeError(message)


def main():
	root = Path(__file__).resolve().parents[1]
	card = json.loads((root/"benchmarks/t7_native_moving_species_evidence.json").read_text())
	require(card.get("schema_version") == 1
		and card.get("kind") == "native_moving_tetra_species_functional_reference_not_production_validation"
		and card.get("external_fem_framework") is False,
		"native moving species evidence is misclassified")
	require(card.get("source_files_sha256") == source_hash(root),
		"native moving species source identity is stale")
	require(card.get("test_command") == "make t7-native-moving-species-test"
		and card.get("two_tetra_test") is True
		and card.get("previous_and_current_mass_matrices") is True
		and card.get("relative_velocity_u_minus_w") is True
		and card.get("inflow_requires_concentration") is True
		and card.get("native_p2_flow_state_adapter") is True
		and card.get("inventory_balance_checked") is True,
		"native moving species functional coverage is incomplete")
	require(card.get("maximum_dense_vertices") == 128
		and card.get("mpi_petsc_scalar_runtime") is True
		and card.get("native_sparse_assembly") is True
		and card.get("native_species_diffusion_spatial_test_command") == "make t7-native-moving-species-spatial-convergence-test"
		and card.get("native_species_diffusion_spatial_mpi_ranks_tested") == [1, 2, 4]
		and card.get("native_species_diffusion_spatial_divisions") == [2, 4, 8]
		and card.get("native_species_diffusion_spatial_reference_is_be_resolvent") is True
		and card.get("native_species_advection_front_reference_is_be_resolvent") is True
		and card.get("native_species_advection_front_zero_diffusivity") is True
		and card.get("native_species_advection_front_monotone_nonnegative_and_conservative") is True
		and card.get("mpi_cell_and_face_ownership") is True
		and card.get("mpi_ranks_tested") == [1, 2, 4]
		and card.get("sparse_test_vertices") == 160
		and card.get("sparse_source_and_balance_checked") is True
		and card.get("multi_step_moving_inventory_checked") is True
		and card.get("labelled_outward_flux_checked") is True
		and card.get("labelled_relative_flow_checked") is True
		and card.get("graph_port_scalar_contract_checked") is True
		and card.get("simulation_graph_edge_contract_checked") is True
		and card.get("zero_flow_donor_concentration_checked") is True
		and card.get("native_ale_hydraulic_port_contract_checked") is True
		and card.get("native_ale_hydraulic_port_petsc_solve_checked") is True
		and card.get("native_ale_flow_to_species_petsc_checked") is True
		and card.get("native_ale_flow_to_species_mpi_ranks_tested") == [1, 2, 4]
		and card.get("native_staged_flow_transport_runtime_checked") is True
		and card.get("native_staged_rollback_mpi_ranks_tested") == [1, 2, 4]
		and card.get("replicated_geometry_and_global_solution") is True
		and card.get("graph_port_integration") is True
		and card.get("one_d_native_tet_one_d_graph_mpi_ranks_tested") == [1, 2, 4]
		and card.get("graph_moving_backflow_rejection_and_retry_checked") is True
		and card.get("graph_rank_local_precommit_failure_agreement_checked") is True
		and card.get("channel_two_step_moving_graph_mpi_ranks_tested") == [1, 2, 4]
		and card.get("channel_two_step_strong_graph_mpi_ranks_tested") == [1, 2, 4]
		and card.get("rcr_species_graph_test_command") == "make t7-native-ale-species-graph-rcr-test"
		and card.get("one_d_native_tet_rcr_explicit_mpi_ranks_tested") == [1, 2, 4]
		and card.get("one_d_native_tet_rcr_fixed_mpi_ranks_tested") == [1, 2, 4]
		and card.get("rcr_rank_zero_precommit_rollback_checked") is True
		and card.get("rcr_distal_reverse_explicit_mpi_ranks_tested") == [1, 2, 4]
		and card.get("rcr_distal_reverse_fixed_mpi_ranks_tested") == [1, 2, 4]
		and card.get("rcr_distal_reverse_full_graph_verified") is True
		and card.get("zero_d_source_graph_test_command") == "make t7-native-ale-species-graph-zero-d-source-test"
		and card.get("zero_d_source_native_tet_rcr_explicit_mpi_ranks_tested") == [1, 2, 4]
		and card.get("zero_d_source_rank_zero_precommit_rollback_checked") is True
		and card.get("zero_d_source_native_tet_rcr_fixed_mpi_ranks_tested") == [1, 2, 4]
		and card.get("zero_d_source_fixed_second_step_iterations") == 12
		and card.get("zero_d_species_paired_file_checkpoint_test_command") == "make coupling-zero-d-flow-species-checkpoint-test"
		and card.get("zero_d_species_paired_file_checkpoint_verified") is True
		and card.get("native_staged_flow_species_paired_restore_mpi_ranks_tested") == [1, 2, 4]
		and card.get("native_staged_paired_next_step_absolute_tolerance") == 1e-12
		and card.get("channel_strong_pressure_relative_tolerance") == 1e-4
		and card.get("one_d_open_loop_hydraulic_rollback_rearm_checked") is True
		and card.get("native_accepted_step_checkpoint") is True
		and card.get("graph_checkpoint_integration") is True
		and card.get("native_source_fem_rcr_graph_checkpoint_test_command") == "make t7-native-ale-species-graph-checkpoint-test"
		and card.get("native_source_fem_rcr_graph_checkpoint_mpi_ranks_tested") == [1, 2, 4]
		and card.get("native_source_fem_rcr_graph_checkpoint_shards") == 5
		and card.get("graph_checkpoint_compatibility_binds_graph_models_controls") is True
		and card.get("graph_checkpoint_changed_model_and_controls_rejected") is True
		and card.get("graph_checkpoint_changed_species_unit_rejected") is True
		and card.get("ftetwild_species_graph_test_command") == "make t7-native-tet-species-ftetwild-test"
		and card.get("ftetwild_species_graph_same_mesh_mpi_test_command") == "make t7-native-tet-species-ftetwild-same-mesh-mpi-test"
		and card.get("ftetwild_species_graph_independent_mpi_ranks_tested") == [1, 2, 4]
		and card.get("ftetwild_species_graph_same_mesh_rank_comparison") is True
		and card.get("ftetwild_species_graph_same_mesh_field_relative_l2_gate") == 1e-10
		and 0 <= card.get("ftetwild_species_graph_same_mesh_field_relative_l2_observed", -1) <= 1e-10
		and 0 <= card.get("ftetwild_species_graph_same_mesh_five_variant_max_relative_l2_observed", -1) <= 1e-8
		and 0 <= card.get("ftetwild_species_graph_same_mesh_six_variant_max_relative_l2_observed", -1) <= 1e-8
		and card.get("ftetwild_species_graph_same_mesh_material_pair") is True
		and card.get("ftetwild_species_graph_same_mesh_four_variants") is True
		and card.get("ftetwild_species_graph_isolated_diffusion_rms_threshold_mol_m3") == 1e-8
		and card.get("ftetwild_species_graph_isolated_source_inventory_checked") is True
		and card.get("ftetwild_species_graph_material_diffusivity_m2_s") == 1e-5
		and card.get("ftetwild_species_graph_material_source_mol_m3_s") == 0.1
		and card.get("ftetwild_species_p1_vtu_output_mpi_ranks_tested") == [1, 2, 4]
		and card.get("ftetwild_species_p1_vtu_accepted_field_checked") is True
		and card.get("ftetwild_species_p1_vtu_reference_current_displacement_checked") is True
		and card.get("ftetwild_species_p1_vtu_moving_step_displacement_x_m") == -1e-5
		and card.get("ftetwild_species_p1_vtu_is_general_production_cli") is False
		and card.get("native_optional_monotone_graph_diffusion_tested") is True
		and card.get("native_optional_monotone_high_peclet_dense_and_mpi_ranks_tested") == [1, 2, 4]
		and card.get("native_optional_monotone_unstabilized_undershoot_mol_m3") < -0.1
		and card.get("native_optional_monotone_nonnegative_and_column_conservative") is True
		and card.get("native_monotone_staged_graph_fTetWild_front_mpi_ranks_tested") == [1, 2, 4]
		and card.get("native_monotone_staged_graph_fTetWild_front_nonnegative") is True
		and card.get("native_monotone_staged_graph_unstabilized_front_rejected") is True
		and card.get("native_monotone_ftetwild_moving_front_mpi_ranks_tested") == [1, 2, 4]
		and card.get("native_monotone_ftetwild_moving_front_steps") == 2
		and card.get("native_monotone_ftetwild_moving_front_displacement_x_m") == -1e-5
		and card.get("native_monotone_ftetwild_moving_front_nonnegative_and_conservative") is True
		and card.get("native_monotone_staged_checkpoint_mode_mismatch_rejected") is True
		and card.get("native_monotone_moving_front_graph_checkpoint_test_command") == "make t7-native-ale-species-monotone-front-checkpoint-test"
		and card.get("native_monotone_moving_front_graph_checkpoint_mpi_ranks_tested") == [1, 2, 4]
		and card.get("native_monotone_moving_front_precommit_rollback_checked") is True
		and card.get("native_monotone_moving_front_restart_fingerprint_fields") == 326
		and card.get("native_monotone_moving_front_two_rank_minimum_mol_m3", -1) >= 0
		and abs(card.get("native_monotone_moving_front_two_rank_global_residual_mol", 1)) <= 1e-8
		and card.get("high_peclet_positivity_verified") is False
		and card.get("physiological_validation_passed") is False,
		"native moving species limitations are missing or overclaimed")
	for mode,minimum_order in (("galerkin",1.5),("monotone",1.2)):
		errors=card.get(f"native_species_diffusion_spatial_{mode}_relative_l2",[])
		orders=card.get(f"native_species_diffusion_spatial_{mode}_orders",[])
		require(len(errors)==3 and len(orders)==2
			and all(isinstance(value,(int,float)) and math.isfinite(value)
				and value>0 for value in errors+orders)
			and all(order>minimum_order for order in orders)
			and all(abs(math.log2(errors[index]/errors[index+1])-orders[index])<1e-10
				for index in range(2)),
			f"native species {mode} spatial convergence evidence is invalid")
	require(all(monotone>galerkin for monotone,galerkin in zip(
		card["native_species_diffusion_spatial_monotone_relative_l2"],
		card["native_species_diffusion_spatial_galerkin_relative_l2"])),
		"native species monotone accuracy tradeoff evidence is missing")
	for mode,minimum_order in (("galerkin",1.4),("monotone",0.4)):
		errors=card.get(f"native_species_advection_front_{mode}_relative_l2",[])
		orders=card.get(f"native_species_advection_front_{mode}_orders",[])
		require(len(errors)==3 and len(orders)==2
			and all(isinstance(value,(int,float)) and math.isfinite(value)
				and value>0 for value in errors+orders)
			and all(order>minimum_order for order in orders)
			and all(abs(math.log2(errors[index]/errors[index+1])-orders[index])<1e-10
				for index in range(2)),
			f"native species {mode} pure-advection convergence evidence is invalid")
	require(card.get("native_species_advection_front_galerkin_coarse_minimum_mol_m3",0)<-0.01
		and len(card.get("native_species_advection_front_monotone_minimum_mol_m3",[]))==3
		and all(value>=0 for value in card["native_species_advection_front_monotone_minimum_mol_m3"])
		and all(monotone>galerkin for monotone,galerkin in zip(
			card["native_species_advection_front_monotone_relative_l2"],
			card["native_species_advection_front_galerkin_relative_l2"])),
		"native species pure-advection positivity/accuracy tradeoff evidence is missing")
	require(card.get("native_species_advection_ale_translation_mpi_ranks_tested")==[1,2,4]
		and card.get("native_species_advection_ale_translation_divisions")==[2,4,8]
		and card.get("native_species_advection_ale_translation_modes")==["galerkin","monotone"]
		and card.get("native_species_advection_ale_translation_displacement_m")==0.02
		and card.get("native_species_advection_ale_translation_gate")==1e-10
		and all(isinstance(card.get(key),(int,float))
			and math.isfinite(card[key]) and 0<=card[key]<=card["native_species_advection_ale_translation_gate"]
			for key in (
				"native_species_advection_ale_translation_max_field_difference_mol_m3",
				"native_species_advection_ale_translation_max_inventory_difference_mol",
				"native_species_advection_ale_translation_max_flux_difference_mol_s",
				"native_species_advection_ale_translation_max_relative_l2_difference")),
		"native species ALE translation-invariance evidence is missing or invalid")
	require(card.get("native_species_affine_ale_expansion_mpi_ranks_tested")==[1,2,4]
		and card.get("native_species_affine_ale_expansion_divisions")==[2,4,8]
		and card.get("native_species_affine_ale_expansion_modes")==["galerkin","monotone"]
		and card.get("native_species_affine_ale_expansion_stretch")==1.2
		and abs(card.get("native_species_affine_ale_expansion_expected_concentration_mol_m3",0)-5./3.)<1e-15
		and card.get("native_species_affine_ale_expansion_expected_inventory_mol")==2.
		and card.get("native_species_affine_ale_expansion_concentration_gate")==1e-10
		and isinstance(card.get("native_species_affine_ale_expansion_max_concentration_error_mol_m3"),(int,float))
		and math.isfinite(card["native_species_affine_ale_expansion_max_concentration_error_mol_m3"])
		and 0<=card["native_species_affine_ale_expansion_max_concentration_error_mol_m3"]<=1e-10
		and card.get("native_species_affine_ale_expansion_inventory_and_flux_checked") is True,
		"native species affine ALE expansion evidence is missing or invalid")
	print("T7 native moving-tetra species evidence: PASS (five-shard functional graph restart; generalized validation/stabilization open)")


if __name__ == "__main__":
	main()
