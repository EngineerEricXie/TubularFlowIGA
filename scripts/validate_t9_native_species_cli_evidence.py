#!/usr/bin/env python3
"""Audit the bounded external-mesh native species CLI evidence."""

import hashlib
import json
from pathlib import Path


SOURCE_FILES = (
	"solvers/cpu/src/native_tet_species_transport.cpp",
	"solvers/cpu/include/NativeTetSpeciesVisualization.hpp",
	"solvers/cpu/include/NativeTetMovingSpeciesCheckpoint.hpp",
	"solvers/cpu/include/NativeTetMovingSpeciesPetscRuntime.hpp",
	"solvers/cpu/include/NativeTetMovingSpeciesTransport.hpp",
	"solvers/cpu/include/NativeTetWallReservoirExchange.hpp",
	"solvers/cpu/tests/test_native_tet_moving_species_transport.cpp",
	"solvers/cpu/tests/test_native_tet_moving_species_petsc_runtime.cpp",
	"solvers/cpu/tests/data/native_tet_hydraulic_star.msh",
	"scripts/test_native_tet_species_cli.py",
	"scripts/test_native_tet_species_workflow.py",
)


def source_hash(root):
	digest = hashlib.sha256()
	for name in SOURCE_FILES:
		contents = (root/name).read_bytes()
		digest.update(name.encode()+b"\0"+str(len(contents)).encode()+b"\0"+contents)
	return digest.hexdigest()


def main():
	root = Path(__file__).resolve().parents[1]
	card = json.loads((root/"benchmarks/t9_native_species_cli_evidence.json").read_text())
	if not (card.get("schema_version") == 1
			and card.get("kind") == "native_tetra_prescribed_velocity_species_cli_functional_only"
			and card.get("external_fem_framework") is False
			and card.get("source_files_sha256") == source_hash(root)
			and card.get("test_command") == "make t9-native-tet-species-cli-test"
			and card.get("external_gmsh41_mesh") is True
			and card.get("newly_generated_ftetwild_pipe_tested") is True
			and card.get("ftetwild_pipe_mpi_ranks") == 2
			and card.get("fixture_mpi_ranks_tested") == [1, 2, 4]
			and card.get("native_p1_assembly_and_petsc_algebra") is True
			and card.get("same_parser_check_input") is True
			and card.get("rigid_ale_relative_velocity_checked") is True
			and card.get("source_response_checked") is True
			and card.get("first_order_decay_implicit_mass_matrix") is True
			and card.get("first_order_decay_analytic_mpi_ranks_tested") == [1, 2, 4]
			and card.get("first_order_decay_moving_restart_ranks") == [2]
			and card.get("reaction_sink_budget_checked") is True
			and card.get("negative_decay_rejected") is True
			and card.get("labelled_wall_exchange_mpi_ranks_tested") == [1, 2, 4]
			and card.get("wall_exchange_bidirectional_and_budget_checked") is True
			and card.get("wall_exchange_matrix_rhs_checked") is True
			and card.get("wall_exchange_restart_ranks") == [2]
			and card.get("wall_exchange_bad_label_parameter_and_flow_rejected") is True
			and card.get("external_tissue_storage_solved") is True
			and card.get("finite_reservoir_mpi_ranks_tested") == [1, 2, 4]
			and card.get("finite_reservoir_paired_restart_ranks") == [2]
			and card.get("finite_reservoir_checkpoint_corruption_rejected") is True
			and card.get("finite_reservoir_workflow_volume_restart") is True
			and card.get("multi_region_reservoir_mpi_ranks_tested") == [1, 2, 4]
			and card.get("multi_region_reservoir_paired_restart_ranks") == [2]
			and card.get("multi_region_reservoir_corrupt_sidecar_rejected") is True
			and card.get("multi_region_workflow_volume_restart") is True
			and card.get("multi_region_ftetwild_surface_workflow") is True
			and card.get("monotone_inflow_front_checked") is True
			and 0.1 < card.get("physical_diffusion_field_change_mol_m3", 0) < 0.12
			and card.get("low_diffusion_masked_by_monotone_numerical_diffusion") is True
			and card.get("missing_inflow_rejected") is True
			and card.get("wrong_label_and_unknown_key_rejected") is True
			and card.get("output_no_overwrite") is True
			and card.get("rank_owned_p1_vtu_all_vertices_cells_checked") is True
			and card.get("summary_case_mesh_sha256_checked") is True
			and card.get("cross_process_restart_mpi_ranks_tested") == [1, 2, 4]
			and card.get("cross_process_restart_field_absolute_gate_mol_m3") == 1e-10
			and card.get("cross_process_checkpoint_bitwise_parity_required") is False
			and card.get("cross_process_moving_ale_field_parity_ranks") == [2]
			and card.get("restart_wrong_case_and_ranks_rejected") is True
			and card.get("restart_missing_corrupt_and_orphan_rejected") is True
			and card.get("accepted_step_vtu_checkpoint_state_publication") is True
			and card.get("concurrent_writer_lock_enabled") is True
			and card.get("solves_fluid_equations") is False
			and card.get("graph_coupling_supported") is False
			and card.get("checkpoint_restart_supported") is True
			and card.get("coupled_graph_checkpoint_restart_supported") is False
			and card.get("physiological_validation_passed") is False):
		raise RuntimeError("native tetra species CLI evidence is stale or overclaimed")
	print("T9 native tetra prescribed-velocity species CLI evidence: PASS (functional only)")


if __name__ == "__main__":
	main()
