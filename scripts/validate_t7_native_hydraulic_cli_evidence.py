#!/usr/bin/env python3
"""Audit the native tetra hydraulic graph CLI functional evidence."""

import hashlib
import json
from pathlib import Path


SOURCE_FILES = (
	"solvers/cpu/src/native_tet_hydraulic_graph.cpp",
	"solvers/cpu/include/NativeTetHydraulicGraphCheckpoint.hpp",
	"solvers/cpu/include/NativeTetAleFlowDomainAdapter.hpp",
	"solvers/cpu/include/NativeTetAleFlowCheckpoint.hpp",
	"solvers/cpu/include/NativeTetAlePetscRuntime.hpp",
	"solvers/cpu/include/NativeTetHydraulicVisualization.hpp",
	"solvers/cpu/tests/test_native_tet_fem.cpp",
	"solvers/cpu/tests/data/native_tet_hydraulic_star.msh",
	"solvers/cpu/tests/data/native_tet_hydraulic_star.json",
	"scripts/test_native_tet_hydraulic_cli.py",
)


def source_hash(root):
	digest = hashlib.sha256()
	for name in SOURCE_FILES:
		contents = (root/name).read_bytes()
		digest.update(name.encode()+b"\0"+str(len(contents)).encode()+b"\0"+contents)
	return digest.hexdigest()


def main():
	root = Path(__file__).resolve().parents[1]
	card = json.loads((root/"benchmarks/t7_native_hydraulic_cli_evidence.json").read_text())
	if not (card.get("schema_version") == 1
		and card.get("kind") == "native_tetra_hydraulic_graph_cli_functional_smoke_only"
		and card.get("external_fem_framework") is False
		and card.get("source_files_sha256") == source_hash(root)
		and card.get("test_command") == "make t7-native-tet-hydraulic-cli-test"
		and card.get("mpi_ranks_tested") == [1, 2, 4]
		and card.get("external_gmsh41_mesh_read") is True
		and card.get("read_only_native_input_preflight_verified") is True
		and card.get("cross_process_cli_restart") is True
		and card.get("quadratic_vtu_full_p2_velocity_output") is True
		and card.get("reference_position_and_displacement_vtu_verified") is True
		and card.get("centroid_fluid_cauchy_stress_vtu_verified") is True
		and card.get("centroid_fluid_cauchy_stress_affine_expected_diagonal_pa") == [3.0, 7.0, -9.0]
		and card.get("restart_cauchy_stress_field_parity") is True
		and card.get("rigid_ale_displacement_mpi_ranks_tested") == [1, 2, 4]
		and card.get("restart_reference_displacement_parity") is True
		and card.get("restart_full_vtu_field_parity") is True
		and card.get("vtk_edge_permutation_verified") is True
		and card.get("published_vtu_no_overwrite_verified") is True
		and card.get("no_slip_channel_tetrahedra") == 24
		and card.get("no_slip_channel_mpi_ranks_tested") == [1, 2, 4]
		and card.get("moving_wall_velocity_checked") is True
		and card.get("no_slip_channel_cross_process_full_field_parity") is True
		and card.get("schema_v1_backward_compatible") is True
		and card.get("two_outlet_schema_v2_mpi_ranks_tested") == [1, 2, 4]
		and card.get("two_outlet_rcr_restart_full_field_parity") is True
		and card.get("two_outlet_fixed_point_mpi_ranks_tested") == [1, 2, 4]
		and card.get("two_outlet_reverse_explicit_mpi_ranks_tested") == [1, 2, 4]
		and card.get("two_outlet_reverse_fixed_point_mpi_ranks_tested") == [1, 2, 4]
		and card.get("fixed_point_failure_edge_diagnostics_verified") is True
		and card.get("real_y_fixed_point_validated") is False
		and card.get("changed_terminal_restart_rejected") is True
		and card.get("duplicate_outlet_label_rejected") is True
		and card.get("relative_flow_balance_gate_enabled") is True
		and card.get("missing_wall_label_rejected") is True
		and card.get("reverse_explicit_mpi_ranks_tested") == [1, 2, 4]
		and card.get("reverse_fixed_point_mpi_ranks_tested") == [1, 2, 4]
		and card.get("changed_case_rejected") is True
		and card.get("wrong_boundary_labels_rejected") is True
		and card.get("inverted_tetra_rejected") is True
		and card.get("generic_0d_species_supported") is False
		and card.get("strict_physiological_validation") is False):
		raise RuntimeError("native tetra hydraulic CLI evidence is stale or overclaimed")
	print("T7 native tetra hydraulic graph CLI evidence: PASS (functional smoke only)")


if __name__ == "__main__":
	main()
