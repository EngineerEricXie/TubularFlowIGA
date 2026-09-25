#!/usr/bin/env python3
"""Audit the native ALE tetra hydraulic source/RCR functional evidence."""

import hashlib
import json
from pathlib import Path


SOURCE_FILES = (
	"solvers/cpu/include/NativeTetAleFlowCheckpoint.hpp",
	"solvers/cpu/include/NativeTetAleFlowDomainAdapter.hpp",
	"solvers/cpu/include/NativeTetHydraulicGraphCheckpoint.hpp",
	"solvers/cpu/tests/test_native_tet_ale_flow_0d_graph.cpp",
	"scripts/run_t7_native_0d_cross_process.py",
	"solvers/cpu/include/NativeTetAleGraphPorts.hpp",
	"include/ZeroDFlowDomain.hpp",
	"include/PressureFlowComponentExecutor.hpp",
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
	card = json.loads((root/"benchmarks/t7_native_0d_graph_evidence.json").read_text())
	require(card.get("schema_version") == 1
		and card.get("kind") == "native_tetra_ale_0d_hydraulic_functional_not_physiological_validation"
		and card.get("external_fem_framework") is False
		and card.get("source_files_sha256") == source_hash(root),
		"native tetra ALE 0D graph source evidence is stale or misclassified")
	require(card.get("test_command") == "make t7-native-ale-flow-0d-graph-test"
		and card.get("source_rcr_explicit_mpi_ranks_tested") == [1, 2, 4]
		and card.get("source_rcr_strong_mpi_ranks_tested") == [1, 2, 4]
		and card.get("explicit_density_kg_m3") == 1000.0
		and card.get("strong_density_kg_m3") == 1.0
		and card.get("strong_pressure_relative_tolerance") == 1e-4
		and card.get("strong_flow_relative_tolerance") == 1e-6
		and card.get("strong_second_step_iterations") == 14
		and card.get("balanced_high_density_explicit_mpi_ranks_tested") == [1, 2, 4]
		and card.get("balanced_high_density_strong_mpi_ranks_tested") == [1, 2, 4]
		and card.get("balanced_high_density_kg_m3") == 1000.0
		and card.get("balanced_source_capacitance_m3_pa") == 1e-9
		and card.get("balanced_source_resistance_pa_s_m3") == 500.0
		and card.get("balanced_source_prescribed_flow_m3_s") == 0.1
		and card.get("balanced_rcr_proximal_resistance_pa_s_m3") == 10.0
		and card.get("balanced_rcr_distal_resistance_pa_s_m3") == 50.0
		and card.get("balanced_rcr_capacitance_m3_pa") == 1e-9
		and card.get("balanced_strong_second_step_iterations") == 19
		and card.get("balanced_strong_outlet_flow_m3_s") == 0.1
		and card.get("moving_second_step_checked") is True
		and card.get("source_and_rcr_volume_balances_checked") is True
		and card.get("rank_local_precommit_rollback_checked") is True
		and card.get("native_accepted_step_checkpoint_checked") is True
		and card.get("three_domain_in_memory_restart_parity_checked") is True
		and card.get("checkpoint_corruption_and_wrong_motion_rejected") is True
		and card.get("native_hydraulic_file_bundle_supported") is True
		and card.get("file_shard_corruption_rejected_without_publication") is True
		and card.get("graph_topology_and_iteration_controls_bound_to_file_identity") is True
		and card.get("cross_process_file_restart_mpi_ranks_tested") == [1, 2, 4]
		and card.get("cross_process_test_command")
			== "make t7-native-ale-flow-0d-cross-process-test"
		and card.get("existing_graph_runner_file_restart_integrated") is False
		and card.get("generic_0d_species_supported") is False
		and card.get("physiological_validation_passed") is False,
		"native tetra ALE 0D graph coverage is incomplete or overclaimed")
	print("T7 native tetra ALE source/RCR hydraulic evidence: PASS (species and physiology open)")


if __name__ == "__main__":
	main()
