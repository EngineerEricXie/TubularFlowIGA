#!/usr/bin/env python3
"""Audit source-bound, functional vessel–0D wall exchange evidence."""

import json
from pathlib import Path


def main():
	root = Path(__file__).resolve().parents[1]
	card = json.loads((root/"benchmarks/t6_native_wall_reservoir_evidence.json").read_text())
	if not (card.get("schema_version") == 1
			and card.get("kind") == "native_tetra_p1_to_fixed_volume_0d_wall_exchange_functional"
			and card.get("test_command") == "make t6-native-wall-reservoir-exchange-test"
			and card.get("external_fem_framework") is False
			and card.get("mpi_ranks_tested") == [1, 2, 4]
			and card.get("outward_positive_pair_balance") is True
			and card.get("forward_and_reverse_transfer") is True
			and card.get("rigid_ale_translation") is True
			and card.get("source_decay_combined_budget") is True
			and card.get("invalid_label_and_flow_rejected") is True
			and card.get("newly_generated_ftetwild_pipe_tested") is True
			and card.get("pipe_transfer_coefficient_sensitivity") is True
			and card.get("solves_tissue_storage") is True
			and card.get("solves_tissue_hydraulics") is False
			and card.get("single_region_cli_paired_file_restart") is True
			and card.get("multi_region_calculation_mpi_ranks_tested") == [1, 2, 4]
			and card.get("multi_region_opposite_direction_and_cross_response") is True
			and card.get("multi_region_cli_paired_file_restart") is True
			and card.get("physiological_validation_passed") is False):
		raise RuntimeError("native wall reservoir evidence is stale or overclaimed")
	print("T6 native vessel–0D wall reservoir exchange evidence: PASS (functional only)")


if __name__ == "__main__":
	main()
