#!/usr/bin/env python3
"""Audit the bounded native tetra Darcy single-physics evidence."""

import hashlib
import json
from pathlib import Path


SOURCE_FILES = (
	"solvers/cpu/include/NativeTetDarcyPetsc.hpp",
	"solvers/cpu/include/NativeTetRt0Flux.hpp",
	"solvers/cpu/include/NativeTetMovingSpeciesTransport.hpp",
	"solvers/cpu/include/NativeTetMovingSpeciesPetscRuntime.hpp",
	"solvers/cpu/include/NativeTetVesselTissueSourceMap.hpp",
	"solvers/cpu/include/NativeTetMatchingInterfaceSource.hpp",
	"solvers/cpu/include/NativeTetBoundaryFlow.hpp",
	"solvers/cpu/include/NativeTetDarcyVisualization.hpp",
	"solvers/cpu/include/NativeTetFem.hpp",
	"solvers/cpu/src/native_tet_darcy.cpp",
	"solvers/cpu/tests/test_native_tet_darcy_petsc.cpp",
	"scripts/test_native_tet_darcy_ftetwild.py",
	"scripts/test_native_tet_darcy_cli.py",
	"scripts/run_native_tet_workflow.py",
	"scripts/test_native_tet_darcy_workflow.py",
)


def source_hash(root):
	digest = hashlib.sha256()
	for name in SOURCE_FILES:
		contents = (root/name).read_bytes()
		digest.update(name.encode()+b"\0"+str(len(contents)).encode()+b"\0"+contents)
	return digest.hexdigest()


def main():
	root = Path(__file__).resolve().parents[1]
	card = json.loads((root/"benchmarks/t6_native_tet_darcy_evidence.json").read_text())
	if not (card.get("schema_version") == 1
			and card.get("kind") == "native_tetra_p1_darcy_single_physics_functional"
			and card.get("source_files_sha256") == source_hash(root)
			and card.get("test_command") == "make t6-native-tet-darcy-test"
			and card.get("external_fem_framework") is False
			and card.get("mpi_ranks_tested") == [1, 2, 4]
			and card.get("affine_pressure_flux_exact") is True
			and card.get("quadratic_source_cell_interior_errors") == [0.0263008, 0.0065752]
			and card.get("outward_neumann_sign_checked") is True
			and card.get("conservative_interior_face_continuity_checked") is True
			and card.get("conservative_cell_balance_with_source_checked") is True
			and card.get("conservative_neumann_and_natural_boundary_checked") is True
			and card.get("conservative_face_flow_vtu_and_summary_checked") is True
			and card.get("rt0_hdiv_postprocessing_face_trace_divergence_checked") is True
			and card.get("mixed_darcy_discretization") is False
			and card.get("darcy_rt0_to_native_species_affine_test") is True
			and card.get("explicit_signed_vessel_tissue_cell_source_map_test") is True
			and card.get("exact_matching_facet_source_map_functional_test") is True
			and card.get("mobility_sensitivity_checked") is True
			and card.get("newly_generated_ftetwild_pipe_tested") is True
			and card.get("invalid_boundary_and_material_rejected") is True
			and card.get("standalone_cli_mpi_ranks_tested") == [1, 2, 4]
			and card.get("cli_rank_owned_pressure_flux_vtu") is True
			and card.get("cli_cell_material_and_source_override_checked") is True
			and card.get("cli_case_mesh_sha256_and_no_overwrite") is True
			and card.get("cli_ftetwild_pipe_tested") is True
			and card.get("standalone_cli_available") is True
			and card.get("workflow_volume_ftetwild_gmsh_routes") is True
			and card.get("workflow_steady_restart_rejected") is True
			and card.get("coupled_vessel_tissue_perfusion") is False
			and card.get("physiological_validation_passed") is False):
		raise RuntimeError("native tetra Darcy evidence is stale or overclaimed")
	print("T6 native tetra Darcy evidence: PASS (Darcy/generic-interface functional; no organ perfusion)")


if __name__ == "__main__":
	main()
