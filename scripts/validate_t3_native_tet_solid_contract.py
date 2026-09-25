#!/usr/bin/env python3
"""Validate the native tetrahedral hyperelastic-solid implementation card."""

import hashlib
import json
from pathlib import Path
import sys

def source_hash(root):
	paths=[root/"solvers/cpu/include/NativeTetHyperelasticSolid.hpp",
		root/"solvers/cpu/include/NativeTetSolidStaticSolver.hpp",
		root/"solvers/cpu/include/NativeTetSolidPrestress.hpp",
		root/"solvers/cpu/include/NativeTetMixedHyperelasticSolid.hpp",
		root/"solvers/cpu/include/NativeTetMixedSolidStaticSolver.hpp",
		root/"solvers/cpu/tests/test_native_tet_hyperelastic_solid.cpp",
		root/"solvers/cpu/tests/test_native_tet_solid_static_solver.cpp",
		root/"solvers/cpu/tests/test_native_tet_solid_prestress.cpp",
		root/"solvers/cpu/tests/test_native_tet_mixed_hyperelastic_solid.cpp",
		root/"solvers/cpu/tests/test_native_tet_mixed_solid_static_solver.cpp",
		root/"solvers/cpu/tests/test_native_tet_mixed_locking.cpp"]
	digest=hashlib.sha256()
	for path in paths:
		label=str(path.relative_to(root));contents=path.read_bytes()
		digest.update(label.encode()+b"\0"+str(len(contents)).encode()+b"\0"+contents)
	return digest.hexdigest()


def main():
	root=Path(__file__).resolve().parents[1];path=root/"benchmarks/t3_native_tet_solid_contract.json"
	try:data=json.loads(path.read_text(encoding="utf-8"))
	except (OSError,json.JSONDecodeError) as error:
		print(f"T3-B contract error: {error}",file=sys.stderr);return 1
	if data.get("schema_version")!=1 or data.get("contract_status")!="frozen" \
			or data.get("backend")!="native_cpp_total_lagrangian_tetrahedral_solid":
		print("T3-B contract error: schema/backend changed",file=sys.stderr);return 1
	if data.get("external_discretization_or_assembly_framework") is not False:
		print("T3-B contract error: external FEM framework became allowed",file=sys.stderr);return 1
	if data.get("source_files_sha256")!=source_hash(root):
		print("T3-B contract error: source hash is stale",file=sys.stderr);return 1
	material=data.get("material",{})
	if "separate stabilized mixed u-p formulation" not in material.get("near_incompressible_policy",""):
		print("T3-B contract error: locking policy is missing",file=sys.stderr);return 1
	mixed=data.get("near_incompressible_mixed_formulation",{})
	if mixed.get("status")!="verified_small_dense_reference" \
			or mixed.get("poisson_ratio_maximum_tested")!=0.4999 \
			or mixed.get("locking_refinement_levels")!=6 \
			or mixed.get("maximum_poisson_sweep_relative_displacement_change",1)>=0.005 \
			or mixed.get("maximum_absolute_volume_change",1)>=2e-4:
		print("T3-B contract error: mixed locking evidence is incomplete",file=sys.stderr);return 1
	reference=data.get("reference_configuration_policy",{})
	if reference.get("loaded_image_geometry_is_not_stress_free_by_default") is not True \
			or reference.get("silent_zero_prestress_for_loaded_images") is not False \
			or reference.get("inverse_elastostatics")!="implemented for known material, dead nodal loads, and explicit zero-displacement supports":
		print("T3-B contract error: reference/prestress policy changed",file=sys.stderr);return 1
	prestress=data.get("inverse_prestress_scope",{})
	if prestress.get("loaded_geometry_reconstruction_status")!="passed" \
			or prestress.get("maximum_reconstruction_error_m",1)>2e-10 \
			or prestress.get("silent_unknown_load_inference") is not False:
		print("T3-B contract error: inverse prestress evidence is incomplete",file=sys.stderr);return 1
	if len(data.get("verification_cards",[]))<11 or len(data.get("deferred_explicit_options",[]))<5:
		print("T3-B contract error: verification or deferred scope is incomplete",file=sys.stderr);return 1
	runtime=data.get("global_runtime_scope",{})
	if runtime.get("assembly")!="repository_owned_dense_reference_assembly" \
			or runtime.get("not_production_scalable") is not True \
			or not isinstance(runtime.get("outputs"),list) or len(runtime["outputs"])<6:
		print("T3-B contract error: global runtime scope is incomplete",file=sys.stderr);return 1
	print("T3-B native tetrahedral solid contract: PASS")
	return 0


if __name__=="__main__":raise SystemExit(main())
