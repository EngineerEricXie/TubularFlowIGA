#!/usr/bin/env python3
import json
import hashlib
from pathlib import Path
import sys

def source_hash(root):
	paths=[
		root/"include/FsiDomainRuntime.hpp",
		root/"solvers/cpu/include/StrongFluidStructureCoupling.hpp",
		root/"solvers/cpu/include/NativeTetMatchingFsiInterface.hpp",
		root/"solvers/cpu/include/NativeTetAleFsiRuntime.hpp",
		root/"solvers/cpu/include/NativeTetAleKinematics.hpp",
		root/"solvers/cpu/include/NativeTetSolidFsiRuntime.hpp",
		root/"solvers/cpu/include/NativeTetFsiCheckpoint.hpp",
		root/"solvers/cpu/tests/test_native_tet_matching_fsi_interface.cpp",
		root/"solvers/cpu/tests/test_native_tet_solid_fsi_runtime.cpp",
		root/"solvers/cpu/tests/test_native_tet_ale_solid_fsi.cpp",
		root/"solvers/cpu/tests/test_native_tet_compliant_channel_fsi.cpp",
	]
	digest=hashlib.sha256()
	for path in paths:
		label=str(path.relative_to(root));contents=path.read_bytes()
		digest.update(label.encode()+b"\0"+str(len(contents)).encode()+b"\0"+contents)
	return digest.hexdigest()

def main():
	root=Path(__file__).resolve().parents[1]
	try:data=json.loads((root/"benchmarks/t5_native_matching_interface_contract.json").read_text(encoding="utf-8"))
	except (OSError,json.JSONDecodeError) as error:print(f"T5 matching interface contract error: {error}",file=sys.stderr);return 1
	if data.get("schema_version")!=1 or data.get("contract_status")!="frozen" or data.get("external_coupling_framework") is not False:
		print("T5 matching interface contract error: schema/ownership changed",file=sys.stderr);return 1
	if data.get("source_files_sha256")!=source_hash(root):
		print("T5 matching interface contract error: source hash is stale; rerun tests and update evidence",file=sys.stderr);return 1
	if data.get("kinematics",{}).get("fluid_velocity_is_not_mesh_velocity_away_from_the_no_slip_interface") is not True:
		print("T5 matching interface contract error: velocity semantics changed",file=sys.stderr);return 1
	gates=data.get("conservation_gates",{})
	if gates.get("total_force_exact_for_constant_triangle_traction") is not True or gates.get("total_moment_exact_for_constant_triangle_traction") is not True or gates.get("maximum_relative_roundoff_error",1)>1e-12:
		print("T5 matching interface contract error: conservation gates changed",file=sys.stderr);return 1
	lifecycle=data.get("strong_coupling_lifecycle",{})
	if lifecycle.get("state_machine")!="shared FsiTrialLifecycle used by existing partitioned coordinator" or lifecycle.get("structure_adapter")!="NativeTetSolidFsiRuntime":
		print("T5 matching interface contract error: lifecycle ownership changed",file=sys.stderr);return 1
	if lifecycle.get("rejected_iteration_changes_committed_state") is not False or lifecycle.get("abort_after_prepare_changes_committed_state") is not False or lifecycle.get("existing_dynamic_aitken_coordinator_instantiated") is not True or lifecycle.get("current_partition_scope")!="single partition":
		print("T5 matching interface contract error: rollback or partition scope changed",file=sys.stderr);return 1
	if lifecycle.get("native_end_to_end_invariant_smoke_is_not_a_physics_benchmark") is not True or lifecycle.get("native_ale_adapter")!="NativeTetAleFsiRuntime dense small-case reference; production PETSc adapter remains separate":
		print("T5 matching interface contract error: native vertical-slice scope changed",file=sys.stderr);return 1
	if len(lifecycle.get("independent_acceptance_gates",[]))!=6 or "infinity" not in lifecycle.get("legacy_gate_compatibility",""):
		print("T5 matching interface contract error: independent gate contract changed",file=sys.stderr);return 1
	nonzero=lifecycle.get("native_nonzero_functional_smoke",{})
	if nonzero.get("classification")!="nonzero coupled operator and transaction smoke, not compliant-channel or added-mass validation" or nonzero.get("coupling_iterations",0)<2 or nonzero.get("committed_solid_displacement_l2_m",0)<=0 or nonzero.get("all_independent_gates_passed") is not True:
		print("T5 matching interface contract error: nonzero functional evidence is missing or overclaimed",file=sys.stderr);return 1
	channel=lifecycle.get("native_compliant_channel_functional_smoke",{})
	if channel.get("classification")!="idealized five-step functional smoke, not physical, added-mass, spatial, or temporal validation" or channel.get("accepted_steps")!=5 or channel.get("final_time_s")!=0.05 or channel.get("coupling_iterations",0)<2 or channel.get("maximum_committed_displacement_m",0)<=0 or channel.get("all_independent_gates_passed") is not True:
		print("T5 matching interface contract error: compliant-channel functional evidence is missing or overclaimed",file=sys.stderr);return 1
	checkpoint=lifecycle.get("paired_checkpoint",{})
	if checkpoint.get("scope")!="single-partition native ALE-fluid plus native solid accepted-step state" or checkpoint.get("checksum_corruption_rejected") is not True or checkpoint.get("uninterrupted_vs_restarted_next_step_identical_committed_identities") is not True or checkpoint.get("uninterrupted_vs_restarted_steps_2_through_5_identical_committed_identities") is not True:
		print("T5 matching interface contract error: paired checkpoint evidence is missing or overclaimed",file=sys.stderr);return 1
	if len(data.get("fail_closed",[]))<6 or len(data.get("deferred",[]))<5:
		print("T5 matching interface contract error: failure/deferred scope incomplete",file=sys.stderr);return 1
	print("T5 native matching interface contract: PASS")
	return 0
if __name__=="__main__":raise SystemExit(main())
