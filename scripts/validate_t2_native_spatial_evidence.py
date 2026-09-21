#!/usr/bin/env python3
import hashlib
import json
import math
from pathlib import Path
import sys

def source_hash(root):
	paths=[root/"solvers/cpu/include/NativeTetFem.hpp",
		root/"solvers/cpu/include/NativeTetAleDenseRuntime.hpp",
		root/"solvers/cpu/tests/test_native_tet_fixed_spatial_convergence.cpp"]
	digest=hashlib.sha256()
	for path in paths:
		label=str(path.relative_to(root));contents=path.read_bytes()
		digest.update(label.encode()+b"\0"+str(len(contents)).encode()+b"\0"+contents)
	return digest.hexdigest()

def require(condition,message):
	if not condition:raise RuntimeError(message)

def main():
	root=Path(__file__).resolve().parents[1]
	data=json.loads((root/"benchmarks/t2_native_spatial_evidence.json").read_text(encoding="utf-8"))
	require(data.get("schema_version")==1,"unsupported T2 spatial evidence schema")
	require(data.get("kind")=="native_fixed_domain_spatial_manufactured_validation",
		"T2 spatial evidence kind changed")
	require(data.get("external_fem_framework") is False,"external FEM cannot count as native evidence")
	require(data.get("source_files_sha256")==source_hash(root),"T2 spatial source hash is stale")
	require(data.get("classification")=="fixed_geometry_spatial_operator_validation_not_poiseuille_or_physiological_validation",
		"T2 spatial evidence is overclaimed")
	convergence=data["convergence"];n=convergence["subdivisions_per_axis"]
	errors=convergence["velocity_relative_l2_errors"];orders=convergence["observed_orders"]
	require(n==[2,3,4] and len(errors)==3 and len(orders)==2 and convergence.get("status")=="passed",
		"T2 spatial refinement levels are incomplete")
	minimum=convergence["minimum_observed_order"]
	for level in range(1,len(n)):
		require(errors[level]<errors[level-1],"T2 spatial error does not decrease")
		computed=math.log(errors[level-1]/errors[level])/math.log(n[level]/n[level-1])
		require(math.isclose(computed,orders[level-1],rel_tol=2e-5),"T2 spatial order is stale")
		require(computed>minimum,"T2 P2 spatial order failed its gate")
	print("T2 native fixed-domain spatial evidence: PASS")
	return 0

if __name__=="__main__":
	try:raise SystemExit(main())
	except (KeyError,OSError,RuntimeError,TypeError,ValueError) as error:
		print(f"T2 native fixed-domain spatial evidence: ERROR: {error}",file=sys.stderr);raise SystemExit(2)
