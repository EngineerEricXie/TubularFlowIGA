#!/usr/bin/env python3
import hashlib
import json
import math
from pathlib import Path
import sys

def source_hash(root):
	paths=[root/"solvers/cpu/include/NativeTetFem.hpp",
		root/"solvers/cpu/include/NativeTetAleTransient.hpp",
		root/"solvers/cpu/include/NativeTetAleDenseRuntime.hpp",
		root/"solvers/cpu/tests/test_native_tet_fixed_temporal_convergence.cpp"]
	digest=hashlib.sha256()
	for path in paths:
		label=str(path.relative_to(root));contents=path.read_bytes()
		digest.update(label.encode()+b"\0"+str(len(contents)).encode()+b"\0"+contents)
	return digest.hexdigest()

def require(condition,message):
	if not condition:raise RuntimeError(message)

def main():
	root=Path(__file__).resolve().parents[1]
	data=json.loads((root/"benchmarks/t2_native_temporal_evidence.json").read_text(encoding="utf-8"))
	require(data.get("schema_version")==1,"unsupported T2 temporal evidence schema")
	require(data.get("kind")=="native_fixed_domain_temporal_manufactured_validation",
		"T2 temporal evidence kind changed")
	require(data.get("external_fem_framework") is False,"external FEM cannot count as native evidence")
	require(data.get("source_files_sha256")==source_hash(root),"T2 temporal source hash is stale")
	require(data.get("classification")=="fixed_geometry_temporal_validation_not_spatial_or_physiological_validation",
		"T2 temporal evidence is overclaimed")
	convergence=data["convergence"];dt=convergence["time_steps_s"]
	errors=convergence["absolute_pressure_gradient_errors_pa_m"];orders=convergence["observed_orders"]
	require(len(dt)==len(errors)==4 and len(orders)==3 and convergence.get("status")=="passed",
		"T2 temporal refinement levels are incomplete")
	lower,upper=convergence["required_order_interval"]
	for level in range(1,len(dt)):
		require(math.isclose(dt[level-1]/dt[level],2.0),"T2 temporal dt is not halved")
		require(errors[level]<errors[level-1],"T2 temporal error does not decrease")
		computed=math.log(errors[level-1]/errors[level],2.0)
		require(math.isclose(computed,orders[level-1],rel_tol=2e-6),"T2 temporal order is stale")
		require(lower<computed<upper,"T2 backward-Euler order failed its gate")
	print("T2 native fixed-domain temporal evidence: PASS")
	return 0

if __name__=="__main__":
	try:raise SystemExit(main())
	except (KeyError,OSError,RuntimeError,TypeError,ValueError) as error:
		print(f"T2 native fixed-domain temporal evidence: ERROR: {error}",file=sys.stderr);raise SystemExit(2)
