#!/usr/bin/env python3
"""End-to-end and fail-closed tests for the initial Gmsh FEM mesh route."""

import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[2]
TOOL = ROOT / "scripts" / "surface_to_fem_volume.py"
SURFACE = ROOT / "examples" / "vascular_flow" / "immersed_aneurysm_chain" / "immersed" / "surface.vtp"
PREFLIGHT = ROOT / "solvers" / "cpu" / "surface_fem_preflight"


def run(*arguments):
	return subprocess.run([sys.executable, str(TOOL), *map(str, arguments)],
		cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def main():
	with tempfile.TemporaryDirectory(prefix="tubularflow-fem-") as directory:
		directory = Path(directory)
		mesh = directory / "volume.msh"
		manifest = directory / "volume.json"
		result = run(SURFACE, mesh, "--manifest", manifest, "--target-size-m", "0.12")
		if result.returncode != 0:
			raise RuntimeError(result.stdout+result.stderr)
		data = json.loads(manifest.read_text(encoding="utf-8"))
		if not mesh.is_file() or not data["gates"]["passed"]:
			raise RuntimeError("healthy surface did not produce a passing mesh")
		preflight = data["surface_preflight"]
		expected_labels = {
			"0": {"triangles": 64, "area_m2": preflight["boundary_labels"]["0"]["area_m2"]},
			"1": {"triangles": 8, "area_m2": preflight["boundary_labels"]["1"]["area_m2"]},
			"2": {"triangles": 8, "area_m2": preflight["boundary_labels"]["2"]["area_m2"]}}
		if preflight["boundary_labels"] != expected_labels:
			raise RuntimeError("input boundary-label inventory changed")
		if data["volume_mesh"]["boundary_labels"] != {
			label: values["triangles"] for label, values in expected_labels.items()}:
			raise RuntimeError("output boundary labels changed")
		if not data["volume_mesh"]["boundary_matches_tetrahedra"]:
			raise RuntimeError("tetrahedral boundary consistency was not checked")
		if data["volume_mesh"]["volume_relative_error"] > 5.0e-3:
			raise RuntimeError("volume changed beyond the configured tolerance")
		if data["volume_mesh"]["boundary_surface_area_relative_error"] != 0.0 \
				or data["volume_mesh"]["maximum_source_surface_distance_m"] != 0.0 \
				or data["volume_mesh"]["boundary_remeshed"] \
				or data["volume_mesh"]["implicit_surface_repair"]:
			raise RuntimeError("exact-boundary Gmsh geometry-change evidence changed")
		for key in ("enclosed_volume_sphere", "surface_area_sphere"):
			if preflight["equivalent_radii_m"][key] <= 0.0:
				raise RuntimeError("equivalent-radius evidence is missing")
		if preflight["bounds_m"]["diagonal_m"] <= 0.0 \
				or preflight["discrete_curvature"]["integrated_absolute_mean_curvature_m"] <= 0.0:
			raise RuntimeError("bounds or discrete-curvature evidence is missing")
		if preflight["geometry_change_policy"] != {"implicit_repair": False,
				"implicit_remesh": False, "welding_is_explicitly_parameterized": True}:
			raise RuntimeError("surface geometry-change policy changed")
		change = data["geometry_change"]
		if change["classification"] != "exact_input_boundary_preservation" \
				or change["implicit_repair"] \
				or any(change[key] != 0.0 for key in (
					"relative_surface_area_change", "area_equivalent_radius_relative_change",
					"integrated_absolute_mean_curvature_relative_change",
					"maximum_source_surface_distance_m")):
			raise RuntimeError("exact-boundary geometry deltas changed")
		contract = data["geometry_contract"]
		if contract["reference_geometry"]["identity_sha256"] != preflight["canonical_sha256"] \
				or contract["current_geometry"]["identity_sha256"] != preflight["canonical_sha256"]:
			raise RuntimeError("reference/current geometry identity changed")
		regions = contract["regions"]
		if contract["length_unit"] != "m" or len(regions) != 1 \
				or {key: regions[0].get(key) for key in
					("dimension", "id", "physical_group", "role")} != {
					"dimension": 3, "id": "fluid", "physical_group": "fluid", "role": "fluid"} \
				or not isinstance(regions[0].get("physical_group_tag"), int) \
				or regions[0]["physical_group_tag"] <= 0:
			raise RuntimeError("geometry units or region role changed")
		if len(data["volume_mesh"]["sha256"]) != 64:
			raise RuntimeError("volume mesh hash is missing")
		if data["mesher"]["wall_s"] <= 0.0 or data["resources"]["pipeline_wall_s"] <= 0.0 \
				or data["resources"]["peak_rss_bytes"] <= 0:
			raise RuntimeError("Gmsh resource evidence is missing")
		canonical_msh = directory / "validated-surface.msh"
		canonical_manifest = directory / "validated-surface.json"
		preflight_result = subprocess.run([str(PREFLIGHT), str(SURFACE), str(canonical_msh),
			"--manifest", str(canonical_manifest), "--max-triangles", "10000000"],
			text=True, capture_output=True)
		if preflight_result.returncode != 0:
			raise RuntimeError(preflight_result.stdout+preflight_result.stderr)
		cached_mesh = directory / "cached-volume.msh"
		cached_manifest = directory / "cached-volume.json"
		cache_options = ("--validated-surface-msh", canonical_msh,
			"--validated-surface-manifest", canonical_manifest)
		cached_result = run(SURFACE, cached_mesh, "--manifest", cached_manifest,
			"--target-size-m", "0.12", *cache_options)
		if cached_result.returncode != 0 \
				or json.loads(cached_manifest.read_text())["surface_preflight_reused"] is not True:
			raise RuntimeError("hash-verified preflight reuse failed: "
				+cached_result.stdout+cached_result.stderr)
		wrong_options = run(SURFACE, directory/"wrong-options.msh",
			"--target-size-m", "0.12", "--default-boundary-id", "7", *cache_options)
		if wrong_options.returncode != 2 or "preflight options" not in wrong_options.stderr:
			raise RuntimeError("stale preflight options did not fail closed")
		unpaired_cache = run(SURFACE, directory/"unpaired.msh",
			"--target-size-m", "0.12", "--validated-surface-msh", canonical_msh)
		if unpaired_cache.returncode != 2 or "supplied together" not in unpaired_cache.stderr:
			raise RuntimeError("unpaired preflight artifacts did not fail closed")
		overwrite_cache = run(SURFACE, canonical_msh, "--target-size-m", "0.12",
			*cache_options)
		if overwrite_cache.returncode != 2 or "must not overwrite" not in overwrite_cache.stderr:
			raise RuntimeError("volume output could overwrite validated surface")
		changed_surface = directory / "changed-source.vtp"
		changed_surface.write_bytes(SURFACE.read_bytes()+b"\n")
		wrong_source = run(changed_surface, directory/"wrong-source.msh",
			"--target-size-m", "0.12", *cache_options)
		if wrong_source.returncode != 2 or "source hash" not in wrong_source.stderr:
			raise RuntimeError("stale preflight source did not fail closed")
		changed_mesh = directory / "changed-canonical.msh"
		shutil.copyfile(canonical_msh, changed_mesh)
		with changed_mesh.open("ab") as stream:
			stream.write(b"\n")
		wrong_mesh = run(SURFACE, directory/"wrong-canonical.msh",
			"--target-size-m", "0.12", "--validated-surface-msh", changed_mesh,
			"--validated-surface-manifest", canonical_manifest)
		if wrong_mesh.returncode != 2 or "canonical mesh hash" not in wrong_mesh.stderr:
			raise RuntimeError("tampered canonical mesh did not fail closed")

		stl = directory / "tetra-mm.stl"
		stl.write_text("""solid tetra
facet normal 0 0 0
outer loop
vertex 0 0 0
vertex 0 1000 0
vertex 1000 0 0
endloop
endfacet
facet normal 0 0 0
outer loop
vertex 0 0 0
vertex 1000 0 0
vertex 0 0 1000
endloop
endfacet
facet normal 0 0 0
outer loop
vertex 0 0 0
vertex 0 0 1000
vertex 0 1000 0
endloop
endfacet
facet normal 0 0 0
outer loop
vertex 1000 0 0
vertex 0 1000 0
vertex 0 0 1000
endloop
endfacet
endsolid tetra
""", encoding="ascii")
		stl_manifest = directory / "tetra.json"
		result = run(stl, directory/"tetra.msh", "--manifest", stl_manifest,
			"--target-size-m", "0.25", "--length-scale-to-m", "0.001",
			"--default-boundary-id", "7")
		if result.returncode != 0:
			raise RuntimeError(result.stdout+result.stderr)
		stl_data = json.loads(stl_manifest.read_text(encoding="utf-8"))
		stl_preflight = stl_data["surface_preflight"]
		if stl_preflight["canonical_vertices"] != 4 or stl_preflight["boundary_labels"]["7"]["triangles"] != 4:
			raise RuntimeError("STL welding or default boundary label changed")
		if abs(stl_preflight["volume_m3"]-1.0/6.0) > 1.0e-14:
			raise RuntimeError("STL millimetre-to-metre scaling changed")

		unsupported = directory / "unsupported.vtp"
		unsupported.write_text("""<?xml version="1.0"?>
<VTKFile type="PolyData" version="1.0" byte_order="LittleEndian"><PolyData><Piece NumberOfPoints="4" NumberOfPolys="4">
<Points><DataArray type="Float64" NumberOfComponents="3" format="binary">AAAA</DataArray></Points>
<Polys><DataArray type="Int32" Name="connectivity" format="ascii">0 1 2</DataArray>
<DataArray type="Int32" Name="offsets" format="ascii">3</DataArray></Polys>
<CellData><DataArray type="Int32" Name="boundary_id" format="ascii">0</DataArray></CellData>
</Piece></PolyData></VTKFile>\n""", encoding="utf-8")
		failure = run(unsupported, directory/"bad.msh", "--target-size-m", "0.1")
		if failure.returncode != 2 or "surface_fem_preflight: ERROR" not in failure.stderr:
			raise RuntimeError("unsupported VTP encoding did not fail closed")

		bad_orientation = directory / "bad-orientation.vtp"
		bad_orientation.write_text("""<?xml version="1.0"?>
<VTKFile type="PolyData" version="1.0" byte_order="LittleEndian"><PolyData><Piece NumberOfPoints="4" NumberOfPolys="4">
<Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">0 0 0 1 0 0 0 1 0 0 0 1</DataArray></Points>
<Polys><DataArray type="Int32" Name="connectivity" format="ascii">0 2 1 0 1 3 0 2 3 1 2 3</DataArray>
<DataArray type="Int32" Name="offsets" format="ascii">3 6 9 12</DataArray></Polys>
<CellData><DataArray type="Int32" Name="boundary_id" format="ascii">0 0 0 0</DataArray></CellData>
</Piece></PolyData></VTKFile>\n""", encoding="utf-8")
		failure = run(bad_orientation, directory/"bad-orientation.msh",
			"--target-size-m", "0.1")
		if failure.returncode != 2 or "opposite edge directions" not in failure.stderr:
			raise RuntimeError("inconsistently oriented VTP did not fail closed: "+failure.stderr)
	print("surface_to_fem_volume_test: PASS")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
