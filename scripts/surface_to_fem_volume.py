#!/usr/bin/env python3
"""Create a labelled first-order tetrahedral Gmsh mesh from a closed surface.

The project C++ surface reader performs the canonical geometry preflight and
exports a labelled surface mesh; this driver only adds and validates the volume.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import tempfile
import time


def reject(message):
	raise ValueError(message)


def file_sha256(path):
	digest = hashlib.sha256()
	with path.open("rb") as stream:
		for block in iter(lambda: stream.read(1024*1024), b""):
			digest.update(block)
	return digest.hexdigest()


def verify_preflight_artifacts(preflight, input_path, canonical_msh, args):
	if not isinstance(preflight, dict):
		reject("surface preflight manifest must be an object")
	checks = preflight.get("checks", {})
	source = preflight.get("source", {})
	if preflight.get("schema_version") != 1 or preflight.get("kind") != "canonical_closed_surface" \
			or not isinstance(checks, dict) or not isinstance(source, dict) \
			or any(checks.get(name) is not True for name in (
				"closed", "oriented", "manifold", "connected",
				"self_intersection_free", "positive_volume")):
		reject("cached surface preflight lacks passing geometry checks")
	if source.get("sha256") != file_sha256(input_path):
		reject("surface preflight source hash does not match input")
	if preflight.get("canonical_msh_sha256") != file_sha256(canonical_msh):
		reject("surface preflight canonical mesh hash does not match cached mesh")
	if (preflight.get("boundary_array") != args.boundary_array
			or preflight.get("default_boundary_id") != args.default_boundary_id
			or preflight.get("max_triangles") != args.max_triangles
			or preflight.get("applied_length_scale_to_m") != args.length_scale_to_m
			or preflight.get("applied_weld_tolerance_m") != args.weld_tolerance_m):
		reject("surface preflight options do not match requested conversion")
	if (not isinstance(preflight.get("boundary_labels"), dict)
			or not preflight["boundary_labels"]
			or not isinstance(preflight.get("triangles"), int)
			or preflight["triangles"] <= 0
			or len(preflight.get("canonical_sha256", "")) != 64):
		reject("surface preflight geometry inventory is incomplete")


def tetra_determinant(a, b, c, d):
	ax, ay, az = (b[index]-a[index] for index in range(3))
	bx, by, bz = (c[index]-a[index] for index in range(3))
	cx, cy, cz = (d[index]-a[index] for index in range(3))
	return ax*(by*cz-bz*cy)-ay*(bx*cz-bz*cx)+az*(bx*cy-by*cx)


def tetra_volume(a, b, c, d):
	return abs(tetra_determinant(a, b, c, d))/6.0


def vector(left, right):
	return [right[index]-left[index] for index in range(3)]


def normalized_jacobian(first, second, third):
	determinant = (first[0]*(second[1]*third[2]-second[2]*third[1])
		-first[1]*(second[0]*third[2]-second[2]*third[0])
		+first[2]*(second[0]*third[1]-second[1]*third[0]))
	denominator = math.sqrt(sum(value*value for value in first)
		*sum(value*value for value in second)*sum(value*value for value in third))
	return abs(determinant)/denominator if denominator > 0.0 else 0.0


def tetra_minimum_scaled_jacobian(a, b, c, d):
	return min(
		normalized_jacobian(vector(a, b), vector(a, c), vector(a, d)),
		normalized_jacobian(vector(b, a), vector(b, d), vector(b, c)),
		normalized_jacobian(vector(c, d), vector(c, a), vector(c, b)),
		normalized_jacobian(vector(d, c), vector(d, b), vector(d, a)))


def main():
	start_wall = time.monotonic()
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("input_surface")
	parser.add_argument("output_msh")
	parser.add_argument("--manifest")
	parser.add_argument("--boundary-array", default="boundary_id")
	parser.add_argument("--surface-preflight", default=str(
		Path(__file__).resolve().parents[1]/"solvers"/"cpu"/"surface_fem_preflight"))
	parser.add_argument("--length-scale-to-m", type=float, default=1.0)
	parser.add_argument("--weld-tolerance-m", type=float, default=0.0)
	parser.add_argument("--default-boundary-id", type=int, default=0)
	parser.add_argument("--max-triangles", type=int, default=10000000)
	parser.add_argument("--validated-surface-msh", type=Path,
		help="canonical surface MSH from a prior successful surface_fem_preflight")
	parser.add_argument("--validated-surface-manifest", type=Path,
		help="matching preflight JSON with source and canonical MSH hashes")
	parser.add_argument("--target-size-m", type=float, required=True)
	parser.add_argument("--minimum-scaled-jacobian", type=float, default=1.0e-3)
	parser.add_argument("--maximum-volume-relative-error", type=float, default=5.0e-3)
	args = parser.parse_args()
	input_argument = Path(args.input_surface)
	input_path = input_argument.resolve()
	output_path = Path(args.output_msh).resolve()
	manifest_path = Path(args.manifest).resolve() if args.manifest else output_path.with_suffix(".json")
	if bool(args.validated_surface_msh) != bool(args.validated_surface_manifest):
		return reject("validated surface MSH and manifest must be supplied together")
	if args.validated_surface_msh and (output_path == args.validated_surface_msh.resolve()
			or manifest_path == args.validated_surface_manifest.resolve()):
		return reject("volume outputs must not overwrite validated surface artifacts")
	if not math.isfinite(args.target_size_m) or args.target_size_m <= 0.0:
		return reject("target size must be finite and positive")
	if not 0.0 < args.minimum_scaled_jacobian <= 1.0:
		return reject("minimum scaled Jacobian must be in (0,1]")
	if not 0.0 <= args.maximum_volume_relative_error < 1.0:
		return reject("maximum volume relative error must be in [0,1)")
	preflight_tool = Path(args.surface_preflight).resolve()
	if not preflight_tool.is_file():
		reject(f"surface preflight tool is missing: {preflight_tool}; run make cpu")

	try:
		import gmsh
	except ImportError as error:
		reject(f"Gmsh Python module is unavailable: {error}")
	gmsh.initialize(sys.argv[:1])
	try:
		gmsh.option.setNumber("General.Terminal", 1)
		gmsh.option.setNumber("General.Verbosity", 2)
		gmsh.option.setNumber("Mesh.MeshSizeMin", args.target_size_m)
		gmsh.option.setNumber("Mesh.MeshSizeMax", args.target_size_m)
		gmsh.option.setNumber("Mesh.ElementOrder", 1)
		# Preserve the validated anatomical surface exactly. Only the empty volume
		# may be meshed; remeshing boundary entities would change the geometry.
		gmsh.option.setNumber("Mesh.MeshOnlyEmpty", 1)
		with tempfile.TemporaryDirectory(prefix="tubularflow-surface-") as directory:
			canonical_msh = Path(directory)/"canonical-surface.msh"
			preflight_path = (args.validated_surface_manifest.resolve()
				if args.validated_surface_manifest else Path(directory)/"surface-preflight.json")
			if args.validated_surface_msh:
				shutil.copyfile(args.validated_surface_msh.resolve(), canonical_msh)
			else:
				command = [str(preflight_tool), str(input_path), str(canonical_msh),
					"--manifest", str(preflight_path), "--boundary-array", args.boundary_array,
					"--length-scale-to-m", str(args.length_scale_to_m),
					"--weld-tolerance-m", str(args.weld_tolerance_m),
					"--default-boundary-id", str(args.default_boundary_id),
					"--max-triangles", str(args.max_triangles)]
				completed = subprocess.run(command, text=True, stdout=subprocess.PIPE,
					stderr=subprocess.PIPE)
				if completed.returncode != 0:
					reject((completed.stderr or completed.stdout).strip())
			try:
				preflight = json.loads(preflight_path.read_text(encoding="utf-8"))
			except (OSError, json.JSONDecodeError) as error:
				reject(f"cannot read surface preflight manifest: {error}")
			verify_preflight_artifacts(preflight, input_path, canonical_msh, args)
			mesher_start = time.monotonic()
			gmsh.open(str(canonical_msh))
		expected_labels = {int(label): value for label, value in preflight["boundary_labels"].items()}
		entity_for_label = {}
		physical_group_for_label = {}
		for dimension, group in gmsh.model.getPhysicalGroups(2):
			name = gmsh.model.getPhysicalName(dimension, group)
			if not name.startswith("boundary_label_"):
				reject(f"unexpected canonical surface physical group {name!r}")
			try:
				label = int(name[len("boundary_label_"):])
			except ValueError:
				reject(f"invalid canonical boundary group {name!r}")
			entities = list(gmsh.model.getEntitiesForPhysicalGroup(2, group))
			if len(entities) != 1 or label in entity_for_label:
				reject(f"boundary label {label} does not map to exactly one surface entity")
			entity_for_label[label] = int(entities[0])
			physical_group_for_label[label] = int(group)
		if set(entity_for_label) != set(expected_labels):
			reject("canonical Gmsh surface groups do not match preflight boundary labels")
		gmsh.model.mesh.createTopology(False, True)
		surfaces = [tag for dimension, tag in gmsh.model.getEntities(2)]
		if not surfaces:
			reject("Gmsh did not retain a surface entity")
		loop = gmsh.model.geo.addSurfaceLoop(surfaces)
		volume_entity = gmsh.model.geo.addVolume([loop])
		gmsh.model.geo.synchronize()
		for label, entity in entity_for_label.items():
			if (2, entity) not in gmsh.model.getEntities(2):
				reject(f"Gmsh topology discarded boundary label {label}")
		volume_group = gmsh.model.addPhysicalGroup(3, [volume_entity])
		gmsh.model.setPhysicalName(3, volume_group, "fluid")
		gmsh.model.mesh.generate(3)
		types, element_tags, element_nodes = gmsh.model.mesh.getElements(3)
		if types != [4] or len(element_tags) != 1:
			reject(f"expected only first-order tetrahedra, got Gmsh element types {types}")
		tags = list(element_tags[0])
		connectivity = list(element_nodes[0])
		if not tags or len(connectivity) != 4*len(tags):
			reject("Gmsh produced no tetrahedra or invalid connectivity")
		all_node_tags, flat_coordinates, _ = gmsh.model.mesh.getNodes()
		coordinates = {int(tag): flat_coordinates[3*index:3*index+3]
			for index, tag in enumerate(all_node_tags)}
		output_label_counts = {}
		surface_facets = set()
		for label, entity in entity_for_label.items():
			surface_types, surface_tags, surface_nodes = gmsh.model.mesh.getElements(2, entity)
			if list(surface_types) != [2] or len(surface_tags) != 1:
				reject(f"boundary label {label} is not represented by first-order triangles")
			count = len(surface_tags[0])
			expected_count = int(expected_labels[label]["triangles"])
			if count != expected_count:
				reject(f"boundary label {label} triangle count changed from "
					f"{expected_count} to {count}")
			output_label_counts[str(label)] = count
			nodes = list(surface_nodes[0])
			for index in range(count):
				facet = tuple(sorted(int(nodes[3*index+corner]) for corner in range(3)))
				if facet in surface_facets:
					reject("generated boundary contains a duplicate triangle")
				surface_facets.add(facet)
		facet_uses = {}
		for index in range(len(tags)):
			tetra = [int(connectivity[4*index+corner]) for corner in range(4)]
			for local in ((0,1,2), (0,1,3), (0,2,3), (1,2,3)):
				facet = tuple(sorted(tetra[corner] for corner in local))
				facet_uses[facet] = facet_uses.get(facet, 0)+1
		tetra_boundary = {facet for facet, uses in facet_uses.items() if uses == 1}
		if surface_facets != tetra_boundary:
			reject("labelled surface triangles do not exactly equal the tetrahedral volume boundary")
		qualities = [tetra_minimum_scaled_jacobian(*(coordinates[int(connectivity[4*index+corner])]
			for corner in range(4))) for index in range(len(tags))]
		determinants = [tetra_determinant(*(coordinates[int(connectivity[4*index+corner])]
			for corner in range(4))) for index in range(len(tags))]
		minimum_determinant = min(determinants)
		if not math.isfinite(minimum_determinant) or minimum_determinant <= 0.0:
			reject(f"tetrahedral mesh has a non-positive Jacobian determinant {minimum_determinant}")
		minimum_quality = min(qualities)
		if not math.isfinite(minimum_quality) or minimum_quality < args.minimum_scaled_jacobian:
			reject(f"minimum tetrahedron scaled Jacobian {minimum_quality} is below "
				f"{args.minimum_scaled_jacobian}")
		volume_m3 = sum(tetra_volume(*(coordinates[int(connectivity[4*index+corner])]
			for corner in range(4))) for index in range(len(tags)))
		surface_volume_m3 = float(preflight["volume_m3"])
		if not math.isfinite(surface_volume_m3) or surface_volume_m3 <= 0.0:
			reject("input surface encloses zero or nonfinite volume")
		volume_relative_error = abs(volume_m3-surface_volume_m3)/surface_volume_m3
		if volume_relative_error > args.maximum_volume_relative_error:
			reject(f"tetrahedral volume relative error {volume_relative_error} exceeds "
				f"{args.maximum_volume_relative_error}")
		output_path.parent.mkdir(parents=True, exist_ok=True)
		gmsh.write(str(output_path))
		mesher_wall_s = time.monotonic()-mesher_start
		volume_mesh_sha256 = hashlib.sha256(output_path.read_bytes()).hexdigest()
		geometry_contract = {
			"schema_version": 1,
			"route": "surface_to_fem_volume",
			"coordinate_system": "cartesian",
			"length_unit": "m",
			"reference_geometry": {
				"identity_sha256": preflight["canonical_sha256"],
				"source_sha256": preflight["source"]["sha256"],
				"source_length_scale_to_m": preflight["applied_length_scale_to_m"]},
			"current_geometry": {"kind": "same_as_reference",
				"identity_sha256": preflight["canonical_sha256"]},
			"stable_ids": {"surface_vertices": "canonical_msh_node_tag",
				"surface_triangles": "canonical_msh_element_tag",
				"volume_nodes": "output_msh_node_tag",
				"volume_elements": "output_msh_element_tag"},
			"regions": [{"id": "fluid", "role": "fluid", "dimension": 3,
				"physical_group": "fluid", "physical_group_tag": volume_group}],
			"boundaries": [{"label": label, "physical_group": f"boundary_label_{label}",
				"physical_group_tag": physical_group_for_label[label],
				"semantic_role": "explicitly_configured_elsewhere"}
				for label in sorted(entity_for_label)]}
		peak_rss_bytes = max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
			resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)*1024
		manifest = {
			"schema_version": 1,
			"route": "surface_to_fem_volume_gmsh",
			"surface_preflight": preflight,
			"surface_preflight_reused": bool(args.validated_surface_msh),
			"geometry_contract": geometry_contract,
			"mesher": {"name": "Gmsh", "version": gmsh.option.getString("General.Version"),
				"target_size_m": args.target_size_m,
				"wall_s": mesher_wall_s},
			"resources": {"pipeline_wall_s": time.monotonic()-start_wall,
				"peak_rss_bytes": peak_rss_bytes},
			"geometry_change": {
				"classification": "exact_input_boundary_preservation",
				"implicit_repair": False,
				"relative_volume_change": volume_relative_error,
				"relative_surface_area_change": 0.0,
				"volume_equivalent_radius_relative_change": abs(
					(volume_m3/surface_volume_m3)**(1.0/3.0)-1.0),
				"area_equivalent_radius_relative_change": 0.0,
				"integrated_absolute_mean_curvature_relative_change": 0.0,
				"maximum_source_surface_distance_m": 0.0},
			"volume_mesh": {"format": "msh4", "element_type": "tetrahedron_p1",
				"sha256": volume_mesh_sha256,
				"nodes": len(coordinates), "elements": len(tags), "volume_m3": volume_m3,
				"volume_relative_error": volume_relative_error,
				"minimum_determinant_m3": minimum_determinant,
				"minimum_scaled_jacobian": minimum_quality,
				"boundary_surface_area_m2": preflight["area_m2"],
				"boundary_surface_area_relative_error": 0.0,
				"maximum_source_surface_distance_m": 0.0,
				"boundary_remeshed": False,
				"implicit_surface_repair": False,
				"boundary_labels": output_label_counts,
				"boundary_physical_group_tags": {str(label): physical_group_for_label[label]
					for label in sorted(physical_group_for_label)},
				"boundary_matches_tetrahedra": True},
			"gates": {"required_minimum_scaled_jacobian": args.minimum_scaled_jacobian,
				"maximum_volume_relative_error": args.maximum_volume_relative_error,
				"passed": True}
		}
		manifest_path.parent.mkdir(parents=True, exist_ok=True)
		manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True)+"\n", encoding="utf-8")
	finally:
		gmsh.finalize()
	print(f"surface_to_fem_volume: PASS elements={len(tags)} minSJ={minimum_quality:.8g} "
		f"volume_error={volume_relative_error:.3g}")
	return 0


if __name__ == "__main__":
	try:
		raise SystemExit(main())
	except (OSError, ValueError, RuntimeError) as error:
		print(f"surface_to_fem_volume: ERROR: {error}", file=sys.stderr)
		raise SystemExit(2)
