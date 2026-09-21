"""Canonical hashes for the frozen T2 contract and its mesh-affecting subset."""

import hashlib
import json


def exact_contract_sha256(contents):
	return hashlib.sha256(contents).hexdigest()


def mesh_contract_sha256(contract):
	projection = {
		"case_id": contract["case_id"],
		"geometry": contract["geometry"],
		"boundary_labels": contract["geometry"]["boundary_labels"],
		"mesh_levels_target_size_m": contract["discretization_contract"][
			"mesh_levels_target_size_m"]}
	contents = json.dumps(projection, sort_keys=True, separators=(",", ":"),
		ensure_ascii=False).encode("utf-8")
	return hashlib.sha256(contents).hexdigest()
