"""Independent signed-flow ledger checks without patient data or MPI."""

import json
from pathlib import Path
import tempfile
import unittest

from scripts.validate_liver_roi_functional_evidence import read_facet_ledger


class FacetLedgerTest(unittest.TestCase):
	def test_signed_flows_and_owner_ids_are_preserved(self):
		ledger = {"schema_version": 1,
			"kind": "native_matching_facet_flow_functional_only",
			"flow_sign": "vessel_outward_positive_tissue_source_positive",
			"facets": [
				{"port_name": "solved_roi_functional", "vessel_cell_id": 7,
					"tissue_cell_id": 11, "vessel_outward_m3_s": 0.3,
					"triangle_m": [[0, 0, 0], [0, 1, 0], [0, 0, 1]]},
				{"port_name": "solved_roi_functional", "vessel_cell_id": 8,
					"tissue_cell_id": 12, "vessel_outward_m3_s": -0.1,
					"triangle_m": [[1, 0, 0], [1, 1, 0], [1, 0, 1]]}]}
		with tempfile.TemporaryDirectory() as directory:
			root = Path(directory)
			path = root/"interface_facets.json"
			path.write_text(json.dumps(ledger), encoding="utf-8")
			parsed = read_facet_ledger(root, 2, 0.2)
			self.assertEqual(sorted(value[2] for value in parsed.values()), [-0.1, 0.3])
			with self.assertRaisesRegex(ValueError, "does not close"):
				read_facet_ledger(root, 2, 0.3)
			ledger["facets"][1]["triangle_m"] = ledger["facets"][0]["triangle_m"]
			path.write_text(json.dumps(ledger), encoding="utf-8")
			with self.assertRaisesRegex(ValueError, "duplicate"):
				read_facet_ledger(root, 2, 0.2)


if __name__ == "__main__":
	unittest.main()
