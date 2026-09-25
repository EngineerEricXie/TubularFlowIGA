"""Fail-closed versioned summary requirements without patient files."""

import unittest

from scripts.validate_liver_roi_functional_evidence import require_v2_summary_bindings


class SummaryBindingsTest(unittest.TestCase):
	def test_v2_requires_audit_and_stage_metrics(self):
		for summary in ({"schema_version": 2},
				{"schema_version": 2, "source_face_audit_sha256": "abc"},
				{"schema_version": 2, "stage_metrics": {}}):
			with self.subTest(summary=summary):
				with self.assertRaisesRegex(ValueError, "requires"):
					require_v2_summary_bindings(summary)
		require_v2_summary_bindings({"schema_version": 2,
			"source_face_audit_sha256": "abc", "stage_metrics": {}})

	def test_v1_legacy_summary_compatible(self):
		require_v2_summary_bindings({"schema_version": 1})


if __name__ == "__main__":
	unittest.main()
