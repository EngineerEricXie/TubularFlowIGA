import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hpc_openmp_matrix import aggregate_openmp, publish


class OpenmpMatrixTests(unittest.TestCase):
    def records(self):
        records = []
        for threads, values in [(1, [1000, 99, 100, 101]), (4, [1000, 49, 50, 51])]:
            for repetition, elapsed in enumerate(values):
                records.append(dict(threads=threads, repetition=repetition, directory=str((threads, repetition)),
                                    status="accepted", field_comparison={"passed": True},
                                    launch={"wall_s": elapsed},
                                    profile={"ranks": 1, "process_wall_s": {"max": elapsed},
                                             "application_elapsed_s": {"max": elapsed},
                                             "peak_rss_bytes": {"max": 100 if threads == 1 else 120},
                                             "phases": {"assembly": {"exclusive_s": {"max": elapsed - 1}}}}))
        return records

    def test_first_excluded_and_thresholds(self):
        modes = aggregate_openmp(self.records(), [1, 4], 3, 1.25)
        self.assertEqual(modes["4"]["observed_speedup"], 2)
        self.assertEqual(modes["4"]["parallel_efficiency"], .5)
        self.assertEqual(modes["4"]["peak_rss_ratio_to_one_thread"], 1.2)
        self.assertTrue(modes["4"]["speedup_claim_allowed"])
        self.assertFalse(modes["1"]["speedup_claim_allowed"])
        records = self.records()
        records[-1]["profile"]["peak_rss_bytes"]["max"] = 126
        self.assertFalse(aggregate_openmp(records, [1, 4], 3, 1.25)["4"]["speedup_claim_allowed"])
        records = self.records()
        for record in records:
            if record["threads"] == 4:
                record["launch"]["wall_s"] = 95
        self.assertFalse(aggregate_openmp(records, [1, 4], 3, 1.25)["4"]["speedup_claim_allowed"])

    def test_missing_duplicate_failed_and_zero_records(self):
        original = self.records()
        cases = [original[:-1], original[:-1] + original[-2:-1]]
        for field in ["status", "field_comparison"]:
            records = copy.deepcopy(original)
            records[-1][field] = "failed" if field == "status" else {"passed": False}
            cases.append(records)
        zero = copy.deepcopy(original)
        for record in zero:
            record["launch"]["wall_s"] = 0
        cases.append(zero)
        for records in cases:
            with self.subTest(records=records), self.assertRaises(ValueError):
                aggregate_openmp(records, [1, 4], 3, 1.25)

    def test_progress_publication_replaces_only_complete_json(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "matrix.json"
            publish(path, {"status": "running"})
            publish(path, {"status": "accepted"})
            self.assertEqual(json.loads(path.read_text()), {"status": "accepted"})
            self.assertFalse(path.with_suffix(".json.tmp").exists())


if __name__ == "__main__":
    unittest.main()
