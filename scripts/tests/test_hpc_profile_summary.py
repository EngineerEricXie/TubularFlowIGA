import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hpc_inventory import digest
from hpc_profile_summary import summarize


class ProfileSummaryTests(unittest.TestCase):
    def fixture(self, root, rank, ranks=2, total=10.0, failed=False, unscoped=3.0):
        directory = root / f"rank-{rank}"
        directory.mkdir()
        profile = {"schema_version": 1, "rank": rank, "ranks": ranks, "status": 0,
                   "elapsed_s": total, "unscoped_s": unscoped,
                   "phases": {"assembly": {"inclusive_s": 7.0, "exclusive_s": 5.0, "calls": 1},
                              "communication": {"inclusive_s": 2.0, "exclusive_s": 2.0, "calls": 1}}}
        log = directory / "stdout.log"
        log.write_text("hpc_profile " + json.dumps(profile) + "\n")
        record = {"schema_version": 1, "kind": "hpc_rank_run", "rank": rank, "ranks": ranks,
                  "returncode": int(failed), "timed_out": False, "status": "process_passed",
                  "wall_s": total+1.0, "resource": {"peak_rss_bytes": 100+rank},
                  "logs": {"stdout.log": digest(log)}}
        (directory / "run.json").write_text(json.dumps(record))

    def test_nested_timers_are_not_added_as_exclusive(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.fixture(root, 0)
            self.fixture(root, 1, total=12.0, unscoped=5.0)
            result = summarize(root)
            self.assertEqual(result["phases"]["assembly"]["exclusive_s"]["per_rank"], [5.0, 5.0])
            self.assertEqual(result["application_elapsed_s"]["max"], 12.0)
            self.assertEqual(result["sum_individual_peak_rss_bytes"], 201)

    def test_incomplete_and_failed_runs_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.fixture(root, 0)
            with self.assertRaises(ValueError):
                summarize(root)
            self.fixture(root, 1, failed=True)
            with self.assertRaises(ValueError):
                summarize(root)

    def test_overlap_and_nonfinite_timing_are_rejected(self):
        for total in (9.0, float("nan")):
            with self.subTest(total=total), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                self.fixture(root, 0, ranks=1, total=total)
                with self.assertRaises(ValueError):
                    summarize(root)

    def test_modified_logs_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.fixture(root, 0, ranks=1)
            with (root / "rank-0/stdout.log").open("a") as output:
                output.write("modified\n")
            with self.assertRaises(ValueError):
                summarize(root)


if __name__ == "__main__":
    unittest.main()
