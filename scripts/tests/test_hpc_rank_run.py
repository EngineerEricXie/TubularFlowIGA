import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "hpc_rank_run.py"


class RankRunTests(unittest.TestCase):
    def launch(self, directory, code, timeout=10, extra_env=None):
        environment = {key: value for key, value in os.environ.items()
                       if not key.startswith(("OMPI_", "PMI_", "SLURM_"))}
        environment.update(extra_env or {})
        return subprocess.run([sys.executable, str(SCRIPT), "--output-dir", directory,
                               "--timeout", str(timeout), "--", sys.executable, "-c", code],
                              env=environment, capture_output=True, text=True, timeout=15)

    def test_success_preserves_literal_arguments_and_records_resources(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.launch(directory, "print('$(literal); `no shell`')")
            self.assertEqual(result.returncode, 0, result.stderr)
            root = Path(directory) / "rank-0"
            self.assertEqual((root / "stdout.log").read_text().strip(), "$(literal); `no shell`")
            report = json.loads((root / "run.json").read_text())
            self.assertEqual(report["status"], "process_passed")
            self.assertGreater(report["resource"]["peak_rss_bytes"], 0)
            self.assertNotEqual(self.launch(directory, "pass").returncode, 0)

    def test_failed_child_is_not_success(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.launch(directory, "raise SystemExit(7)")
            self.assertEqual(result.returncode, 7)
            report = json.loads((Path(directory) / "rank-0/run.json").read_text())
            self.assertEqual(report["status"], "failed")

    def test_timeout_is_recorded(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.launch(directory, "import time; time.sleep(20)", 0.1)
            self.assertEqual(result.returncode, 124)
            report = json.loads((Path(directory) / "rank-0/run.json").read_text())
            self.assertTrue(report["timed_out"])

    def test_mismatched_launcher_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.launch(directory, "pass", extra_env={"OMPI_COMM_WORLD_RANK": "0", "OMPI_COMM_WORLD_SIZE": "2"})
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse((Path(directory) / "rank-0").exists())


if __name__ == "__main__":
    unittest.main()
