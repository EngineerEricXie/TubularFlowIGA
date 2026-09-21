"""Bounded stage timing and immutable-output checks, without patient data."""

from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_liver_roi_from_raw import run


class RawRoiStageTest(unittest.TestCase):
	def test_stage_metrics_and_log_no_overwrite(self):
		with tempfile.TemporaryDirectory() as directory:
			root = Path(directory)
			metrics = run("tiny", [sys.executable, "-c", "pass"], root)
			self.assertGreaterEqual(metrics["wall_seconds"], 0)
			self.assertGreater(metrics["gnu_time_maximum_rss_kb"], 0)
			with self.assertRaisesRegex(ValueError, "already exists"):
				run("tiny", [sys.executable, "-c", "pass"], root)


if __name__ == "__main__":
	unittest.main()
