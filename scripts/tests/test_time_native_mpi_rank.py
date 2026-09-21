"""Rank-scoped GNU time records, without MPI or patient data."""

import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1]/"time_native_mpi_rank.py"


class RankTimeTest(unittest.TestCase):
	def test_record_and_no_overwrite(self):
		with tempfile.TemporaryDirectory() as directory:
			env = dict(os.environ, OMPI_COMM_WORLD_RANK="0")
			command = [sys.executable, str(SCRIPT), "--metrics-dir", directory,
				"--", sys.executable, "-c", "pass"]
			self.assertEqual(subprocess.run(command, env=env, capture_output=True).returncode, 0)
			values = (Path(directory)/"rank0.time.txt").read_text().split()
			self.assertEqual(len(values), 2)
			self.assertGreater(int(values[1]), 0)
			self.assertNotEqual(subprocess.run(command, env=env, capture_output=True).returncode, 0)


if __name__ == "__main__":
	unittest.main()
