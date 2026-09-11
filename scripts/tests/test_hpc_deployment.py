import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPTS = Path(__file__).resolve().parents[1]
ROOT = SCRIPTS.parent
sys.path.insert(0, str(SCRIPTS))


def load(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


build = load("hpc_build_manifest")
scaling = load("hpc_cross_node_scaling")
tiers = load("hpc_test_tiers")


class DeploymentTests(unittest.TestCase):
    def test_petsc_macros_and_default_cuda_architectures_are_recorded(self):
        with tempfile.TemporaryDirectory() as directory:
            header = Path(directory) / "petscconf.h"
            header.write_text("#define PETSC_HAVE_MUMPS 1\n#define IGNORED 1\n")
            self.assertEqual(build.macros(header)["values"], {"PETSC_HAVE_MUMPS": "1"})
        architectures = build.cuda_architectures()
        self.assertIn(architectures["source"], ("CUDA_ARCHS", "solvers/cuda/Makefile default"))
        self.assertTrue(architectures["values"])

    def test_tier_timeout_is_failed_and_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            result = tiers.run([sys.executable, "-c", "import time; time.sleep(5)"],
                               Path(directory), os.environ.copy(), 0.05, "timeout")
            self.assertEqual(result["returncode"], 124)
            self.assertTrue(result["timed_out"])
            self.assertEqual(result["status"], "failed")

    def test_scheduler_skip_has_machine_readable_reason(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "scheduled"
            environment = {key: value for key, value in os.environ.items()
                           if not key.startswith("SLURM_")}
            result = subprocess.run([sys.executable, str(SCRIPTS / "hpc_test_tiers.py"),
                                     "--tier", "scheduled", "--output-dir", str(output)],
                                    cwd=ROOT, env=environment, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            report = json.loads((output / "result.json").read_text())
            self.assertEqual(report["status"], "skipped")
            self.assertIn("two nodes", report["skip_reason"])

    def test_scaling_statistics_reject_incomplete_samples(self):
        with self.assertRaises(ValueError):
            scaling.samples([1.0, 2.0])
        result = scaling.samples([3.0, 1.0, 2.0])
        self.assertEqual(result["median"], 2.0)

    def test_slurm_wrappers_parse(self):
        for path in (ROOT / "solvers/cpu/slurm/multinode_graph.sbatch",
                     ROOT / "solvers/cpu/slurm/cross_node_scaling.sbatch",
                     ROOT / "solvers/cpu/slurm/cross_node_fsi.sbatch"):
            result = subprocess.run(["bash", "-n", str(path)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)

    def test_hash_tree_is_nonempty_and_immutable(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tree = root / "tree"
            tree.mkdir()
            (tree / "value").write_text("accepted")
            output = root / "hashes.json"
            result = subprocess.run([sys.executable, str(SCRIPTS / "hpc_hash_tree.py"),
                                     "--root", str(tree), "--output", str(output)],
                                    cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            hashes = json.loads(output.read_text())
            self.assertEqual(list(hashes), [str((tree / "value").resolve())])
            self.assertNotEqual(subprocess.run([sys.executable,
                str(SCRIPTS / "hpc_hash_tree.py"), "--root", str(tree),
                "--output", str(output)], capture_output=True, text=True).returncode, 0)


if __name__ == "__main__":
    unittest.main()
