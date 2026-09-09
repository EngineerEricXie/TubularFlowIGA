import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "hpc_inventory.py"
SPEC = importlib.util.spec_from_file_location("hpc_inventory", SCRIPT)
inventory = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(inventory)


class InventoryTests(unittest.TestCase):
    def test_content_changes_and_missing_files_are_visible(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "input").write_bytes(b"original")
            first = inventory.file_records(root, ["input"])
            (root / "input").write_bytes(b"modified")
            second = inventory.file_records(root, ["input"])
            self.assertNotEqual(inventory.records_digest(first), inventory.records_digest(second))
            (root / "input").unlink()
            self.assertTrue(inventory.file_records(root, ["input"])[0]["missing"])

    def test_input_cannot_escape_repository(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "repo"
            root.mkdir()
            (root.parent / "outside").write_text("not a fixture")
            with self.assertRaises(ValueError):
                inventory.input_records(root, ["../outside"])

    def test_probe_failure_is_retained(self):
        result = inventory.probe(["/nonexistent/tubularflow-hpc-probe"])
        self.assertEqual(result["status"], "unavailable")
        self.assertIn("error", result)

    def test_petsc_prefix_must_identify_configured_build(self):
        with tempfile.TemporaryDirectory() as directory:
            self.assertEqual(inventory.petsc_configuration(directory)["status"], "unresolved")
            include = Path(directory) / "include"
            include.mkdir()
            (include / "petscconf.h").write_text("#define PETSC_HAVE_MUMPS 1\n")
            result = inventory.petsc_configuration(directory)
            self.assertEqual(result["features"]["PETSC_HAVE_MUMPS"], "1")
            self.assertIsNone(result["features"]["PETSC_USE_COMPLEX"])

    def test_catalog_references_exist(self):
        root = SCRIPT.parents[1]
        catalog = json.loads((root / "benchmarks/hpc_baselines.json").read_text())
        self.assertEqual(len({case["id"] for case in catalog["cases"]}), 4)
        for case in catalog["cases"]:
            self.assertTrue(inventory.input_records(root, case["source_inputs"]))
            for path in case["native_gate_sources"] + case["reference_evidence"]:
                self.assertTrue((root / path).is_file(), path)


if __name__ == "__main__":
    unittest.main()
