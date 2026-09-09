import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hpc_inventory import digest
from hpc_reference_states import POLICY, compare_states, state_manifest, execution_configuration, assembly_diagnostics
from types import SimpleNamespace
from unittest.mock import patch


class ReferenceStateTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.policy = json.loads(POLICY.read_text())
        self.reference = self.make_state("reference")
        self.candidate = self.make_state("candidate")

    def make_state(self, name, fixture="fsi"):
        path = self.root / name
        path.mkdir()
        schema = self.policy["fixtures"][fixture]
        manifest = {"schema_version": 1, "kind": "fixture_reference_state", "fixture": fixture,
                    "native_gates_passed": True, "time_s": schema["time_s"], "step": schema["step"], "fields": []}
        for field, metadata in schema["fields"].items():
            # Absolute coordinates remain O(1), displacement O(1e-8). They
            # must never share a reference norm that masks a motion defect.
            value = 1e-8 if "displacement" in field else 1.0
            text = "".join(str(node) + " " + " ".join([format(value, ".17g")] * metadata["columns"]) + "\n"
                           for node in [9007199254740993, 9007199254740995])
            file = path / (field + ".txt")
            file.write_text(text)
            manifest["fields"].append({"name": field, "file": file.name, "rows": 2,
                                       **metadata, "sha256": digest(file)})
        (path / "manifest.json").write_text(json.dumps(manifest))
        return path

    def edit_manifest(self, action):
        file = self.candidate / "manifest.json"
        manifest = json.loads(file.read_text())
        action(manifest)
        file.write_text(json.dumps(manifest))

    def edit_field(self, name, action, directory=None, refresh_hash=True):
        directory = directory or self.candidate
        file = directory / (name + ".txt")
        file.write_text(action(file.read_text()))
        if refresh_hash:
            manifest_file = directory / "manifest.json"
            manifest = json.loads(manifest_file.read_text())
            next(f for f in manifest["fields"] if f["name"] == name)["sha256"] = digest(file)
            manifest_file.write_text(json.dumps(manifest))

    def test_complete_roundtrip_and_integer_ids(self):
        result = compare_states(self.reference, self.candidate, self.policy)
        self.assertTrue(result["passed"])
        self.assertEqual(len(result["fields"]), 9)
        self.assertTrue(all(f["absolute_l2"] == 0 for f in result["fields"].values()))

    def test_motion_not_hidden_by_coordinate_scale(self):
        self.edit_field("material_displacement", lambda _: "9007199254740993 1.1e-8 1e-8 1e-8\n9007199254740995 1e-8 1e-8 1e-8\n")
        result = compare_states(self.reference, self.candidate, self.policy)
        self.assertFalse(result["passed"])
        self.assertFalse(result["fields"]["material_displacement"]["passed"])
        self.assertTrue(result["fields"]["material_position"]["passed"])

    def test_zero_reference_uses_absolute_gate(self):
        for directory in [self.reference, self.candidate]:
            self.edit_field("fluid_pressure", lambda _: "9007199254740993 0\n9007199254740995 0\n", directory)
        self.assertTrue(compare_states(self.reference, self.candidate, self.policy)["passed"])
        self.edit_field("fluid_pressure", lambda _: "9007199254740993 2e-12\n9007199254740995 0\n")
        self.assertFalse(compare_states(self.reference, self.candidate, self.policy)["passed"])

    def test_id_mismatch_above_double_precision(self):
        self.edit_field("fluid_pressure", lambda text: text.replace("9007199254740993", "9007199254740992"))
        with self.assertRaisesRegex(ValueError, "IDs must match"):
            compare_states(self.reference, self.candidate, self.policy)

    def test_tampered_field_hash(self):
        self.edit_field("surface_force", lambda text: text + "99 0 0 0\n", refresh_hash=False)
        with self.assertRaisesRegex(ValueError, "hash"):
            state_manifest(self.candidate, self.policy)

    def test_truncated_or_nonfinite_field(self):
        self.edit_field("fluid_pressure", lambda _: "9007199254740993 1\n")
        with self.assertRaisesRegex(ValueError, "dimensions"):
            state_manifest(self.candidate, self.policy)
        self.edit_field("fluid_pressure", lambda _: "9007199254740993 nan\n9007199254740995 1\n")
        with self.assertRaisesRegex(ValueError, "nonfinite"):
            state_manifest(self.candidate, self.policy)

    def test_required_field_cannot_be_omitted(self):
        self.edit_manifest(lambda manifest: manifest["fields"].pop())
        with self.assertRaisesRegex(ValueError, "fields"):
            state_manifest(self.candidate, self.policy)

    def test_wrong_units_or_epoch(self):
        self.edit_manifest(lambda manifest: manifest["fields"][0].update(units="mm/s"))
        with self.assertRaisesRegex(ValueError, "schema"):
            state_manifest(self.candidate, self.policy)
        self.edit_manifest(lambda manifest: manifest.update(time_s=2.0))
        with self.assertRaisesRegex(ValueError, "epoch"):
            state_manifest(self.candidate, self.policy)

    def test_missing_completion_manifest(self):
        (self.candidate / "manifest.json").unlink()
        with self.assertRaises(FileNotFoundError):
            state_manifest(self.candidate, self.policy)

    def test_path_escape_and_false_native_gate(self):
        self.edit_manifest(lambda manifest: manifest["fields"][0].update(file="../field.txt"))
        with self.assertRaisesRegex(ValueError, "schema"):
            state_manifest(self.candidate, self.policy)
        self.edit_manifest(lambda manifest: manifest.update(native_gates_passed=False))
        with self.assertRaisesRegex(ValueError, "incomplete"):
            state_manifest(self.candidate, self.policy)

    def test_immersed_requires_controller(self):
        reference = self.make_state("immersed1", "immersed")
        candidate = self.make_state("immersed2", "immersed")
        self.assertTrue(compare_states(reference, candidate, self.policy)["passed"])
        with self.assertRaisesRegex(ValueError, "different fixture"):
            compare_states(reference, self.candidate, self.policy)





class ExecutionConfigurationTests(unittest.TestCase):
    def options(self, **changes):
        values = dict(cpu=0, cpus=None, threads=1, timeout=60)
        values.update(changes)
        return SimpleNamespace(**values)

    @patch("hpc_reference_states.os.sched_getaffinity", return_value={0, 2, 4, 6})
    def test_serial_and_multicore(self, affinity):
        cpus, env = execution_configuration(self.options())
        self.assertEqual(cpus, [0])
        self.assertEqual(env["IGA_ASSEMBLY_THREADS"], "1")
        cpus, env = execution_configuration(self.options(cpus=[0, 2, 4, 6], threads=4))
        self.assertEqual(cpus, [0, 2, 4, 6])
        self.assertEqual(env["OMP_THREAD_LIMIT"], "4")
        self.assertEqual(env["IGA_ASSEMBLY_BATCH_SIZE"], "4")
        self.assertEqual(env["OPENBLAS_NUM_THREADS"], "1")

    @patch("hpc_reference_states.os.sched_getaffinity", return_value={0, 2, 4, 6})
    def test_reject_unavailable_duplicate_and_oversubscribed(self, affinity):
        for changes in [dict(cpus=[1]), dict(cpus=[0, 0]), dict(threads=2),
                        dict(threads=0), dict(cpus=[]), dict(timeout=float("inf"))]:
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                execution_configuration(self.options(**changes))


class AssemblyDiagnosticsTests(unittest.TestCase):
    def test_observed_team_required(self):
        self.assertEqual(assembly_diagnostics("", 1), [])
        with self.assertRaisesRegex(ValueError, "lacks actual"):
            assembly_diagnostics("", 4)
        line = "element_assembly threads_requested=4 team_size=4 cells=27 batches=7 maximum_resident_items=4"
        self.assertEqual(assembly_diagnostics(line, 4)[0]["cells"], 27)
        for bad in [line.replace("team_size=4", "team_size=2"),
                    line.replace("maximum_resident_items=4", "maximum_resident_items=8"),
                    line + " junk"]:
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                assembly_diagnostics(bad, 4)


if __name__ == "__main__":
    unittest.main()
