from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hpc_compare_fields import compare


class FieldComparisonTests(unittest.TestCase):
    def check(self, left, right, node_ids=False):
        with tempfile.TemporaryDirectory() as directory:
            reference, candidate = Path(directory)/"ref", Path(directory)/"candidate"
            reference.write_text(left)
            candidate.write_text(right)
            return compare(reference, candidate, 1e-6, 1e-12, node_ids)

    def test_nonzero_reference_uses_relative_norm(self):
        self.assertTrue(self.check("3 4\n", "3.000001 4\n")["passed"])
        self.assertFalse(self.check("3 4\n", "3.1 4\n")["passed"])

    def test_zero_reference_uses_explicit_absolute_tolerance(self):
        self.assertTrue(self.check("0 0\n", "0 1e-13\n")["passed"])
        self.assertFalse(self.check("0 0\n", "0 1e-10\n")["passed"])

    def test_invalid_fields_are_rejected(self):
        for left, right in [("1\n", "nan\n"), ("1\n", "1 2\n"),
                            ("1\n2\n", "1\n"), ("", ""), ("1\n", "inf\n")]:
            with self.assertRaises(ValueError):
                self.check(left, right)

    def test_node_ids_are_checked_not_used_as_field_values(self):
        self.assertTrue(self.check("0 3\n1 4\n", "0 3\n1 4\n", True)["passed"])
        for text in ["1 3\n2 4\n", "0 3\n0 4\n"]:
            with self.assertRaises(ValueError):
                self.check("0 3\n1 4\n", text, True)


if __name__ == "__main__":
    unittest.main()
