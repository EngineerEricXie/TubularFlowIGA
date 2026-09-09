import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hpc_inventory import digest
from hpc_profile_summary import summarize
from hpc_serial_matrix import aggregate_serial, cpu_reference, device_allocation, native_diagnostics


class SerialMatrixTests(unittest.TestCase):
    def test_immersed_solve_only_duplicate_and_nonfinite_fd_rejected(self):
        fd = [f"aneurysm zero-state centered-FD {direction} output-block {block} relative-defect=1e-12"
              for direction in ("velocity", "pressure", "controller") for block in range(3)]
        final = ["aneurysm static solve final-residual=1e-15", "aneurysm conservation open-normalized-balance=1e-12"]
        native_diagnostics("immersed", "\n".join(fd + final))
        for lines in (final, fd[:-1] + [fd[0]] + final, fd + fd[:1] + final,
                      [fd[0].replace("1e-12", "nan")] + fd[1:] + final):
            with self.subTest(lines=lines), self.assertRaises(ValueError):
                native_diagnostics("immersed", "\n".join(lines))

    def test_fsi_missing_convergence_or_iteration_is_rejected(self):
        history = "fsi_iteration=0 residual_rms_m=1e-8 threshold_m=1e-7\n"
        final = "compliant_channel_fsi converged=true iterations=1 center_displacement_m=1e-5"
        native_diagnostics("fsi", history + final)
        for log in (final, history, history + final + "\n" + final,
                    history.replace("fsi_iteration=0", "fsi_iteration=1") + final,
                    history + final.replace("1e-5", "inf")):
            with self.subTest(log=log), self.assertRaises(ValueError):
                native_diagnostics("fsi", log)

    def test_allocation_requires_scope_peak_and_no_live_buffers(self):
        line = "cuda_allocations scope=project_device_buffers requested_peak_bytes=1024 requested_live_bytes=0"
        self.assertEqual(device_allocation(line)["requested_peak_bytes"], 1024)
        for log in ("", line + "\n" + line, line.replace("bytes=1024", "bytes=0"),
                    line.replace("live_bytes=0", "live_bytes=4"), line.replace("project_device_buffers", "whole_device")):
            with self.subTest(log=log), self.assertRaises(ValueError):
                device_allocation(log)

    def records(self):
        return [{"repetition": i, "directory": str(i), "status": "accepted", "launch": {"wall_s": v + 1},
                 "profile": {"ranks": 1, "process_wall_s": {"max": v}, "application_elapsed_s": {"max": v - 1},
                             "peak_rss_bytes": {"max": 100}, "phases": {"assembly": {"exclusive_s": {"max": v - 2}}}},
                 "device_allocation": {"requested_peak_bytes": 1024}}
                for i, v in enumerate([1000, 9, 10, 11])]

    def test_first_run_excluded_from_host_and_device_statistics(self):
        result = aggregate_serial(self.records(), 3, True)
        self.assertEqual(result["process_wall_s"]["median"], 10)
        self.assertEqual(result["project_device_buffer_peak_bytes"]["count"], 3)
        self.assertEqual(result["phases_exclusive_s"]["assembly"]["values"], [7, 8, 9])

    def test_failed_missing_duplicate_and_parallel_records_rejected(self):
        records = self.records()
        failed = copy.deepcopy(records)
        failed[-1]["status"] = "failed"
        parallel = copy.deepcopy(records)
        parallel[-1]["profile"]["ranks"] = 2
        for group in (records[:-1], records[:-1] + records[:1], failed, parallel):
            with self.subTest(group=group), self.assertRaises(ValueError):
                aggregate_serial(group, 3, False)

    def reference(self, root):
        first = root / "first"
        rank = first / "rank-0"
        rank.mkdir(parents=True)
        profile = {"schema_version": 1, "rank": 0, "ranks": 1, "status": 0,
                   "elapsed_s": 1.0, "unscoped_s": 1.0, "phases": {}}
        (rank / "stdout.log").write_text("hpc_profile " + json.dumps(profile) + "\n")
        record = {"kind": "hpc_rank_run", "schema_version": 1, "rank": 0, "ranks": 1,
                  "returncode": 0, "timed_out": False, "status": "process_passed", "wall_s": 2,
                  "resource": {"peak_rss_bytes": 100}, "logs": {"stdout.log": digest(rank / "stdout.log")}}
        (rank / "run.json").write_text(json.dumps(record))
        fields = {}
        for name in ("field.txt", "field.txt.pressure"):
            (first / name).write_text("1 2 3\n")
            fields[name] = {"candidate_sha256": digest(first / name)}
        artifact = root / "input"
        artifact.write_text("recorded input")
        metadata = {"kind": "hpc_cpu_matrix", "system": "flow", "input_and_binary_sha256": {str(artifact): digest(artifact)}}
        result = {"status": "field_comparisons_passed", "observations": [
            {"ranks": 1, "repetition": 0, "status": "field_comparison_passed", "directory": str(first),
             "profile": summarize(first), "fields": fields}]}
        (root / "matrix.json").write_text(json.dumps(metadata))
        (root / "summary.json").write_text(json.dumps(result))
        return first, artifact

    def test_changed_reference_fields_inputs_and_failed_matrix_rejected(self):
        for changed in ("field", "input", "status"):
            with self.subTest(changed=changed), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                first, artifact = self.reference(root)
                self.assertEqual(cpu_reference(root, "flow")[1], first)
                if changed == "field":
                    (first / "field.txt").write_text("9 9 9\n")
                elif changed == "input":
                    artifact.write_text("different input")
                else:
                    (root / "summary.json").write_text('{"status": "failed"}')
                with self.assertRaises(ValueError):
                    cpu_reference(root, "flow")


if __name__ == "__main__":
    unittest.main()
