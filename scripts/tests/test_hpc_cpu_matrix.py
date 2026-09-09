import copy
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hpc_cpu_matrix import aggregate, hybrid_execution, launch, samples, thread_plan


class CpuMatrixTests(unittest.TestCase):
    def test_fixed_core_plan_includes_pure_mpi_and_rejects_oversubscription(self):
        self.assertEqual(thread_plan([1, 2, 4], 4), {1: 4, 2: 2, 4: 1})
        self.assertEqual(thread_plan([1, 2, 4], None), {1: 1, 2: 1, 4: 1})
        for ranks, cores in (([1, 2], 4), ([1, 3, 4], 4), ([1, 2, 8], 4), ([0, 4], 4),
                             ([1, 1], None), ([], 4), ([1], 0)):
            with self.subTest(ranks=ranks, cores=cores), self.assertRaises(ValueError):
                thread_plan(ranks, cores)

    def hybrid_fixture(self, root):
        topology = root / "topology"
        for cpu in range(8):
            local = topology / f"cpu{cpu}/topology"
            local.mkdir(parents=True)
            (local / "physical_package_id").write_text("0")
            (local / "core_id").write_text(str(cpu // 2))
        for rank in range(2):
            local = root / f"rank-{rank}"
            local.mkdir()
            (local / "run.json").write_text(json.dumps(dict(
                rank=rank, ranks=2, returncode=0, timed_out=False, status="process_passed",
                affinity_cpus=list(range(4 * rank, 4 * rank + 4)))))
            (local / "stdout.log").write_text(
                f"body_fitted_element_assembly rank={rank} threads_requested=2 team_size=2 "
                "elements=27 batches=14 maximum_resident_items=2\n")
        return topology

    def test_actual_teams_and_physical_core_bindings_ignore_smt_siblings(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            topology = self.hybrid_fixture(root)
            result = hybrid_execution(root, 2, 2, topology)
            self.assertEqual(result[0]["physical_cores"], [(0, 0), (0, 1)])
            self.assertEqual(result[1]["physical_cores"], [(0, 2), (0, 3)])
            self.assertEqual(len(result[0]["assembly_reports"]), 1)

    def test_overlapping_rank_cores_and_actual_worker_mismatch_are_rejected(self):
        for failure in ("overlap", "one_core", "missing_report", "wrong_team", "oversized_batch", "timeout"):
            with self.subTest(failure=failure), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                topology = self.hybrid_fixture(root)
                local = root / "rank-1"
                run = json.loads((local / "run.json").read_text())
                log = (local / "stdout.log").read_text()
                if failure == "overlap":
                    run["affinity_cpus"] = [2, 3, 4, 5]
                elif failure == "one_core":
                    run["affinity_cpus"] = [4, 5]
                elif failure == "timeout":
                    run["timed_out"] = True
                elif failure == "missing_report":
                    log = ""
                elif failure == "wrong_team":
                    log = log.replace("team_size=2", "team_size=1")
                else:
                    log = log.replace("maximum_resident_items=2", "maximum_resident_items=8")
                (local / "run.json").write_text(json.dumps(run))
                (local / "stdout.log").write_text(log)
                with self.assertRaises(ValueError):
                    hybrid_execution(root, 2, 2, topology)

    def records(self):
        records = []
        for n in (1, 2):
            for repetition in range(4):
                duration = (1000 if repetition == 0 else 8 + repetition) / n
                records.append({"ranks": n, "repetition": repetition,
                                "status": "field_comparison_passed", "directory": f"np{n}-{repetition}",
                                "launch": {"wall_s": duration},
                                "profile": {"process_wall_s": {"max": duration - 1},
                                            "application_elapsed_s": {"max": duration - 2},
                                            "peak_rss_bytes": {"max": 100},
                                            "sum_individual_peak_rss_bytes": 100 * n,
                                            "phases": {"assembly": {"exclusive_s": {"max": 2.0 / n}}}}})
        return records

    def test_first_run_excluded_and_memory_scope_retained(self):
        result = aggregate(self.records(), [1, 2], 3)
        self.assertEqual(result["1"]["launcher_wall_s"]["values"], [9, 10, 11])
        self.assertEqual(result["2"]["observed_speedup"], 2)
        self.assertEqual(result["2"]["observed_parallel_efficiency"], 1)
        self.assertEqual(result["2"]["sum_peak_rss_ratio_to_one_rank"], 2)
        self.assertEqual(result["2"]["max_rank_peak_rss_bytes"]["median"], 100)

    def test_failed_missing_duplicate_and_unknown_observations_rejected(self):
        original = self.records()
        failed = copy.deepcopy(original)
        failed[-1]["status"] = "failed"
        unknown = copy.deepcopy(original)
        unknown[-1]["ranks"] = 4
        for records in (failed, original[:-1], original[:-1] + [original[0]], unknown):
            with self.subTest(records=records[-1]), self.assertRaises(ValueError):
                aggregate(records, [1, 2], 3)

    def test_phase_schema_change_rejected(self):
        records = self.records()
        records[-1]["profile"]["phases"] = {}
        with self.assertRaisesRegex(ValueError, "phase schemas"):
            aggregate(records, [1, 2], 3)

    def test_nonfinite_negative_and_insufficient_samples_rejected(self):
        for values in ([1, 2], [1, 2, float("nan")], [1, 2, float("inf")], [1, 2, -1]):
            with self.subTest(values=values), self.assertRaises(ValueError):
                samples(values)

    def test_failed_launcher_is_preserved_and_output_is_exclusive(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            result = launch([sys.executable, "-c", "import sys; print('failure'); sys.exit(7)"],
                            directory, os.environ.copy(), 10)
            self.assertEqual(result["returncode"], 7)
            self.assertFalse(result["timed_out"])
            self.assertEqual(json.loads((directory / "launch.json").read_text()), result)
            self.assertEqual((directory / "launcher.stdout").read_text(), "failure\n")
            with self.assertRaises(FileExistsError):
                launch([sys.executable, "-c", "pass"], directory, os.environ.copy(), 10)

    def test_timeout_is_a_failed_observation(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = launch([sys.executable, "-c", "import time; time.sleep(30)"],
                            Path(temporary), os.environ.copy(), 0.1)
            self.assertEqual(result["returncode"], 124)
            self.assertTrue(result["timed_out"])

    def test_literal_arguments_do_not_invoke_a_shell(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            literal = "$(exit 8); `exit 9` with spaces"
            result = launch([sys.executable, "-c", "import sys; print(sys.argv[1])", literal],
                            directory, os.environ.copy(), 10)
            self.assertEqual(result["returncode"], 0)
            self.assertEqual((directory / "launcher.stdout").read_text().strip(), literal)


if __name__ == "__main__":
    unittest.main()
