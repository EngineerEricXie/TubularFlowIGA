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
prepare = load("hpc_prepare_scaling_cases")
tiers = load("hpc_test_tiers")
finalize = load("hpc_finalize_cross_node")


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value) + "\n")


def build_manifest(revision):
    return {"kind": "hpc_build_manifest", "dependency_status": "passed",
            "source": {"commit": {"returncode": 0, "output": revision},
                       "status_porcelain": {"returncode": 0, "output": ""}},
            "binaries": [{"exists": True}]}


def rank_run(hostname, rank=0, ranks=2):
    return {"kind": "hpc_rank_run", "status": "process_passed", "returncode": 0,
            "timed_out": False, "hostname": hostname, "rank": rank, "ranks": ranks}


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

    def test_weak_dimensions_hold_work_per_rank(self):
        for ranks in (1, 64, 128, 256):
            transverse, axial = prepare.weak_dimensions(ranks, 256, 8)
            self.assertLessEqual(max(transverse, axial), 128)
            self.assertEqual(transverse * transverse * axial // ranks, 256)

    def test_slurm_wrappers_parse(self):
        for path in (ROOT / "solvers/cpu/slurm/multinode_graph.sbatch",
                     ROOT / "solvers/cpu/slurm/cross_node_scaling.sbatch",
                     ROOT / "solvers/cpu/slurm/cross_node_fsi.sbatch"):
            result = subprocess.run(["bash", "-n", str(path)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)

    def test_spooled_wrappers_resolve_submission_checkout(self):
        for name, boundary in (("multinode_graph", "source_case="),
                               ("cross_node_fsi", "output="),
                               ("cross_node_scaling", "source_root=")):
            script = ROOT / "solvers/cpu/slurm" / (name + ".sbatch")
            with tempfile.TemporaryDirectory() as directory:
                spool = Path(directory) / "slurm_script"
                spool.write_text("module() { :; }\n" +
                                 script.read_text().split(boundary, 1)[0] +
                                 'printf "%s" "$repo_root"\n')
                environment = dict(os.environ, SLURM_SUBMIT_DIR=str(ROOT))
                environment.pop("IGA_REPO_ROOT", None)
                result = subprocess.run(["bash", str(spool)], cwd=directory,
                                        env=environment, capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stdout, str(ROOT))

    def test_spooled_wrapper_preflight_failure_preserves_scheduler_record(self):
        for name in ("multinode_graph", "cross_node_fsi", "cross_node_scaling"):
            with tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                spool = root / "slurm_script"
                spool.write_text("module() { :; }\n" +
                    (ROOT / "solvers/cpu/slurm" / (name + ".sbatch")).read_text())
                environment = {key: value for key, value in os.environ.items()
                               if not key.startswith(("SLURM_", "IGA_"))}
                environment.update(SLURM_SUBMIT_DIR=str(ROOT), PETSC_DIR=str(root),
                    IGA_GRAPH_CASE=str(root / "case"), IGA_OUTPUT_ROOT=str(root / "graph"),
                    IGA_CHECKPOINT_ROOT=str(root / "checkpoint"),
                    IGA_FSI_OUTPUT=str(root / "fsi"),
                    IGA_SCALING_CASE_ROOT=str(root / "cases"),
                    IGA_SCALING_OUTPUT=str(root / "scaling"))
                result = subprocess.run(["bash", str(spool)], cwd=directory,
                                        env=environment, capture_output=True, text=True)
                self.assertEqual(result.returncode, 2, result.stderr)
                record = {"multinode_graph": root / "graph/attempt-0/scheduler.json",
                          "cross_node_fsi": root / "fsi/scheduler.json",
                          "cross_node_scaling": root / "scaling.scheduler.json"}[name]
                report = json.loads(record.read_text())
                self.assertEqual(report["status"], "failed")
                self.assertEqual(report["returncode"], 2)

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

    def test_cross_node_finalizer_accepts_complete_same_revision_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            revision = "a" * 40
            graph, fsi, scaling_root = root / "graph", root / "fsi", root / "scaling"
            for attempt, status in ((0, "checkpointed"), (1, "passed")):
                location = graph / f"attempt-{attempt}"
                write_json(location / "scheduler.json", {
                    "kind": "hpc_scheduler_attempt", "status": status, "returncode": 0,
                    "checkpoint_requested": attempt == 0,
                    "command": ["solver"] + (["--restart-dir", "checkpoint"] if attempt else []),
                    "slurm": {"SLURM_JOB_NUM_NODES": "2", "SLURM_JOB_ID": "91",
                              "SLURM_RESTART_COUNT": str(attempt)}})
                write_json(location / "build.json", build_manifest(revision))
                (location / "nodes.txt").write_text("node-a\nnode-b\n")
            write_json(graph / "attempt-1/results/graph_binding_manifest.json", {})

            write_json(fsi / "scheduler.json", {"kind": "hpc_scheduler_attempt",
                "status": "passed", "returncode": 0,
                "slurm": {"SLURM_JOB_NUM_NODES": "2", "SLURM_JOB_ID": "92"}})
            write_json(fsi / "build.json", build_manifest(revision))
            write_json(fsi / "scheduled/result.json", {"kind": "hpc_test_tier",
                "tier": "scheduled", "status": "passed", "commands": [{}],
                "source_commit": revision})
            write_json(fsi / "strong-fsi.json", {"status": "passed"})
            write_json(fsi / "paired-restart.json", {"status": "passed",
                "source_ranks": 4, "target_ranks": 2})
            for run in ("writer-4", "reader-2"):
                write_json(fsi / run / "rank-0/run.json", rank_run("node-a", 0))
                write_json(fsi / run / "rank-1/run.json", rank_run("node-b", 1))

            ranks, repetitions = [1, 256], 3
            inputs, runs, summary = {}, [], {"strong": {}, "weak": {}}
            for mode in ("strong", "weak"):
                for count in ranks:
                    elements = 16384 if mode == "strong" else 256 * count
                    inputs[f"{mode}-{count}"] = {
                        "elements": elements, "elements_per_rank": elements / count}
                    summary[mode][str(count)] = {"elements_per_rank": elements / count,
                        "max_rank_process_wall_s": {}, "max_rank_peak_rss_bytes": {},
                        "solver_iterations": [2] * repetitions, "speedup_vs_one_rank": 1,
                        "parallel_efficiency": 1,
                        "phases_max_rank_exclusive_s": {"communication": {}}}
                    for repetition in range(repetitions):
                        run_dir = scaling_root / f"{mode}-np{count}-repeat{repetition}"
                        if count == 256:
                            write_json(run_dir / "rank-0/run.json", rank_run("node-a", 0, count))
                            write_json(run_dir / "rank-1/run.json", rank_run("node-b", 1, count))
                        runs.append({"mode": mode, "ranks": count, "repetition": repetition,
                            "directory": str(run_dir), "status": "passed",
                            "physical_validation": {"returncode": 0},
                            "field_comparison": {"velocity": {"passed": True},
                                "pressure": {"passed": True}} if mode == "strong" else {}})
            write_json(scaling_root / "summary.json", {"kind": "hpc_cross_node_scaling",
                "status": "passed", "source_status": [], "source_commit": revision,
                "allocation": {"SLURM_JOB_NUM_NODES": "2", "SLURM_JOB_ID": "93"},
                "ranks": ranks, "repetitions": repetitions, "inputs": inputs,
                "runs": runs, "summary": summary})
            write_json(root / "scaling.build.json", build_manifest(revision))
            acceptance = root / "acceptance.json"
            result = subprocess.run([sys.executable, str(SCRIPTS / "hpc_finalize_cross_node.py"),
                "--graph-output", str(graph), "--fsi-output", str(fsi),
                "--scaling-output", str(scaling_root), "--output", str(acceptance)],
                cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            report = json.loads(acceptance.read_text())
            self.assertEqual(report["status"], "passed")
            self.assertEqual(len(report["completed_items"]), 7)

    def test_cross_node_finalizer_rejects_dirty_build(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "build.json"
            report = build_manifest("b" * 40)
            report["source"]["status_porcelain"]["output"] = " M solver.cpp"
            write_json(path, report)
            with self.assertRaisesRegex(ValueError, "not clean"):
                finalize.build_revision(path)


if __name__ == "__main__":
    unittest.main()
