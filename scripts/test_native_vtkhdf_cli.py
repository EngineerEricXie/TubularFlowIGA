#!/usr/bin/env python3
"""Bounded MPI integration fixtures for native temporal VTKHDF and resume.

Keep outputs in a fresh --work directory, then validate comparisons.json with
pvpython --no-mpi scripts/validate_native_vtkhdf_reader.py.
"""
import argparse
import copy
import json
import os
import shutil
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=False)
    binaries = args.bin_dir.resolve()
    source = Path(__file__).resolve().parents[1]
    mesh = source/"solvers/cpu/tests/data/native_tet_hydraulic_star.msh"
    comparisons, executions = [], []
    environment = {**os.environ, "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                   "OPENBLAS_NUM_THREADS": "1"}

    def run(name, binary, ranks, case, *options, success=True):
        command = ["mpiexec", "-np", str(ranks), str(binaries/binary), str(case), *map(str, options)]
        result = subprocess.run(command, capture_output=True, text=True, env=environment, timeout=90, check=False)
        (work/(name+".log")).write_text(result.stdout+result.stderr)
        executions.append({"name": name, "command": command, "returncode": result.returncode,
                           "expected_success": success})
        if (result.returncode == 0) != success:
            raise RuntimeError(f"unexpected result: {name}: {result.stdout}\n{result.stderr}")

    def save(name, case):
        path = work/(name+".json")
        path.write_text(json.dumps(case, indent=2)+"\n")
        return path

    def compare(hdf, legacy, prefix, steps, dt, first_step=1):
        comparisons.append({"hdf": str(hdf), "snapshots": [
            {"time_s": dt*step, "file": str(legacy/f"{prefix}{step}"/"snapshot.pvtu")}
            for step in range(first_step, steps+1)]})

    def check_pvd(directory, steps, dt):
        entries = ET.parse(directory/"flow.pvd").getroot().find("Collection").findall("DataSet")
        assert len(entries) == len(steps)
        for entry, step in zip(entries, steps):
            assert abs(float(entry.attrib["timestep"])-dt*step) < 1e-12
            assert entry.attrib["file"] == f"step_{step}/snapshot.pvtu"
            assert (directory/entry.attrib["file"]).is_file()

    base = json.loads((source/"solvers/cpu/tests/data/native_tet_hydraulic_star.json").read_text())
    base["mesh_file"] = str(mesh)
    for moving in (False, True):
        case = copy.deepcopy(base)
        if not moving:
            case["motion"] = {"kind": "fixed", "speed_x_m_s": 0.}
        case["species"] = [{"id": "tracer", "initial_concentration_mol_m3": 2.,
                            "diffusivity_m2_s": .01, "source_mol_m3_s": 0.,
                            "inflow_concentration_by_label_mol_m3": {"1": 2., "2": 2.}}]
        path = save(f"hydraulic_{moving}", case)
        for ranks in (1, 2, 8):
            label = f"hydraulic_{moving}_{ranks}"
            hdf, legacy = work/(label+"_hdf"), work/(label+"_pvtu")
            run(label+"_hdf", "native_tet_hydraulic_graph", ranks, path, "--output-dir", hdf)
            if ranks in (1, 2):
                run(label+"_pvtu", "native_tet_hydraulic_graph", ranks, path,
                    "--output-dir", legacy, "--visualization-format", "pvtu")
            else:
                legacy = work/f"hydraulic_{moving}_2_pvtu"
            compare(hdf/"flow.vtkhdf", legacy, "step_", 2, .05)
            compare(hdf/"species_tracer.vtkhdf", legacy, "species_tracer_step_", 2, .05)
            assert len(list(hdf.rglob("*"))) == 2
        # A root-only file creation failure must terminate all ranks.
        run(f"overwrite_{moving}", "native_tet_hydraulic_graph", 2, path,
            "--output-dir", work/f"hydraulic_{moving}_2_hdf", success=False)

    for moving in (False, True):
        velocity = [.01, 0., 0.] if moving else [0., 0., 0.]
        case = {"schema_version": 1, "mesh_file": str(mesh), "species_id": "tracer",
                "initial_concentration_mol_m3": 2., "diffusivity_m2_s": .01,
                "source_mol_m3_s": .1, "fluid_velocity_m_s": velocity,
                "mesh_velocity_m_s": velocity, "inflow_concentration_by_label_mol_m3": {},
                "monotone": False, "time": {"dt_s": .1, "steps": 3}}
        path = save(f"species_{moving}", case)
        for ranks in (1, 2):
            label = f"species_{moving}_{ranks}"
            full, resumed, legacy = (work/(label+suffix) for suffix in ("_full", "_resumed", "_pvtu"))
            run(label+"_full", "native_tet_species_transport", ranks, path, "--output-dir", full)
            run(label+"_stop", "native_tet_species_transport", ranks, path,
                "--output-dir", resumed, "--stop-after-step", 1)
            run(label+"_format_mismatch", "native_tet_species_transport", ranks, path,
                "--resume", resumed, "--visualization-format", "pvtu", success=False)
            run(label+"_resume", "native_tet_species_transport", ranks, path, "--resume", resumed)
            run(label+"_pvtu", "native_tet_species_transport", ranks, path,
                "--output-dir", legacy, "--visualization-format", "pvtu")
            compare(full/"species.vtkhdf", legacy, "step_", 3, .1)
            compare(resumed/"species.vtkhdf", legacy, "step_", 3, .1)
            assert (full/"checkpoint_step_3.bin").read_bytes() == (resumed/"checkpoint_step_3.bin").read_bytes()
            assert len(list(resumed.glob("*.vtkhdf"))) == 1
            assert not list(resumed.glob("step_*"))
    # Hydraulic checkpoint and temporal file must advance together.
    path = save("flow_restart", base)
    checkpoint, resumed = work/"flow_checkpoints", work/"flow_resumed"
    run("flow_stop", "native_tet_hydraulic_graph", 2, path,
        "--checkpoint-dir", checkpoint, "--output-dir", resumed, "--stop-after-step", 1)
    run("flow_resume", "native_tet_hydraulic_graph", 2, path,
        "--checkpoint-dir", checkpoint, "--restart-dir", checkpoint, "--output-dir", resumed)
    compare(resumed/"flow.vtkhdf", work/"hydraulic_True_2_pvtu", "step_", 2, .05)

    # The repository workflow also defaults to HDF and resumes the same file.
    for backend, binary, case in (
            ("native_tet_p2p1_ale_hydraulic", "native_tet_hydraulic_graph", path),
            ("native_tet_p1_prescribed_species", "native_tet_species_transport", work/"species_False.json")):
        configuration = save("workflow_"+backend, {
            "schema_version": 1, "backend": backend, "input_route": "volume",
            "input_file": str(mesh), "case_file": str(case), "mpi_ranks": 2,
            "output_directory": str(work/("workflow_"+backend))})
        for phase, options in (("partial", ["--stop-after-step", "1"]), ("resume", ["--resume"])):
            command = [os.sys.executable, str(source/"scripts/run_native_tet_workflow.py"),
                       str(configuration), "--solver", str(binaries/binary), *options]
            result = subprocess.run(command, capture_output=True, text=True, env=environment, timeout=120, check=False)
            (work/f"workflow_{backend}_{phase}.log").write_text(result.stdout+result.stderr)
            executions.append({"name": f"workflow_{backend}_{phase}", "command": command,
                               "returncode": result.returncode, "expected_success": True})
            if result.returncode:
                raise RuntimeError(result.stdout+result.stderr)
        folder = work/("workflow_"+backend)/"fields"
        if "hydraulic" in backend:
            compare(folder/"flow.vtkhdf", work/"hydraulic_True_2_pvtu", "step_", 2, .05)
        else:
            compare(folder/"species.vtkhdf", work/"species_False_2_pvtu", "step_", 3, .1)
    # Repeat restarts both in a fresh series and in the original directory.
    # An older checkpoint must still be rejected after a fresh series advances.
    for moving in (False, True):
        case = copy.deepcopy(base)
        case["time"]["steps"] = 3
        if not moving:
            case["motion"] = {"kind": "fixed", "speed_x_m_s": 0.}
        path = save(f"repeated_restart_{moving}", case)
        for ranks in (1, 2):
            outputs = {}
            for fmt in ("vtkhdf", "pvtu"):
                label = f"repeated_{moving}_{ranks}_{fmt}"
                original, fresh = work/(label+"_original"), work/(label+"_fresh")
                checkpoint, saved = work/(label+"_checkpoints"), work/(label+"_step1")
                common = ("--visualization-format", fmt)
                run(label+"_stop", "native_tet_hydraulic_graph", ranks, path, *common,
                    "--checkpoint-dir", checkpoint, "--output-dir", original, "--stop-after-step", 1)
                shutil.copytree(checkpoint, saved)
                run(label+"_fresh", "native_tet_hydraulic_graph", ranks, path, *common,
                    "--checkpoint-dir", checkpoint, "--restart-dir", checkpoint,
                    "--output-dir", fresh, "--stop-after-step", 2)
                run(label+"_stale", "native_tet_hydraulic_graph", ranks, path, *common,
                    "--restart-dir", saved, "--output-dir", fresh, success=False)
                run(label+"_again", "native_tet_hydraulic_graph", ranks, path, *common,
                    "--checkpoint-dir", checkpoint, "--restart-dir", checkpoint, "--output-dir", fresh)
                run(label+"_original", "native_tet_hydraulic_graph", ranks, path, *common,
                    "--checkpoint-dir", saved, "--restart-dir", saved, "--output-dir", original)
                outputs[fmt] = (original, fresh)
                if fmt == "pvtu":
                    check_pvd(original, [1, 2, 3], .05)
                    check_pvd(fresh, [2, 3], .05)
            compare(outputs["vtkhdf"][0]/"flow.vtkhdf", outputs["pvtu"][0], "step_", 3, .05)
            compare(outputs["vtkhdf"][1]/"flow.vtkhdf", outputs["pvtu"][0], "step_", 3, .05, first_step=2)
    save("comparisons", {"comparisons": comparisons, "executions": executions})
    print(json.dumps({"passed": True, "executions": len(executions), "series": len(comparisons)}))


if __name__ == "__main__":
    main()
