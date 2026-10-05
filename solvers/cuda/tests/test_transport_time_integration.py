#!/usr/bin/env python3
"""Exercise CUDA time-scheme dispatch with an existing packed control mesh.

Usage: python3 test_transport_time_integration.py IGA_CUDA DATABASE CASE_DIR
Requires a CUDA device and CASE_DIR/controlmesh.vtk from the packed database.
All generated configurations and results stay in a temporary directory.
"""

import json
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile


def main():
    binary, database, source = (Path(arg).resolve(strict=True) for arg in sys.argv[1:])
    mesh = (source / "controlmesh.vtk").read_text().split()
    nodes = int(mesh[mesh.index("POINTS") + 1])
    first_label = mesh.index("LOOKUP_TABLE") + 2
    labels = {int(value) for value in mesh[first_label:first_label + nodes]}
    configuration = {
        "schema_version": 2,
        "fields": [{"name": "c", "kind": "scalar", "initial_value": 1.0}],
        "time": {"dt": 0.1, "steps": 2},
        "equation_systems": [{
            "name": "decay", "kind": "linear_transport", "unknowns": ["c"],
            "terms": [{"operator": operator, "equation": "c", "coefficient": 1.0}
                      for operator in ("time_derivative", "linear_coupling")],
        }],
        "boundaries": [{"label": label, "name": f"boundary_{label}", "conditions": []}
                       for label in sorted(labels) if label >= 0],
    }
    with tempfile.TemporaryDirectory(prefix="cuda-transport-scheme-") as directory:
        root = Path(directory)
        for scheme in (None, "backward_euler", "bdf2", "bdf2_restart"):
            case = root / (scheme or "default")
            case.mkdir()
            shutil.copyfile(source / "controlmesh.vtk", case / "controlmesh.vtk")
            (case / "initial_velocityfield.txt").write_text("0 0 0\n" * nodes)
            system = configuration["equation_systems"][0]
            if scheme:
                system["time_integration"] = "bdf2" if scheme == "bdf2" else "backward_euler"
            (case / "simulation_config.json").write_text(json.dumps(configuration))
            output = case / "result.txt"
            command = [str(binary), "solve", str(database), str(case),
                       "--output", str(output), "--visualization-format", "vtu"]
            if scheme == "bdf2_restart":
                metadata = {"schema_version": 2, "nodes": nodes, "fields": ["c"],
                            "system": "decay", "velocity_source": "prescribed",
                            "completed_step": 1, "physical_time": 0.1, "dt": 0.1,
                            "state_file": "restart.state", "state_format": "raw_float64",
                            "time_integration": "bdf2", "history_file": "restart.history"}
                (case / "restart.json").write_text(json.dumps(metadata))
                (case / "restart.state").write_bytes(struct.pack(f"={nodes}d", *([1 / 1.1] * nodes)))
                (case / "restart.history").write_bytes(struct.pack(f"={nodes}d", *([1.0] * nodes)))
                command += ["--restart", str(case / "restart")]
            result = subprocess.run(command, capture_output=True, text=True, timeout=60)
            diagnostic = result.stdout + result.stderr
            if scheme in ("bdf2", "bdf2_restart"):
                expected = ("CUDA transport supports backward_euler only" if scheme == "bdf2"
                            else "CUDA transport restart supports backward_euler checkpoints only")
                if result.returncode == 0 or expected not in diagnostic or output.exists():
                    raise RuntimeError(f"{scheme} was not rejected before output:\n{diagnostic}")
            else:
                if result.returncode != 0:
                    raise RuntimeError(f"{scheme or 'default'} failed:\n{diagnostic}")
                rows = [line.split() for line in output.read_text().splitlines()]
                if len(rows) != nodes or any(len(row) != 2 or int(row[0]) != i
                        or not abs(float(row[1]) - 1 / 1.1**2) < 1e-5 for i, row in enumerate(rows)):
                    raise RuntimeError("backward-Euler decay differs from the analytic recurrence")
            print(f"CUDA transport {scheme or 'default'} passed")


if __name__ == "__main__":
    main()
