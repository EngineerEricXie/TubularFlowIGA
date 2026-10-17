"""Regression for native_tet_transient_transport on a small Y bifurcation.

Checks flow conservation, species budgets, the inlet concentration schedule,
the BDF2 start-up (backward Euler on the first step), output sizes, and that a
run stopped after step 2 and restarted reproduces the uninterrupted checkpoint.

    python3 tests/check_native_tet_transient_transport.py BINARY CASE [RANKS]
"""
import json
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

binary, case = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
ranks = sys.argv[3] if len(sys.argv) > 3 else "1"


def run(output, *extra):
    command = ["mpiexec", "-np", ranks, str(binary), str(case), "--output-dir", str(output), *extra]
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL, timeout=900)


def table(points, t):
    if t <= points[0][0]:
        return points[0][1]
    for (t0, v0), (t1, v1) in zip(points, points[1:]):
        if t <= t1:
            return v0 + (v1 - v0) * (t - t0) / (t1 - t0)
    return points[-1][1]


spec = json.loads(case.read_text())
with tempfile.TemporaryDirectory() as root:
    full, split = Path(root) / "full", Path(root) / "split"
    run(full)
    steps = [json.loads(line) for line in (full / "steps.jsonl").read_text().splitlines()]
    assert len(steps) == spec["time"]["steps"]
    tracer = spec["species"][0]
    for row in steps:
        assert row["relative_mass_imbalance"] <= 1e-5 and row["relative_inlet_error"] <= 1e-8, row
        for name, item in row["species"].items():
            scale = max(abs(item["outward_transfer_mol_s"]), abs(item["source_mol_s"]), 1e-30)
            assert abs(item["balance_defect_mol_s"]) <= 1e-6 * scale or abs(item["balance_defect_mol_s"]) <= 1e-15, (name, item)
        expected = table(tracer["inflow_concentration_mol_m3"], row["time_s"])
        assert abs(row["species"]["tracer"]["inlet_value"] - expected) <= 1e-12
        assert row["species"]["tracer"]["scheme"] == ("backward_euler" if row["step"] == 1 else "bdf2")
        assert row["species"]["age"]["scheme"] == "backward_euler"
    assert steps[-1]["species"]["tracer"]["maximum"] >= 0.99 * table(tracer["inflow_concentration_mol_m3"], steps[-1]["time_s"])
    assert steps[-1]["species"]["age"]["inventory_mol"] > 0
    series = json.loads((full / "series.json").read_text())
    dofs = 3 * series["p2_nodes"] + series["vertices"]
    for step in (2, 4):
        assert (full / f"flow_{step:06d}.f32").stat().st_size == 4 * dofs
        assert (full / f"species_tracer_{step:06d}.f32").stat().st_size == 4 * series["vertices"]
    assert (full / "flow.pvd").exists() and (full / "completed.json").exists()
    run(split, "--stop-after-step", "2")
    assert not (split / "completed.json").exists()
    run(split, "--restart")
    a, b = (full / "checkpoint.f64").read_bytes(), (split / "checkpoint.f64").read_bytes()
    assert len(a) == len(b)
    values = struct.unpack(f"{len(a) // 8}d", a), struct.unpack(f"{len(b) // 8}d", b)
    worst = max(abs(x - y) / max(1.0, abs(x)) for x, y in zip(*values))
    assert worst <= 1e-10, worst
print(f"native transient transport ({ranks} rank(s)): PASS")
