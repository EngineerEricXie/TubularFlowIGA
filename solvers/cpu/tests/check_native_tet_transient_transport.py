"""Regression for native_tet_transient_transport on a small Y bifurcation.

Checks flow conservation, species budgets, the inlet concentration schedule,
the BDF2 start-up (backward Euler on the first step), output sizes, and that a
run stopped after step 2 and restarted reproduces the uninterrupted checkpoint,
step log and snapshot series. A nonzero outlet pressure must appear once, not
twice, in the VTK pressure, and a changed mesh must be rejected on restart.

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


def run(output, *extra, case_path=None, check=True):
    command = ["mpiexec", "-np", ranks, str(binary), str(case_path or case), "--output-dir", str(output), *extra]
    return subprocess.run(command, check=check, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        timeout=900).returncode


def vtk_pressure(snapshot):
    import xml.etree.ElementTree as ET
    values = []
    for piece in snapshot.glob("*.vtu"):
        node = ET.parse(piece).getroot().find(".//Piece/PointData/DataArray[@Name='pressure_pa']")
        values += [float(value) for value in node.text.split()]
    return values


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
    # Snapshots on both sides of the restart, so flow.pvd must keep both.
    mesh = Path(root) / "mesh.msh"
    mesh.write_text((case.parent / spec["mesh_file"]).read_text())
    spec["mesh_file"] = str(mesh)
    spec["output"]["visualization_steps"] = [2, 4]
    case = Path(root) / "case.json"
    case.write_text(json.dumps(spec))
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
    # Same case text, different mesh contents: the restart identity rejects it.
    original = mesh.read_text()
    lines = original.splitlines()
    coordinate = lines.index("$Nodes") + 3 + int(lines[lines.index("$Nodes") + 2].split()[3])
    moved = [float(value) for value in lines[coordinate].split()]
    moved[0] += 1e-9
    lines[coordinate] = " ".join(repr(value) for value in moved)
    mesh.write_text("\n".join(lines) + "\n")
    assert mesh.read_text() != original
    assert run(split, "--restart", check=False) != 0, "a changed mesh was accepted on restart"
    mesh.write_text(original)
    run(split, "--restart")
    a, b = (full / "checkpoint.f64").read_bytes(), (split / "checkpoint.f64").read_bytes()
    assert len(a) == len(b)
    values = struct.unpack(f"{len(a) // 8}d", a), struct.unpack(f"{len(b) // 8}d", b)
    worst = max(abs(x - y) / max(1.0, abs(x)) for x, y in zip(*values))
    assert worst <= 1e-10, worst
    resumed = [json.loads(line)["step"] for line in (split / "steps.jsonl").read_text().splitlines()]
    assert resumed == list(range(1, spec["time"]["steps"] + 1)), resumed
    series_text = (split / "flow.pvd").read_text()
    assert "vtk_step_000002" in series_text and "vtk_step_000004" in series_text, series_text
    # The outlet gauge pressure is applied once: the VTK pressure range equals
    # the P1 range of the flow file (midside values interpolate vertices).
    pressure = dict(spec, outlets={"pressure_pa": 100.0}, time={"dt_s": spec["time"]["dt_s"], "steps": 1},
        output={"every_steps": 1, "precision": "float64", "visualization_steps": [1]})
    pressure_case = Path(root) / "pressure.json"
    pressure_case.write_text(json.dumps(pressure))
    gauge = Path(root) / "gauge"
    run(gauge, case_path=pressure_case)
    raw = (gauge / "flow_000001.f64").read_bytes()
    flow = struct.unpack(f"{len(raw) // 8}d", raw)
    vertex_pressure = flow[3 * series["p2_nodes"]:]
    shown = vtk_pressure(gauge / "vtk_step_000001")
    outlets = json.loads((gauge / "steps.jsonl").read_text().splitlines()[0])["boundaries"]
    assert all(abs(outlets[str(label)]["average_pressure_pa"] - 100.0) < 5.0
        for label in spec["boundaries"]["outlet_labels"]), outlets
    assert min(vertex_pressure) - 1e-9 <= min(shown) and max(shown) <= max(vertex_pressure) + 1e-9, (
        min(shown), max(shown), min(vertex_pressure), max(vertex_pressure))
print(f"native transient transport ({ranks} rank(s)): PASS")
