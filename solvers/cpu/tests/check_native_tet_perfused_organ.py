"""Regression for the native_tet_perfused_organ CLI on a small synthetic organ.

Builds a Kuhn-split box of tissue, an inflow and an outflow SWC tree and
nearest-terminal territories with the standard library only, then checks:
input checks and the geometry file, flow conservation and solve diagnostics,
the opt-in acid inhibition of glycolysis and its validation,
that an advance split in two through state_out/state_in reproduces one
uninterrupted advance bit for bit, per-substep advance diagnostics, rejection
of invalid restored states and plant parameters, and the --serve protocol
(errors answered per line, only {"command": "exit"} stops the server).

    python3 tests/check_native_tet_perfused_organ.py BINARY
"""
import json
import math
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

binary = Path(sys.argv[1]).resolve()
SPECIES = ["o2_total_mol_m3", "glucose_mol_m3", "lactate_mol_m3", "buffer_base_mol_m3", "co2_total_mol_m3",
    "insulin_mU_L"]
KUHN = ((0, 1, 3, 7), (0, 3, 2, 7), (0, 2, 6, 7), (0, 6, 4, 7), (0, 4, 5, 7), (0, 5, 1, 7))


def write_box_mesh(path, counts, edge_mm):
    """Gmsh 4.1 box of Kuhn tetrahedra, boundary label 1, volume named fluid."""
    nx, ny, nz = counts
    index = lambda i, j, k: (k*(ny+1)+j)*(nx+1)+i
    points = [(i*edge_mm, j*edge_mm, k*edge_mm) for k in range(nz+1) for j in range(ny+1) for i in range(nx+1)]
    tets = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                corners = [index(i+dx, j+dy, k+dz) for dz in (0, 1) for dy in (0, 1) for dx in (0, 1)]
                tets += [tuple(corners[c] for c in pattern) for pattern in KUHN]
    uses = {}
    for tet in tets:
        for opposite in range(4):
            face = tuple(sorted(tet[n] for n in range(4) if n != opposite))
            uses[face] = uses.get(face, 0)+1
    triangles = [face for face, count in uses.items() if count == 1]
    box = f"0 0 0 {nx*edge_mm} {ny*edge_mm} {nz*edge_mm}"
    lines = ["$MeshFormat", "4.1 0 8", "$EndMeshFormat", "$PhysicalNames", "2", '2 1 "boundary_label_1"',
        '3 2 "fluid"', "$EndPhysicalNames", "$Entities", "0 0 1 1", f"1 {box} 1 1 0", f"1 {box} 1 2 1 1",
        "$EndEntities", "$Nodes", f"1 {len(points)} 1 {len(points)}", f"3 1 0 {len(points)}"]
    lines += [str(n+1) for n in range(len(points))]
    lines += [f"{x} {y} {z}" for x, y, z in points]
    total = len(triangles)+len(tets)
    lines += ["$EndNodes", "$Elements", f"2 {total} 1 {total}", f"2 1 2 {len(triangles)}"]
    lines += [f"{n} {a+1} {b+1} {c+1}" for n, (a, b, c) in enumerate(triangles, start=1)]
    lines.append(f"3 1 4 {len(tets)}")
    first = len(triangles)+1
    lines += [f"{n} {a+1} {b+1} {c+1} {d+1}" for n, (a, b, c, d) in enumerate(tets, start=first)]
    lines.append("$EndElements")
    path.write_text("\n".join(lines)+"\n")
    ids = list(range(first, first+len(tets)))
    cells = []
    for tet in tets:
        p = [points[n] for n in tet]
        u, v, w = ([p[m][a]-p[0][a] for a in range(3)] for m in (1, 2, 3))
        volume = abs(u[0]*(v[1]*w[2]-v[2]*w[1])-u[1]*(v[0]*w[2]-v[2]*w[0])+u[2]*(v[0]*w[1]-v[1]*w[0]))/6.0
        cells.append(([sum(q[a] for q in p)/4.0 for a in range(3)], volume))
    return ids, cells


def write_tree(path, rows):
    path.write_text("".join(" ".join(str(value) for value in row)+"\n" for row in rows))
    parents = {row[6] for row in rows}
    return {row[0]: row[2:5] for row in rows if row[0] not in parents and row[6] != -1}


def territories(ids, cells, terminals):
    groups = {}
    for cell_id, (centroid, volume) in zip(ids, cells):
        nearest = min(terminals, key=lambda t: sum((centroid[a]-terminals[t][a])**2 for a in range(3)))
        groups.setdefault(nearest, []).append((cell_id, volume))
    result = {}
    for terminal, members in groups.items():
        total = sum(volume for _, volume in members)
        result[str(terminal)] = {"cell_ids": [cell for cell, _ in members],
            "weights": [volume/total for _, volume in members]}
    return result


def run(root, *arguments, expect_success=True):
    result = subprocess.run([str(binary), str(root/"case.json"), *map(str, arguments)], cwd=root,
        capture_output=True, text=True, timeout=600)
    if expect_success and result.returncode:
        raise SystemExit(f"{arguments} failed: {result.stderr}")
    if not expect_success and result.returncode == 0:
        raise SystemExit(f"{arguments} unexpectedly succeeded")
    return result


def request(root, name, body, expect_success=True):
    (root/f"{name}.json").write_text(json.dumps(body))
    result = run(root, "--request", f"{name}.json", "--output", f"{name}_out.json", expect_success=expect_success)
    return json.loads((root/f"{name}_out.json").read_text()) if expect_success else result.stderr


def main():
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        ids, cells = write_box_mesh(root/"tissue.msh", (4, 3, 3), 5.0)
        portal = write_tree(root/"portal.swc", [(1, 1, -4, 7.5, 7.5, 1.5, -1), (2, 3, 4, 7.5, 7.5, 1.2, 1),
            (3, 3, 10, 4, 4, 0.8, 2), (4, 3, 10, 11, 11, 0.8, 2), (5, 3, 16, 4, 11, 0.8, 2), (6, 3, 16, 11, 4, 0.8, 2)])
        hepatic = write_tree(root/"hepatic.swc", [(1, 1, 24, 7.5, 7.5, 1.5, -1), (2, 3, 16, 7.5, 7.5, 1.2, 1),
            (3, 3, 5, 4, 7.5, 0.8, 2), (4, 3, 5, 11, 7.5, 0.8, 2), (5, 3, 12, 7.5, 3, 0.8, 2)])
        (root/"territories.json").write_text(json.dumps({"schema_version": 1, "trees": {
            "portal": territories(ids, cells, portal), "hepatic": territories(ids, cells, hepatic)}}))
        (root/"case.json").write_text(json.dumps({"schema_version": 1, "kind": "perfused_organ",
            "tissue_mesh_file": "tissue.msh", "mesh_length_scale_to_m": 1e-3, "trees": [
                {"id": "portal", "role": "inflow", "skeleton_file": "portal.swc", "length_scale_to_m": 1e-3},
                {"id": "hepatic", "role": "outflow", "skeleton_file": "hepatic.swc", "length_scale_to_m": 1e-3}],
            "territories_file": "territories.json", "reaction_model": "hepatic_v1", "lobular_zones": 3}))
        checked = json.loads(run(root, "--check-input").stdout)
        assert checked["ok"] and checked["cells"] == len(cells) and checked["trees"] == 2, checked
        described = json.loads(run(root, "--describe", "geometry.f64").stdout)
        assert (root/"geometry.f64").stat().st_size == 5*8*len(cells)
        volume = sum(v for _, v in cells)*1e-9
        assert math.isclose(described["volume_m3"], volume, rel_tol=1e-12), described["volume_m3"]

        ph = 6.1+math.log10(24.0/(0.0307*40.0))
        buffer_base = 24.0+17.7+0.42*4.6+(1.43*4.6+7.7)*(ph-7.4)
        inlet = dict(zip(SPECIES, [3.0, 8.0, 1.0, buffer_base, 24.0+0.0307*40.0, 100.0]))
        fractions = dict(zip(SPECIES, [0.2, 0.7, 0.7, 0.3, 0.5, 0.25]))
        kinetics = {"vmax_o2_mol_m3_s": 0.02, "glycolysis_capacity_mol_m3_s": 0.02, "lactate_vmax_mol_m3_s": 0.004,
            "glycogenolysis_vmax_mol_m3_s": 0.003, "glucose_uptake_vmax_mol_m3_s": 0.003,
            "insulin_clearance_s_inv": 1e-3, "temperature_c": 37.0}
        conductance = {"portal": 3e-5, "hepatic": 9e-5}
        solve = {"command": "solve", "dt_s": 1e9, "viscosity_pa_s": 3e-3,
            "root_pressure_pa": {"portal": 1000.0, "hepatic": 0.0}, "specific_conductance_inv_pa_s": conductance,
            "mobility_m2_pa_s": 1e-8, "inlet": inlet, "distribution_volume_fraction": fractions,
            "chemistry": {"hemoglobin_g_dl": 4.6}, "kinetics": kinetics, "glycogen_mol_m3": 100.0,
            "field_output": {"file": "fields.f64", "names": ["po2_mmHg", "tissue_pressure_pa"]}}
        responses = {}
        for resolution in ("cell", "portal_territory", "organ"):
            solved = responses[resolution] = request(root, "solve_"+resolution, dict(solve, resolution=resolution))
            trees = solved["trees"]
            inflow = trees["portal"]["root_inflow_m3_s"]
            assert inflow > 0 and math.isclose(-trees["hepatic"]["root_inflow_m3_s"], inflow, rel_tol=1e-9), trees
            diagnostics = solved["diagnostics"]
            assert diagnostics["newton_failures"] == 0 and diagnostics["gauss_seidel_unconverged"] == 0, diagnostics
            for name, totals in solved["species_totals"].items():
                assert abs(totals["balance_defect"]) <= 1e-9*abs(totals["inflow"]), (resolution, name, totals)
        assert (root/"fields.f64").stat().st_size == 8*(3*len(cells)+len(cells))
        inhibited = request(root, "solve_inhibited", dict(solve, kinetics=dict(kinetics, glycolysis_ph50=6.8,
            glycolysis_ph_hill=4.0)))
        assert inhibited["diagnostics"]["newton_failures"] == 0, inhibited["diagnostics"]
        lactate = lambda response: response["species_totals"]["lactate_mol_m3"]["reaction"]
        assert lactate(inhibited) <= lactate(responses["cell"]), (lactate(inhibited), lactate(responses["cell"]))
        error = request(root, "solve_bad_ph50", dict(solve, kinetics=dict(kinetics, glycolysis_ph50=5.0)),
            expect_success=False)
        assert "pH50" in error, error

        parameters = {"viscosity_pa_s": 3e-3, "outflow_pressure_pa": {"hepatic": 0.0},
            "specific_conductance_inv_pa_s": conductance, "mobility_m2_pa_s": 1e-8,
            "distribution_volume_fraction": fractions, "kinetics": kinetics, "chemistry": {}}
        controls = {"inlet_pressure_pa": 1000.0, "fio2": 0.6, "bolus": {"glucose_mol_m3": 1e-5},
            "vasodilator_bolus_mg": 2.0}
        plain = {"inlet_pressure_pa": 1000.0, "fio2": 0.6}
        initial = {"reservoir_volume_m3": 2e-3, "hematocrit_percent": 18.0, "reservoir": inlet, "glycogen_mol_m3": 100.0}
        advance = lambda duration, substeps, ctl, state, out: dict({"command": "advance", "duration_s": duration,
            "substeps": substeps, "parameters": parameters, "controls": ctl,
            "state_out": {"json": out+".json", "fields": out+".f64"}}, **state)
        full = request(root, "full", advance(1200.0, 6, controls, {"initial_state": initial}, "full_state"))
        first = request(root, "first", advance(600.0, 3, controls, {"initial_state": initial}, "half_state"))
        second = request(root, "second", advance(600.0, 3, plain,
            {"state_in": {"json": "half_state.json", "fields": "half_state.f64"}}, "second_state"))
        assert full["records"] == first["records"]+second["records"], "split advance differs from one advance"
        assert (root/"full_state.f64").read_bytes() == (root/"second_state.f64").read_bytes()
        assert json.loads((root/"full_state.json").read_text()) == json.loads((root/"second_state.json").read_text())
        assert (root/"full_state.f64").stat().st_size == 8*((len(SPECIES)+2)*3*len(cells)+len(cells))
        assert full["newton_failures"] == 0 and full["gauss_seidel_unconverged"] == 0, full
        for record in full["records"]:
            diagnostics = record["diagnostics"]
            assert diagnostics["newton_failures"] == 0 and diagnostics["gauss_seidel_unconverged"] == 0, diagnostics
            assert diagnostics["flow_residual_norm"] < 1e-15, diagnostics
            for name in SPECIES:
                scale = record["flow_m3_s"]*inlet[name]
                assert abs(diagnostics["balance_defect"][name]) <= 1e-9*scale, (name, diagnostics)

        # A restored state outside what an advance can produce is rejected.
        header = json.loads((root/"half_state.json").read_text())
        fields = bytearray((root/"half_state.f64").read_bytes())
        nodes = 3*len(cells)
        struct.pack_into("<d", fields, 8*len(SPECIES)*nodes, 2.0)  # viability of node 0
        (root/"bad_state.f64").write_bytes(bytes(fields))
        (root/"bad_state.json").write_text(json.dumps(header))
        error = request(root, "reject_state", advance(60.0, 1, plain,
            {"state_in": {"json": "bad_state.json", "fields": "bad_state.f64"}}, "unused"), expect_success=False)
        assert "viability" in error, error
        assert not (root/"unused.json").exists()
        error = request(root, "bad_parameters", dict(advance(60.0, 1, plain, {"initial_state": initial}, "unused"),
            parameters=dict(parameters, oedema_maximum=1.0)), expect_success=False)
        assert "oedema maximum" in error, error
        error = request(root, "bad_cost", dict(advance(60.0, 1, plain, {"initial_state": initial}, "unused"),
            parameters=dict(parameters, shear_injury_rate_inv_pa_s=-1.0)), expect_success=False)
        assert "control-cost" in error, error

        # An independent sweep-gas CO2 fixes arterial pCO2 whatever the FiO2.
        pco2 = lambda fio2: request(root, "gas_co2_%g" % fio2, advance(60.0, 1, dict(plain, fio2=fio2,
            gas_co2_fraction=0.05), {"initial_state": initial}, "unused_gas"))["records"][0]["arterial"]["pco2_mmHg"]
        assert math.isclose(pco2(0.3), pco2(0.9), rel_tol=1e-9), (pco2(0.3), pco2(0.9))
        error = request(root, "bad_gas_co2", advance(60.0, 1, dict(plain, gas_co2_fraction=0.5),
            {"initial_state": initial}, "unused"), expect_success=False)
        assert "gas CO2" in error, error

        lines = [json.dumps(dict(solve, tag="a")), "{", json.dumps(dict(advance(60.0, 1, plain,
            {"initial_state": initial}, "served"), substeps=0)), json.dumps({"tag": "exit"}),
            json.dumps({"command": "exit"}), json.dumps(dict(solve, tag="after"))]
        served = subprocess.run([str(binary), str(root/"case.json"), "--serve"], cwd=root, input="\n".join(lines)+"\n",
            capture_output=True, text=True, timeout=600, check=True)
        answers = [json.loads(line) for line in served.stdout.splitlines()]
        assert answers[0].get("ready") and [a["ok"] for a in answers[1:]] == [True, False, False, False], answers
        assert answers[1]["tag"] == "a" and "substeps" in answers[3]["error"], answers
    print("native_tet_perfused_organ CLI checks passed")


if __name__ == "__main__":
    main()
