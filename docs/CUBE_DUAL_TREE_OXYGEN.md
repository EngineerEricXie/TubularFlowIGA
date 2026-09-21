# Idealized cube dual-tree passive tracer

This case advances a passive concentration tracer through the arterial lumen,
fixed Darcy region, and venous lumen of the
[refined dual-tree FSI case](CUBE_DUAL_TREE_FSI.md). The converged post-FSI
geometry, fluid velocity, and Darcy RT0 face fluxes remain fixed during
transport. Vessel walls are impermeable; transfer occurs only at the four
arterial and four venous terminal interfaces.

The project P1 tetrahedral transport solver uses backward Euler and monotone
graph diffusion in each region. Solved arterial mass outflow is deposited into
adjacent tissue cells as a volumetric source. Tissue outflow divided by the
corresponding water flux defines each venous inlet donor concentration. Every
step checks solver convergence, nonnegative concentration, regional
inventories, both interface transfers, and the global mass balance.

The case is a numerical transport test, not an oxygen-physiology model. It uses
porosity one and zero consumption and does not include hemoglobin binding,
partial pressure, or patient material data.

## Configuration and run

[`cases/idealized_cube_dual_tree_oxygen.json`](../cases/idealized_cube_dual_tree_oxygen.json)
defines:

| Parameter | Value |
|---|---:|
| Arterial inlet concentration | `1 mol/m³` |
| Diffusion coefficient | `10⁻⁶ m²/s` |
| Time step | `5 s` |
| Step range | 10 to 60 steps |
| Stop condition | Mean venous outlet concentration reaches 1% of inlet |

First produce the refined flow run described in the FSI guide, then execute:

```bash
make -C solvers/cpu native_tet_cube_oxygen \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 scripts/run_idealized_cube_dual_tree_oxygen.py \
  /tmp/idealized-cube-refined-fsi-r8 \
  /tmp/idealized-cube-oxygen \
  --case cases/idealized_cube_dual_tree_oxygen.json --ranks 8
```

The runner rebuilds and partitions the geometry, verifies its manifest against
the flow run, and extracts post-FSI velocity, displacement, and Darcy flux by
global ID. If a matching partitioned mesh already exists, pass
`--split-directory PATH` to reuse it. The transport output still receives its
own mesh copy. Output directories must not exist.

## Output and visualization

Open `transport/tissue_oxygen.pvd` in ParaView and color by
`concentration_mol_m3`. Use Clip or Slice with **Surface With Edges** to see the
interior tetrahedra. `transport/artery_oxygen.pvd` and
`transport/vein_oxygen.pvd` share the same time axis. Fix a common color range
when comparing regions or time steps.

The tissue coordinates do not move during the sequence, and both vessel
regions use their fixed post-FSI geometry. The animation shows concentration
transport, not additional wall motion. `ledger.csv` records inventories,
interface transfers, outlet concentration, and balance residuals for every
step. `provenance.json` binds the case, flow summary, mesh, extracted fields,
and solver by SHA-256.

## Recorded validation

Independent two- and eight-rank WSL runs completed 12 steps (`t=0` through
`60 s`, 13 states total). Both reached a mean venous outlet concentration of
approximately `0.01071475877 mol/m³`, or `1.0715%` of the inlet value.

| Quantity | Result |
|---|---:|
| Tissue tetrahedra per state | `101,716` |
| Maximum global balance defect, 8 ranks | `5.84×10⁻²⁰ mol/s` |
| Final tissue nodal concentration range | `0.00155–0.342 mol/m³` |
| Maximum two-/eight-rank field difference | `9.9067×10⁻¹² mol/m³` |

All 39 field arrays were compared by global node ID. A negative test limited
to two steps with a 90% breakthrough target exited with status 2 and did not
publish a successful PVD time index.

```bash
python3 scripts/compare_idealized_cube_oxygen.py \
  /tmp/idealized-cube-oxygen-verified-r2 \
  /tmp/idealized-cube-oxygen-verified-r8
python3 -m unittest scripts.tests.test_idealized_cube_oxygen
make -C solvers/cpu native-tet-moving-species-transport-test
```

The reported breakthrough time belongs to this synthetic parameter set. A
physiological oxygen model would require an oxygen-state definition,
dissolved/bound relationships, wall exchange, tissue porosity, consumption,
material data, and corresponding experimental validation.
