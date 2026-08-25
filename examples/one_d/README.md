# Native one-dimensional examples

These cases use schema version 3 and run directly from their committed SWC and
`simulation_config.json`; they do not require the 3D mesh-generation or spline
pipeline, Python, FEniCS, or HexSim.

```bash
make one-d-petsc

./solvers/one_d/iga_1d examples/one_d/rigid_straight --check
./solvers/one_d/iga_1d examples/one_d/rigid_straight \
  --output-dir /tmp/tubularflowiga-1d-rigid

./solvers/one_d/iga_1d examples/one_d/compliant_bifurcation \
  --output-dir /tmp/tubularflowiga-1d-compliant

./solvers/one_d/iga_1d examples/one_d/vca_transport \
  --output-dir /tmp/tubularflowiga-1d-vca
```

The compliant bifurcation uses the explicit A/Q solver and RCR terminal beds.
The VCA-style case uses the PETSc pressure-network formulation, six transported
species, wall exchange, metabolism, oxygen-derived fields, and vasodilation.
Open `profile_1d.pvd` in ParaView for the time series.

| Case | Scheme | Main check |
|---|---|---|
| `rigid_straight` | rigid `steady_poiseuille` | Poiseuille pressure drop and equal inlet/outlet flow |
| `compliant_bifurcation` | compliant `explicit_rusanov` | pulsatile A/Q update, junction loss, and two RCR states |
| `vca_transport` | PETSc `pressure_network` | six conservative species, Robin oxygen exchange, metabolism, blood-gas arrays, and vasodilation |

Generated CSV, VTP/PVD, summary, and checkpoint files belong in a separate
output directory and are ignored by Git. The full schema, SI units, PETSc/MPI
usage, restart contract, and manual Hex concept map are in
[`docs/ONE_D.md`](../../docs/ONE_D.md).
