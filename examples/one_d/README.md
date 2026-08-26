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

./solvers/one_d/iga_1d examples/one_d/multispecies_physiology \
  --output-dir /tmp/tubularflowiga-1d-multispecies
```

Run the native VCA PFC closed-loop smoke case with:

~~~bash
./solvers/one_d/iga_1d examples/one_d/vca_pfc_closed_loop \
  --output-dir /tmp/tubularflowiga-1d-vca
~~~

It records the arterial history, aggregated two-outlet venous return, and
external-circuit balance in coupling_manifest.json.

The compliant bifurcation uses the explicit A/Q solver and RCR terminal beds.
The multispecies physiology case uses a pulsatile PETSc pressure network, six
config-selected transported species, wall exchange, metabolism,
oxygen-derived fields, and vasodilation.
The example runs for `0.24 s` with mean inlet flow `4e-7 m3/s`, giving a mean
cell velocity near `0.095 m/s`. This is long enough for the oxygen, glucose,
and lactate inlet fronts to pass the junction and enter both outlet branches.
Open `profile_1d.pvd` in ParaView for the time series.

## Multispecies result

![Animated 1D velocity, oxygen, glucose, and lactate transport](../../docs/images/multispecies-1d-pulse.gif)

At `t=0`, the network starts at oxygen `0.12`, glucose `5.0`, and lactate
`1.0`; the inlet prescribes `0.14`, `6.0`, and `0.5`, respectively. At
`t=0.24`, the terminal cells in both branches reached approximately oxygen
`0.132`, glucose `5.59`, and lactate `0.705`. The common color ranges across
all GIF frames make that inlet-to-branch transport visible without per-frame
rescaling.

Reproduce the README animation from the repository root with:

```bash
pvbatch scripts/render_multiphysics_gif.py --dimension 1d \
  --transport /tmp/tubularflowiga-1d-multispecies/profile_1d.pvd \
  --output docs/images/multispecies-1d-pulse.gif \
  --title "1D pulse multispecies physiology" --frames 12
```

| Case | Scheme | Main check |
|---|---|---|
| `rigid_straight` | rigid `steady_poiseuille` | Poiseuille pressure drop and equal inlet/outlet flow |
| `compliant_bifurcation` | compliant `explicit_rusanov` | pulsatile A/Q update, junction loss, and two RCR states |
| `multispecies_physiology` | PETSc `pressure_network` | pulse-driven six-species transport, Robin oxygen exchange, metabolism, blood-gas arrays, and vasodilation |
| `liver_vein_obj_segment` | rigid `steady_poiseuille` | radius-annotated line OBJ parsing on vertices 1–12 of the liver input |

Generated normalized SWC, skeleton VTP, simulation CSV, VTP/PVD, summary, and
checkpoint files belong in a separate output directory and are ignored by
Git. The full schema, SI units, PETSc/MPI usage, restart contract, and manual
Hex concept map are in
[`docs/ONE_D.md`](../../docs/ONE_D.md).
