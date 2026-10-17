# Perfused-organ step solver

`solvers/cpu/native_tet_perfused_organ` advances a machine-perfused organ by
one control interval. Rigid 1D vascular trees supply and drain a 3D Darcy
tissue at their terminals, and six perfusate species are carried through the
tissue with organ kinetics. The caller keeps the slow tissue state (viability,
glycogen, previous concentrations) and passes it in each time, so every call can
change pressures, perfusate composition and parameters. Nothing is restarted
from a checkpoint.

**Status: functional.** The equations are verified against analytic limits
(see [Verification](#verification)). The kinetic constants are
literature-order defaults that the caller supplies; they are not calibrated
against perfusion data, and the solver makes no clinical claim.

## Model

### Hemodynamics (steady)

Each tree is a rooted SWC network of Poiseuille segments with conductance
`π r⁴ / (8 μ L)`. Every terminal `k` owns a territory of tissue cells `i` with
weights `w_ki` that sum to one. The unresolved microvasculature between the
terminal and cell `i` carries

```text
q_ki = g_tree · V_k · w_ki · e_i · (p_k − p̄_i)      [m³/s, vessel → tissue]
```

where:

- `V_k` is the territory volume;
- `p̄_i` is the cell's mean nodal pressure;
- `e_i > 0` is an optional per-cell exchange multiplier for local
  microvascular state.

The tissue obeys P1 Darcy flow, `−∇·(K∇p) = Σ_k q_ki / V_i`, with zero flux on
its boundary. Tree nodes, tissue nodes and exchange terms form one symmetric
positive-definite system. Its only Dirichlet data are the root pressures. By
default MUMPS Cholesky solves it; `-organ_flow_ksp_type` and
`-organ_flow_pc_type` override that. Conservative RT0 face flows are recovered
with the existing Darcy routine.

### Transport (backward Euler over the interval)

Every cell is split into `lobular_zones` well-mixed zones in series. This is a
tanks-in-series model of the sinusoid; with 3 zones it mirrors hepatic zones
1–3.

- Exchange inflow carries the inlet concentration into zone 1.
- Inflow from neighbouring cells also enters zone 1.
- Exchange outflow and outflow to neighbours leave from the last zone.

Each species `s` in each zone node satisfies

```text
φ_s V (c − c_old)/Δt + (F_out + X_out) c − Σ_j F_ji c_j − X_in c_inlet = V R_s(c, state)
```

When `Δt` far exceeds the residence time, this reaches the quasi-steady state.
Unperfused nodes relax over `Δt` instead of having no steady state.

The upwind graph is ordered by strongly connected components. Each node is
solved by a local Newton iteration on the six species with a nonnegativity line
search. Cyclic components fall back to Gauss–Seidel.

### Species and acid–base

| Species | Unit | Meaning |
|---|---|---|
| `o2_total_mol_m3` | mol/m³ | dissolved + haemoglobin-bound O₂ |
| `glucose_mol_m3` | mol/m³ | |
| `lactate_mol_m3` | mol/m³ | |
| `buffer_base_mol_m3` | mol/m³ | HCO₃⁻ + non-bicarbonate buffer anions |
| `co2_total_mol_m3` | mol/m³ | HCO₃⁻ + dissolved CO₂ |
| `insulin_mU_L` | mU/L | |

O₂ content reuses the shared `OxygenCapacity` model:
`C = α pO₂ + 1.34·Hb·S(pO₂)`, where `S` is a Hill curve. The P50 is shifted by
the Bohr effect, and pO₂ is recovered by safeguarded Newton inversion.

pH comes from solving `HCO₃⁻ (1 + 10^(pK−pH)) = TCO₂`, with:

- `HCO₃⁻ = BB − A0 − β (pH − 7.4)`;
- `β = 1.43 Hb + 7.7` mmol/(L·pH) (Van Slyke);
- `A0 = 17.7 + 0.42 Hb` mmol/L, both by default (Hb in g/dL).

Strong acid (lactate) lowers buffer base one-for-one. CO₂ production raises
total CO₂, and most of that becomes HCO₃⁻ through the non-bicarbonate buffer.

### Kinetics `hepatic_v1`

Rates are per tissue volume. Each scales with the local viable fraction `v`
and with `Q10^((T−37)/10)`.

- **O₂ uptake** is supply-limited: `v Vmax pO₂/(Km + pO₂)`.
- **Glycolysis (Pasteur effect).** The aerobic ATP shortfall is covered by
  glycolysis, up to `v · glycolysis_capacity`, at one ATP per lactate. Glycogen
  is used first, then perfusate glucose. The part still unmet is reported per
  node as `energy_deficit`, an input to the caller's injury model.
- **Lactate clearance** needs O₂: `v Vmax_L f_O₂ L/(Km_L + L)`. A fraction of the
  cleared lactate becomes glucose.
- **Glycogen flux.** Glycogenolysis releases glucose and is suppressed by
  insulin. Insulin promotes glucose uptake into glycogen. `glycogen_rate` is
  returned per node for the caller to integrate.
- **Insulin** is cleared first-order.
- **CO₂** is produced at the respiratory quotient.

## Case file (static)

```json
{
  "schema_version": 1,
  "kind": "perfused_organ",
  "tissue_mesh_file": "liver_h5.msh",
  "mesh_length_scale_to_m": 0.001,
  "trees": [
    {"id": "portal", "role": "inflow", "skeleton_file": "portal.swc", "length_scale_to_m": 0.001},
    {"id": "hepatic", "role": "outflow", "skeleton_file": "hepatic.swc", "length_scale_to_m": 0.001}
  ],
  "territories_file": "territories_h5.json",
  "reaction_model": "hepatic_v1",
  "lobular_zones": 3
}
```

- **Tissue mesh.** It is labelled ASCII Gmsh 4.1. Every boundary face must be
  labelled, and the volume physical name is `fluid`.
- **Territories file.**
  `{"schema_version": 1, "trees": {"<tree id>": {"<SWC terminal id>": {"cell_ids": [...], "weights": [...]}}}}`.
  The cell ids are Gmsh element ids, and the weights of each terminal sum to
  one.
- **Trees.** At least one inflow tree and one outflow tree are required.

## Request (one per call)

Required keys:

- **Interval and flow:** `dt_s`, `viscosity_pa_s`, `mobility_m2_pa_s`.
- **Per tree:** `root_pressure_pa {tree: Pa}` and
  `specific_conductance_inv_pa_s {tree: 1/(Pa s)}`.
- **Per species:** `inlet {species: value}` and
  `distribution_volume_fraction {species: value}`.
- `chemistry {hemoglobin_g_dl, ...}`. These keys follow the shared `OxygenCapacityParameters`.
- `kinetics {vmax_o2_mol_m3_s, glycolysis_capacity_mol_m3_s,
  lactate_vmax_mol_m3_s, glycogenolysis_vmax_mol_m3_s,
  glucose_uptake_vmax_mol_m3_s, insulin_clearance_s_inv, temperature_c, ...}`.

Optional keys:

- `initial`, plus the scalar defaults `viability` and `glycogen_mol_m3`.
- `resolution`: `cell` (default), `portal_territory`, or `organ`.
- `field_input {file, names}` and `field_output {file, names}`.
- `vtk_file`, `hypoxic_po2_mmhg` (default 5), `tag`.

Field files are raw native-endian float64, one array per name, in order:

| Field kind | Names | Length |
|---|---|---|
| Per cell | `mobility_scale`, `exchange_scale`, `tissue_pressure_pa`, `exchange_in_m3_s`, `exchange_out_m3_s`, `perfusion_s_inv`, `darcy_flux_m_s` | cells |
| Per zone node | species, `viability`, `glycogen_mol_m3`, `po2_mmHg`, `ph`, `pco2_mmHg`, `bicarbonate_mol_m3`, `o2_uptake_mol_m3_s`, `lactate_production_mol_m3_s`, `lactate_clearance_mol_m3_s`, `glycogen_rate_mol_m3_s`, `energy_deficit` | cells × zones, ordered by cell then zone (`node = cell·zones + zone`) |

The response is one JSON object. It contains:

- **Per tree:** root and terminal flows, terminal pressures, and the mixed
  concentration of the flow from tissue into that tree. For the outflow tree
  this is the organ's venous effluent.
- **Per species totals:** inflow, outflow, storage rate, reaction, and balance
  defect.
- **Tissue quantiles**, volume-weighted:
  - pO₂, pH, pCO₂, energy deficit, lactate, glucose, perfusion, and tissue
    pressure;
  - the last zone's pO₂ and energy deficit;
  - the hypoxic volume fraction.
- **Diagnostics**, including timings.

Modes:

```bash
native_tet_perfused_organ CASE.json --check-input
native_tet_perfused_organ CASE.json --describe geometry.f64   # centroids, volumes, cell ids
native_tet_perfused_organ CASE.json --request request.json --output response.json
native_tet_perfused_organ CASE.json --serve                   # one request JSON per stdin line
```

## Plant advance (`"command": "advance"`)

The full perfusion plant reuses the shared `VcaExternalCircuit` (reservoir,
infusion, dialyzer, bolus) with the backward-compatible extensions described
below. One request advances the plant by `duration_s` in `substeps` organ
solves.

Each substep runs in this order:

1. The reservoir perfusate passes the gas exchanger. The gas is air blended with
   carbogen (95 % O₂, 5 % CO₂), so the gas pCO₂ rises with FiO₂, matching
   Model41. Arterial O₂ and total CO₂ are efficiency-weighted blends of the gas
   equilibrium and the venous values.
2. The organ is solved at the controlled portal pressure. Oedema and
   vasodilator scale the per-cell conductances.
3. `AdvanceWithUptake` integrates the reservoir. It linearises organ uptake as
   first-order clearance and treats production as a constant source, then
   solves each species' linear ODE exactly, including volume change from
   dialysis. The update is therefore stable for any step.
4. The slow tissue states advance:
   - viability falls while the energy deficit exceeds a threshold;
   - glycogen integrates its rate;
   - sinusoidal oedema grows with the energy deficit and with tissue pressure
     above a threshold, resolves first-order, and compresses conductance by
     `(1 − E)^n`;
   - the vasodilator is eliminated first-order.

Each request carries:

- `parameters`: every physical constant, so it is hidden in the truth and
  guessed by the tools;
- `controls`: `inlet_pressure_pa`, `fio2`, dialysis in/out flows and
  dialysate, `infusion_rates`, `bolus`, `vasodilator_bolus_mg`,
  `temperature_c`;
- the state: either `initial_state`, or `state_in` with `json` and `fields`
  paths;
- optional `state_out` and `heterogeneity` (per-cell `mobility_scale` and
  `exchange_scale`);
- `resolution`.

The state field file holds, in order:

1. the six tissue species, per node;
2. viability, per node;
3. glycogen, per node;
4. oedema, per cell.

The response holds one record per substep, each with:

- flow, pressure, reservoir volume and haematocrit;
- the arterial and venous O₂, saturation, pH, pCO₂, HCO₃⁻ and species;
- the hypoxic volume fraction, mean viability, mean oedema, zone-3 pO₂ and
  energy-deficit quantiles;
- the organ uptake per species.

`VcaExternalCircuit` extensions (defaults keep previous behaviour):

- a `pressure_control` pump;
- `ApplyBolus`;
- `DialyzerDefinition::inflow_m3_s/outflow_m3_s`;
- `AdvanceWithUptake`.

`Advance` rejects the new modes rather than integrating them explicitly.

## Geometry from a DICOM SEG

```bash
python3 scripts/organ/seg_to_masks.py seg.dcm masks.npz --segment Liver --segment Portal \
    --segment Hepatic --segment Tumor_1
python3 scripts/organ/mask_to_tissue_mesh.py masks.npz liver_h5.msh --mask Liver \
    --fill Tumor_1 --size-mm 5 --manifest liver_h5.json
python3 scripts/organ/vessel_mask_to_swc.py masks.npz portal.swc --mask Portal \
    --spacing-mm 2 --min-radius-mm 1.0     # needs scikit-image and networkx
python3 scripts/organ/build_territories.py liver_h5.msh territories_h5.json \
    --tree portal=portal.swc --tree hepatic=hepatic.swc
```

- **Tissue mesh.** It block-averages the mask onto cubes and keeps the largest
  face-connected set. Each cube is split into six Kuhn tetrahedra, which
  conforms across neighbours and keeps every Jacobian positive.
- **Territories.** Each cell goes to its nearest terminal, which works for any
  organ. Check the reported distances for strongly non-convex organs.
- **Patient data.** Keep patient-derived outputs out of the repository.

## Verification

`make -C solvers/cpu perfused-organ-test` runs `perfused_organ_test`, which
checks the following:

- **Flow.**
  - Organ inflow equals outflow within 1e-9.
  - Terminal flows sum to the root flow.
  - The Darcy cell balance holds.
  - With near-infinite mobility, the network reduces to the analytic
    series/parallel resistance.
  - A uniform exchange multiplier equals scaling the conductances.
- **Chemistry.**
  - The O₂ content curve inverts correctly.
  - Henderson–Hasselbalch holds at the reference state.
  - The buffer response to CO₂ and to strong acid is correct.
- **Kinetics.** The anoxic glycolysis cap, the energy deficit and the zero rates
  of dead tissue match their limits.
- **Transport.**
  - A single perfused cell matches a scalar root solve for O₂.
  - A stagnant cell relaxes as backward Euler predicts.
  - Three zones reproduce tanks in series exactly.
  - Inert species stay at the inlet value.
  - Species balance holds for cell, zone, and pooled graphs.

## Limitations

- **No sub-mesh vessel geometry.** Exchange happens only at tree terminals, and
  the microvascular conductances `g_tree` are calibration parameters.
- **Zonation is only a sub-grid model.** The mesoscale Darcy flow redistributes
  between cells, but lobular zonation comes only from the tanks-in-series zones.
- **Simplified tissue physics.** There is one tissue compartment. There is no
  oxygen diffusion between zones and no tissue deformation.
- **Temperature is a rate scale only.** No heat transport is solved.
- **Uncalibrated parameters.** The kinetic constants must be supplied, and
  their provenance must be recorded by the caller.
