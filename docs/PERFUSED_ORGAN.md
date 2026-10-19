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
search. Cyclic components fall back to Gauss–Seidel, capped at 500 sweeps.

A node counts in `newton_failures` when its Newton iteration ends without
reaching the residual tolerance. That includes a step blocked by the
nonnegativity clip.

The inner pO₂ and acid–base solves leave a round-off floor near 1e-10 in the
scaled residual. So when the iteration stalls, through no descent or a
vanishing update, it counts as converged only below 1e-8. A cyclic component that reaches the sweep cap is counted in
`gauss_seidel_unconverged`. Both solutions are still returned, so callers must
check the counts.

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

- **O₂ uptake** is supply-limited: `v Vmax f(pO₂)`, where
  `f = (pO₂/Km)^n / (1 + (pO₂/Km)^n)`. With `n = 1` (the default) this is
  Michaelis–Menten. With `n > 1`, respiration stays nearly O₂-independent
  above a critical pO₂, so a well-oxygenated organ shows no aerobic shortfall.
- **Glycolysis (Pasteur effect).** The aerobic ATP shortfall is covered by
  glycolysis, up to `v · glycolysis_capacity`, at one ATP per lactate. Glycogen
  is used first, then perfusate glucose. The part still unmet is reported per
  node as `energy_deficit`, an input to the caller's injury model.
- **Acid inhibition of glycolysis (opt-in).** Setting `glycolysis_ph50` scales
  lactate production by `1/(1 + 10^(n (pH50 − pH)))`, where
  `n = glycolysis_ph_hill` (default 1).
  - The inhibition is off by default, which leaves `hepatic_v1` unchanged.
  - `glycolysis_ph50` must lie in the acid–base range [5.5, 9], and `n` must be
    positive.
  - Production unmet by inhibition is added to `energy_deficit`.
  - Without the inhibition, anoxic tissue can exhaust buffer base (see
    Limitations).
  - No values are built in; choose them and record their provenance.
- **Lactate clearance** needs O₂: `v Vmax_L f_O₂ L/(Km_L + L)`. A fraction of the
  cleared lactate becomes glucose.
- **Glycogen flux.** Glycogenolysis releases glucose and is suppressed by
  insulin. Insulin promotes glucose uptake into glycogen. `glycogen_rate` is
  returned per node for the caller to integrate.
- **Aerobic glucose oxidation.** An optional `aerobic_glucose_fraction` of
  the O₂ uptake oxidises glucose, at 6 O₂ per glucose and RQ 1. It draws on
  glycogen first and then on perfusate glucose. The remaining uptake oxidises
  other fuel at `respiratory_quotient`.
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
  The optional keys include `glycolysis_ph50` and `glycolysis_ph_hill`.

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
- **Diagnostics**, including timings: the flow residual and exchange
  imbalance, the transport component counts, `gauss_seidel_sweeps`,
  `gauss_seidel_unconverged` and `newton_failures`.

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

1. The organ flow is solved at the controlled portal pressure. Oedema and
   vasodilator scale the per-cell conductances.
2. The reservoir perfusate passes the gas exchanger. The gas is air blended with
   carbogen (95 % O₂, 5 % CO₂), so the gas pCO₂ rises with FiO₂, matching
   Model41. Arterial O₂ and total CO₂ are efficiency-weighted blends of the
   equilibrium with the exhaust gas and the venous values.
   - Without `sweep_gas_flow_m3_s` the sweep is infinite: the exhaust equals
     the inlet gas.
   - With it, the gas side is well mixed. The exhaust pCO₂ and pO₂ close the
     gas-side balances `pCO₂_out = pCO₂_in + D·q·(C_v − C_a)/G` and
     `pO₂_out = pO₂_in − D·q·(C_a − C_v)/G`, where `D` is the dry gas
     pressure, `q` the organ flow and `G` the molar sweep flow at
     `sweep_gas_temperature_k`. Both are solved by bisection, CO₂ first.
   - CO₂ removal then rises with sweep and saturates at the membrane
     efficiency, while O₂ uptake depends little on sweep. This matches the
     sweep-driven CO₂ removal measured in membrane oxygenators (de Villiers
     Hugo et al., J Extra Corpor Technol 2017).
3. The organ transport is solved with that arterial inflow.
4. `AdvanceWithUptake` integrates the reservoir. It linearises organ uptake as
   first-order clearance and treats production as a constant source, then
   solves each species' linear ODE exactly, including volume change from
   dialysis. The update is therefore stable for any step.
5. The slow tissue states advance:
   - viability falls while the energy deficit exceeds a threshold. Two
     optional modifiers scale that rate. The pH paradox factor
     `(1 + 10^(n(pH50 − 7.4)))/(1 + 10^(n(pH50 − pH)))`, at the node pH, lets
     acidosis protect energy-depleted hepatocytes (Gores et al., Am J Physiol
     1988; Currin et al., FASEB J 1991). The factor `1 + gain·(glucose −
     threshold)₊` lets hyperglycaemia amplify ischaemic injury (Behrends et
     al., J Gastrointest Surg 2010);
   - glycogen integrates its rate;
   - sinusoidal oedema grows with the energy deficit and with tissue pressure
     above a threshold, resolves first-order, and compresses conductance by
     `(1 − E)^n`;
   - the vasodilator is eliminated first-order. With
     `vasodilator_respiration_ic50_mg_l` finite, it also scales the
     mitochondrial respiratory capacity (O₂ uptake and the O₂-dependent lactate
     clearance, not the ATP demand) by `1 − Imax·C/(C + IC50)`, as papaverine
     inhibits complex I (Benej et al., PNAS 2018). `vasodilator_respiration_imax`
     (default 1) is the share of respiration that passes complex I. The unmet
     demand appears as glycolysis and energy deficit.

Each request carries:

- `parameters`: every physical constant, so it is hidden in the truth and
  guessed by the tools;
- `controls`: `inlet_pressure_pa`, `fio2`, dialysis in/out flows and
  dialysate, `infusion_rates`, `bolus`, `vasodilator_bolus_mg`,
  `temperature_c`, and the optional `gas_co2_fraction` and
  `sweep_gas_flow_m3_s`;
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
- the hypoxic volume fraction (pO₂ < 5 mmHg), the deficit volume fraction
  (energy deficit above `injury_deficit_threshold`, which also counts tissue
  whose respiration is inhibited), mean viability, mean oedema, zone-3 pO₂
  and energy-deficit quantiles;
- the organ uptake per species;
- `diagnostics` for that substep's organ solve: `newton_failures`,
  `gauss_seidel_sweeps`, `gauss_seidel_unconverged`, `flow_residual_norm`,
  `exchange_imbalance_m3_s`, and the species `balance_defect`.

The response also totals `newton_failures` and `gauss_seidel_unconverged` over
the substeps. A substep is committed, and written to `state_out`, even when
these counts are nonzero.

The advance rejects a plant state that it could not have produced itself:
- non-finite or negative tissue concentrations;
- viability outside [0, 1];
- negative glycogen;
- oedema outside [0, 1);
- a negative vasodilator concentration.

It also rejects plant parameters that would otherwise be clamped silently or
become NaN:
- gas fractions outside `0 < air O₂ < carbogen O₂ ≤ 1`;
- an oxygenator efficiency or injury threshold outside [0, 1];
- an oedema maximum outside [0, 1);
- negative rates or exponents, including the control-cost rates below;
- non-finite control-cost thresholds;
- a non-positive vasodilator EC50.

`VcaExternalCircuit` extensions (defaults keep previous behaviour):

- a `pressure_control` pump;
- `ApplyBolus`;
- `DialyzerDefinition::inflow_m3_s/outflow_m3_s`;
- `AdvanceWithUptake`.

`Advance` rejects the new modes rather than integrating them explicitly.

Three optional costs of aggressive control are off by default:

| Cost | Trigger | Effect |
|---|---|---|
| Hyperoxic vasoconstriction | arterial pO₂ above `hyperoxic_constriction_threshold_mmhg` | exchange conductance divided by `1 + c·excess/100 mmHg` |
| Oxidative injury | node pO₂ above `oxidative_injury_threshold_mmhg` | viability decays |
| Shear injury | tissue pressure above `shear_injury_threshold_pa` | viability decays |

Without them, maximal oxygenation and pressure carry no penalty. The shear
injury is driven by mean cell tissue pressure, a proxy for sinusoidal shear;
no wall shear stress is computed.

The injury modifiers are also off by default: `injury_ph50` is NaN, and
`injury_glucose_gain_per_mol_m3` is 0.

The optional control `gas_co2_fraction`, in [0, 0.2], sets the sweep-gas CO₂
independently of FiO₂. FiO₂ plus the fraction must not exceed 1. Without it,
the Model41 air/carbogen blend is used, in which gas CO₂ follows FiO₂.

Each record also reports `oxygenator`: the exhaust pO₂ and pCO₂, and
`co2_removal_mol_s`.

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
  - A cell whose buffer base is pinned by the clip counts as a Newton failure.
  - With acid inhibition, the same cell converges with positive buffer base.
  - A strongly recirculating loop that reaches the sweep cap counts as unconverged.
- **Kinetics.** Acid inhibition scales lactate production by exactly its
  factor, leaves the rates unchanged when disabled, and rejects invalid
  parameters.
- **Plant.**
  - The anoxic unit-cube plant reports Newton failures in every substep.
  - With acid inhibition (pH50 6.8, n = 4; and pH50 6.6, n = 6) the same cube
    keeps buffer base positive. Every substep converges and conserves every
    species to 1e-8 of its inflow at every resolution.
  - An oxygenated variant converges and conserves every species to 1e-8 of
    its inflow at every resolution.
  - On that variant the control costs lower flow and viability without
    Newton failures.
  - The sweep exchanger closes both gas-side balances to 1e-9.
    - CO₂ removal rises with sweep.
    - A very large sweep reproduces the infinite-sweep exchanger.
    - In the plant, a lower sweep raises arterial and exhaust pCO₂.
  - With acid inhibition, the pH-paradox factor raises viability and the
    glucose factor lowers it.
  - Respiratory inhibition halves O₂ uptake and lactate clearance at a scale
    of 0.5 and raises lactate production. A vasodilator bolus with a finite
    IC50 lowers the plant's O₂ uptake, and a nonpositive IC50 is rejected.
  - Invalid parameters and states are rejected.

`cpu.native_tet_perfused_organ_cli_test`
(`tests/check_native_tet_perfused_organ.py`) runs the CLI on a synthetic box
organ that it builds with the standard library. It checks:
- the input check and the geometry file;
- solve conservation and diagnostics at every resolution;
- the acid-inhibition keys and the rejection of an out-of-range `glycolysis_ph50`;
- that an advance split through `state_out`/`state_in` reproduces one advance
  bit for bit;
- the per-substep diagnostics;
- rejection of invalid states and parameters, including a negative
  control-cost rate;
- that `gas_co2_fraction` fixes arterial pCO₂ whatever the FiO₂, and that an
  out-of-range fraction is rejected;
- that a lower `sweep_gas_flow_m3_s` raises arterial and exhaust pCO₂, that a
  nonpositive sweep is rejected, and that the injury-modifier keys are
  accepted;
- the `--serve` exit protocol.

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
- **Buffer base can be exhausted.** By default `hepatic_v1` releases one H⁺ per
  lactate regardless of pH. Anoxic tissue can therefore consume buffer base
  faster than perfusion supplies it.
  - Once buffer base reaches zero, the node has no nonnegative solution. The
    clip holds buffer base at zero, the node counts as a Newton failure, and
    the species balance no longer closes.
  - `glycolysis_ph50` lets production fall with pH first. It lowers production
    sharply as pH falls but never stops it, so tissue with almost no
    buffer-base supply can still run out.
  - Treat results with nonzero `newton_failures` as outside the model's
    validity.
