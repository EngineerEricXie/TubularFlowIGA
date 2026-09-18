# Native 0D R/RC/RLC flow networks

The native 0D solver solves lumped hydraulic circuits on an SWC or
radius-annotated line-OBJ tree. Model names follow the circuit equations,
not the dimension of the input centerline. It uses the electrical analogy

| Circuit | Flow network |
|---|---|
| voltage | pressure, Pa |
| current | volume flow, m3/s |
| resistance | viscous resistance, Pa s/m3 |
| capacitance | vascular compliance, m3/Pa |
| inductance | fluid inertance, Pa s2/m3 |

Four canonical 0D formulations are available:

| Formulation | States and terms | Time model |
|---|---|---|
| `steady_r` | segment R, nodal pressure | steady |
| `transient_rc` | segment R, nodal C and pressure | backward Euler |
| `transient_rlc` | segment R/L and flow, nodal C and pressure | linearized backward Euler |
| `nonlinear_rlc` | pressure-dependent R/L, flow, nodal C and pressure | nonlinear SNES |

`steady_r` is the canonical replacement for the old
`dimension: "1d"`, `rigid + steady_poiseuille` spelling. The equations and
Poiseuille resistances are unchanged; only their physical classification is
corrected. `transient_rlc` and `nonlinear_rlc` similarly replace the old
`linearized_aq` and `nonlinear_aq` names. Deprecated spellings remain readable
and produce a migration warning.

| Deprecated input | Canonical input |
|---|---|
| `dimension: "1d"`, `rigid`, `steady_poiseuille` | `dimension: "0d"`, `lumped`, `petsc`, `steady_r` |
| `pressure_network` or `lumped_rc` | `transient_rc` |
| `linearized_aq` | `transient_rlc` |
| `nonlinear_aq` | `nonlinear_rlc` |

RCR denotes resistor-capacitor-resistor and does not contain an inductance.
Segment inertance is present only in the RLC formulations.

In `transient_rc`, each segment is a resistor. Its total compliance is split
equally between its two endpoint nodes. At every non-Dirichlet node, backward
Euler solves

```text
C_i (P_i^(n+1)-P_i^n)/dt
  + sum_j (P_i^(n+1)-P_j^(n+1))/R_ij
  + Q_terminal^(n+1) = Q_source^(n+1).
```

The resulting matrix is sparse, symmetric for pressure-terminal cases, and
assembled in O(nodes + segments) work. PETSc distributes its rows over MPI
ranks with exact diagonal/off-diagonal nonzero preallocation. The topology,
matrix, vectors, KSP, and preconditioner workspace are retained across time
steps. OpenMP support used by the shared 1D runtime remains available.

The RLC formulations additionally retain one flow state per segment and solve

```text
L_ij (Q_ij^(n+1)-Q_ij^n)/dt + R_ij Q_ij^(n+1)
  = P_i^(n+1)-P_j^(n+1),
```

with geometric inertance `L_ij = density * length / area`. They remain 0D
because a segment has one lumped flow state; no axial field is resolved inside
the segment. `implicit_1d_pde` becomes 1D by expanding every segment into its
configured number of spatial cells.

## Build and run

```bash
make zero-d-petsc
make zero-d-test

./solvers/one_d/iga_0d examples/zero_d/rc_rcr_bifurcation --check
./solvers/one_d/iga_0d examples/zero_d/rc_rcr_bifurcation

mpiexec -np 2 ./solvers/one_d/iga_0d \
  examples/zero_d/rc_rcr_bifurcation \
  -ksp_type cg -pc_type jacobi -ksp_rtol 1e-10
```

`iga_0d` and `iga_1d` share the validated geometry, boundary, output, and
checkpoint implementation. A 0D case selects the lumped family explicitly
with `dimension: "0d"`, `kind: "network_flow_0d"`, and `model: "lumped"`.
The numerical scheme is named `petsc`; `implicit_petsc` is reserved for the
spatially distributed 1D PDE.

## Flow-system configuration

```json
{
  "name": "lumped_circuit",
  "kind": "network_flow_0d",
  "unknowns": ["area", "flow_rate", "pressure"],
  "model": "lumped",
  "scheme": "petsc",
  "formulation": "transient_rc",
  "dynamic_viscosity": 0.004,
  "density": 1060.0,
  "discretization": {"cells_per_segment": 1},
  "lumped_parameters": {
    "resistance_scale": 1.0,
    "compliance_scale": 0.0,
    "segment_resistance": {"2": 1.0e8},
    "segment_compliance": {"2": 1.0e-10}
  }
}
```

For a steady Poiseuille resistance network, use the same structure with:

```json
"scheme": "petsc",
"formulation": "steady_r"
```

`steady_r` ignores vessel compliance and permits `pressure` or `resistance`
outlets. RC/RCR outlets are transient and therefore require
`transient_rc`, `transient_rlc`, or `nonlinear_rlc`.

The keys in `segment_resistance` and `segment_compliance` are child node IDs;
each identifies the segment joining that child to its parent. Resistance must
be positive and compliance may be zero. Missing resistance values use
Poiseuille resistance multiplied by `resistance_scale`. Missing compliance
values use wall-law volume compliance multiplied by `compliance_scale`.
Setting `compliance_scale` to zero and specifying selected segment
compliances creates a fully explicit lumped circuit.

## Terminal models

The following outlet closures are supported:

- `pressure`: fixed terminal pressure;
- `resistance`: `P = P_reference + R Q`;
- `windkessel_rc`: two-element parallel RC model;
- `windkessel_rcr`: proximal resistance in series with a parallel distal
  resistance and capacitor.

An RC terminal is configured as:

```json
{
  "field": "pressure",
  "type": "windkessel_rc",
  "resistance": 1.0e9,
  "capacitance": 1.0e-10,
  "reference_pressure": 0.0,
  "initial_pressure": 0.0
}
```

An RCR terminal uses `proximal_resistance`, `distal_resistance`, and
`capacitance`. For both models, backward Euler advances the capacitor pressure
without an artificial update at time zero. Checkpoints preserve terminal
pressure, capacitor pressure, flow, and every network pressure/flow state.

## Outputs and validation

In addition to the shared branch, profile, VTKHDF, summary, and checkpoint
files, `node_timeseries.csv` records every circuit-node pressure and
`outlet_timeseries.csv` records terminal inflow, terminal pressure,
capacitor pressure, distal flow, and capacitor storage rate. The visualization
file is named `profile_0d.vtkhdf`; spatially distributed 1D cases use
`profile_1d.vtkhdf`. Use `--visualization-format vtp` for the legacy PVD/VTP
collection. The flow summary reports both network balance

```text
inlet - terminal inflow - vascular storage
```

and total-circuit balance

```text
inlet - distal outflow - vascular storage - terminal capacitor storage.
```

The committed mixed RC/RCR bifurcation regression reaches relative balance
residuals below `1e-13`, and checkpoint/restart reproduces the uninterrupted
final branch and outlet states byte-for-byte.
