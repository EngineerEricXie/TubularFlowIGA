# Native tetra transient transport CLI

`solvers/cpu/native_tet_transient_transport` runs transient P2/P1
Taylor–Hood Navier–Stokes on a rigid labelled tetrahedral mesh, driven by a
prescribed inflow waveform, and transports one or more passive species with
each step's P2 velocity (one-way coupling). It uses the project-owned native
FEM. PETSc supplies only algebra and MPI.

It complements the hydraulic graph CLI:

| | `native_tet_hydraulic_graph` | `native_tet_transient_transport` |
|---|---|---|
| Inflow | 0D source reservoir | Prescribed Q(t) waveform |
| Outlets | Terminal RCRs | Prescribed outlet pressure |
| Species | Backward Euler, weak inflow, constant inlet concentration | Backward Euler or BDF2, weak or strong Dirichlet inlet, time-dependent c_in(t) |

The waveform form suits dataset generation and verification against other
solvers with the same boundary conditions.

```bash
make -C solvers/cpu native_tet_transient_transport PETSC_DIR=/path/to/petsc PETSC_ARCH=arch
mpiexec -np 2 solvers/cpu/native_tet_transient_transport case.json --check-input
mpiexec -np 8 solvers/cpu/native_tet_transient_transport case.json --output-dir run_001
# Continue an interrupted run (same case file and rank count):
mpiexec -np 8 solvers/cpu/native_tet_transient_transport case.json --output-dir run_001 --restart
```

## Case (schema version 1)

```json
{
  "schema_version": 1,
  "mesh_file": "volume.msh",
  "boundaries": {"inlet_label": 1, "outlet_labels": [2, 3], "wall_labels": [0]},
  "fluid": {"density_kg_m3": 1060.0, "dynamic_viscosity_pa_s": 0.0035,
            "nonlinear_tolerance": 1e-9, "maximum_iterations": 20},
  "inflow": {"profile": "parabolic", "flow_m3_s": [[0.0, 1.5e-8], [0.3, 2.0e-8]],
             "startup_ramp_s": 1.0},
  "outlets": {"pressure_pa": 0.0},
  "time": {"dt_s": 0.1, "steps": 300},
  "species": [
    {"id": "tracer", "initial_concentration_mol_m3": 0.0, "diffusivity_m2_s": 1e-5,
     "time_integration": "bdf2", "inlet_condition": "dirichlet",
     "inflow_concentration_mol_m3": [[3.0, 0.0], [3.5, 1.0], [13.5, 1.0], [14.0, 0.0]]}
  ],
  "output": {"every_steps": 1, "precision": "float32", "visualization_steps": [100, 300]}
}
```

### Boundaries and flow

- **Label assignment.** Every mesh boundary label must be assigned exactly
  once, as the inlet, an outlet or a wall.
- **Inlet cap.** The inlet cap must be planar. Its area-weighted centroid and
  normal, and its largest in-plane vertex distance (the rim radius), define
  the parabolic profile `2 (1 − r²/R²)` on the P2 inlet nodes. Wall labels
  override shared rim nodes.
- **Inflow scaling.** The profile is normalized by its discrete inward flux,
  so the prescribed flow is matched exactly. It is then scaled by `Q(t)` and
  by the optional cosine start-up ramp `½(1 − cos(π t / t_ramp))`.
- **Tables.** `flow_m3_s` and every species `inflow_concentration_mol_m3`
  accept a constant or `[time_s, value]` pairs. Values are interpolated
  linearly and held constant outside the samples.
- **Outlets.** Outlets receive the prescribed gauge pressure as natural
  traction. The solver works in gauge pressure. The optional
  `backflow_stabilization_beta` enables kinetic-energy backflow
  stabilization; it is off by default.

### Species options

| Key | Meaning |
|---|---|
| `time_integration` | `backward_euler` (default) or `bdf2`; BDF2 takes the first step with backward Euler |
| `inlet_condition` | `inflow_flux` (default, weak advective inflow) or `dirichlet` (strong; see `docs/validation/T7_NATIVE_MOVING_SPECIES.md`) |
| `outlet_backflow_concentration_mol_m3` | Concentration carried by any reverse flow through an outlet face (default 0) |
| `source_mol_m3_s`, `first_order_decay_rate_s_inv` | Constant volumetric source and linear decay; a unit source with zero inflow gives the mean age |
| `wall_exchange_by_label` | Linear wall exchange (Robin) per wall label: `transfer_coefficient_m_s`, `external_concentration_mol_m3` |
| `monotone` | Backward-Euler graph diffusion. Without it, Galerkin undershoot is allowed and reported |

## Acceptance gates and output

**Gates.** Every step checks two things and stops otherwise:

- **Flow.** Boundary fluxes are re-evaluated by independent quadrature, and
  the step stops if the net relative imbalance exceeds `1e-5` or the inlet
  flow differs from `Q(t)` by more than `1e-8` relative.
- **Species.** The species budget defect must be within `1e-6` of the
  largest budget term.

**Output** in `--output-dir`:

- `steps.jsonl`: one line per step, with the inflow, the per-label outward
  flow and mean pressure, Newton and Krylov counts, and per species the
  inventory, outward transfer, source, wall exchange, balance defect,
  minimum and maximum, inlet value and scheme;
- `flow_NNNNNN.f32|f64` every `every_steps` steps: interleaved P2 velocity xyz
  in native vertex then edge order, followed by P1 gauge pressure;
- `species_<id>_NNNNNN.f32|f64`: P1 concentration;
- `checkpoint.f64` and `checkpoint.json`, written every step in double
  precision and including the BDF2 history. The checkpoint identity is the
  SHA-256 of the case text and the mesh file contents, plus the rank count;
- PVTU snapshots at `visualization_steps`, with `flow.pvd`. Their
  `pressure_pa` is the same gauge pressure as the flow files, so outlets show
  `pressure_pa` once;
- `series.json` (layout, inlet geometry, precision) and `completed.json`.

`--stop-after-step N` stops after an accepted step for staged execution.
`--restart` continues from `checkpoint.json` and rejects a changed case, mesh
or rank count. It keeps `steps.jsonl` to one line per accepted step, dropping
lines logged after the last checkpoint, and keeps earlier snapshots in
`flow.pvd`.

## Regression

`native_tet_transient_transport_test` (1 and 2 ranks) runs
`tests/check_native_tet_transient_transport.py` on a 5,973-tetra Y bifurcation
generated by `preprocessing/tet/skeleton_to_tet.py` from
`examples/vascular_flow/y_bifurcation`. It has two species:
- a BDF2 tracer with a strong Dirichlet inlet ramp;
- a backward-Euler mean-age field with a unit source.

The test checks the flow gates and the species budgets at every step, the
inlet schedule, the BDF2 start-up, the output sizes and the PVTU series. It
also checks that a run stopped after step 2 and restarted reproduces the
uninterrupted checkpoint, to `1e-10` relative, as well as its step log and
snapshot series. A changed mesh under the same case is rejected on restart,
and a 100 Pa outlet pressure appears once in the VTK pressure.

This is a functional regression, not a convergence study. Mesh and time
refinement are case-specific responsibilities of the user.
