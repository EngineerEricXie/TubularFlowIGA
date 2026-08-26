# Native C++ 1D flow and transport

The `iga_1d` executable solves rooted vascular networks directly from SWC. It
is a C++17/PETSc implementation: Python, FEniCS, petsc4py, HexSim, 3D mesh
generation, and `.ntiga` packing are not runtime dependencies.

Use 1D when network-scale wave propagation, terminal-bed response, flow split,
or many transported species matter more than a resolved 3D velocity field. Use
the 3D IGA path for local three-component velocity, pressure, and wall-resolved
geometry. The current release does not couple the two models and does not
include FSI or a 1D CUDA backend.

## Build and run

```bash
./scripts/check_dependencies.sh one-d
make one-d-petsc
make one-d-test

./solvers/one_d/iga_1d CASE_DIR --check
./solvers/one_d/iga_1d CASE_DIR --system blood_flow_1d
./solvers/one_d/iga_1d CASE_DIR --output-dir OUTPUT_DIR
```

`--check` parses the complete schema, reads and validates the SWC tree, resolves
topological boundaries, and checks referenced node IDs without advancing time.
Without `--system`, a case must have exactly one 1D flow system; the solver then
runs every transport system whose `flow_system` names it.

PETSc options remain command-line options rather than JSON keys:

```bash
mpiexec -np 2 ./solvers/one_d/iga_1d CASE_DIR \
  -snes_rtol 1e-9 -ksp_rtol 1e-10
```

The implicit nonlinear solver uses LU on one rank and distributed MUMPS LU on
multiple ranks by default. A PETSc installation used for multi-rank nonlinear
solves must therefore provide MUMPS, or the defaults must be replaced with a
working site-specific `-ksp_*`/`-pc_*` configuration. Explicit flow and
multi-species updates use OpenMP where the problem size benefits; one MPI rank
is normally sufficient for the small explicit examples.

## Skeleton and SI contract

An SWC geometry block is:

```json
"geometry": {
  "kind": "swc_network",
  "file": "skeleton_initial.swc",
  "length_scale_to_m": 1.0e-6
}
```

Each non-comment SWC row has `id type x y z radius parent_id`. Exactly one row
must have parent `-1`; every other parent must exist, and the graph must be an
acyclic tree connected to that root. IDs are preserved but need not be
contiguous. Coordinates, segment lengths, and radii must be finite and
positive after multiplication by `length_scale_to_m`. Zero-length segments,
cycles, duplicate IDs, missing parents, and invalid configured boundary IDs are
rejected before simulation.

Radius-annotated line OBJ is selected with `kind: "obj_network"`. Its `v`
records carry `x y z radius auxiliary auxiliary`, while `l` records define the
undirected tree. `root_node_id` optionally selects a 1-based OBJ vertex;
otherwise the largest-radius terminal is the root. See the strict
[skeleton-format contract](SKELETON_FORMATS.md).

The 1D schema uses SI for flow:

| Quantity | Unit |
|---|---|
| coordinates, radius, wall thickness | m |
| area | m² |
| time | s |
| flow rate | m³/s |
| centerline velocity | m/s |
| pressure | Pa |
| density | kg/m³ |
| dynamic viscosity | Pa·s |
| resistance | Pa·s/m³ |
| RCR capacitance | m³/Pa |
| scalar diffusivity | m²/s |

Scalar concentration and physiology values use the units chosen for that
species, but initial, inlet, source, reaction, wall-exchange, metabolism, and
derived-field parameters must be mutually consistent. The blood-gas helper
uses mmHg for `p50_mmhg`, `pO2`, and `pCO2`, and the configured conventional
oxygen solubility/Hb capacity values for dissolved, bound, and total oxygen.

## Schema v3

Schema v2 remains the 3D case format. Native 1D cases set:

```json
{
  "schema_version": 3,
  "dimension": "1d",
  "geometry": {"kind": "swc_network", "file": "skeleton_initial.swc", "length_scale_to_m": 1.0},
  "fields": [],
  "time": {"dt": 0.001, "steps": 100, "output_every": 10},
  "temporal_functions": [],
  "equation_systems": [],
  "boundaries": []
}
```

Unknown keys are errors. `fields`, named temporal functions, time controls, and
named equation systems follow the same design as the 3D configuration, while
the equation and boundary kinds below are specific to a topological network.

### Flow schemes

| Model and scheme | Formulation | Purpose |
|---|---|---|
| `rigid` + `steady_poiseuille` | no `formulation` key | Segment Poiseuille resistance, downstream reduction, pressure, and flow split; a time-varying inlet is solved quasi-statically |
| `compliant` + `explicit_rusanov` | no `formulation` key | Conservative finite-volume A/Q equations, Rusanov flux, friction, internal CFL substeps, and linear or Olufsen wall law |
| `compliant` + `implicit_petsc` | `pressure_network` | Lumped compliant pressure network |
| `compliant` + `implicit_petsc` | `linearized_aq` | Linearized nodal-pressure/branch-flow system |
| `compliant` + `implicit_petsc` | `nonlinear_aq` | SNES nonlinear nodal-pressure/branch-flow system with analytic sparse Jacobian |
| `compliant` + `implicit_petsc` | `implicit_1d_pde` | Multi-cell nonlinear implicit A/Q network |

`rigid` cannot be paired with a compliant scheme, and `formulation` is accepted
only for `implicit_petsc`. Wall laws are `linear` and `olufsen`.
`cells_per_segment`, `cfl`, momentum correction `alpha`, and
`min_area_fraction` live in `discretization`. The explicit solver terminates on
a non-finite state or an area below the configured physical bound; it does not
silently clamp the solution.

### Inlets, outlets, and junctions

Boundaries use `role: inlet|outlet|wall` and optional SWC `node_ids`, not 3D
mesh-face labels. Omitting outlet IDs applies the condition to all leaf nodes.
An inlet Dirichlet condition selects `quantity: flow_rate` or
`centerline_velocity` and names a temporal function. Outlet types are:

- `pressure`: prescribed terminal pressure;
- `resistance`: terminal resistance and reference pressure;
- `windkessel_rcr`: proximal resistance, distal resistance, capacitance,
  reference pressure, and initial capacitor pressure.

Junctions conserve flow and select `pressure_balance: static|total`.
`loss_model` is `none`, `constant`, `table`, `angle_sin2`, or
`mynard_valen_sendstad`. A constant coefficient can be overridden by SWC node
ID in `node_coefficients`. A `table` may instead contain sorted interpolation
data such as:

```json
"angle_table": [
  {"angle_degrees": 0.0, "coefficient": 0.0},
  {"angle_degrees": 90.0, "coefficient": 0.5}
]
```

`reference_velocity` chooses the parent or child dynamic-pressure scale.

### Transport and physiology

A `network_transport_1d` system names its upstream `flow_system` and lists any
number of scalar species. Each species has diffusivity, nonnegative first-order
reaction rate, and signed volume source. The conservative unknown is `A*C`.
Inlet concentration is Dirichlet; wall exchange is `no_flux`,
`constant_flux`, or `robin` with coefficient and exterior value.

The optional `physiology` block supplies signed metabolism source/sink rates,
oxygen-capacity constants, and a transported vasodilator feedback. Radius
feedback relaxes toward a bounded target derived from the original SWC radius,
then affects the next flow step. Supported requested derived arrays are `pO2`,
`pCO2`, `pH`, `SaO2`, `SvO2`, `dissolved_oxygen`, `bound_oxygen`,
`total_oxygen`, and `hematocrit`. Explicitly requested fields with missing
transport dependencies are configuration errors; the output manifest records
all solved and derived fields.

The `unknowns` and `species` arrays are the config controls for what is
simulated. There is no hard-coded oxygen, glucose, or neuron species set. The
`multispecies_physiology` example selects six fields; another case may select a
smaller or different scalar set and enable physiology helpers only when their
named dependencies are present.

### VCA vascular coupling

Native 1D can expose a stable SI port between the vascular network and an
external VCA circuit. The top-level simulation_scope.mode selects flow_only,
vascular_open_loop, vca_replay, or vca_closed_loop. Replay reads a relative
JSON or CSV inlet history. Closed-loop runs advance a well-mixed reservoir from
the conservative, signed aggregation of all outlet flows and species fluxes.
The currently supported native boundary mode is a positive flow-controlled pump
with explicit_staggered ordering.

The perfusate block selects rbc or pfc and the transported oxygen state
(dissolved_oxygen or total_oxygen). A PFC perfusate must use zero hematocrit
and hemoglobin. Closed-loop runs require every reservoir species, including the
selected oxygen state, to be a transported 1D field. The optional oxygenator,
dialyzer, and infusion blocks operate outside the vascular mesh; their source
terms and reservoir mass changes are included in coupling_manifest.json.

Run the committed two-outlet PFC smoke case with:

~~~bash
./solvers/one_d/iga_1d examples/one_d/vca_pfc_closed_loop \
  --output-dir /tmp/tubularflowiga-1d-vca
~~~

The manifest records the exact arterial history, aggregated venous return,
reservoir/device balance, total vascular mass, source integrals, and residuals.
Closed-loop checkpoint/restart is intentionally rejected until reservoir state
is included in the checkpoint format.

## Checkpoint and restart

```bash
./solvers/one_d/iga_1d CASE_DIR \
  --checkpoint OUTPUT/checkpoint --checkpoint-every 50 --stop-after-step 100

./solvers/one_d/iga_1d CASE_DIR \
  --restart OUTPUT/checkpoint --output-dir OUTPUT/resumed
```

The PETSc binary state and JSON metadata include A/Q/pressure, node pressure,
outlet and RCR capacitor state, all species, dynamic vessel radii, completed
step, physical time, `dt`, and explicit substep count. Restart rejects a changed
config, changed SWC geometry, mismatched species ordering, inconsistent time,
truncated state, and corrupt metadata using config/network fingerprints and
strict size checks.

## Outputs

All generated files go to `--output-dir` (default:
`CASE_DIR/results/one_d/SYSTEM`) and should not be committed:

- `skeleton_normalized.swc`: validated, explicitly rooted skeleton in the
  source coordinate unit;
- `skeleton.vtp`: static skeleton in solver SI coordinates, including radius,
  topology, role, segment, and branch arrays;
- `flow_timeseries.csv`: inlet/outlet flow, sampled storage/continuity
  diagnostic, pressure drop, and area bounds;
- `branch_timeseries.csv`: segment flow, endpoint pressure, and resistance;
- `profile_1d.csv`: cell A/Q/pressure/velocity;
- `species_profile_1d.csv` and `derived_profile_1d.csv`;
- `profile_1d_*.vtp` and `profile_1d.pvd` for ParaView;
- `summary.json` with model, timing, completion, sampled conservation, and peak
  RSS;
- `physiology_fields.json` with solved/derived/skipped status.

Open `skeleton.vtp` for the static network or `profile_1d.pvd` for simulated
fields. Choose a point-data array and use **Tube** with `radius` as an absolute
scalar when a finite-width skeleton rendering is desired.

## Hex/FEniCS concept map

This table is a manual migration aid, not an accepted legacy schema. `iga_1d`
does not read Hex `sim_config.json`, invoke Hex code, or promise matching file
names or step-by-step floating-point values.

| Hex/FEniCS concept | Native schema-v3 location |
|---|---|
| centerline/network file and coordinate scale | `geometry.file`, `geometry.length_scale_to_m` |
| 1D blood-flow model | `equation_systems[].kind: network_flow_1d`, `model`, `scheme`, `formulation` |
| blood density and viscosity | flow system `density`, `dynamic_viscosity` |
| elastic/Olufsen vessel parameters | flow system `wall` |
| axial resolution and CFL | flow system `discretization` |
| inlet flow or centerline-velocity waveform | inlet condition `quantity` + named `temporal_functions` entry |
| terminal pressure/resistance/Windkessel | outlet `pressure`, `resistance`, or `windkessel_rcr` condition |
| junction pressure/loss selection | flow system `junctions` |
| transported biochemical variables | declared scalar `fields` + `network_transport_1d.species` |
| inlet concentration and wall permeability | inlet/wall boundary conditions |
| reaction and volume production/consumption | species `reaction_rate`, `volume_source` |
| metabolism and blood-gas derived values | top-level `physiology` |
| checkpoint/restart | native CLI `--checkpoint`, `--checkpoint-every`, `--restart` |

Start from one of the three source-only cases in [`examples/one_d`](../examples/one_d/README.md)
and transfer values deliberately with explicit SI conversion.
