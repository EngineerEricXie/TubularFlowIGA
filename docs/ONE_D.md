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

For the shared multi-case entry point, set `CASE_ROOT` and `RANKS` in
`execution.conf` and run `./scripts/run_cases.sh CASE_NAME`. The runner reads
the validated execution plan, bypasses the 3D mesh/database pipeline, invokes
`iga_1d`, and writes the completed run below `CASE/generated/` (or the selected
`OUTPUT_ROOT`). `BACKEND=cuda` is rejected for 1D before launch.

`--check` parses the complete schema, reads and validates the network, resolves
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
undirected network. `root_node_id` optionally selects a 1-based OBJ vertex;
otherwise the largest-radius terminal is the root. Connected cyclic OBJ graphs
are supported only by `implicit_petsc` + `implicit_1d_pde`, currently with no
transport system and `junctions.loss_model: "none"`. Explicit 1D, 0D, SWC,
and the 3D mesh pipeline remain tree-only. See the strict
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
| `compliant` + `explicit_rusanov` | no `formulation` key | Conservative finite-volume A/Q equations, Rusanov flux, friction, internal CFL substeps, and linear or Olufsen wall law |
| `compliant` + `implicit_petsc` | `implicit_1d_pde` | Multi-cell nonlinear implicit A/Q network |

Only formulations that resolve A/Q along centerline cells are classified as
1D. Steady Poiseuille R networks, compliant pressure networks, and one-state-
per-segment RLC networks are 0D circuits and are documented in
[ZERO_D.md](ZERO_D.md). Their former 1D spellings remain deprecated input
aliases and print a migration warning.

Wall laws are `linear` and `olufsen`.
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
- `windkessel_rc`: two-element resistance/capacitance terminal;
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

### Inlet waveforms

The 1D inlet evaluates its named temporal function at every physical time step.
Unlike the 3D inlet, where a temporal function normally multiplies a spatial
velocity profile, a 1D waveform is the absolute inlet value: use `m3/s` for
`quantity: "flow_rate"` and `m/s` for `quantity: "centerline_velocity"`.
The latter is converted to mean flow assuming a parabolic profile,
`Q = 0.5 A u_centerline`.

The supported temporal-function kinds are:

- `constant`: requires `value`;
- `sinusoid`: requires `mean`, `amplitude`, and positive `period`, with optional
  phase in radians;
- `periodic_table`: requires a positive `period`, relative CSV `file`, and
  `interpolation: "linear"`;
- `fourier`: requires `mean`, positive `period`, and equal-length `cosine` and
  `sine` coefficient arrays, with optional phase in radians.

For example, a measured periodic inflow is configured as:

```json
"temporal_functions": [
  {
    "name": "measured_inlet_flow",
    "kind": "periodic_table",
    "units": "m3/s",
    "period": 0.8,
    "file": "inlet_flow.csv",
    "interpolation": "linear"
  }
]
```

with a case-local CSV such as:

```csv
time,value
0.0,8.0e-7
0.2,1.4e-6
0.5,7.0e-7
```

Sample times must be strictly increasing, start at zero, and lie in
`[0, period)`. Linear interpolation also joins the last sample to the first
sample at the period boundary. The solver resolves the CSV inside the case
directory and caches it on first use, so it is not reread at every time step.

For Fourier input, coefficient index 0 is the first harmonic. With
`theta = 2*pi*(t mod period)/period + phase`, the evaluated value is
`mean + sum(cosine[n]*cos((n+1)*theta) + sine[n]*sin((n+1)*theta))`.

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
- `outlet_timeseries.csv`: terminal inflow/pressure, capacitor pressure,
  distal flow, capacitor storage rate, and RC/RCR parameters;
- `branch_timeseries.csv`: segment flow, endpoint pressure, and resistance;
- `profile_1d.csv`: cell A/Q/pressure/velocity;
- `species_profile_1d.csv` and `derived_profile_1d.csv`;
- `profile_1d.vtkhdf`: the default compressed ParaView time series; geometry
  and topology are stored once, while A/Q/pressure/velocity, transported
  species, and derived fields are appended at each output time;
- `summary.json` with model, timing, completion, sampled conservation, and peak
  RSS;
- `physiology_fields.json` with solved/derived/skipped status.

Open `skeleton.vtp` for the static network or `profile_1d.vtkhdf` for simulated
fields. The temporal file also includes static segment/node IDs, cell-local
position, and reference radius. Choose a point-data array and use **Tube** when
a finite-width rendering is desired.

VTKHDF is the default because one file contains every output time without
repeating the network geometry. Use `--visualization-format vtp` only when a
legacy ParaView workflow requires `profile_1d_*.vtp` plus `profile_1d.pvd`.
Checkpoint and CSV files remain independent of the visualization selection.
On restart, an existing compatible VTKHDF file is safely appended; the restart
time replaces an identical final time instead of creating a duplicate step.

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
