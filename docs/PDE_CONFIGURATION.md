# Configurable PDE systems

`simulation_config.json` is the canonical case configuration. Schema version 2
describes fields, equation systems, time integration, and named boundary
conditions in one file. Field names such as `tracer`, `oxygen`, or `Nplus` have
no built-in numerical meaning.

Schema version 3 adds required `dimension: "1d"|"3d"`. A v3/3d document uses
the same systems as v2; v2 remains supported unchanged. A v3/1d document uses
the native SWC-network flow, transport, physiology, and topological-boundary
contract documented in the [native 1D guide](ONE_D.md).

## Fields and equation systems

Each field has a unique `name` and a `kind`: `scalar`, `vector3`, or `pressure`.
An equation system lists its unknown fields. The CPU `linear_transport` backend
accepts any positive number of scalar unknowns and composes their weak form from:

- `time_derivative`
- `diffusion`
- `advection` with the `prescribed` velocity field
- `linear_coupling` between any two unknowns
- `volume_source`

Each term names its test equation with `equation`; `trial` defaults to that same
field. Coefficients are numeric and already include any products of physical
parameters. SUPG can be enabled per equation. Unsupported operator names and
unknown field references are rejected before assembly.

Navier–Stokes is selected with `kind: "navier_stokes"`. Its two unknowns are a
`vector3` velocity followed by pressure. `viscosity` is the positive dynamic
viscosity and is required; `density` is positive and defaults to `1`. The
optional `time_integration` is either `steady` (the default) or
`backward_euler`. Backward Euler uses the top-level positive `time.dt` and
`time.steps` values.

The CPU finite-element kernel and solver assemble the physical inertial term
`density * (du/dt + u dot grad(u))`, the Cauchy viscous term using dynamic
viscosity, and pressure. Its stabilization includes the backward-Euler temporal
scale. The CPU solver stores the previous converged state, advances
`time.steps` physical steps, resets nonlinear convergence checks per step, and
evaluates velocity Dirichlet or pressure-traction waveforms at the new-step time. Velocity profiles
are spatial data multiplied by a `scale` and, when configured, the temporal
waveform value.

Schema version 2 optionally accepts named `temporal_functions` of kind
`constant`, `sinusoid`, `periodic_table`, or `fourier`. Each declares
`units`; periodic functions declare a positive `period`. Tables name a
strict `time,value` CSV with samples in `[0, period)` and linear,
periodically wrapped interpolation. Fourier `cosine` and `sine` arrays use
one-based harmonic order. A Dirichlet or `pressure_traction` condition may name
one with `waveform`.
Parsing and dependency-free evaluation are implemented, including negative-time
and period wrapping. CPU and CUDA configured transport materialize Dirichlet
values at each physical step; a constant waveform reproduces constant-boundary results.
CPU backward-Euler Navier–Stokes executes waveform-backed velocity Dirichlet
and pressure-traction conditions. Steady flow still rejects them because a
time-varying
boundary has no meaning without physical time integration. The CUDA transient
kernel, physical-time loop, raw-state checkpoint/restart, and time-indexed
output passed multi-cycle Womersley and restart GPU validation. See the
[CPU](../solvers/cpu/VALIDATION.md) and
[CUDA](../solvers/cuda/VALIDATION.md) validation reports.

## Boundary conditions

Every mesh label is assigned a user-selected name and a list of field-specific
conditions. Supported schema values are `dirichlet`, `no_flux`, `flux`, `robin`,
`advective_outflow`, `pressure_traction`, `resistance`, `windkessel_rc`, and
`windkessel_rcr`. The configured transport assembler executes scalar
Dirichlet conditions and treats `no_flux` and `advective_outflow` as natural
conditions. CPU and CUDA also assemble matching `flux` and `robin` surface
terms from packed element-face labels.
A velocity Dirichlet condition may provide a three-component `value`, or set
`"profile": "initial_velocityfield.txt"` plus a numeric `scale`. The latter
preserves the case-generated vessel inlet profile, as in the HexSim workflow.
For Navier–Stokes, `pressure_traction` supplies the scalar pressure in the
natural boundary force `-p n`; unlike pressure Dirichlet, it retains all
continuity rows. R/RC/RCR outlet models update this same traction through their
3D/0D fixed-point coupling.

### 3D VCA port contract

Schema-v3 3D cases may declare the physical VCA ports explicitly. This avoids
assuming that a particular mesh label is an inlet or that every other label is
an outlet:

```json
"coupling": {
  "three_d_ports": {"inlet_label": 1, "outlet_labels": [2, 3]}
}
```

The inlet must have an `initial_velocityfield.txt` velocity Dirichlet profile;
the bridge scales that profile to the signed SI pump flow. The listed outlets
are the only faces eligible for conservative return aggregation.

`iga_navier_stokes` supports this contract for CPU backward-Euler
`vca_closed_loop`. It advances the pump-controlled flow step, then advances
exactly one in-memory `linear_transport` system with the new flow velocity.
All reservoir species must be transported fields and must have inlet Dirichlet
conditions on the declared VCA inlet. Outlet flow, area-averaged pressure, and
outward species fluxes are MPI-reduced before the external circuit advances.

The runner writes `coupling_manifest.json` next to `--output`, or to
`CASE_DIR/results/vca_flow/` without an output path. The initial implementation
requires `velocity_source: "prescribed"` because it consumes the in-memory
flow state rather than a velocity-file series. Checkpoint/restart remains
rejected until flow, transport, outlet, and reservoir state share one
checkpoint. CUDA and 3D VCA replay/open-loop modes remain unsupported.
Native `iga_1d` remains the more complete VCA backend while those limitations
are removed.

## Time-resolved velocity sources

CPU configured transport can replace the built-in `prescribed` velocity with a
named snapshot series:

```json
{
  "velocity_sources": [
    {
      "name": "flow_snapshots",
      "kind": "snapshot_series",
      "manifest": "velocity_series.csv",
      "interpolation": "linear",
      "out_of_range": "error"
    }
  ]
}
```

Advection and SUPG entries then set `"velocity": "flow_snapshots"`. The
manifest is relative to the case directory and contains strictly increasing
`time,file` records. Snapshot files use the same three-column nodal velocity
format written by `iga_navier_stokes`:

```csv
time,file
0.1,flow/velocity.step000001.txt
0.2,flow/velocity.step000002.txt
```

At each transport step the CPU solver loads the bracketing snapshots, linearly
interpolates every nodal component when the time grids differ, and rebuilds the
advection/SUPG operators. `out_of_range` is either `error` or `hold`. The
original positional velocity file remains the reproducible offline path for
systems that use `prescribed`. The CUDA snapshot-series assembly path passed
CPU parity and checkpoint/restart gates on both straight and branching cases.

The runnable vascular-flow configuration is
[`examples/vascular_flow/straight_tube/simulation_config.json`](../examples/vascular_flow/straight_tube/simulation_config.json).
The runnable
[`examples/neuron_transport/straight_neurite/simulation_config.json`](../examples/neuron_transport/straight_neurite/simulation_config.json)
expresses the previous two-field neuron equations using the same generic
operators. Application-specific examples deliberately do not mix unrelated
flow and transport systems in one configuration.

The `vascular_flow/multispecies_pulse` showcase intentionally combines them.
Its `linear_transport.unknowns` array is the config selection of species to
simulate; each selected field has matching weak-form terms and boundary
conditions. The optional top-level `physiology` block adds signed metabolism
rates to those species and requests output-time blood-gas arrays. Dependencies
are validated before assembly. Three-dimensional vasodilation is rejected
because the current rigid mesh has no FSI; radius feedback remains available
in the 1D solver.

## CPU use

After preparing the mesh database and loading the documented PETSc/MPI setup:

```bash
mpiexec -np 2 ./solvers/cpu/iga_solve \
  case.ntiga CASE_DIR tracer_transport CASE_DIR/tracer.txt
```

Omit the system name only when the case has one linear transport system. The
output contains `node_id` followed by fields in configured order; a neighboring
`.fields` file records their names.

Every requested 3D result also has a ParaView XML companion. A final solve
writes `NAME.vtu`; time-indexed output writes `NAME.stepNNNNNN.vtu` and a
`NAME.pvd` collection. Flow point data contain `velocity` and `pressure`.
Transport point data contain every configured unknown plus requested derived
physiology arrays. Legacy text, `.pressure`, and `.fields` outputs remain for
validation and CPU/CUDA comparison.

`iga_navier_stokes DATABASE CASE_DIR` reads the `navier_stokes` system and its
velocity/pressure boundary values from the same file. `steady` performs one
nonlinear solve; `backward_euler` performs `time.steps` solves with step size
`time.dt` and writes the final state to the optional output path. For transition only, it
still accepts cases containing `simulation_parameter.txt` plus
`case_config.json`; that input is converted to the shared representation rather
than selecting a legacy numerical model.

For transient output and restart, use `--output PATH --output-every N`,
`--checkpoint PREFIX --checkpoint-every N`, and `--restart PREFIX`. The
checkpoint JSON records completed-step and physical-time metadata; the PETSc
state file stores the distributed four-field vector.
`--stop-after-step N` intentionally ends an otherwise unchanged full run at
step `N` and writes its final checkpoint, which makes restart equivalence tests
independent of configuration edits.

Configured transport accepts the same output/restart pattern:

```bash
mpiexec -np 2 ./solvers/cpu/iga_solve DATABASE.ntiga CASE_DIR \
  --system tracer_transport --output tracer.txt --output-every 10 \
  --checkpoint tracer-checkpoint --checkpoint-every 10 --stop-after-step 10

mpiexec -np 2 ./solvers/cpu/iga_solve DATABASE.ntiga CASE_DIR \
  --system tracer_transport --output tracer-resumed.txt \
  --restart tracer-checkpoint
```

Transport metadata validates node count, ordered field names, system name,
velocity-source name, completed step, physical time, and `dt`. This ensures a
snapshot-series restart resumes with the correct interpolation time rather than
replaying the velocity sequence from zero.

## CUDA use

The CUDA backend lowers the same compiled term list and supports coupled scalar
block widths from one through eight:

```bash
./solvers/cuda/iga_cuda solve \
  case.ntiga CASE_DIR tracer_transport CASE_DIR/tracer-cuda.txt
```

The configured output and `.fields` ordering match CPU. The older CUDA
`transport` command remains only for old two-file input compatibility.
Configured CUDA transport implements matching raw-state checkpoint/restart and
time-indexed output. A 14,565-node branching snapshot-series case matched CPU
to `3.74673e-12` and its uninterrupted CUDA restart reference to
`6.64217e-20` relative L2.

## Extension boundary

JSON composes registered weak-form operators; it is not a raw mathematical
expression evaluator. A new differential operator therefore needs a reviewed
C++ operator implementation and CPU/CUDA lowering. New field names, field
counts, coefficients, couplings, and supported boundary types do not.
