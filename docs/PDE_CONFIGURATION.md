# Configurable PDE systems

`simulation_config.json` is the canonical case configuration. Schema version 2
describes fields, equation systems, time integration, and named boundary
conditions in one file. Field names such as `tracer`, `oxygen`, or `Nplus` have
no built-in numerical meaning.

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
evaluates velocity or pressure waveforms at the new-step time. Velocity profiles
are spatial data multiplied by a `scale` and, when configured, the temporal
waveform value.

Schema version 2 optionally accepts named `temporal_functions` of kind
`constant`, `sinusoid`, `periodic_table`, or `fourier`. Each declares
`units`; periodic functions declare a positive `period`. Tables name a
strict `time,value` CSV with samples in `[0, period)` and linear,
periodically wrapped interpolation. Fourier `cosine` and `sine` arrays use
one-based harmonic order. A Dirichlet condition may name one with `waveform`.
Parsing and dependency-free evaluation are implemented, including negative-time
and period wrapping. CPU and CUDA configured transport materialize Dirichlet values at
each physical step; a constant waveform reproduces constant-boundary results.
CPU backward-Euler Navier–Stokes executes waveform-backed velocity and pressure
Dirichlet conditions. Steady flow still rejects them because a time-varying
boundary has no meaning without physical time integration. CUDA transient flow,
restart, and time-indexed output remain on the [development roadmap](ROADMAP.md).

## Boundary conditions

Every mesh label is assigned a user-selected name and a list of field-specific
conditions. Supported schema values are `dirichlet`, `no_flux`, `flux`, `robin`,
and `advective_outflow`. The configured transport assembler executes scalar
Dirichlet conditions and treats `no_flux` and `advective_outflow` as natural
conditions. CPU and CUDA also assemble matching `flux` and `robin` surface
terms from packed element-face labels.
A velocity Dirichlet condition may provide a three-component `value`, or set
`"profile": "initial_velocityfield.txt"` plus a numeric `scale`. The latter
preserves the case-generated vessel inlet profile, as in the HexSim workflow.

The public smoke configuration is
[`examples/smoke/input/simulation_config.json`](../examples/smoke/input/simulation_config.json).
It solves a field called `tracer`; there are no NeuronTransport names in that
equation system. [`examples/neuron_transport/simulation_config.json`](../examples/neuron_transport/simulation_config.json)
expresses the previous two-field neuron equations using the same operators.

## CPU use

After preparing the mesh database and loading the documented PETSc/MPI setup:

```bash
mpiexec -np 2 ./solvers/cpu/iga_solve \
  case.ntiga CASE_DIR tracer_transport CASE_DIR/tracer.txt
```

Omit the system name only when the case has one linear transport system. The
output contains `node_id` followed by fields in configured order; a neighboring
`.fields` file records their names.

`iga_navier_stokes DATABASE CASE_DIR` reads the `navier_stokes` system and its
velocity/pressure boundary values from the same file. `steady` performs one
nonlinear solve; `backward_euler` performs `time.steps` solves with step size
`time.dt` and writes the final state to the optional output path. For transition only, it
still accepts cases containing `simulation_parameter.txt` plus
`case_config.json`; that input is converted to the shared representation rather
than selecting a legacy numerical model.

## CUDA use

The CUDA backend lowers the same compiled term list and supports coupled scalar
block widths from one through eight:

```bash
./solvers/cuda/iga_cuda solve \
  case.ntiga CASE_DIR tracer_transport CASE_DIR/tracer-cuda.txt
```

The configured output and `.fields` ordering match CPU. The older CUDA
`transport` command remains only for old two-file input compatibility.

## Extension boundary

JSON composes registered weak-form operators; it is not a raw mathematical
expression evaluator. A new differential operator therefore needs a reviewed
C++ operator implementation and CPU/CUDA lowering. New field names, field
counts, coefficients, couplings, and supported boundary types do not.
