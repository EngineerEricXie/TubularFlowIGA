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
`vector3` velocity followed by pressure, and `viscosity` is required in the
configuration; it is no longer a solver constant.

The current Navier–Stokes backend is steady. The top-level `time` object drives
configured transport only; it does not make flow transient. Velocity profiles
are spatial data multiplied by a constant `scale`.

Schema version 2 optionally accepts named `temporal_functions` of kind
`constant`, `sinusoid`, `periodic_table`, or `fourier`. Each declares
`units`; periodic functions declare a positive `period`. Tables name a
strict `time,value` CSV with samples in `[0, period)` and linear,
periodically wrapped interpolation. Fourier `cosine` and `sine` arrays use
one-based harmonic order. A Dirichlet condition may name one with `waveform`.
Parsing and dependency-free evaluation are implemented, including negative-time
and period wrapping, but CPU/CUDA boundary resolvers still reject waveform
execution until transient solver integration is complete. Pulsatile flow and
previous-step velocity state therefore remain unavailable. See the
[development roadmap](ROADMAP.md).

## Boundary conditions

Every mesh label is assigned a user-selected name and a list of field-specific
conditions. Supported schema values are `dirichlet`, `no_flux`, `flux`, `robin`,
and `advective_outflow`. The configured transport assembler executes scalar
Dirichlet conditions and treats `no_flux` and `advective_outflow` as natural
conditions. CPU also assembles `flux` and `robin` surface terms from packed
element-face labels. CUDA currently rejects those two conditions pending parity;
it never silently treats them as executed numerical terms.
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
velocity/pressure boundary values from the same file. For transition only, it
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
