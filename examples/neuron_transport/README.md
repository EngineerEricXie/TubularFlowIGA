# Neuron-transport examples

These cases solve the legacy two-field axonal material-transport model through
the generic scalar PDE backend. They do not solve membrane voltage, action
potentials, synapses, or network electrophysiology.

## Cases

| Case | Geometry | Boundary labels |
|---|---|---|
| `straight_neurite` | Short constant-radius neurite | wall 0, inlet 1, outlet 2 |
| `branched_neurite` | One parent splitting into two arms | wall 0, inlet 1, outlets 2 and 3 |

Both configurations solve `N0` and `Nplus`. `N0` diffuses, `Nplus` advects with
the generated `initial_velocityfield.txt`, and linear coupling transfers mass
between the fields. Walls are no-flux, label 1 prescribes both inlet values,
and terminal labels use advective outflow.

## Prepare and run on CPU

From the repository root:

```bash
NEURON_WORK="$(mktemp -d /tmp/tubularflowiga-neuron.XXXXXX)"
RANKS=2 ./scripts/prepare_example.sh \
  neuron_transport/straight_neurite "$NEURON_WORK"
NEURON_DB="$NEURON_WORK/straight_neurite-2.ntiga"

mpiexec -np 2 ./solvers/cpu/iga_solve \
  "$NEURON_DB" "$NEURON_WORK" \
  --system neuron_transport \
  --output "$NEURON_WORK/neuron-cpu.txt"
```

Build the MPI/PETSc solver first if `solvers/cpu/iga_solve` is absent:

```bash
make cpu-petsc PETSC_DIR="$PETSC_DIR" PETSC_ARCH="${PETSC_ARCH:-}"
```

The result has one row per database node: `node_id N0 Nplus`. The neighboring
`neuron-cpu.txt.fields` file records the field order.

## Run on CUDA

```bash
./solvers/cuda/iga_cuda solve \
  "$NEURON_DB" "$NEURON_WORK" \
  --system neuron_transport \
  --output "$NEURON_WORK/neuron-cuda.txt"
```

Build with `make cuda CUDA_ARCHS=YOUR_GPU_ARCH` and verify
`./solvers/cuda/iga_cuda device-info` before running. CPU and CUDA use the same
case configuration and `.ntiga` database.

## Customize

- Change geometry and radius in `skeleton_initial.swc`.
- Change discretization controls in `mesh_parameter.txt`.
- Change diffusivity, coupling coefficients, time step, inlet values, or number
  of steps in `simulation_config.json`.

Field names and coefficients are explicit configuration values; the solver
does not silently apply additional neuron biology. See the
[PDE configuration guide](../../docs/PDE_CONFIGURATION.md) before adding terms.
