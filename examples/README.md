# Example catalog

The committed examples contain source inputs only. Meshes, packed databases,
partitions, checkpoints, and simulation results are generated in a separate
work directory so that `examples/` stays small and reproducible.

## Runnable cases

| Application | Case | Purpose | Equation system |
|---|---|---|---|
| Neuron transport | `neuron_transport/straight_neurite` | Smallest two-field transport case | `neuron_transport` |
| Neuron transport | `neuron_transport/branched_neurite` | Transport through a parent and two arms | `neuron_transport` |
| Vascular flow | `vascular_flow/straight_tube` | Smallest rigid 3D flow case | `blood_flow` |
| Vascular flow | `vascular_flow/bent_tube` | Curved rigid vessel | `blood_flow` |
| Vascular flow | `vascular_flow/y_bifurcation` | One inlet and two independent outlets | `blood_flow` |
| Vascular flow | `vascular_flow/iga_wordmark` | Connected IGA-letter showcase pipe | `blood_flow` |

[`validation/womersley`](validation/womersley/) contains analytical validation
inputs and is not a standalone geometry case.

## Common input contract

Every runnable case directory contains exactly:

- `skeleton_initial.swc`: rooted centerline coordinates, radius, and parent;
- `mesh_parameter.txt`: five control-mesh generation parameters;
- `simulation_config.json`: schema-v2 fields, equations, time, and boundaries.

The preparation pipeline generates all downstream files, including
`controlmesh.vtk`, `initial_velocityfield.txt`, `bzmeshinfo.txt`, the spline
cache, METIS partitions, and the `.ntiga` database. Do not commit these
generated artifacts.

## Prepare a case

Run from the repository root and provide an empty work directory when a stable
path is useful:

```bash
EXAMPLE_WORK="$(mktemp -d /tmp/tubularflowiga-example.XXXXXX)"
RANKS=2 ./scripts/prepare_example.sh \
  vascular_flow/straight_tube "$EXAMPLE_WORK"
```

The first argument must be the application-qualified case name. The script
checks required dependencies, builds preprocessing tools, prepares and packs
the geometry, validates the configuration and boundary labels, and prints only
the solver commands appropriate for that application.

The partition count controls the generated database filename and must equal the
MPI process count used by the CPU solver. A database packed for multiple CPU
ranks can also be read by the single-GPU CUDA backend.

## Run and inspect

Follow the CPU or CUDA command printed at the end of preparation. Before a CPU
simulation, build the PETSc executables with `make cpu-petsc`; before a CUDA
simulation, build with `make cuda` and confirm `iga_cuda device-info` succeeds.

A successful preparation reports:

- the control-mesh Jacobian metrics and successful packed-element validation;
- node, element, rank, and boundary-face counts;
- resolved boundary labels and field conditions;
- the exact database and case paths.

After a vascular run, use `iga_flow_validate` to check inlet/outlet mass balance
and the divergence theorem. After a transport run, confirm the result row count
matches the database node count and inspect the neighboring `.fields` file.

## Render with ParaView

ParaView's `pvbatch` can attach CPU coefficient outputs to `controlmesh.vtk`,
take a center-plane slice, add a color bar, and save a PNG without opening the
GUI. For vascular flow:

```bash
pvbatch --force-offscreen-rendering scripts/render_example.py \
  --kind flow \
  --mesh "$VASCULAR_WORK/controlmesh.vtk" \
  --result "$VASCULAR_WORK/velocity-cpu.txt" \
  --output "$VASCULAR_WORK/vascular-flow.png" \
  --title "Vascular flow" \
  --legend "Velocity magnitude"
```

For neuron transport:

```bash
pvbatch --force-offscreen-rendering scripts/render_example.py \
  --kind transport \
  --mesh "$NEURON_WORK/controlmesh.vtk" \
  --result "$NEURON_WORK/neuron-cpu.txt" \
  --array Nplus \
  --output "$NEURON_WORK/neuron-transport.png" \
  --title "Neuron transport" \
  --legend "Nplus coefficient"
```

The render is a piecewise-linear display of IGA coefficients on the control
mesh. Use solver validation utilities and IGA quadrature for quantitative
claims. The committed gallery images were produced with ParaView 5.13.1.

The README hero uses a specified `z=0` slice because its connector arcs pass
behind the visible lettering:

```bash
pvbatch --force-offscreen-rendering scripts/render_example.py \
  --kind flow \
  --mesh "$IGA_WORK/controlmesh.vtk" \
  --result "$IGA_WORK/velocity-cpu.txt" \
  --output "$IGA_WORK/iga-wordmark-flow.png" \
  --slice-z 0 \
  --title "TubularFlowIGA" \
  --legend "Velocity magnitude"
```

## Keep or remove generated data

The work directory is fully reproducible from the three committed inputs. Keep
it while inspecting results; remove the entire work directory when those
results are no longer needed. Never point a cleanup command at the repository
root or at a directory containing unrelated case data.

Application-specific instructions:

- [Neuron transport](neuron_transport/README.md)
- [Vascular flow](vascular_flow/README.md)
- [Womersley validation](validation/womersley/README.md)
- [Latest public-example validation](VALIDATION.md)
