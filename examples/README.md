# Example catalog

The committed examples contain source inputs only. Meshes, packed databases,
partitions, checkpoints, and simulation results are generated in a separate
work directory so that `examples/` stays small and reproducible.

## Runnable cases

| Application | Case | Purpose | Equation system |
|---|---|---|---|
| Neuron transport | `neuron_transport/straight_neurite` | Smallest two-field transport case | `neuron_transport` |
| Neuron transport | `neuron_transport/branched_neurite` | Transport through a parent and two arms | `neuron_transport` |
| Neuron transport | `neuron_transport/nmo_06840_bifurcation` | Large real-morphology transport regression | `neuron_transport` |
| Vascular flow | `vascular_flow/straight_tube` | Smallest rigid 3D flow case | `blood_flow` |
| Vascular flow | `vascular_flow/bent_tube` | Curved rigid vessel | `blood_flow` |
| Vascular flow | `vascular_flow/y_bifurcation` | One inlet and two independent outlets | `blood_flow` |
| Vascular flow | `vascular_flow/iga_wordmark` | Connected IGA-letter showcase pipe | `blood_flow` |
| Vascular flow | `vascular_flow/multispecies_pulse` | Pulsatile 3D Navier--Stokes and six-species physiology | `blood_flow` + `multispecies_physiology_3d` |
| Vascular flow | `vascular_flow/vca_bifurcation` | CPU 3D VCA closed-loop bifurcation and restart smoke case | `blood_flow` + one VCA transport system |
| Native 1D | `one_d/rigid_straight` | Analytic rigid Poiseuille pressure drop | `blood_flow_1d` |
| Native 1D | `one_d/compliant_bifurcation` | Pulsatile compliant A/Q flow, junction loss, and RCR | `blood_flow_1d` |
| Native 1D | `one_d/multispecies_physiology` | Pulsatile six-species transport, physiology, and vasodilation | `blood_flow_1d` + `multispecies_physiology_1d` |
| Native 1D | `one_d/liver_vein_obj_segment` | Radius-annotated OBJ skeleton input | `blood_flow_1d` |
| Vascular flow | `vascular_flow/liver_vein_obj_segment` | 3D mesh and flow from the same OBJ excerpt | `blood_flow` |

[`validation/womersley`](validation/womersley/) contains analytical validation
inputs and is not a standalone geometry case.

## Common input contract

Every 3D runnable case directory requires:

- the `.swc` or `.obj` file named by `geometry.file`: centerline coordinates and radius;
- `simulation_config.json`: schema-v4 geometry, adaptive mesh controls, quality
  gates, fields, equations, time, and boundaries.

A case may also contain a `README.md` with provenance, resource requirements,
case-specific commands, and validation evidence. It is documentation rather
than a generated solver input.

The preparation pipeline generates all downstream files, including
`skeleton_normalized.swc`, ParaView-ready `skeleton.vtp`, geometry diagnostics,
`controlmesh.vtk`, `mesh_quality.json`,
`initial_velocityfield.txt`, `bzmeshinfo.txt`, the spline cache, METIS
partitions, and the `.ntiga` database. Do not commit these generated artifacts.

Native 1D directories instead contain only `skeleton_initial.swc` or
`skeleton_initial.obj` and a schema-v3 `simulation_config.json`. Run them
directly with `iga_1d`; they do not use the 3D mesh block, the preparation
pipeline, or a packed database. See the [1D example guide](one_d/README.md).

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

The work directory is fully reproducible from the committed geometry and
configuration. Keep it while inspecting results; remove the entire work
directory when those
results are no longer needed. Never point a cleanup command at the repository
root or at a directory containing unrelated case data.

Application-specific instructions:

- [Neuron transport](neuron_transport/README.md)
- [Large NMO_06840 neuron transport](neuron_transport/nmo_06840_bifurcation/README.md)
- [Vascular flow](vascular_flow/README.md)
- [3D VCA bifurcation validation](vascular_flow/vca_bifurcation/README.md)
- [Native 1D flow and transport](one_d/README.md)
- [Womersley validation](validation/womersley/README.md)
- [Latest public-example validation](VALIDATION.md)
