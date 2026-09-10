# Visualization output

CPU and CUDA production solvers share the default visualization formats.
The CPU flow CLI additionally supports the partitioned mode below. Checkpoint
outputs are independent of this selection.

## Moving immersed snapshots

Moving immersed flow snapshots are published as recoverable `<stem>.epoch…`
directories. Each complete epoch contains `fields.vtu`, a typed unstructured
grid of independent `VTK_VERTEX` cells at the volume-quadrature support points,
and `metrics.json`; `<stem>.pvd` is rebuilt from complete epochs. These are
point clouds supporting volume-field inspection, not boundary-conforming volume
meshes. A completed epoch that survives an interruption before the PVD update
is incorporated by the next rebuild. Static geometry continues to use the
Bézier VTKHDF path described below.

The standalone `phase8_compliant_channel_fsi_paraview` exporter requires one
MPI rank and rejects larger communicators before creating or modifying output.
It runs the single-process compliant-channel FSI fixture; distributed FSI is
tracked in [HPC-07](WORKSTATION_HPC_TODO.md#hpc-07分散式-fsi).

## Default format

`--visualization-format auto` is the default:

| Simulation | Visualization result |
| --- | --- |
| Transient flow or transport | Temporal `<output-stem>.vtkhdf` |
| Steady flow | Final `<output-stem>.vtu` |

Use `--visualization-format vtkhdf` or `--visualization-format vtu` to override
the automatic choice. `--output-every N` controls transient snapshots. The
initialized state is timestep zero.

## Partitioned CPU flow output

`iga_navier_stokes ... --output results/flow.txt --visualization-format pvtu`
writes cubic Bézier cells to `flow.stepNNNNNN/rankR.vtu`, one piece per MPI
rank, with `snapshot.pvtu` and a top-level `flow.pvd` time index. Open the PVD
in ParaView. Shared points carry partition-independent Int64 global IDs;
empty partitions produce valid empty pieces. Geometry is repeated each frame.

This mode extracts only the control rows required by each rank's owned elements
and exchanges shared point values without gathering the complete solution to
root. Initialization still builds and validates the complete Bézier geometry
on root, including overlap checks, then releases that mesh. This remaining
geometry peak and actual large-case RSS have not yet been eliminated or measured.

`--output-every N` includes the initialized/restarted state and every Nth step;
the final state is always emitted once. Without that option, only the final
state is emitted. Checkpoints retain their existing format and frequency.
PVTU mode produces no text velocity/pressure fields or velocity-series CSV;
use the existing formats when those transport/pipeline interfaces are needed.
This option currently applies only to `iga_navier_stokes`.

Use a fresh output prefix, including after `--restart`: an existing PVD or
snapshot directory is rejected. Each completed snapshot is published before
an atomic replacement of the PVD index. A failed later snapshot leaves the
previous index readable; incomplete pieces remain for diagnosis. A completed
snapshot whose index update fails is not automatically recovered. This protocol
assumes one producer and does not provide crash durability or MPI process-loss
recovery.

## Temporal Bézier VTKHDF

The VTKHDF file contains cubic `VTK_BEZIER_HEXAHEDRON` cells. Points,
connectivity, cell types, element IDs, partition owners, and higher-order
degrees are stored once. Each timestep appends only compressed point arrays and
step metadata, so geometry is not repeated across time.

The point arrays are obtained by applying each extraction column to the control
point solution. A global extraction-signature registry gives shared Bézier
points one visualization point ID. Coordinate equality alone never merges
points with different solution-space identities. `controlmesh.vtk` is not
rewritten and remains the canonical preprocessing/pipeline interface.

CPU and CUDA flow/transport CLIs explicitly close VTKHDF before printing their
final success summaries. CPU close errors are shared across ranks. Configured CPU
transport also coordinates root visualization initialization failures, including
nonstandard exceptions, before peers continue. When embedding
`TemporalVtkHdfWriter`, call `Close()` before reporting success: it checks flush,
root-group close, remaining local file objects, and file close. Repeated successful
`Close()` calls are harmless; `Append()` is rejected once closing has started.
Destructors provide best-effort cleanup during unwinding. These checks do not
provide atomic publication or crash durability. See the
[finalization validation report](progress/HPC_01C_IO_FINALIZATION_PROGRESS.md).

## Geometry report

Before VTKHDF creation, `<output-stem>.bezier_geometry.json` records:

- element-local and unique point counts;
- shared extraction-signature references;
- coordinate coincidences with different signatures;
- small coordinate repairs caused by legacy `%.6g` cache quantization;
- collapsed element points and sampled non-positive Jacobians;
- bounding-box overlap candidates and volume overlaps confirmed by inverse
  mapping.

Jacobian and overlap checks use the final merged visualization coordinates.
Invalid geometry prevents VTKHDF creation. Candidate counts are retained
separately so bounding-box contact is not reported as confirmed volume overlap.

## Restart behavior

Checkpoint restart resumes an existing VTKHDF only when its static geometry
hash and point-array schema match the current database and configuration.
Interrupted uncommitted rows are truncated to the recorded step count. Writing
the same final physical time replaces that timestep rather than duplicating it.

The 2026-09-10 cubic Bezier ordering correction fixes the interior point order
on two faces. Existing VTKHDF files written with the old order have a different
geometry hash and cannot be resumed by the corrected writer; regenerate the
visualization in a new output file. Solver checkpoint state is unchanged.
The ParaView regression now checks interpolation inside the cell, as well as
stored arrays. See [partitioned output validation](progress/HPC_06B_PARTITIONED_VTK_PROGRESS.md).

## Validation

The HDF5 schema regression is part of `make cpu-test`. With ParaView `pvpython`
installed, run the actual reader test:

```bash
make -C solvers/cpu vtkhdf-paraview-test
```

On Linux with a linker supporting GNU `--wrap`, exercise returned HDF5 close
errors and retry behavior with `make -C solvers/cpu vtkhdf-close-failure-test`.
