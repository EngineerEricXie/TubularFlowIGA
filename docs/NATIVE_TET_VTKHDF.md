# Native tetrahedral temporal VTKHDF

The CPU `native_tet_hydraulic_graph` and `native_tet_species_transport` CLIs
default to `--visualization-format auto`, which selects VTKHDF. Set
`--visualization-format pvtu` for the previous partitioned XML snapshots.
The workflow JSON accepts the same optional `visualization_format` key.
Existing IGA visualization and solver checkpoint formats are unchanged.

| Runner | Files in the output directory |
| --- | --- |
| Hydraulic graph | `flow.vtkhdf` |
| Hydraulic graph with species, including T7 | `flow.vtkhdf`, `species_<id>.vtkhdf` |
| Prescribed-velocity species | `species.vtkhdf`, existing checkpoint/summary files |

Each file contains every accepted output time. File count is independent of
MPI rank count and timestep count. Fixed coordinates, connectivity, cell types,
and Int64 global point/cell IDs are stored once. ALE writes current coordinates
at each time while retaining connectivity. Point and cell fields use lossless
FP64 storage with chunked DEFLATE level 4. Flow retains quadratic tetrahedra,
P2 velocity, P1 pressure interpolated to midside nodes, reference coordinates,
displacement and the nine-component cell Cauchy stress. Species retains linear
tetrahedra and concentration. No averaging across cells or rank partitions is
introduced. Array names and physical units match the PVTU path.

The implementation follows the VTKHDF 2.1
[temporal unstructured-grid specification](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/vtkhdf_specifications.html#temporal-data).
Open the `.vtkhdf` directly in a compatible ParaView and use the time controls.
This format needs an HDF5 development library at build time, discovered through
the existing Makefile settings `HDF5_ALL_CFLAGS` and `HDF5_LIBS`.

## MPI and memory

The native Tet runtime already replicates mesh and state. Only rank zero builds
the complete visualization grid and writes serial HDF5; no solution gather is
added. Initialization, append and explicit close errors propagate collectively.
This removes shared-node duplication and the per-rank/per-time XML files. It
does **not** provide distributed parallel HDF5 I/O: rank zero needs memory for
one complete visualization frame. For a future nonreplicated or much larger
runtime, use explicit PVTU until distributed HDF5 output is implemented.

## Restart and failure behavior

New outputs use exclusive file creation and reject overwrites. Same-directory
hydraulic checkpoint restart appends only when topology, IDs, coordinate mode,
field schemas, exact published step count and final physical time agree. A
hydraulic restart can also use a fresh output directory, starting a new series.
That series records its initial accepted-step offset so later restarts can
append to it while still checking its frame count against the checkpoint.
Files without offset metadata retain the original zero-offset interpretation.
PVTU restart restores the existing `flow.pvd` entries before appending; its last
time must match the checkpoint. A fresh directory starts a new PVD collection.
The prescribed-species `--resume DIR` continues its single file in place.
Use the same visualization format when resuming a workflow.

Each append flushes data before updating `NSteps`. Resume checks dataset row
counts and rejects incomplete or extra rows, including an output step ahead of
the solver checkpoint. This is strict rejection, not automatic crash repair.
HDF5 flush/close checks do not promise filesystem crash durability or MPI
process-loss recovery. On failure, preserve the files for diagnosis. Native
checkpoint binaries and per-step checkpoint retention remain independent of
visualization selection. Historical PVTU runs must explicitly select `pvtu`.

## Export existing native flow states without solving again

```bash
make -C solvers/cpu native_tet_flow_export
solvers/cpu/native_tet_flow_export series.json flow.vtkhdf
```

The input is a fixed labelled Gmsh 4.1 mesh and little-endian IEEE float64
states: interleaved P2 velocity xyz in native vertex/edge order, followed by P1
vertex pressure. The manifest declares the pressure reference explicitly:

```json
{
  "mesh_file": "mesh.volume.msh",
  "dynamic_viscosity_pa_s": 0.0035,
  "pressure_offset_pa": 400.0,
  "snapshots": [
    {"time_s": 0.1, "state_file": "state_1.f64"},
    {"time_s": 0.2, "state_file": "state_2.f64"}
  ]
}
```

All input paths resolve relative to the manifest. Use `pressure_offset_pa: 0`
when the stored pressure is already physical. The exporter validates state
size, finiteness and strictly increasing times, and rejects an existing output.
It does not initialize MPI/PETSc or execute a PDE solve. It cannot infer units,
mesh identity or ordering from a raw array; retain source hashes in the calling
campaign's provenance record.

## Related native mesh and matrix fixes

Native ALE PETSc matrices now preallocate the union of all 34-DOF cell stencils
on each owned row, including off-rank element contributions and the flow-rate
controller rows/columns. Structural zeros are established before element-only
assembly so PETSc cannot discard controller capacity. New allocations are an
error. Element integration, residuals, boundary conditions and Newton controls
are unchanged. Symbolic adjacency construction uses temporary host memory.

The skeleton mesher removes unreferenced surface nodes before transferring the
surface to the volume model, preserves retained node tags, rejects empty required
surfaces and rejects volume output containing unused nodes. Its manifest records
removed surface nodes and zero unused volume nodes. Optional
`--centerline-search exact_kdtree` requires NumPy/SciPy and uses a conservative
segment bound with original-order tie breaking; `linear` remains the default.
Experimental OCC union/scaling changes and portal-specific pruning/coarsening
parameters are not adopted by this change.

## Reproducible bounded validation

Use fresh output directories and an MPI-capable compute host:

```bash
make mesh-test
make -C solvers/cpu temporal_unstructured_vtkhdf_test
solvers/cpu/temporal_unstructured_vtkhdf_test /tmp/tet-writer-fixture
make -C solvers/cpu native_tet_hydraulic_graph native_tet_species_transport \
  native_tet_ale_petsc_runtime_test PETSC_DIR=/path/to/petsc
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -np 2 \
  solvers/cpu/native_tet_ale_petsc_runtime_test
python3 scripts/integration/test_native_vtkhdf_cli.py --bin-dir solvers/cpu \
  --work /tmp/tet-vtkhdf-mpi
pvpython --no-mpi scripts/validate_native_vtkhdf_reader.py \
  /tmp/tet-vtkhdf-mpi/comparisons.json
```

The writer fixture tests mixed point/cell fields, Int64 IDs above 2^53, fixed and
moving geometry, exclusive creation, resume, invalid times/topology/schema,
append-after-close and injected HDF5 flush failure. MPI fixtures cover 1/2/8
ranks (including ranks with no owned cells), flow plus species, standalone
species restart, repeated hydraulic restarts in original and fresh directories,
PVD history preservation, stale checkpoints, workflow restart and collective output
failure. The independent ParaView reader compares every field and connectivity
against explicit PVTU output. These are bounded format and numerical regression
checks, not a mesh-convergence or physiological-validation study.
