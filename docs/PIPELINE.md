# Pipeline and Case Preparation

## 1. Generate the control mesh

A case directory starts with exactly one `skeleton_initial.swc` or
radius-annotated `skeleton_initial.obj`, plus `mesh_parameter.txt`. The
parameter file contains, in order:

1. noise-smoothing iterations;
2. bifurcation-node smoothing ratio;
3. noise-smoothing ratio;
4. target Bezier segment length;
5. bifurcation refinement ratio.

Build and test the standalone C++ generator from the repository root:

```bash
make mesh
make mesh-test
./preprocessing/mesh/tubular_mesh pipeline \
  "$CASE_DIR" meshgeneration/template
```

After strict parsing and topology validation, it writes
`skeleton_normalized.swc` and `skeleton.vtp`. OBJ is thereby converted to an
explicitly rooted SWC before smoothing continues. Valid SWC follows the same
normalization path. The original input is never overwritten.

The pipeline then writes `skeleton_smooth.swc`, `controlmesh.vtk`, and
`initial_velocityfield.txt`. The implementation evaluates its own cubic
B-splines, constructs tube and bifurcation hexahedra, and has no MATLAB, TREES,
Eigen, or VTK-library dependency. A 3D run requires a tree whose nodes have at
most two children after root orientation. It stops before producing a mesh and
reports the offending node ID and child count if this constraint is violated.

The smoothed SWC is intentionally written to eight decimal places and read
back before meshing. This preserves the legacy file-interface behavior at
layer-count boundaries.

MATLAB remains as an optional reference workflow. Install TREES separately,
add both TREES and this repository recursively to the MATLAB path, set
`io_path` near the top of `TreeSmooth.m` and `Hexmesh_main.m`, then run:

```matlab
cd('/path/to/TubularFlowIGA/meshgeneration')
addpath(genpath('/path/to/TREES'))
TreeSmooth
Hexmesh_main
```

Both implementations write the same file interfaces. The C++ generator rejects
invalid lengths and diameters, limits unsafe point motion, checks element
corners and 4x4x4 quadrature samples, and
stops before writing a mesh below its scaled-Jacobian floor. Inspect the result
in ParaView as a second gate. A positive determinant is mandatory; a minimum
scaled Jacobian above 0.1 is a practical production target.

## 2. Extract the IGA representation

Build with `make spline`. The argument is a directory prefix, so include its
trailing slash:

```bash
OMP_NUM_THREADS=8 ./preprocessing/spline/spline \
  "$CASE_DIR/" --no-legacy-text
```

Set the thread count to the allocated CPU cores. The extractor processes
elements in bounded-memory chunks and formats element records in parallel while
preserving deterministic legacy text output.

It reads `controlmesh.vtk` and writes:

- `bzmeshinfo.txt`: Bezier element connectivity;
- `spline_cache.igacache`: versioned sparse coefficients and Bezier points;
- `bzmesh.vtk`: visualization output.

### Coordinate normalization

Before extraction, the spline code subtracts the minimum coordinate on each
axis and divides every coordinate by the smallest domain-axis extent. The
Bezier mesh and packed `.ntiga` geometry therefore use translated, normalized
coordinates rather than the original SWC coordinate units. Record the scale
for every physical study and transform flow, time, material, pressure, and
transport parameters consistently. The public examples use internally
consistent numerical values and are not presented as patient-specific SI
calibrations.

Omit `--no-legacy-text` to additionally reproduce `cmat.txt` and `bzpt.txt`.
The binary cache records a control-mesh content hash; `iga_pack` rejects stale,
truncated, malformed, or version-incompatible caches.

## 3. Partition and pack

CPU rank count, METIS partition suffix, and `mpiexec -np` must agree.

```bash
RANKS=8
mpmetis "$CASE_DIR/bzmeshinfo.txt" "$RANKS"
./solvers/cpu/iga_pack "$CASE_DIR" "$RANKS" "$DATABASE"
./solvers/cpu/iga_inspect "$DATABASE"
```

The packer prefers `spline_cache.igacache` and falls back to legacy text when
the cache is absent. Pass `--legacy-text` after the output path to force that
fallback for regression testing. It creates a binary database with direct
element offsets, sparse extraction rows, adjacency, per-rank touching-element
indices, and six boundary-face labels per element. Version 4 writers derive
face labels from control-mesh topology and point labels; readers remain
compatible with version 3 databases, whose face labels are unavailable. Repack
when changing the CPU rank count. CUDA ignores ownership records and may reuse
any valid packed database.

## 4. Validate and solve

Use a workstation only for cases that fit its resources. On a shared cluster,
request a compute allocation through its scheduler. Validate with
`iga_mesh_check` or `iga_cuda mesh-check` before timing either physics
stage. CPU runs use `mpiexec`; CUDA execution is one process on one GPU.

Platform-specific examples:

- [PSC Bridges-2](BRIDGES2.md)

Record the case snapshot, compiler, MPI/PETSc or CUDA versions, rank count or
GPU, scheduler job ID when applicable, solver tolerances, numerical norms,
stage timers, host RSS, and CUDA peak allocation. Compare output fields, not
runtime alone.
