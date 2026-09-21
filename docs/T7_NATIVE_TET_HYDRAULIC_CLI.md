# Native tetra hydraulic graph CLI

`solvers/cpu/native_tet_hydraulic_graph` reads a labelled ASCII Gmsh 4.1
tetrahedral mesh and a versioned JSON case, then runs the project-owned P2/P1
ALE Navier–Stokes FEM between an existing 0D source reservoir and one or more
terminal RCRs (version 2 for multiple outlets). PETSc supplies algebra/MPI, not FEM basis, quadrature, assembly, or
coupling physics. A mesh produced by the fTetWild adapter can be used after its
canonical Gmsh conversion and geometry gates; fTetWild itself is only a
preprocessing mesher.

Build and run on an allocated compute resource:

```bash
make -C solvers/cpu native_tet_hydraulic_graph PETSC_DIR=/path/to/petsc
mpiexec -np 2 ./solvers/cpu/native_tet_hydraulic_graph case.json --check-input
mpiexec -np 2 ./solvers/cpu/native_tet_hydraulic_graph case.json \
  --checkpoint-dir /path/to/new/checkpoints \
  --output-dir /path/to/new/fields --stop-after-step 1
mpiexec -np 2 ./solvers/cpu/native_tet_hydraulic_graph case.json \
  --restart-dir /path/to/new/checkpoints --output-dir /path/to/new/fields
```

`case.json` follows the concrete
[`native_tet_hydraulic_star.json`](../solvers/cpu/tests/data/native_tet_hydraulic_star.json)
example. `mesh_file` is relative to the JSON file unless absolute. Every
boundary label in the mesh must be assigned exactly once as inlet, outlet, or
wall. The inlet receives the source's flow; each outlet receives its RCR's
pressure. 0D and 3D graph flows use SI units and the existing outward-positive
sign. `motion.kind` is `fixed` or `rigid_translation_x`; the latter moves the
whole mesh by `speed_x_m_s × dt_s` on every accepted step. `coupling.method` is
`explicit` or fixed-point `fixed`. The JSON parser rejects unknown fields,
unsupported method names, overlapping/missing labels, and nonpositive
model/time parameters. The read-only `--check-input` path additionally reads
the referenced native mesh, checks label equality and positive tetrahedral
Jacobians, and builds P2 topology without advancing time or writing solver
output. It cannot be combined with output, restart, or stop options.

Schema version 2 retains one source/inlet but supports distinct RCR loads on
multiple outlet labels. Omit `boundaries.outlet_label` and `terminal_rcr`, and
use `terminal_rcrs` instead; each entry has a `boundary_label` plus the same
five RCR parameters as the version-1 `terminal_rcr` object. For example:

```json
"boundaries": {"inlet_label": 1, "wall_labels": [0]},
"terminal_rcrs": [
  {"boundary_label": 2, "proximal_resistance_pa_s_m3": 1.0,
   "distal_resistance_pa_s_m3": 10.0, "capacitance_m3_pa": 0.01,
   "distal_pressure_pa": 0.0, "initial_pressure_pa": 0.0},
  {"boundary_label": 3, "proximal_resistance_pa_s_m3": 1.0,
   "distal_resistance_pa_s_m3": 20.0, "capacitance_m3_pa": 0.01,
   "distal_pressure_pa": 0.0, "initial_pressure_pa": 0.0}
]
```

Each outlet has its own graph edge, pressure state, storage balance, and
checkpoint shard. Log keys are `outlet_<label>_m3_s` and
`terminal_<label>_pa`. Version-1 files and their old log keys remain valid.
The accepted graph checks edge and aggregate 3D flow imbalance relative to
actual flow scale, with a `1e-12 m³/s` absolute floor; small-flow cases no
longer pass solely because their flows are below `1e-6 m³/s`.

For restart, the case and mesh contents, graph/model settings, executable,
PETSc version/options, and MPI rank count are bound to the checkpoint identity.
The dedicated bundle atomically publishes each accepted source–FEM–RCR epoch.
The CLI initializes hydraulic edge pressure guesses to zero on every step,
including after restart; this keeps its current stateless coupling policy
reproducible. It rejects edge/3D flow imbalance or 0D storage imbalance before
writing an epoch. With `--output-dir`, each accepted step creates an immutable
`step_N/snapshot.pvtu` with rank-owned quadratic tetra VTU pieces. Point data
contains every P2 velocity DOF and pressure interpolated from P1 vertices at
the midside nodes. It also contains `reference_position_m` and
`displacement_m` for every P2 node; current coordinates equal reference
coordinates plus displacement, including after restart. The six native edge
nodes are reordered to the
[VTK quadratic tetra convention](https://vtk.org/doc/nightly/release/8.1/html/classvtkQuadraticTetra.html).
Cell data additionally contains row-major `cauchy_stress_pa` (nine components),
evaluated at each tetrahedron centroid from the native P2 velocity gradient,
P1 pressure and case viscosity as `−pI+μ(∇u+∇uᵀ)`. This is fluid Cauchy stress
in Pa, not solid wall stress or a surface-averaged wall-shear field.
An existing `step_N` directory is rejected without overwriting its published
snapshot. This is full nodal velocity/pressure visualization, not graph-history
output.

Validation:

```bash
make t7-native-tet-hydraulic-cli-test PETSC_DIR=/path/to/petsc \
  T7_SPECIES_MPI_RANKS=2
```

The test uses an actual Gmsh 4.1 mesh, runs uninterrupted and split-process
checkpoint/restart executions, and compares second-step hydraulic outputs,
complete VTU fields, and quadratic cell topology. It rejects a changed case,
wrong boundary label, and inverted tetrahedron; the latter two are also
rejected by read-only input preflight. It also generates a 24-tetra channel
with separate inlet, outlet, and no-slip wall labels, verifies
throughflow and moving-wall P2 velocity, and compares its entire second-step
field, including reference/displacement, after cross-process restart at 1, 2,
and 4 ranks. Both steps verify the prescribed rigid translation. Removing the wall
label is rejected. With a negative source flow, the same channel reverses
both port flow signs in explicit and fixed-point coupling at all three rank
counts. These reverse-flow results cover this idealized, low-inertia parameter
set only.
The 24-tetra regression also splits its outlet into labels 2 and 4, assigns
different RCR impedances, and checks both branch flows and full-field
cross-process restart at 1, 2, and 4 ranks. Duplicate outlet labels and a
changed terminal model on restart are rejected. On this controlled mesh, both
explicit and fixed-point coupling also pass positive and reversed source flow
at 1, 2, and 4 ranks. A deliberately one-iteration fixed-point run is rejected
with per-edge pressure and flow residual diagnostics; it is not silently
accepted as an explicit step. A separately scaled, real two-outlet Y
functional case also passes fixed-point coupling on Gmsh and fTetWild meshes
at one and two ranks (see [T9](T9_NATIVE_TET_WORKFLOW.md)); an earlier 0D load failed.
Neither result establishes fixed-point robustness for arbitrary loads or Y
meshes.
An optional local end-to-end smoke uses an installed upstream fTetWild binary
to generate a labelled 0.03 m circular pipe, then solves one native source–FEM–RCR
step and checks the full P2 VTK field and flow signs. It is separate from the
dependency-free CLI regression because the fTetWild executable is not vendored:

```bash
make t7-native-tet-ftetwild-smoke PETSC_DIR=/path/to/petsc \
  FTETWILD_BIN=/path/to/FloatTetwild_bin T7_SPECIES_MPI_RANKS=2
```

The local WSL run passed at 1 and 2 MPI ranks with upstream commit `d7d99bb`.
Repeated mesher invocations yielded 607–738 tetrahedra; the test checks mesh
gates, flow signs, and complete quadratic
field structure rather than exact mesh or QoI equality. These values are
functional evidence only, not a mesh-convergence or physiological result.

These are functional smoke tests at 1, 2, and 4 ranks—not a physiological
parameter calibration, strict spatial/temporal accuracy claim, generic 0D
species path, large-mesh scaling result, or integration with the older IGA
multidomain runner. The current metadata payload is limited to 16 MiB, so
large native flow fields need a streamed/sharded checkpoint extension.
