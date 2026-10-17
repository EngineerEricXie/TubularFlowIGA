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
1 to 1024 outlet labels. Omit `boundaries.outlet_label` and `terminal_rcr`, and
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
By default the CLI initializes hydraulic edge pressure guesses to zero on
every step, including after restart; with `coupling.predictor` set to
`implicit_0d` (below) the guesses come from a monolithic solve of the
committed states instead. Both policies are stateless, so restarts stay
reproducible. It rejects edge/3D flow imbalance or 0D storage imbalance before
writing an epoch. With `--output-dir`, accepted steps append to one `flow.vtkhdf` by default.
`--visualization-format pvtu` selects immutable `step_N/snapshot.pvtu`
with rank-owned quadratic tetra VTU pieces. See [Native Tet VTKHDF](NATIVE_TET_VTKHDF.md)
for compression, MPI memory and strict restart checks. Point data
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

## Implicit 0D predictor (`coupling.predictor`)

The default coupling sends one pressure per edge to the native solver, reads
the flow back, and relaxes the pressure (`fixed`, `aitken`). That outer loop
is a Dirichlet–Neumann iteration whose gain is roughly the terminal impedance
divided by the 3D resistance. With liver-scale terminal resistances
(`Rp+Rd ~ 1e9 Pa·s/m³`) and a 3D tree that drops a few pascals at
`3e-6 m³/s` (`~1e6 Pa·s/m³`), each relaxed iteration multiplies the
difference between outlet pressures by hundreds. The second iteration then
drives a non-physical high-Reynolds flow and the native Newton solve fails.
That is how the patient portal case below failed after 60 minutes.

Opt in with

```json
"coupling": {"method": "fixed", "maximum_iterations": 4,
  "pressure_relative_tolerance": 1e-4, "flow_relative_tolerance": 1e-6,
  "relaxation_factor": 0.5, "predictor": "implicit_0d"}
```

At the start of each step the CLI probes every 0D neighbor through the same
backward-Euler trial evaluation the executor uses. Over one step a terminal
RCR is exactly affine, `p = R_eff·Q + b`, and the source is affine in the
inlet pressure. The native system is then solved once with these relations
as bordered unknowns: one resistance control per outlet and a
pressure-dependent flow control at the inlet. The predicted pressures seed
the unchanged executor. Its first iteration re-solves the native trial with
those pressures, which needs no linear solve when the prediction is
consistent. Acceptance still requires the executor's pressure and flow
residuals and every existing balance gate. A wrong prediction can only cost
iterations; it cannot change the accepted state. The predictor needs a 0D
source and no species. Predictor runs also use the scaled Newton
line-search merit described under solver options.

`"fluid": {"steady": true}` drops the 3D momentum time derivative; it
requires `motion.kind` `fixed` and the same 0D source restriction. The 0D
models still step with `dt_s`; choose small capacitances when they should
also be quasi-steady.

## Controlled restart with boundary updates

A restart normally rejects any change to the case. A case may instead
declare boundary values that are allowed to change between a checkpoint and
its restart:

```json
"boundary_updates": {"mutable": ["source.prescribed_flow_m3_s",
  "terminal_rcrs[].proximal_resistance_pa_s_m3",
  "terminal_rcrs[].distal_resistance_pa_s_m3",
  "terminal_rcrs[].distal_pressure_pa"]}
```

Only these four paths (and the schema-1 `terminal_rcr.` forms) are accepted.
With the declaration:

- The checkpoint binds a canonical copy of the case with exactly those
  values masked, plus the mesh, the 0D structural identities (role,
  capacitance and every undeclared value), the native model, the executor
  controls, the PETSc options, the binary and the rank count. The
  declaration itself, labels, terminal order, fluid, numerics and time
  horizon therefore stay identity-checked.
- A new checkpoint shard stores the exact 0D models of the latest accepted
  step and the boundary-value history: contiguous step ranges, each with the
  hash of every 0D model. Each 0D state is validated against the model that
  produced it.
- If the restored values differ from the case, the restart is rejected
  unless `--accept-boundary-update` is given. With the flag, the restored 0D
  states keep their stored pressures and the next step uses the new values.
  Capacitances and roles can never change.
- Epoch ids include the boundary-value hash, so branching several decisions
  from one checkpoint into separate `--checkpoint-dir` roots stays
  unambiguous.

`--summary-json` records the values, their hash, which domains changed at
the restart, and the checkpoint history. Typical use, once per decision with
a large `time.steps` horizon:

```bash
mpiexec -np 8 native_tet_hydraulic_graph decision_k.json \
  --restart-dir ck_{k-1} --accept-boundary-update \
  --checkpoint-dir ck_k --stop-after-step K --summary-json summary_k.json
```

Species in this CLI remain independent passive scalars with constant source,
first-order decay and linear wall exchange. Coupled or state-dependent
(nonlinear) reactions are not available, and the predictor, steady and
boundary-update options apply to hydraulic-only cases.

## Summary JSON and solver statistics

`--summary-json FILE` is rewritten atomically after each accepted step. It
holds the case, execution and model hashes, the rank, tetrahedron and DOF
counts, and the restart record. Each step has:

- `inlet`: outward and total inflow, area-mean 3D pressure, source pressure.
- Per outlet label: outward flow, 3D face mean pressure, terminal port
  pressure, stored pressure.
- Coupling iterations, plus predictor and trial Newton and Krylov counts.
- Assembly and linear-solve seconds, wall time, peak RSS (maximum and sum
  over ranks).
- The boundary values in force.

It covers the steps run by that process. Every step also prints a numeric
`native_hydraulic_solver` line. `--summary-json` and `-native_ale_monitor`
(per-Newton-iteration log) are excluded from the execution identity.

## Solver options for large meshes

PETSc options reach this CLI through `PETSC_OPTIONS` (including
`-options_file`); they are part of the execution identity. For meshes of
roughly 0.5–1 million P2/P1 unknowns the measured configuration is MUMPS
inside GMRES with a lagged factorization
([`native_tet_hydraulic_mumps_lagged.opts`](../solvers/cpu/options/native_tet_hydraulic_mumps_lagged.opts)):

```
-ksp_type gmres -ksp_gmres_restart 100 -ksp_rtol 1e-6 -ksp_max_it 200
-native_ale_pc_reuse_max_iterations 30
-pc_type lu -pc_factor_mat_solver_type mumps
-mat_mumps_icntl_7 5 -mat_mumps_icntl_14 50 -mat_mumps_icntl_35 2 -mat_mumps_cntl_7 1e-8
```

`-native_ale_pc_reuse_max_iterations N` keeps the factorization across
Newton iterations while GMRES converges within `N` iterations. Otherwise it
rebuilds the factorization and repeats that solve. Block low-rank
compression (`icntl_35`) is safe here because the factorization only
preconditions GMRES.

`-native_ale_scaled_line_search`, which is on for predictor runs, makes the
Newton backtracking compare residuals weighted by a characteristic force and
flow. The raw 2-norm mixes momentum rows in N with flow rows in m³/s. When
the first residual is only an inlet flow constraint (3e-6 m³/s), it rejects
a good Stokes step whose convective momentum residual is 2.5e-5 N. The
convergence test and the linear system are unchanged.

Schur-complement fieldsplit preconditioning with a pressure-mass,
`selfp`/BoomerAMG or LSC Schur approximation is far slower on these bordered,
long-tube systems. It did not finish one Newton iteration in about 10
minutes on the 2 mm portal mesh. The replicated runtime also places every
pressure row on the last rank.

## Patient portal-vein functional case

This is a numerical and functional result only, with no physiological claim.
It uses TCIA CRLM-CT-1072 `Portal`
([provenance](LIVER_GEOMETRY_CANDIDATES.md)), 15 grouped terminal RCRs,
source flow `3e-6 m³/s`, density 1060, viscosity 0.0035, `dt=0.1 s`, and one
step. Plain fixed-point coupling failed (see above). With
`coupling.predictor=implicit_0d` and the lagged MUMPS options, each step is
accepted after one coupling iteration and four Newton iterations.

- **Meshes:** remeshed meshes from
  [`cap_vessel_ports.py`](NATIVE_TET_WORKFLOW.md#labelled-surfaces-from-vessel-segmentations),
  with 0.58–0.86 million unknowns.
- **Cost at 8 ranks:** 87–121 s per process and at most 3.1 GiB per rank.
- **Per decision:** a steady restart after changed resistances takes 68 s
  at 8 ranks on 2 mm.
- **Downstream mesh:** the original sliver mesh (minimum scaled Jacobian
  `4e-7`) still stalls near a residual of `1.5e-6`, even with an exact LU.

Tables and hashes are in
[the liver geometry audit](LIVER_GEOMETRY_CANDIDATES.md#portal-vein-truth-model-runs).

Validation:

```bash
make t7-native-tet-hydraulic-cli-test PETSC_DIR=/path/to/petsc \
  T7_SPECIES_MPI_RANKS=2
```

The same target also runs `scripts/integration/test_native_tet_hydraulic_control.py`.
It checks that the predictor needs one coupling iteration and reproduces a
tightly converged fixed point to `1e-8`, and that a steady run succeeds while
a steady moving mesh is rejected. It checks the summary schema against the
step lines, and a 25-outlet grid with checkpoint/restart. For the controlled
restart it checks:

- An unchanged restart matches the uninterrupted run.
- A changed restart is rejected without `--accept-boundary-update`.
- With the flag, the new proximal resistance holds exactly at every terminal.
- The history is `[1..1, 2..3]`.
- Changed capacitance, viscosity, undeclared values, policy or outlet order
  are rejected.
- The flag without a declaration, or without a restart, is rejected.

The runtime unit test checks that a resistance-controlled solve equals a
prescribed-pressure solve at the returned pressures. The 0D checkpoint test
covers restores under a previous model.

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
at one and two ranks (see [T9](NATIVE_TET_WORKFLOW.md)); an earlier 0D load failed.
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
