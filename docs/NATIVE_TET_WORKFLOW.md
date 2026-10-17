# Explicit native tetra hydraulic workflow

`scripts/run_native_tet_workflow.py` is a small user-facing entry point from
labelled surface or tetrahedral volume data to the project-owned P2/P1 ALE
Navier–Stokes source→FEM→one-or-more-RCR runner. PETSc supplies sparse algebra and MPI;
fTetWild or Gmsh may supply a volume mesh but never a FEM solver. This is a
hydraulic functional route, not a generic organ model, coupled FSI route,
or physiological validation. The same entry point also accepts the separate
`native_tet_p1_prescribed_species` backend described in
[the species guide](NATIVE_TET_SPECIES_CLI.md); it does not couple species
to this hydraulic run. A third explicit
`native_tet_p1_darcy_steady` backend uses the same surface/volume meshing
and provenance path for an independent, single-compartment
[Darcy solve](validation/T6_NATIVE_TET_DARCY.md); it is not vessel–tissue coupling.

Transient backends default to one compressed VTKHDF series per field family.
The optional workflow key `visualization_format` accepts `auto`, `vtkhdf`, or
`pvtu`; select `pvtu` explicitly for legacy XML output. Keep the selection
unchanged across resume. See [Native Tet VTKHDF](NATIVE_TET_VTKHDF.md).

First build the native executable:

```bash
make -C solvers/cpu native_tet_hydraulic_graph surface_fem_preflight \
  PETSC_DIR=/path/to/petsc
```

The version-1 workflow JSON explicitly selects the backend and input route.
For a labelled VTP/STL surface, use one of these `mesher` objects:

```json
{
  "schema_version": 1,
  "backend": "native_tet_p2p1_ale_hydraulic",
  "input_route": "surface",
  "input_file": "pipe.vtp",
  "case_file": "pipe_hydraulic.json",
  "output_directory": "pipe_run_001",
  "mpi_ranks": 2,
  "mesher": {
    "kind": "ftetwild",
    "target_size_m": 0.007,
    "envelope_m": 0.0001,
    "executable": "/path/to/FloatTetwild_bin"
  }
}
```

For Gmsh, replace `mesher` with
`{"kind":"gmsh","target_size_m":0.007}`. For an already labelled ASCII
Gmsh 4.1 tetrahedral volume, set `input_route` to `"volume"`, point
`input_file` at the `.msh`, and omit `mesher`. All paths are relative to the
workflow JSON unless absolute. `case_file` follows
[`native_tet_hydraulic_star.json`](../solvers/cpu/tests/data/native_tet_hydraulic_star.json);
the workflow replaces its `mesh_file` field with the generated/copied mesh.
For a two-outlet bifurcation, use the version-2 `terminal_rcrs` case form in
[T7](NATIVE_TET_HYDRAULIC_CLI.md); each outlet label keeps an independent
native RCR state.
The case must still assign every inlet, outlet, and wall label correctly and
provide the source, outlet RCR(s), fluid, motion, coupling, and time data. Unsupported
backend and input routes fail before creating an output directory; there is no
organ-name inference or external FEM fallback.

Run on an allocated compute resource:

```bash
python3 scripts/run_native_tet_workflow.py workflow.json
```

### Labelled surfaces from vessel segmentations

A segmentation isosurface, such as the output of
`scripts/dicom_seg_to_surface.py`, is closed but has blind vessel ends and no
port labels. `scripts/cap_vessel_ports.py` turns it into a labelled closed
surface for the surface route above. It takes the surface and a ports JSON
listing each port's label, a point, an outward tangent and a radius; several
ports may share a label:

```bash
python3 scripts/cap_vessel_ports.py lumen.stl ports.json capped.vtp \
  --length-scale-to-m 0.001 --smoothing-iterations 30 \
  --target-size 1.0 --centerline-swc centerline.swc
python3 scripts/surface_to_fem_volume.py capped.vtp volume.msh \
  --length-scale-to-m 0.001 --target-size-m 0.001
```

The tool works in these stages:

1. **Smoothing (optional).** A windowed-sinc filter removes voxel
   staircasing. The per-vertex displacement is recorded and can be bounded
   with `--max-smoothing-displacement`.
2. **Cutting.** Each port removes the blind end beyond a plane, set inward
   from the port point by the first of `--cut-offsets-radii` that yields one
   simple local loop. Vertices within `--snap-fraction` of their shortest
   edge are snapped onto the plane, so the cut adds no slivers.
3. **Cap checks.** A plane through a rounded or tapering tip is rejected when
   its cap is much smaller than the port radius or than the section one
   radius deeper. If no offset passes, the valid cut with the largest cap is
   used and marked `fallback_largest_valid_cap` in the manifest.
4. **Remeshing (optional).** Gmsh remeshes the wall. The element size is
   bounded by `--target-size`, by `--radius-factor` times the distance to
   the centerline, and by `--thickness-factor` times the local wall-to-wall
   distance. Each cap is meshed as a planar surface on its remeshed rim, so
   labels never depend on feature angles.

The tool fails closed in these cases:

- a port has no valid cut;
- a remesh folds back onto itself after Gmsh's other 2D algorithms are
  tried;
- the surface is not a closed two-manifold.

The C++ surface preflight still runs inside `surface_to_fem_volume.py`.

On CRLM-CT-1072 `Portal`, the plane-cut plus fan-cap surface previously used
downstream passed only with `--minimum-scaled-jacobian 1e-9` (minimum
scaled Jacobian `4.05e-7`). This route instead passes the default `1e-3`
volume gate; see the
[liver geometry audit](LIVER_GEOMETRY_CANDIDATES.md#portal-vein-truth-model-runs).
`preprocessing/tet/skeleton_to_tet.py` remains a template-free alternative
for sparse trees. Its OCC fragment of every capsule fails or stalls on dense
real trees (2,076 primitives on the pruned portal tree). `--timeout-s` and
the stage reports bound and attribute that cost.

For a planned pause after an accepted step, use
`--stop-after-step 1`. The status becomes `partial`, and a subsequent
`--resume` with the same workflow JSON, solver binary, MPI rank count, case
template, source input, and mesher binary (if used) finishes the remaining
steps without remeshing or overwriting accepted fields:

```bash
python3 scripts/run_native_tet_workflow.py workflow.json --stop-after-step 1
python3 scripts/run_native_tet_workflow.py workflow.json --resume
```

These pause/restart options apply only to transient hydraulic and species
backends. The steady Darcy backend rejects them before creating output and
records one completed solve as `accepted_steps=1` in the shared status schema.

The same resume path accepts an interrupted `running` or `failed` simulation
only when it finds a complete checkpoint and all corresponding VTK steps. It
rejects altered input/effective case/mesh/solver identity, missing steps,
orphan visualization steps, and a concurrently running workflow. A failure
before the first completed checkpoint cannot be resumed. The underlying CLI
rechecks checkpoint checksums and model/execution identity before restoring.

The output directory must not already exist. It retains the input-derived
`volume.msh`, effective `case.json`, full P2 velocity/P1 pressure VTK fields,
P2 reference coordinates and displacement in metres, centroid fluid Cauchy
stress in Pa on each tetrahedral cell, accepted-step checkpoints,
logs, and `workflow_status.json`. Before stepping, the workflow runs the
project-owned solver's read-only `--check-input` path on the effective case and
mesh. It records `input_preflight.log` and its timing, and rejects missing or
overlapping labels and non-positive tetrahedral Jacobians before producing a
`simulation.log`; a failed preflight remains inspectable but cannot be resumed
as a simulation. Surface routes
also retain `mesh_manifest.json` with geometry and quality gates. Status
contains source, template, configuration, solver, mesh, and optional fTetWild
binary SHA-256 hashes, stage timings, and an explicit `functional_only` label.
Its `evidence_classification` separates three questions: a completed native
run gets `numerical_acceptance=functional_checks_only`, never a strict
numerical-validation claim; `physical_acceptance=not_established` for every
run; and `demonstration_visualization=available` only after an accepted step.
Partial and failed runs keep `numerical_acceptance=not_established`. These
labels do not replace the underlying mesh, solver, conservation, or QoI tests.
On failure, partial output and its `failed` status are left for inspection;
nothing is silently overwritten. `simulation_resume_N.log` files retain each
resume attempt. The underlying native CLI and its checkpoint rules are
described in [T7](NATIVE_TET_HYDRAULIC_CLI.md).

Local regression:

```bash
make t9-native-tet-workflow-test PETSC_DIR=/path/to/petsc
make t9-native-tet-workflow-test PETSC_DIR=/path/to/petsc \
  FTETWILD_BIN=/path/to/FloatTetwild_bin
make t9-native-y-rcr-smoke PETSC_DIR=/path/to/petsc \
  FTETWILD_BIN=/path/to/FloatTetwild_bin T7_SPECIES_MPI_RANKS=2
make t9-native-y-rcr-strong-smoke PETSC_DIR=/path/to/petsc \
  FTETWILD_BIN=/path/to/FloatTetwild_bin T7_SPECIES_MPI_RANKS=2
```

The first command covers the volume route, interruption/resume field parity,
changed-identity rejection, and fail-closed configuration. The
second additionally generates a real labelled pipe through fTetWild and Gmsh,
then runs the native solver and checks complete quadratic VTK topology. Local
WSL runs passed with upstream fTetWild commit `d7d99bb` and Gmsh 4.8.0.
The Y command constructs a four-label bifurcation surface, meshes it with
Gmsh and optionally fTetWild, and solves two independently loaded native RCR
outlets. Outlet label 2 uses `Rp=1000`, `Rd=10000 Pa·s/m³`; label 3 uses
`Rp=2000`, `Rd=20000 Pa·s/m³`, each with `C=1e-6 m³/Pa`. The regression checks
the two terminal pressure states are distinct, both outlet flows are positive,
and the global mass balance closes. In the local
two-rank smoke, Gmsh produced 307 tetrahedra while fTetWild counts varied;
both completed two steps with positive branch outlet flows and relative mass
closure under the native gate. fTetWild's stochastic mesh counts are not
frozen assertions. Both remain functional, low-inertia examples rather than
validated physiological bifurcation solutions.
With source capacitance `1e-11 m³/Pa`, source resistance `1e8 Pa·s/m³`, and
the above asymmetric outlet loads, the separate fixed-point Y regression passed two steps on
both meshers at one and two MPI ranks. The second step took 12 outer iterations
in each local run, with independently positive outlet flows and relative mass
closure under the native gate. The regression requires more than one outer
iteration and retains the native nonlinear tolerance `1e-8` and coupling
flow-relative tolerance `1e-6`; no conservation gate was relaxed.

An earlier, differently scaled 0D load failed in the native Newton line
search on Gmsh (`residual=0.000236`) and fTetWild (`0.000025`). Relaxing the
inner tolerance alone left a `0.817569` normalized inlet-edge flow residual
at iteration 30. These are negative results for that load, not evidence that
the current fixed-point example is generally robust. A further local probe
with label 3 set to `Rp=5e6`, `Rd=5e7 Pa·s/m³` and `C=1e-10 m³/Pa` also
failed the inner Newton solve (`residual=0.001364`); it is not a passing case.
Re-run the supported functional example with:

```bash
python3 scripts/integration/test_native_tet_hydraulic_y.py \
  --solver solvers/cpu/native_tet_hydraulic_graph --ranks 2 --coupling fixed
```

The failing high-impedance load is now a regression for the
[implicit 0D predictor](NATIVE_TET_HYDRAULIC_CLI.md#implicit-0d-predictor-couplingpredictor).
Without the predictor, label 3 (`Rp=5e6`, `Rd=5e7 Pa·s/m³`, `C=1e-10 m³/Pa`)
still stops at `residual=0.001364`. With `coupling.predictor=implicit_0d`
each step is accepted after one coupling iteration, and the test checks that
the same case without the predictor still fails:

```bash
make t9-native-y-rcr-implicit-smoke PETSC_DIR=/path/to/petsc T7_SPECIES_MPI_RANKS=2
```

The workflow passes `coupling.predictor`, `fluid.steady` and
`boundary_updates` through its case file. Restarts that change declared
boundary values use the CLI directly with `--accept-boundary-update`; the
workflow `--resume` path does not change the case.
This proves functional data routing only; spatial/temporal accuracy and
physiological QoIs remain separate validation tasks.

A separate fixed-case source→native tetra FEM→RCR species graph driver supports
P1 concentration VTU snapshots on a labelled fTetWild pipe, including
reference/current mesh and ALE displacement. See
[T7 native moving species](validation/T7_NATIVE_MOVING_SPECIES.md). The generic
prescribed-velocity species workflow backend is separate from that coupled
graph driver.
The separate [prescribed-velocity native species CLI](NATIVE_TET_SPECIES_CLI.md)
now accepts a versioned case and external labelled tetra mesh, with its own
scalar/ALE checkpoint and restart. It does not turn this hydraulic workflow
into a coupled species graph.
