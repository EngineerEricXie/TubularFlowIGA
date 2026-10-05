# VascularFM mesher, flow, and transport patches

These four commits import the existing VascularFM local patches onto
TubularFlowIGA `43b93de8eebd5abfa492176fe46ffc1229a490e7`. They retain the
existing file-format interfaces and defaults except for the documented
residual-only convergence assembly. They do not modify the tetrahedral route.

## Verification status

The tests and workstation benchmark below are historical evidence recorded
with the original patches on 2026-09-29. During upstream PR preparation on
2026-10-05, the patches were applied in order, Git whitespace checks passed,
and all 22 changed source files were compared byte-for-byte with the existing
VascularFM source copy. No build, solver, model training, or held-out evaluation
was launched for this transfer. Existing recorded tests were not rerun.

The benchmark used 8 local CPU MPI ranks. The original notes do not record
CPU model or PETSc/MPI versions; this transfer supplies no new timing or
hardware validation. A small native rebuild and targeted test run remain
review gates before merging this draft PR.

The 835 unused vertices previously diagnosed in a tetrahedral fixture are
outside these patches and remain unfixed. Full-tree simulation and scaling
are not validated by this transfer. BDF2 support is limited to the standalone
`iga_solve` linear-transport path; coupled transport explicitly rejects it.

The original patch files and SHA bindings remain in the VascularFM repository
under `docs/provenance/tubularflowiga_vendored_local_patches/`. The local
VascularFM revision marker is intentionally not installed into this upstream
Git checkout.

## 0001 Mesher: port dataset fixes; diameter-scaled port extensions

- Ports the VascularFM dataset mesher fixes to the current mesher: synthetic
  tangent port extensions (`--synthetic-port-extensions`), the spline
  `--preserve-port-rims` option, persisted B-spline layer tangents, nonuniform
  tangent/curvature estimates, and diagnostic mesh export.
- Junction clearance is the larger of the upstream pairwise arc bound and the
  local chord criterion `|x_layer - x_junction| >= factor * max(D_junction,
  D_layer)`, which arc length alone does not guarantee on curved or widening
  sections.
- Synthetic extension length is `max(3 D_max, 2 min(target_spacing,
  max_spacing_over_diameter D_max))`. Previously `4 target_spacing` lengthened
  the synthetic domain when axial spacing was coarsened, so raising
  `target_spacing` above about 1 mm increased the element count.
- Tests: `mesh_core_test` (upstream and ported cases).

## 0002 Mesher: template radius compensation

- Optional `--compensate-template-radius` for `tubular_mesh pipeline`. Cubic
  spline walls through N control points on the unit circle pass their knots at
  radius `(4 + 2 cos(2 pi / N)) / 6`; all cross-section, merge and branch
  template points are scaled by the inverse so the wall matches the prescribed
  radius. Default off keeps existing geometry.
- Measured lumen area error at a port: circle `target_size` 0.5 (N = 16)
  -5.0 % without, below 0.1 % with; `target_size` 1.0 (N = 8) -18.7 % without,
  -0.1 % with.
- Tests: exact factor for N = 8 and uniform scaling of generated templates.

## 0003 Flow: residual-only convergence checks and opt-in Jacobian reuse

- The body-fitted Newton loop now assembles the Jacobian only when a linear
  solve follows. The convergence check after an update assembles the residual
  alone (`NavierStokesAssemblyRequest::ResidualOnly`, already used by the
  immersed runtime), so the accepted residual is identical and the unused
  tangent is skipped. No option; results are unchanged.
- `--jacobian-reuse-steps N` (transient only, default 0): modified Newton that
  reuses the factorized Jacobian for up to N accepted steps while every reused
  iterate at least halves the residual; otherwise it refreshes within the same
  step. Convergence is always tested on the exact residual. Reuse is dropped on
  rollback, abort, failed assembly and checkpoint restore; a reused matrix is
  never flushed during failure recovery. A reused-Jacobian update that raises the
  residual is undone (state restored from a backup vector) and the iteration
  continues as exact Newton from the accepted iterate; this fixed a start-up
  divergence seen on 3 of 276 portal cases in the first transient step.
- Tests: `navier_stokes_test`, `body_fitted_accepted_checkpoint_test`,
  `collective_failure_test`, `body_fitted_checkpoint_bundle_test`
  (codec/save/resume).
- Benchmark: see the table below.

## 0004 Transport: BDF2 for linear_transport in iga_solve

- `equation_systems[].time_integration` is accepted for `linear_transport`:
  `backward_euler` (default) or `bdf2`. BDF2 takes one backward-Euler step, then
  solves `(3/2 M + dt K) c^(n+1) = M (2 c^n - c^(n-1)/2) + dt f`; SUPG and all
  other terms are unchanged.
- Checkpoints: backward Euler still writes schema 1 byte-identically. BDF2
  writes schema 2 with `time_integration` and `history_file` (state one step
  back), so a restart continues with the same two-level history.
- The coupled `TransientTransportRuntime` rejects `bdf2` explicitly until it
  carries a history vector through its trial/commit transactions. Programmatic
  definitions that keep the generic `steady` default compile to
  `backward_euler`, so existing coupled callers are unaffected (verified with
  `body_fitted_accepted_checkpoint_test` and `collective_failure_test`).
- Tests: `generic_transport_test` (leading-mass scaling, configuration
  accept/reject), `transport_checkpoint_test` (schema 1/2 round trip and
  rejections). End-to-end `iga_solve` decay `c' = -c` on a portal mesh:
  observed order 0.97/0.99 (backward Euler) and 2.07/2.03 (BDF2) for
  dt = 1/10, 1/20, 1/40.
- Not covered: transport budget tooling that assumes a backward-Euler time
  derivative should be reviewed before using it with BDF2 output.

## Benchmark (0003)

Portal subtree, 2,166 control points (circle template 1.0), backward Euler,
dt = 0.1 s with a pulsatile inflow, 10 steps, 8 MPI ranks on one workstation,
FGMRES + ASM/ILU(1), no other load:

| Variant | Solver time | Per step | Newton iterates | Fresh / reused Jacobians | Final state vs upstream |
|---|---|---|---|---|---|
| upstream 43b93de | 333.4 s | 33.3 s | 31 | 31 / 0 (plus unused check tangents) | reference |
| residual-only checks | 252.0 s | 25.2 s | 31 | 21 / 0 | identical (0 difference) |
| `--jacobian-reuse-steps 10` | 72.7 s | 7.3 s | 35 | 3 / 22 | velocity 4.5e-5, pressure 5.0e-6 relative L2 |
