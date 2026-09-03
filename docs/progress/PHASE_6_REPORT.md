# Phase 6 Progress

## PR6.1 — Backend-neutral immersed 3D domain kind

Status: complete.

The coupling graph now distinguishes `three_d_immersed_flow` from the existing
body-fitted 3D backend while treating both as semantically three-dimensional.
Schema-v5 accepts immersed domains without a body-fitted `.ntiga` database;
their case directory and boundary-label ports remain the narrow configuration
contract for a later runtime factory.

Immersed ports use canonical nonnegative integer `boundary_label` locators,
outward orientation `+1`, and the current pressure/flow/mean-normal-traction
interface quantities.  Total-pressure and species extensions are not part of
PR6.1.

The approved first Phase 6 integration target is quasi-static immersed 3D.
This milestone does not claim backward-Euler equivalence.

Validated by `simulation_graph_test`, `multidomain_config_test`, and
`pressure_flow_executor_test`: schema parsing/canonical serialization,
dimension-aware graph validation, exact runtime-kind matching, and
backend-neutral pressure/flow executor behavior.

## PR6.2 — Static immersed open ports, loads, and controllers

Status: complete (validated 2026-09-02).

`ImmersedFlowPort.hpp` adds catalog-bound open-patch definitions for positive
surface labels.  Port labels are unique, source-area-audited, and disjoint
from selected Nitsche wall labels.  Pressure and mean-normal-traction inputs
use the established body-fitted convention `sigma*n=-p_b*n`, with
`p_b=-t_n` for a normal-traction value; total pressure is rejected because no
approved averaging definition exists.  Measurements use the retained physical
surface rule and its outward normals: area, outward flow, mean pressure, and
mean normal traction.

The static runtime fixes port modes and scalar topology at construction, while
`SetPortControlValue(id,value)` permits a finite value refresh before a later
assembly/solve lifecycle.  A flow port owns one multiplier and contributes
`R_u += lambda c`, `R_lambda=Q-Qtarget`; the stored negative residual contains
`-lambda c` and `Qtarget-Q`, and both Jacobian off-diagonal blocks are `+c`.
Pressure-like ports remove the mean-zero gauge; no-port and all-flow cases
retain it, with scale-aware all-flow net-target compatibility checked before
assembly.  Positive cut cells always require ghost coverage, including
cap-only open-port cells; wall diagnostics count only scattered wall terms.

Required real-PETSc validation from the repository root (system PETSc 3.15,
`PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real`):

- `make -C solvers/cpu immersed_flow_port_test PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real && ./solvers/cpu/immersed_flow_port_test`
- `make -C solvers/cpu immersed_static_flow_test PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real && ./solvers/cpu/immersed_static_flow_test`
- `make -C solvers/cpu phase5_closure_test PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real && ./solvers/cpu/phase5_closure_test --manufactured-level 3`
- `./solvers/cpu/phase5_closure_test --sliver-m 8`

Observed results: `make -B immersed_flow_port_test immersed_static_flow_test
phase5_closure_test PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real`
completed successfully.  `ldd immersed_flow_port_test` resolved
`libpetsc_real.so.3.15`.  `immersed_flow_port_test` completed successfully on
two consecutive invocations; `immersed_static_flow_test`,
`phase5_closure_test --manufactured-level 3`, and
`phase5_closure_test --sliver-m 8` each completed successfully once.

The focused port test covers cap-only selection/ghost coverage, independent
cap area/load auditing, pressure and normal-traction signs, isolated
controller residual/Jacobian blocks, tiny-target rejection, affine traction
measurement, overflow rejection, scalar-row gauge cases, and
compact/expanded determinism.  No volume residual path was changed.

## PR6.3 — Transactional immersed flow-domain runtime

Status: complete (validated 2026-09-02).

`ThreeDImmersedFlowDomain` is the final `CoupledDomainRuntime` adapter for
the static immersed backend.  It reports exactly
`DomainKind::ThreeDImmersedFlow`, validates an immutable runtime port catalog
against canonical boundary-label metadata (including count, ID, parsed label,
and control mode), and uses the coupling step only to advance interface time;
the backend remains a `dt=0` steady solve.  Each configured pressure,
mean-normal-traction, or flow controller receives exactly one finite hydraulic
input at the end of an active step.  Total pressure and species remain
rejected.

Trials always restart from the backend committed vector.  Inputs, controls,
measurements, and committed time/index are staged separately.  A rollback is
valid only for a solved trial; cancelling an unsolved or prepared step uses
`AbortStep`, which restores the backend image and leaves no provisional port
result visible.
The backend now owns a preallocated prepared `Vec`: `PrepareCommit()` copies
the solved trial into that buffer, while `FinalizeCommit() noexcept` only
swaps already-owned vector handles and updates counters.  The adapter stages
validated port states before that publication and finalizes its own maps and
clock with non-allocating swaps.  The compatibility `Commit()` API is retained
as prepare followed by finalize.

Focused real-PETSc validation completed with system PETSc 3.15
(`PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real`):

- `make -B -C solvers/cpu immersed_flow_domain_test PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real`
- `./solvers/cpu/immersed_flow_domain_test`
- `./solvers/cpu/immersed_flow_domain_test`
- `git diff --check`
- `make clean`

The domain test uses a nonzero flow-controlled inlet, pressure-controlled
outlet, and wall label 0.  It checks exact vector/serialized-port replay and
deterministic Newton/KSP work after rollback; changed-input abort recovery;
prepare-failure recovery; exact adapter/backend rollback, prepare, finalize,
and abort counter deltas; stable idempotent finalization; invalid transitions;
traction-target forwarding, nontrivial traction output, and gauge suppression;
output-subset visibility; and registry-kind matching.

## PR6.4 — Production immersed-case factory integration

Status: complete (validated 2026-09-02).

`ImmersedFlowCase` is itself the registry-owned `CoupledDomainRuntime` and
owns the production serial immersed dependency chain in its required lifetime
order: classified closed VTP surface, compact cut-cell volume quadrature,
surface quadrature, ghost catalog, static runtime, then the internal
`ThreeDImmersedFlowDomain` adapter.  Thus the adapter never borrows a runtime
owned by a separate container.  It reads a contained
`simulation_config.json` and `immersed_geometry.json`, rejects non-contained
assets, requires MPI size one, and accepts only `steady` backend
configurations for flow-only 3D Navier--Stokes.  Coupling still advances
quasi-static 3D load samples; this is not a `quasi_static` backend config
value.  Transport/species, body-fitted
mesh/profile/gauge assumptions, total pressure, and outlet models are
rejected.  The loader requires the VTP reader's closed, connected surface and
canonical hash, exact graph/runtime port ID-label-control agreement, and an
exact label partition between Nitsche walls and graph ports.

All JSON-to-integral conversions are finite, integral, and target-range
checked before conversion (including grid, quadrature/ghost limits,
`max_depth`, and PETSc nonlinear/KSP iteration caps).  Runtime density and
viscosity have one canonical source: the Navier--Stokes system configuration;
geometry-side physics overrides are rejected.  Floating-point manifest fields
are emitted with precision 17.

The generic multidomain runner now creates immersed cases through this factory
while preserving the existing body-fitted preflight and adapter path.  The
coupling executor is unchanged.  Balance checks derive 3D participation from
the graph kind/dimension and accepted `PortState` values rather than
body-fitted-native maps.  The generic binding manifest records immersed kind,
surface hash, Cartesian grid, and compact volume/surface/ghost catalog audit.
Schema-v6 transport remains an explicit rejected gate for immersed domains.

Focused validation from `solvers/coupling` used system PETSc 3.15
(`PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real`):

- `make PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real immersed_case_factory_test multidomain_config_test pressure_flow_executor_test iga_multidomain_flow multidomain_flow_smoke_test` — succeeded; all five targets were current.
- `./immersed_case_factory_test && ./multidomain_config_test && ./pressure_flow_executor_test && ./multidomain_flow_smoke_test` — succeeded (exit status 0).  The final smoke command invokes `iga_multidomain_flow` from its executable directory, exercising the generic immersed graph there.
- `git diff --check` — succeeded with no whitespace errors.

The factory fixture uses a tiny closed labelled VTP and covers registry-owned
lifetime binding, contained paths, topology, transport/species, outlet,
profile/gauge, total-pressure, grid/quadrature/ghost/runtime range, physics
override, and exact port ID/label/control negatives.  The generic-runner smoke
also executes an immersed fixture and audits kind, hash, grid, and catalogs;
the existing body-fitted generic smoke remains in that target.

## PR6.5 — Quasi-static immersed aneurysm-chain closure benchmark

Status: complete.

`examples/vascular_flow/immersed_aneurysm_chain/` contains the deterministic
steady low-Reynolds-number 1D--immersed-3D--1D case: a gradual octagonal
surface-of-revolution bulge, planar labelled inlet/outlet caps (1/2), and wall
label 0. It is deliberately a quasi-static steady backend example and makes
no transient or backward-Euler-equivalence claim.

On 2026-09-03, the composed production closure command passed on local serial
PETSc 3.15 in 1302.15 s real time (1300.81 s user, 0.20 s system):

```bash
make -C solvers/coupling phase6-aneurysm-closure-test PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
```

Hardware: Intel(R) Core(TM) i9-14900KF, one socket, 8 cores, 16 logical CPUs.
The target ran three explicit and three strong-fixed quasi-static 3D load
samples, rather than claiming transient equivalence. Its emitted coupled
maxima were `1.49078e-15` edge normalized flow residual,
`3.92737e-09` strong-fixed normalized pressure residual, 5 strong-fixed
iterations, `5.43931e-12` immersed 3D normalized balance, and
`5.4408e-12` net external 1D normalized balance. All manifests were audited;
the injected precommit failure published no completion marker; and retry output
was byte-for-byte replayed.

The direct immersed runtime evidence emitted all controller, open-port, and
wall metrics for depth-3 flow samples of `5e-05`, `7.5e-05`, and `1e-04 m^3/s`,
plus the depth-2 `1e-04 m^3/s` comparison. The largest reported controller
absolute residual was `2.71051e-20 m^3/s`, the largest normalized open balance
was `3.57677e-15`, and the largest normalized wall leakage was `1.573e-05`.
All samples converged in two nonlinear and two KSP iterations; the maximum
independently measured true-linear relative residual was `1.49106e-13`.
Cap area and normal errors were both exactly zero (8 triangles per cap).
Depth-3 resistance was `418.089 Pa s/m^3`, depth-2 resistance was
`431.908 Pa s/m^3`, the lubrication estimate was `483.139 Pa s/m^3`, the
lubrication difference was `0.134641`, and the depth difference was
`0.0319953`. The transaction replay recorded exact committed state, committed
port state, and nonlinear/KSP work replay; the retry counters were one commit,
prepare, and finalize, plus one rollback and one abort.

Already-observed phase-close regressions were not rerun: `make -C
solvers/coupling test` passed; `make -C solvers/coupling multidomain-test
PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real` passed with species
maxima `6.78e-20` edge, `9.56e-15` domain, and `8.89e-15` global; and the
root-supervised `make -C solvers/cpu petsc-test
PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real` passed after
test-local fixes, including the VCA two-rank `rtol=1e-12` check.

## PR6.6 — Conservative resolved mixed form

The immersed volume path now explicitly selects the conservative resolved
mixed form, while every other `BuildNavierStokesElementFromPoints` caller keeps
the legacy body-fitted default.  Its resolved terms are

\[
R_{u,p}^{\Omega}=\int_\Omega N_a\nabla p,\qquad
R_{p,u}^{\Omega}=-\int_\Omega\nabla N_a\cdot u.
\]

Every retained physical surface point on a positive cut cell supplies the
all-label completion

\[
R_{u,p}^{\Gamma}=-\int_\Gamma N_a p n,\qquad
R_{p,u}^{\Gamma}=+\int_\Gamma N_a u\cdot n.
\]

Thus the stored negative residual adds `+N_a p n` and
`-N_a(u\cdot n)`, and the two trace Jacobian blocks are exact negatives of one
another under transpose.  VMS and PSPG residuals and derivatives, including
the committed `25fd551` tangent, are unchanged.  The trace is assembled when
volume assembly is enabled, independently of Nitsche-wall selection.

The focused algebraic coverage uses a cubic body-fitted completion tolerance
of `2e-14` and a trace-block `J_{p,u}+J_{u,p}^T` tolerance of `2e-12`.
The multi-label runtime fixture checks that the summed physical pressure
residual is the open velocity trace plus a nonzero prescribed wall-normal
trace, rather than volume divergence minus measured wall flow; it also checks
total/open/wall diagnostic and per-label sums.

On local serial PETSc 3.15, the focused conservative-form depth-2 aneurysm
Jacobian/conservation regression completed in `35.18 s`. Its first accepted
candidate residual was `1.09978e-8`; the final nonlinear residual was
`2.50713e-17`; and the first true linear relative residual was
`1.49106e-13`. The inlet controller error was `2.71051e-20 m^3/s`. Production
quadrature reported open, wall, total, and volume-divergence flows of
`-2.08177e-19`, `1.57300e-9`, `1.57300e-9`, and `1.43606e-9 m^3/s`,
respectively. The open normalized balance and wall leakage normalized by
throughflow both pass the focused `1e-3` gates. The volume-divergence value is
retained as a diagnostic only. The focused PR6.6 evidence is now complemented
by the completed depth-3 coupled-chain closure in PR6.5.

## Phase 6

Status: complete. The closure claim is limited to deterministic quasi-static
3D load samples and does not assert transient or backward-Euler equivalence.
