# Phase 7 Closure Report

Status: **complete**. Validation date: **2026-09-04**.

Phase 7 closes prescribed moving anatomy on a fixed Eulerian Cartesian
background. It is explicitly not ALE, FSI, patient-specific/clinical, or a
tracer/residence-time model.

## Gate disposition

| User-defined gate | Closure evidence |
| --- | --- |
| Geometry epochs | Every evaluated surface time owns an immutable full rebuild: surface, classification, cut volume, surface quadrature, and ghost chain. Geometry/publication identities bind the epoch and transition; both layout directions ran. |
| State/history transfer | Old velocity transfers only to the exact target Cartesian layout. Committed anchors remain bitwise exact; extension must cover the target and passes deterministic map/operator/history identities plus fail-closed coverage, pivot, and resource checks. |
| Moving-wall Nitsche | Target material provenance evaluates prescribed wall velocity at retained quadrature points; moving-wall Nitsche diagnostics and the wall-relative gate passed. |
| Moving mass/geometric conservation | Outward, target-endpoint diagnostics report volume-rate, material/fluid flux, Reynolds, moving mass, continuity, and wall-relative leakage. |
| Nontrivial temporal/refinement benchmark | The 16-step ED/ES/ED benchmark passed and the one-step same-field depth audit reports the U1/U2/U3 and divergence results below. |
| Idealized LV end-to-end | The deterministic 12-sector/four-ring prolate LV, with basal inlet/outlet, completed its prescribed cycle and captured Q/enstrophy/turnover/stagnation proxies. |
| Fresh Sol numerical closure approval | Fresh Sol review found no blocker, high, or medium issue; gates 1--7 are approved. No full-cycle rerun was needed: post-cycle changes were documentary, test-fixture, edge-semantics, or post-assembly diagnostic aggregation, and the current one-step run exercised the changed path. |

## Formulation and transactional state

The normal convention is outward. At the target endpoint,

\[
G_{BE}=(V^{n+1}-V^n)/dt,\qquad R_{Reynolds}=G_{BE}-Q_w,
\qquad R_{moving}=G_{BE}+Q_u-Q_w.
\]

The endpoint continuity record and wall-relative leakage are

\[
R_{cont}=Q_{port}+Q_{w,wall},\qquad Q_{u,wall}-Q_{w,wall}.
\]

`R_div = integral(div u) - Q_u` remains a volume/surface quadrature
consistency diagnostic, not a moving-wall continuity gate or hidden GCL
correction. Phase 7 introduces no ALE convection, mesh Jacobian,
swept-volume/space-time integration, or FSI.

The symmetric Nitsche form retains its traction, adjoint-consistency,
pressure-gap, and velocity-pressure blocks. Its penalty is

\[
\eta=\eta_\mu+\eta_t,\qquad
\eta_\mu=16\gamma_\mu c_\alpha\mu/h_n,\qquad
\eta_t=16\gamma_t\rho h_n/dt.
\]

The inertial term is zero for `dt=0`; it is impedance stabilization, not a new
wall model.

Geometry, layout, global state, port records, velocity history, and
conservation record have separate ownership. A step is `idle -> trial ->
prepared -> finalized`: `PrepareCommit()` creates an unpublished replacement,
while `FinalizeCommit()` makes the non-throwing owner exchange. Rollback
restores the frozen seed/history and exact replay identities; abort discards an
unpublished trial/prepared owner. Transfer requires exact map and immutable
source/target geometry identities, complete target history, and bitwise
preservation of committed anchors.

## Geometry, quadrature, and output contract

The LV benchmark uses a fixed `6x6x7` grid, depth-2 cut quadrature, and
empty-rule rescue through depth 9. The rescue rule fails closed (the same probe
rejects a cap of 8), including the tested rescue-9 cells. Expanded and compact
catalogs preserve their supported polynomial moments; snapshot parity is
required to `1e-12` for volume, integrated/mean enstrophy, and mean Q.
Q-positive and stagnant fractions are finite sampling proxies, not required
compact/expanded pointwise matches.

The sufficient Phase 7 output contract is schema-v3 `fields.vtu` and
`metrics.json` per committed epoch plus the PVD collection. [Post-Phase-7
publication I/O hardening](../POST_PHASE_7_IO_HARDENING.md) is explicitly
deferred and is not a Phase 7 gate.

## Numerical evidence

### Accepted full cycle

The earlier accepted full-cycle run completed 16 ED/ES/ED steps with depth
2/rescue 9, `rho=1050`, `mu=0.012`, period `0.8 s`, nominal `dt=0.05 s`, in
wall-clock `4144.761158443 s`. ES was at `0.4 s`; final ED was at `0.8 s`. It exercised
both layout directions, wrote 17 snapshots/PVD entries, and verified
rollback/retry plus Q, enstrophy, turnover, and stagnation.

| Metric | Accepted maximum/value |
| --- | ---: |
| true linear residual | `1.5491925338140507e-12` |
| controller residual | `3.8116482626443515e-21` |
| wall-relative ratio | `0.07218720443602622` |
| moving-wall continuity | `3.7581576788760377e-13` |
| divergence quadrature (informational) | `0.011798596501757136` |
| Reynolds defect | `0.010851693521875436` |
| moving-mass defect | `0.010852383873688134` |
| wall leakage | `0.0002429819483780136` |

### Current one-step closure path

The current `--one-step` evidence exited 0: `558.551811729 s` internal,
`9:20.27` external, `243500 KiB` maximum RSS, and 99% CPU. Solve took
`273.728010662 s`, rollback replay `223.055530644 s`; it converged in Newton
3/KSP 3.

| Metric | Current value |
| --- | ---: |
| true linear residual | `1.1734064164541041e-13` |
| controller residual | `1.6940658945086007e-21` |
| wall-relative ratio | `0.04898440928839063` |
| moving-wall continuity | `2.1368688405307986e-13` |
| divergence quadrature (informational) | `0.0041981219914684238` |
| Reynolds defect | `0.010037561945652409` |
| moving-mass defect | `0.0099243580409942422` |
| wall leakage | `0.00011320390445221881` |

The same solved endpoint field was reintegrated, not re-solved, for the
refinement audit. `U1/U2/U3` were `5.7567301587301505e-5`,
`2.9320444444444379e-5`, and `1.459561904761912e-5`; `|Rdiv|` was
`2.77346899600701e-7`, `1.4241138895048124e-7`, and
`5.9829120961597169e-8`. Endpoint differences refined accordingly. This is
same-field quadrature evidence, not PDE spatial or temporal convergence.

The focused Newton regression verifies the exact solve in 7 updates rather
than the obsolete hard-coded 8. Its final residual was
`5.722864349262003e-14 <= 1.0945723152615691e-9`; true linear residual was
`1.4736950373295145e-15`. Robust formulation gating and exact replay remain.
This is test hardening, not numerical weakening.

Latest aggregate components passed through extension. Focused moving
transient/capture and publisher tests passed, as did the current one-step and
the earlier full 16-step run. The monolithic immersed-transient full rerun
reached the now-corrected brittle assertion; it is not claimed as a complete
latest rerun. The focused exact-solve evidence above replaced that obsolete
assertion.

## Targets

Use the configured PETSc location:

```bash
make phase7-focused-test PETSC_DIR=/path/to/petsc PETSC_ARCH=your-arch
make phase7-lv-closure-test PETSC_DIR=/path/to/petsc PETSC_ARCH=your-arch
make -C solvers/cpu phase7_lv_closure_test PETSC_DIR=/path/to/petsc PETSC_ARCH=your-arch
./solvers/cpu/phase7_lv_closure_test --one-step
./solvers/cpu/phase7_lv_closure_test /path/to/retained-output/lv_cycle
```

`phase7-focused-test` covers moving-cut/runtime/snapshot functionality and
cut, ghost, Nitsche, port, history, and publisher prerequisites. The LV target
is full-cycle mode; an output directory retains the schema-v3 artifacts rather
than temporary output.

## Performance and Phase 8 carry-forward

This workstation has 16 logical CPUs, 8 physical cores, 47 GiB RAM, and one
NUMA node. The one-step run used 99% of one CPU, identifying serial local
element, cut-quadrature, Nitsche, and ghost work. Phase 8 should benchmark
1/2/4/8/16 MPI ranks only for truly distributed representative solves and
pursue safe OpenMP/local parallelism and cut-cell load balancing; it should
not run duplicate rank-local cases.

## Known limits

There is no PDE temporal-order convergence claim; refinement is same-field
quadrature. There is no FSI/ALE/tracer residence calculation. Q, stagnation,
turnover, and washout are proxies. Advanced I/O remains deferred as linked
above.
