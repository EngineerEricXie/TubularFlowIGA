# Phase 8 Closure Report

Status: **complete**.  Fresh Sol review: **APPROVE**.  Validation date:
**2026-09-05**.

Phase 8 closes the roadmap's first two-way FSI vertical slice: a controlled,
small-displacement compliant membrane over a rectangular immersed-flow
channel.  It is foundational evidence only.  It does not claim a compliant
tube, aneurysm wall, thin-shell valve, valve opening/closing, leaflet contact,
patient-specific valve, ALE/remeshing, nonmatching transfer, monolithic FSI,
or distributed FSI execution.

## Architecture and contracts

The slice composes the typed `FsiCouplingEdge` and exact surface contracts;
the patch-authoritative `MaterialSurfacePatchMap` and closed-material
composition; conservative fluid traction extraction/projection; the local
pre-tensioned membrane runtime; the moving immersed-flow runtime; and the
strong Dirichlet--Neumann coordinator with area-weighted dynamic Aitken
relaxation.  Fluid-on-structure traction is

\[
t_{\mathrm{on\ structure}}=-\sigma_f n_f.
\]

The interface carries structural displacement and velocity to the fluid and
fluid traction to the structure.  Material/topology, layout/partition, field
stamp, composition, and producer-state identities are exact contracts rather
than caller assertions.  The patch fixture is label 7, a 3x3 top-wall patch;
its outward normal is `+z`.  The accepted response therefore requires positive
center normal displacement/velocity and a positive `z` component of the
nodal traction resultant.

State remains transactional.  Trial input is accepted once, solve output is
unpublished until prepare, rejection/abort preserves the last committed epoch,
and finalization uses prevalidated non-throwing ownership exchanges.  The
coordinator rejects both trial owners before applying Aitken and commits the
fluid composition and structural state as the same FSI epoch.

## Numerical evidence

The manufactured membrane check measured spatial rates `1.93863` and
`1.97457`, and temporal rates `.98151` and `.990864`.

The controlled real compliant-channel execution converged in four coupling
iterations.  The accepted residual sequence (m) was
`[2.24949e-05, 9.19914e-06, 4.99426e-07, 2.71104e-08]`; the final threshold
was `4.80669e-08`.  Aitken factors were `.5`, `.8`, and `.8`; each fluid trial
reported Newton `2` and KSP `2` iterations.

| Metric | Value |
| --- | ---: |
| Center displacement | `3.80669e-05 m` |
| Center velocity | `7.61337e-04 m/s` |
| Open outward flow | `-6.84227e-05 m3/s` |
| Normalized moving mass | `.0155613` |
| Wall-relative leakage | `.0155613` |
| Discrete moving-wall continuity | `1.32707e-14` |
| Top-patch traction resultant magnitude | `.0171412 N` |

The final accepted values satisfy the hardened gates: moving mass `< .03`,
wall-relative leakage `< .03`, and discrete moving-wall continuity `< 1e-8`.
No additional coupled rerun is needed after this assertion-only tightening:
the captured result already satisfies those gates, and neither production code
nor the formulation changed.

## Runtime and performance context

The retained command measured `11:58.33` wall time, `716.23 s` user time,
`1.66 s` system time, `1136620 KiB` peak RSS (about `1.084 GiB`), and 99% CPU.
It used PETSc at `/home/tsungyeh/petsc` on a host with 16 logical and 8
physical CPUs.  The current runtime is sequential `PETSC_COMM_SELF` with one
partition, so an MPI rank sweep is inapplicable; distributed scalability and
rank/partition evidence carry into Phase 9.

## Focused coverage and lineage

Focused coverage includes the surface contracts and lifecycle, patch
composition/kinematics, conservative fluid traction projection, pre-tensioned
membrane and its runtime, moving immersed-flow FSI runtime, deterministic
strong-coordinator checks, and the controlled compliant-channel target.

Closure lineage is `d648af1`, `f895e04`, `93e59d0`, `b3ab9e4`, `c0a4790`,
`ab4842b`, `6547745`, `166bb92`, `d04b8a8`, `eaee2e3`, and `f7f80c7`, plus
this uncommitted closure slice.

## Deferred work

Phase 9 must make collective/distributed execution real before claiming MPI
scalability, including partition agreement and representative rank sweeps.
The roadmap's later compliant-tube, aneurysm-wall, and valve milestones remain
separate follow-ons, with nonmatching transfer, contact, ALE/remeshing, and
monolithic coupling still out of scope.
