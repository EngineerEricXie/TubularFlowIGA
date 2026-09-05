# Foundational Roadmap Final Report

Status: **complete**. The foundational Phases 0--9 are closed; this report
consolidates their authoritative evidence. It does not represent a
patient-specific clinical workflow or a fully production-hardened digital
twin.

## Final architecture

The delivered system has dependency-free coupling contracts over native,
transactional 1D, body-fitted 3D, immersed 3D, moving-domain, FSI, and 0D
runtimes. `SimulationGraph` owns validated heterogeneous topology;
`DomainRuntimeRegistry` owns exactly one matching runtime per domain; component
executors own deterministic graph traversal, strong coupling, and all-domain
rollback/prepare/finalize. Interfaces are SI and outward-positive, so pressure,
flow, and species conservation do not depend on domain ordering.

Model routing is explicit rather than inferred from names: schema v3/v4 retain
standalone 1D/3D behavior; schema v5 routes executable flow-only heterogeneous
graphs; schema v6 adds transactional 1D/3D species routing. Schema v5 also
routes one-port source/RCR 0D flow models, whereas schema v6 rejects 0D until
species semantics exist. Body-fitted and immersed 3D domains share coupling
semantics but retain backend-appropriate materialization.

## Phase closure record

| Phase | Delivered evidence |
| --- | --- |
| 0 | Common ports and rollback-safe native 1D/3D lifecycle, with legacy parity preserved. [Report](PHASE_0_REPORT.md) |
| 1 | Explicit/strong 1D--3D--1D straight-vessel coupling, spatial and temporal verification. [Report](PHASE_1_REPORT.md) |
| 2 | Immutable heterogeneous graph, schema-v5 runtime binding, branches, and multiple 3D regions. [Report](PHASE_2_REPORT.md) |
| 3 | Reversal-safe, conservative schema-v6 multi-species routing and MPI balance evidence. [Report](PHASE_3_REPORT.md) |
| 4 | Explicit immutable quadrature migration for CPU flow and transport kernels. [Report](PHASE_4_REPORT.md) |
| 5 | Closed-surface immersed Cartesian cubic B-spline flow with cut integration, Nitsche/ghost stabilization, and manufactured/sliver gates. [Report](PHASE_5_REPORT.md) |
| 6 | Transactional quasi-static immersed 3D graph runtime and aneurysm-chain closure. [Report](PHASE_6_REPORT.md) |
| 7 | Prescribed moving immersed anatomy, exact history/retry, moving-wall diagnostics, and idealized LV cycle. [Report](PHASE_7_REPORT.md) |
| 8 | Controlled two-way immersed FSI with membrane, traction transfer, and strong coordinator. [Report](PHASE_8_REPORT.md) |
| 9 | Transactional source/RCR 0D integration and real five-domain strong multiscale closure. [Report](PHASE_9_REPORT.md) |

## Major numerical and design decisions

- All interfaces use SI units and a local outward-flow convention. A pressure/flow edge requires conservative opposed outward flows; species donor choice is measured-flow based and reversal-safe.
- Physical state is committed/trial state, never implicit solver scratch. Rejected trials roll back; all-domain commit is prepare-all followed by no-throw finalization.
- 0D source and RCR capacitors use backward Euler with explicit signed storage, pump/sink, and graph-port amount accounting. Phase 9 proves component conservation after pairwise edge cancellation.
- CPU integration is explicit by rule type. Immersed flow uses a closed triangulated surface, Cartesian cubic B-spline background, adaptive cut-cell volume/surface quadrature, Nitsche wall terms, and ghost stabilization.
- Moving anatomy is fixed-Eulerian prescribed motion, not ALE. The Phase 8 FSI slice is partitioned Dirichlet--Neumann with dynamic Aitken, not monolithic or contact-capable.

## Completion-definition traceability

| # | Completion item | Authoritative evidence |
| ---: | --- | --- |
| 1 | Preserved standalone 1D functionality | [Phase 0](PHASE_0_REPORT.md) documented and checked legacy 1D cases. |
| 2 | Preserved standalone body-fitted 3D functionality | [Phase 0](PHASE_0_REPORT.md) compatibility/parity gates retain the validated 3D behavior. |
| 3 | Verified strongly coupled 1D–3D pulsatile flow | [Phase 1](PHASE_1_REPORT.md) provides the strong-coupling and pulsatile reference gates. |
| 4 | Multidomain heterogeneous simulation graph | [Phase 2](PHASE_2_REPORT.md) provides validated topology, runtime binding, and component execution. |
| 5 | Conservative cross-dimensional species transport | [Phase 3](PHASE_3_REPORT.md) verifies staged 1D/3D multi-species transfer and balances. |
| 6 | Generic quadrature abstraction | [Phase 4](PHASE_4_REPORT.md) records the immutable quadrature API migration. |
| 7 | Static arbitrary-surface immersed IGA flow | [Phase 5](PHASE_5_REPORT.md) provides closed-surface immersed-flow and manufactured-closure evidence. |
| 8 | 1D + immersed-3D multiscale coupling | [Phase 6](PHASE_6_REPORT.md) records the transactional immersed graph runtime and aneurysm-chain closure. |
| 9 | Prescribed moving-anatomy flow | [Phase 7](PHASE_7_REPORT.md) closes fixed-Eulerian prescribed-motion flow. |
| 10 | Verified foundational FSI capability | [Phase 8](PHASE_8_REPORT.md) verifies the bounded two-way immersed FSI vertical slice. |
| 11 | Integration of 0D/1D/3D components into the multiscale architecture | [Phase 9](PHASE_9_REPORT.md) records the real five-domain strong multiscale closure. |
| 12 | Persistent numerical diagnostics and documentation | [Phase 0](PHASE_0_REPORT.md), [Phase 1](PHASE_1_REPORT.md), [Phase 2](PHASE_2_REPORT.md), [Phase 3](PHASE_3_REPORT.md), [Phase 4](PHASE_4_REPORT.md), [Phase 5](PHASE_5_REPORT.md), [Phase 6](PHASE_6_REPORT.md), [Phase 7](PHASE_7_REPORT.md), [Phase 8](PHASE_8_REPORT.md), and [Phase 9](PHASE_9_REPORT.md) retain phase numerical evidence. |
| 13 | Relevant backward compatibility | [Phase 0](PHASE_0_REPORT.md) compatibility matrix and parity gates; [Phase 4](PHASE_4_REPORT.md) preserves body-fitted quadrature behavior. |
| 14 | Phase-specific verification evidence | [Phase 0](PHASE_0_REPORT.md), [Phase 1](PHASE_1_REPORT.md), [Phase 2](PHASE_2_REPORT.md), [Phase 3](PHASE_3_REPORT.md), [Phase 4](PHASE_4_REPORT.md), [Phase 5](PHASE_5_REPORT.md), [Phase 6](PHASE_6_REPORT.md), [Phase 7](PHASE_7_REPORT.md), [Phase 8](PHASE_8_REPORT.md), and [Phase 9](PHASE_9_REPORT.md). |

The reports above carry the numerical evidence: Phase 1's spatial/temporal
coupling verification; Phase 3's MPI species balances; Phase 5's manufactured
convergence and sliver spectra; Phase 6's real immersed closure; Phase 7's
accepted 16-step LV cycle; Phase 8's controlled FSI convergence; and Phase
9's analytic, temporal, conservation, retry, and serial/two-rank parity gates.

## Compatibility and final closing checks

The final Phase 9 source/harness review received fresh Sol **APPROVE**. The
sequential closing matrix passed `make coupling-test` and the focused real
PETSc Phase 6 aneurysm regression. CPU configuration/schema dispatch and its
initial unit gates passed; the intentionally expensive full CPU target was not
claimed as a whole-target pass after it entered cut-volume work. The 1D serial
core, coupling, runtime, and PETSc gates passed. Its two-rank test was blocked
only because this local PETSc lacks MUMPS LU support for `mpiaij`, after MPI
itself launched correctly. Details and measured values are in the [Phase 9
report](PHASE_9_REPORT.md).

## Known limits and recommended next steps

The foundational scope deliberately excludes patient-specific segmentation or
calibration, advanced valves/leaflet contact, a full closed-loop 0D heart,
0D species, ALE/remeshing, monolithic/nonmatching FSI, and broad production
I/O/restart hardening. Phase 9's two-rank root case is a parity check, not a
scaling study; four ranks are unsupported because the C2 root has only two
owned elements. Immersed FSI is still serial `PETSC_COMM_SELF`.

Recommended next work, in order:

1. Implement one-solve distributed immersed ownership, cut-cell load balance,
   and local OpenMP assembly; then establish meaningful rank sweeps.
2. Add robust multirate coupled clocks and restart/checkpoint ownership across
   0D/1D/3D graph components.
3. Extend 0D circulation to a closed-loop heart and species, with corresponding
   conservation and reversal contracts.
4. Pursue advanced valves, leaflet contact, and clinically grounded
   patient-specific inputs only as separately verified programs.
