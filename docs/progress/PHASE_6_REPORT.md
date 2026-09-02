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
