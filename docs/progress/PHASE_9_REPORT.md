# Phase 9 Closure Report

Status: **complete**. Fresh Sol review: **APPROVE**. Validation date:
**2026-09-05**.

Phase 9 closes the foundational roadmap with a bounded, real multiscale
circulation case. It is not a claim of a patient-specific digital twin, a
closed-loop heart model, advanced valves/contact, distributed immersed flow,
or production-hardened restart/I/O.

## Lineage and architecture

The accepted PR9 lineage is `25a3acd` (dependency-free 0D contracts),
`3bc3346` (transactional 0D runtime), `f9d7332` (generic schema-v5 runner
integration), and `dd6b4a1` (bounded closure harness).

The production `iga_multidomain_flow` runner executes an acyclic five-domain,
four-edge tree:

```text
source 0D -> body-fitted 3D root -> rigid 1D tree -> rcr_a 0D
                                                   -> rcr_b 0D
```

The generic pressure/flow component executor, not a bespoke closure solver,
uses strong fixed-point iteration with relaxation `.5`. The 1D horizon equals
the 3D macro grid, so each accepted macro step advances the native 1D runtime
once. Every domain is registry-owned behind `CoupledDomainRuntime`; graph
topology is immutable and the component owns the ordered iteration state.

## 0D formulation, signs, and state ownership

All 0D quantities are SI and `outward_flow_m3_s > 0` means flow leaves that
0D subsystem. A source accepts interface pressure `p_Gamma` and its
backward-Euler state is

\[
p_s^{n+1}=\frac{(C_s/dt)p_s^n+Q_p+p_\Gamma/R_s}
{C_s/dt+1/R_s}, \qquad Q_{s,out}=\frac{p_s^{n+1}-p_\Gamma}{R_s}.
\]

For a terminal, the port convention is `Q_port=-Q_in`; reverse flow is valid
and is never clamped. Its RCR update and reported port pressure are

\[
p_c^{n+1}=\frac{(C_t/dt)p_c^n+Q_{in}+p_v/R_d}
{C_t/dt+1/R_d}, \qquad p_\Gamma=p_c^{n+1}+R_pQ_{in}.
\]

`ZeroDFlowDomainRuntime` owns immutable model/port metadata and one committed
pressure/time/index image. `BeginStep` takes only the exact next clock state;
each `SolveTrial` is evaluated from that frozen image. Rollback or abort
discard every trial input, port record, and balance. Prepare validates a staged
state and accounting record; idempotent no-throw finalization publishes the
pressure and advances the clock exactly once. SHA-256 identities bind model,
state, and per-step accounting.

The terminal models deliberately differ: `rcr_a` has `Rp=10`, `Rd=100`,
`C=5e-4`; `rcr_b` has `Rp=20`, `Rd=80`, `C=6.25e-4` (SI), with both capacitor
pressures initialized to zero. Strict graph/manifest checks retain the exact
`tree.terminal_a -> rcr_a.port` and `tree.terminal_b -> rcr_b.port` bindings,
distinct model identities, one legal `zero_d_port` per 0D domain, contained
relative model files, one source start, terminal-RCR leaves, and no 0D--0D
edge. Schema v6 rejects a 0D domain because a 0D species contract does not yet
exist.

## Closure evidence

The final source/harness review approved the following production-gate
evidence. The harness checks the actual schema-v5 runner outputs, including
every accepted step, port, edge, iteration, 0D-history, and 3D-balance row.

- Strong coupling is fixed `.5`, with a configured cap of 50 iterations and
  an observed maximum of 24. The final accepted maxima were pressure residual
  `9.96463e-7`, flow residual `2.88598e-16`, 0D accounting residual
  `4.74338e-20 m3`, and normalized 3D balance `2.55346e-14`.
- The independent analytic settled oracle combines source resistance, the C2
  square-duct 3D resistance, rigid branch resistances
  `8 mu L/(pi r^4)`, and each RCR branch/parallel combination. It verifies
  separately the source, both capacitor pressures, both terminal port
  pressures, and both flows; the intentional A/B ordering prevents a copied
  terminal model from passing.
- Per-step conservation independently recomputes
  `sum(delta V_0D) - V_pump + V_sink,a + V_sink,b`; graph-port amounts are
  checked against edge flows and cancel pairwise. Every edge also checks
  outward-flow cancellation and `p_measured = p_applied`; the existing 3D
  boundary mass audit remains active.
- The clock is exact by construction: each step starts at the preceding
  `EndTime`, rather than `step*dt`, and all emitted timestamps are checked.
  Same-end-time `dt`, `dt/2`, and `dt/4` runs require first-order refinement
  ratios near three when `dt/4` is the reference.
- Serial and two-rank results are compared numerically for accepted pressure,
  flow, 0D state, 3D balance, and complete per-edge iteration diagnostics.
  The harness also injects a precommit failure, requires no completion marker,
  and requires the clean retry to reproduce accepted CSV history exactly.

Earlier timing measurements (`19.01 s` / `44440 KiB` serial and `15.12 s` /
`41872 KiB` at two ranks) predate the asymmetric terminal configuration and
are therefore not presented as final performance evidence. The validation host
has 16 logical / 8 physical CPUs. Four ranks are deliberately unsupported for
this fixture: its C2 3D root owns only two elements, so additional ranks would
be empty rather than a meaningful scaling point. The dominant work remains
sequential strong-trial replay. Immersed FSI remains `PETSC_COMM_SELF`; its
distributed scalability is deferred.

## Compatibility and fresh closing matrix

The phase preserves standalone schema-v3/v4 paths, `.ntiga` format, native
1D behavior, and schema-v5 flow graphs. Schema-v6 species routing remains
supported for its existing 1D/3D scope but explicitly rejects 0D models.

The proportional closing matrix on this revision was sequential:

| Command | Result | Scope |
| --- | --- | --- |
| `make coupling-test` | pass | Eight dependency-free coupling gates, including graph, pressure/flow, species, schema/config, 0D, surface, and FSI contracts. |
| `make cpu-test` | compatibility dispatch and its initial CPU unit gates passed before the deliberately expensive cut-volume test was stopped; not recorded as a whole-target pass | CPU configuration/schema dispatch and base geometry path. |
| `make one-d-test PETSC_DIR=/home/tsungyeh/petsc` | serial core, coupling, runtime, and PETSc gates passed; two-rank subtest blocked by PETSc configuration | MPI reached successfully, then requested unavailable MUMPS LU for `mpiaij`; not a production numerical failure. |
| `make -C solvers/cpu phase6-aneurysm-depth2-regression PETSC_DIR=/home/tsungyeh/petsc` | pass | Real focused immersed Jacobian, solve, conservation, and port regression. |

The focused immersed gate measured FD block defects through `1.41286e-11`,
final nonlinear residual `2.49118e-17`, inlet-target error `2.71051e-20`,
open normalized balance `2.58187e-15`, and wall-leakage normalized
`1.573e-5`. No Phase 7 LV, Phase 8 compliant-channel, or Phase 9 full
multiscale rerun was required: their authoritative reports and final Phase 9
source/harness review cover paths unaffected by subsequent production changes.

## Deferred work

Next work should make immersed assembly/ownership distributed in one solve,
add cut-cell load balance and local OpenMP assembly, then address multirate and
restart for coupled 0D/1D/3D graphs. Advanced valve/contact, a full closed-loop
0D heart, and 0D species remain out of scope.
