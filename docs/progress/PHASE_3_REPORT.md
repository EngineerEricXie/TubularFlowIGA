# Phase 3 report

Status: **in progress**. PR 3.1 provides the schema and dependency-free
species contract. Native conservative transfer remains open for PR 3.2, and
the two-species production gate remains open for PR 3.3.

## Objective

Conserve arbitrary transported quantities across 1D/3D interfaces using
outward-positive total species flux and report independently normalized
per-species global balances.

## PR 3.1: general species port data

Schema v6 extends the flow graph without changing standalone schema v3/v4 or
flow-only schema v5. It declares logical species and concentration units at
graph scope, maps each logical ID to a native scalar field in every domain,
and lists sorted species IDs on every participating port and pressure-flow
edge. Flux units are explicitly derived as concentration unit times `m^3/s`.
Arbitrary logical and native names are supported; no species is privileged.

Both endpoints must provide and accept concentration and total outward flux.
This is a deliberate reversal-safe capability contract, independent of the
hydraulic pressure/flow receiver roles. Validation rejects empty, duplicate,
or undeclared species, missing domain bindings, reused native fields,
endpoint/edge set differences, and incomplete species capabilities.

`ResolveSpeciesDonor` uses only measured outward flows. It requires a
conservative opposite-sign pair, supports forward and reversed flow, and uses
the last committed donor as the near-zero hysteresis owner. Same-sign flow,
excess flow residual, nonfinite input, and near-zero flow without committed
ownership fail explicitly. Existing native runners reject schema v6 until
PR 3.2 makes transport transactional.

Focused validation:

```text
make -C solvers/coupling test
make -C solvers/cpu coupling_port_test config-check-dispatch-test
./solvers/cpu/coupling_port_test
./solvers/cpu/vca_3d_runtime_test
```

The focused legacy VCA runtime regression passed with the local PETSc/MPI
toolchain. The PR-closing Sol review accepted schema ownership, unit and sign
semantics, reversal/hysteresis routing, graph binding invariants, and v3/v4/v5
compatibility.

## PR 3.2 numerical gate

### Transactional 3D transport foundation

`TransientTransportRuntime` now has the same trial lifecycle needed by a
composite flow/transport adapter: `BeginStep`, deterministic `SolveTrial`,
`RollbackTrial`, idempotent `AbortStep`, two-stage prepare/finalize, and the
existing `Advance` compatibility wrapper. The committed PETSc concentration
vector and logical step count are snapshotted and restored together. Failed or
replayed trials cannot increment physical time state, and finalization is
nonthrowing. The legacy VCA runtime regression covers exact rollback, exact
replay, one commit, open-step abort, and `Advance` parity.
The focused Sol review accepted physical-state ownership, failed-solve
recovery, prepared-step abort, nonthrowing finalization, checkpoint refresh,
and legacy wrapper semantics.

The remaining numerical gate is:

Flow and species will share one graph transaction: converge flow, hold the
accepted velocity, perform rollback-safe species trials, validate interface
and global balances, then prepare all domains before any finalization. The 3D
port flux is

```text
integral((c*u - D*grad(c)) dot n dA).
```

The 1D runtime must report the exact finite-volume face flux used by its
conservative `A*C` update, including diffusion. Global balance must use
time-integrated boundary transfer and sources, not only end-step flux. A
species failure aborts the still-open flow transaction.

The first 1D flux slice now centralizes the numerical face formula
`Q*C_upwind - A*D*(C_right-C_left)/dx` for the conservative `A*C` update. The
runtime stores and reports the actual final-substep numerical boundary flux;
an immediate root bifurcation is aggregated across every root face before the
root outward sign is applied. Terminal reporting uses its stored numerical
face flux with the currently supported zero-gradient exterior trace. Focused
core and runtime regressions pass, and focused numerical re-review accepted
the aggregation, publication timing, and orientation. Time integration of
those face transfers and reversal-safe outlet boundary inputs remain part of
the PR 3.2 gate.
