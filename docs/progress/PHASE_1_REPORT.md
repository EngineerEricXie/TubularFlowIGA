# Phase 1 report

Status: PR 1.2 complete; Phase 1 remains in progress pending Aitken relaxation,
subcycling, and the full convergence benchmark. Restart and a
general multidomain graph remain deferred.

## PR 1.1 explicit 1D--3D--1D smoke evidence

The PETSc driver `iga_1d_3d_explicit` advances a flow-only upstream 1D case,
one body-fitted 3D case, and a flow-only downstream 1D case in that explicit
order. It records all six external/interface port pressure, outward-flow, and
area values; two conservative edge residuals; total 3D boundary mass residual;
wall contribution; full external pressure drop; and trial linear iterations.
It writes a rank-zero CSV and manifest only after all requested steps commit.

`solvers/coupling/tests/test_explicit_coupling_smoke.cpp` creates a unique
temporary version-5 `.ntiga` database plus 1 m straight one-element 3D and
one-segment rigid 1D cases at run time. No generated input, result, or database
is committed. It runs the driver in one and two ranks and checks fieldwise
history agreement using `1e-10 * max(1, |a|, |b|)` hybrid absolute/relative
tolerance, all outward-flow signs,
interface residuals, positive/matched areas, total 3D relative mass imbalance
below `1e-3`, finite pressure diagnostics, and rejected-step output
suppression. `TUBULARFLOWIGA_KEEP_TEST_OUTPUT=1` retains its temporary files
for diagnosis.

The low-Re constant-flow fixture uses `Q=1e-3 m^3/s`, `rho=1e-6 kg/m^3`,
`mu=1 Pa s`, `dt=0.01 s`, and three steps. Its predeclared analytic comparison
is two circular 1D Poiseuille segments plus one unit-square straight-duct
resistance: `0.0787196354574 Pa`. The intentionally single cubic 3D element
measured `0.0718643347844 Pa`, passing the predeclared relative tolerance of
`0.5`. The all-circular three-segment reference is `0.0753982236862 Pa`.
Maximum cap-flow, wall-flow, total-3D-mass, and external-rigid-1D relative
imbalances were respectively `1.95e-15`, `0`, `1.95e-15`, and `2.39e-15`
(the two-rank maxima were no larger than `1.09e-15`). Pressure exchange remains
one-step lagged throughout. This is a coarse execution smoke gate only, not
the final straight-vessel discretization-convergence benchmark.

Validated locally with PETSc 3.15:

```text
make coupling-test
make -C solvers/coupling petsc-test PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
git diff --check
```

The MPI smoke passed one rank, two ranks, and injected pre-commit failure. It
uses `-ksp_type preonly -pc_type lu` internally so trial linear-iteration
diagnostics as well as physical history are rank-invariant for this tiny test.

## PR 1.2 strong fixed-relaxation evidence

`--coupling-mode strong-fixed` retains the PR 1.1 straight-chain scope and
starts all three runtime trials once per physical step. It iterates the two
applied pressures with fixed `omega=0.5`; the two fixed-point residuals are
normalized by the required positive `Pref`, while conservative flow-transfer
residuals remain a hard gate. Each rejected sweep rolls back downstream, 3D,
then upstream before reapplying every input. Only the converged sweep commits.
The strong-only output is rank-zero `strong_coupling_history.csv`,
`strong_coupling_iterations.csv`, and `strong_coupling_manifest.json`; no
strong artifact is written for max-iteration or injected pre-commit failure.

The same generated low-Re fixture passed one and two ranks with `Pref =
0.0787196354574 Pa`, pressure tolerance `1e-6`, flow tolerance `1e-10`, and
maximum 50 iterations. The first step took 21 sweeps with final normalized
residuals `5.6614468717e-7` and `8.8146926286e-17`, accepting 21 3D KSP
iterations from 441 total (420 rejected). The two constant subsequent steps
converged in one sweep (6 and 0 accepted/all KSP iterations respectively), as
expected from seeding with the accepted fixed point. The smoke also checks the
exact `x + 0.5(G-x)` update, rank-parity of both physical and iteration
histories, flow/mass/area gates, a forced one-iteration nonconvergence, and
the existing injection hook in strong mode. This remains a coupling-lifecycle
smoke, not numerical convergence evidence.

The strong iteration CSV records physical step and sweep index; applied,
measured, and relaxed-next pressure guesses; signed Pa and normalized pressure
residuals; all four constituent outward interface flows; separate signed/raw
and normalized flow residuals; wall and total 3D mass diagnostics;
per-attempt/cumulative 3D KSP work; and convergence. The
step CSV records final signed/normalized pressure residuals and accepted/all/
rejected KSP work. The manifest records the fixed-point and outward-flow
formulas, physical inputs, port orientations, Newton controls, initial
guesses, and post-success-only output/work semantics. Rejected 1D work is not
reported, while the 1D internal-substep counter remains rollback-owned.

Architecture Gates A, B, and C pass for PR 1.2. The serialized signed pressure
residual vector and applied/measured/next pressure states provide the history
required by PR 1.3 without changing the runtime lifecycle or port contract.
