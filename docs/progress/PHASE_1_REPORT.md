# Phase 1 report

Status: PR 1.4 complete; Phase 1 remains in progress pending the full
straight-vessel convergence benchmark. Restart and a
general multidomain graph remain deferred.

The 3D nonlinear budget remains backward-compatible at 30 iterations and is
now configurable in every coupling mode with `--three-d-max-newton N` (positive
integer). The selected value is recorded as
`newton_controls.maximum_iterations` in the corresponding manifest. A
preliminary H2 steady-ladder run exceeded the default budget and passed with a
budget of 60 on the current fixture; it is evidence for the configurable
budget only, not a completed convergence gate.

## PR 1.4 time-subcycling status

The 1D runtime now executes an integer number of configured-time substeps
inside one rollback-safe macro lifecycle transaction, sampling configured
upstream waveforms at endpoints and holding interface data across the macro
interval. The dependency-free preflight validates ratio, reconstructed time,
step count, and horizon consistency; fast tests cover `N=1`, `N=2`, `N=4` and
invalid ratios/horizons. Runtime tests establish rigid `N=4` equivalence with
four committed `N=1` steps plus exact rollback replay. This remains
zero-order-held coupling, not a subcycling convergence benchmark.

The generated PETSc smoke additionally runs independent upstream/downstream
`dt=0.0025 s`, `steps=12` cases beneath the unchanged 3D `dt=0.01 s`,
`steps=3` macro grid for explicit, strong-fixed, and strong-Aitken modes on
one and two MPI ranks. It verifies endpoint histories against the unsplit
rigid reference, rank parity, and per-domain work identities: an explicit step
has accepted/all/rejected configured work `4/4/0`; strong attempts have four
configured substeps each, so accepted work is four and all/rejected work is
the actual sum over attempts and its accepted-work difference. For the
constant-work fixture these reduce to `4*iteration_count` and
`4*(iteration_count-1)`. A synthetic variable-work audit proves that the
implementation does not infer totals from the final attempt. The N=4 fixture
retained the existing deterministic 21-sweep fixed and 3-sweep Aitken first
macro step.
The runtime fast suite additionally verifies N=4 RCR capacitor/pressure state
against four committed N=1 steps and an injected implicit failure at configured
substep three: diagnostics report planned/attempted/completed `4/3/2`, rollback
restores the initial state, and retry executes four `0.0025 s` callbacks.
Configured open-loop diagnostics record the exact endpoint sequence
`.0025,.005,.0075,.01 s` for the N=4 macro trial. Manifests now record both
plans, ratio tolerance, ZOH/endpoint/CFL semantics, and summed per-domain 1D
accepted/all/rejected configured and CFL work.

The runtime suite also compares every transported species field and dynamic
vasodilation state after one N=4 macro trial with four committed N=1 trials.
An explicit Rusanov failure after one completed configured substep preserves
the partial CFL-work delta, restores the committed state, and replays the
frozen inlet schedule exactly on retry. The MPI smoke independently sums the
iteration CSV into step accepted/all/rejected totals and then cross-checks the
manifest for all three coupling modes on one and two ranks. Invalid ratios and
horizons fail before producing output, and injected or maximum-iteration
failures leave no artifacts.

Architecture Gates A, B, and C pass for PR 1.4. The full PETSc/MPI suite and
the dependency-free 1D and coupling suites passed on the final snapshot.

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

## PR 1.3 Aitken relaxation evidence

`strong-aitken` applies a rollback-safe scalar Aitken proposal to the existing
two-pressure strong vector, resetting at every physical step. The generated
one- and two-rank smoke fixture completed its first strong step in 3 sweeps,
versus 21 with fixed `omega=0.5`; subsequent seeded constant steps require one
sweep. The smoke independently reconstructs every serialized Aitken proposal
from `r=G-x`, including the scaled-before-subtract dot products, clamp/fallback
status, next guess, reset-first-row state, and the final unused converged
proposal. It also verifies the per-step summary: updates equal sweeps minus
one, the last applied factor equals the penultimate proposal (or explicit zero
for a one-sweep step), the final proposed factor equals the final row, and the
legacy step factor is the last applied value or the configured initial factor.

The default Aitken bound is `1`; the deterministic MPI fixture uses the
documented `0.999` maximum only to avoid a last-bit rank-dependent exact-bound
classification while auditing status codes. This does not weaken the default
fast test, which separately proves the default `[2,0] -> [1.5,0]` proposal
clamps an unclamped `2` to `1` with `ClampedMaximum` status. The generated
manifest is checked for the signed recurrence, `omega_hat=-omega_previous`,
scale-before-subtract rule, fallback retention, reset, and final-unused
semantics; the fixed manifest is checked to contain no Aitken block.

The current one-rank physical-history comparison reports a maximum pressure
difference of `2.1200915592545222e-08 Pa` (`external_pressure_drop_pa`) versus
the predeclared `3.1487954182974675e-07 Pa` tolerance
`4*Pref*pressure_tol + 1e-12`; the maximum nonpressure difference is
`1.0842021724855044e-18 m3/s` (`downstream_root_outward_flow_m3_s`) versus
the `1e-10*max(1,|a|,|b|)` hybrid tolerance. The smoke checks full strong
physics/work/conservation audit for fixed and Aitken, one/two-rank parity,
negative rejection of Aitken-only minimum/maximum options in explicit and
fixed modes, and no artifacts on Aitken max-iteration or injected failures.
Those failure paths also assert the precision-17 complete current-step
iteration diagnostics, including `G`, normalized pressure and flow residuals,
omega/status, and KSP totals. This is an execution smoke only, not a
convergence benchmark.

Architecture Gates A, B, and C pass for PR 1.3. The dependency-free vector
Aitken implementation is reusable by later partitioned coupling work without
changing the 1D or 3D runtime lifecycle.

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
