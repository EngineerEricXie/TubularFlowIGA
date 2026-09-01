# Phase 1 report

Status: PR 1.4 and the straight-vessel steady spatial gate are complete. The
pulsatile harness establishes rank-one, fixed-spatial-mesh temporal
self-consistency and coupling-tolerance sensitivity, but is not yet the
physically equivalent all-1D wave-reference gate required to close Phase 1.
Restart and a general multidomain graph remain deferred.

## Pulsatile straight-vessel temporal and tolerance verification

`make coupling-temporal-convergence-test` is deliberately separate from the
fast smoke and steady spatial targets.  It creates every input at runtime and
does not commit a database, VTK/control mesh, velocity field, case, or result.
The fixture is the existing open-uniform cubic C2 unit-square duct with two
transverse and two axial elements.  This is a fixed, modest spatial mesh: it
keeps the test practical and is intentionally *not* a replacement for the
completed spatial-convergence gate above. The temporal report therefore labels
coupled-versus-reference pressure differences as a combination of fixed
spatial error, temporal/coupling effects, and the rigid algebraic 1D versus
unsteady 3D model mismatch. Only comparisons made at the same spatial mesh and
between coupled cases isolate temporal or coupling changes.

The generated upstream 1D input is a sinusoid with mean `1e-3 m3/s`, 20%
amplitude, and period `0.16 s`.  Density is `1 kg/m3` and viscosity is
`1 Pa s`, retaining low Reynolds number while making the unsteady 3D response
measurable over the pulse. The simulation advances four periods and measures
the final period, separating the convergence metrics from the initial
adjustment of the steady mean-flow profile; an explicit cycle-to-cycle audit
also requires the third and fourth pressure and flow waveforms to agree within
`1e-10` relative L-infinity. The three-dimensional macro grids
are `dt = 0.01, 0.005, 0.0025 s`, giving respectively 16, 32, and 64 endpoint
samples per period. Both 1D domains use `N=2` configured substeps per macro
step, so their open-loop sinusoid is sampled at substep endpoints while the
interface data are held zero-order over the macro interval, exactly as
documented by the coupling architecture. Strong Aitken is used with a tight
pressure reference tolerance `1e-10`, flow tolerance `1e-10`, and a
60-Newton-iteration 3D budget.

For every run the harness parses the rank-zero strong history, iteration
history, and manifest.  It checks macro time alignment, finite P/Q data,
conservative interface transfer, 3D mass and wall-flow diagnostics, committed
iteration counts, Aitken factors, accepted/all/rejected 1D work identities,
and the manifest's endpoint/ZOH/subcycling semantics.  It calculates L2 and
L-infinity waveform differences on nested macro grids, first-harmonic pressure
drop amplitude and phase, and checks that the medium--fine waveform difference
does not exceed the coarse--medium difference apart from floating-point
roundoff.  This is a temporal self-convergence comparison at fixed spatial
resolution; it makes no claim of a spatially converged pulsatile field.

The same generated geometry also supplies an independently advanced all-1D
reference: a single rigid three-segment runtime with the area-matched circular
radius and the square duct's hydraulic-equivalent middle length.  It advances
on the 1D substep grid, independently of the coupling driver.  The harness
requires external Q(t) agreement at macro endpoints and reports pressure-drop
waveform, amplitude, and phase differences.  The all-1D pressure waveform is
a hydraulic reference, not a claim that a circle and a square are geometrically
identical; it is most directly meaningful for the external-flow transfer and
the steady hydraulic component.  The fixed square-duct spatial error and any
resolved 3D inertial response remain visible rather than being hidden.

The fine `dt=0.0025 s` case is repeated at pressure tolerances `1e-4`, `1e-6`,
and `1e-10` (the last is the tight reference).  The reported pressure-drop and
upstream-root waveform distance to the tight result must be non-increasing as
the configured tolerance tightens.  This is an ordering check with only a
machine-roundoff allowance, not an arbitrary absolute waveform threshold.

The dependency-free generation/audit command is:

```text
make coupling-temporal-fixture-test
```

For a quick PETSc diagnostic of the first 16-sample-per-period macro grid,
without claiming the complete ladder, run:

```text
TUBULARFLOWIGA_TEMPORAL_SINGLE_CASE=1 \
  ./solvers/coupling/pulsatile_straight_vessel_temporal_test
```

The full command must be run where MPI local sockets and the full runtime are
permitted:

```text
make coupling-temporal-convergence-test
```

Set `TUBULARFLOWIGA_KEEP_TEST_OUTPUT=1` to retain the generated cases and CSV
artifacts for review.

The full target passed locally with PETSc 3.15. Across the final period, the
maximum normalized pressure L2 difference decreased from `0.00955211` on the
coarse--medium comparison to `0.00538689` on medium--fine; the corresponding
L-infinity difference decreased from `0.0157750` to `0.00911057`. Every
pressure field passed separately. In particular, the 3D-outlet pressure L2
difference decreased `0.00353307 -> 0.00304419` and its L-infinity difference
decreased `0.00422934 -> 0.00416772`; this per-field gate prevents the aggregate
external drop from hiding a nonconvergent interface quantity. Interface-flow
differences remained at roundoff.

The coupled external-drop first-harmonic amplitudes were
`0.0202688, 0.0195659, 0.0191344 Pa`, so consecutive differences decreased
`0.000702961 -> 0.000431453 Pa`. Endpoint phases referred to the common period
origin were `-1.08715, -1.05078, -1.03304 rad`, with consecutive differences
decreasing `0.0363708 -> 0.0177448 rad`. The independently advanced all-1D
reference had amplitude `0.0157439 Pa`; at the fine grid the coupled/reference normalized
pressure-drop L2 difference was `0.0884427` and phase difference was
`0.537761 rad`, while endpoint flow agreed to `4.27e-14` relative L-infinity.
Those nonzero pressure and phase differences are consistent with the fixed
square-duct discretization and its resolved inertial response; they are
reported, not interpreted as temporal error or hidden by calibration. Because
the current rigid `steady_poiseuille` 1D scheme is algebraic, however, this
reference has no fluid inertia, compliance, or pulse propagation. It therefore
cannot by itself validate physical pressure-wave phase or pulse transit, and
Phase 1 remains open pending a transient reference or manufactured comparison
that supplies those dynamics.

At the fine timestep, waveform distances to the `1e-10` pressure-tolerance
reference decreased from `7.37420e-5` at `1e-4` to `1.01375e-7` at `1e-6` and
zero by definition for the tight reference. Every run also passed time alignment,
finite P/Q, conservation, wall flow, residual, iteration/Aitken, subcycling
work, and artifact/manifest audits. The maximum cycle-3-to-cycle-4 pressure
relative L-infinity difference was `1.45499e-11`; the maximum flow relative
L-infinity difference was `2.40693e-14`. The corrected production target
exited zero in `10:19.12` with `44884 kB` peak RSS.

## Straight-vessel steady spatial verification

The production opt-in target `make coupling-convergence-test` constructs its
`.ntiga`, VTK control mesh, velocity profile, and 1D/3D cases at run time. The
3D section is a unit-square straight duct with an independent analytic
Poiseuille coefficient `28.45415376956191 Pa s/m3` for unit length and unit
viscosity. The all-1D middle section uses the area-matched circular radius and
a hydraulic-equivalent length, giving the independently evaluated total drop
`0.07871963622699861 Pa` at `Q=1e-3 m3/s`; the native 1D runtime reproduced it
to roundoff.

The fixture uses open-uniform cubic C2 tensor-product splines. Dependency-free
audits cover direct/extracted basis agreement, partition and linear precision,
C2 derivative continuity, nonsymmetric tensor orientation, sparse serialized
extraction, affine Bézier coordinates, positive Jacobians, boundary labels,
wall traces, shared connectivity, one/two-rank ownership, profile symmetry and
walls, L2 convergence, energy resistance, and a discrete Ritz weak residual.
`-UNDEBUG` is explicit on the fixture target so these assertions remain active
under caller-supplied release flags.

Each refinement `n=2,4,8` solves lengths `L=1,1.5,2 m` with axial counts
`nx=n*L`. Consecutive half-length slopes make the assumed length-independent
end correction testable. The worst combined 3D/external slope disagreement
decreased `5.8373% -> 0.21937% -> 0.028604%`; the outlet static/traction gap
range decreased `3.9580e-4 -> 3.0233e-5 -> 3.2919e-7 Pa`. End-cancelled 3D
bulk coefficients were `20.7369671`, `27.3152612`, and `28.3726195 Pa s/m3`,
with relative errors `27.1215%`, `4.00255%`, and `0.286546%` and observed
orders `2.76044`, `3.80408`. External-path bulk errors were `25.7305%`,
`3.89630%`, and `0.285389%`, with orders `2.72330`, `3.77110`. The fine
Richardson order was `2.63725`, extrapolated coefficient `28.5751226`, and
fine GCI `0.892159%` (gate: at most 1%).

All nine cases had zero final two-step resistance change, seven total Aitken
sweeps across three physical steps, roundoff interface and 3D mass residuals,
and zero wall flow. Fine direct all-1D pressure-drop errors were
`1.7129%`, `1.4965%`, and `1.3339%` for increasing length (gate: 2%). The
upstream fixed-point pressure jump was zero; the fine outlet static/traction
jump was about `-3.42e-4 Pa` and remained within 0.5% of the absolute reference
drop. The final production run passed in `54:54.62`, used `815288 kB` peak RSS,
and exited zero on the local PETSc 3.15 installation.

An earlier raw cap-to-cap ladder was rejected because its length-independent
end pressure offset produced only first-order-looking convergence. A numerical
review also rejected accepting two-length cancellation alone: the production
gate now uses three lengths, enforces converging/fine length linearity, retains
absolute external-drop and interface-jump bounds, and treats neither static
cap pressure nor the natural traction parameter as interchangeable.

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
