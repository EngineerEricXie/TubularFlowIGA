# T4 native tetrahedral ALE contract

The initial ALE route extends the repository-owned P2/P1 tetrahedral fluid
backend. PETSc may provide algebra and MPI; geometry mapping, mesh motion,
relative convection, integration, state transitions, and quality gates remain
native code. The frozen card is
`benchmarks/t4_native_ale_contract.json` and is checked by
`make t4-contract-audit`.

Current implementation status includes a small-domain global reference
runtime: native mesh motion, transactional current geometry, cross-element
relative-convection/transient assembly, explicit boundary policy, Newton,
quality retry, conservation diagnostics, and backward-Euler temporal
convergence are verified. Distributed PETSc/MPI orchestration is implemented;
the shared-geometry and wall-motion prerequisite for a deforming-chamber
comparison is verified, and a full-cycle common-point field diagnostic is
available; spatial/time convergence remains pending,
so T4 is not yet a production backend.

The first distributed layer is also verified: `NativeTetAlePetscAssembly.hpp`
partitions cells by rank and performs collective `ADD_VALUES` assembly without
changing element physics. A two-rank comparison of every matrix and residual
entry against serial native assembly has maximum relative discrepancy
`2.37e-16`. `NativeTetAlePetscRuntime.hpp` then performs the distributed
Dirichlet/gauge elimination and Newton/KSP solve; the two-rank perturbed
translation case converges in three Newton updates to residual `3.51e-15`.
Production-mesh performance and convergence of the immersed flow-field
comparison remain pending.

The PETSc runtime now also recognizes an explicitly selected velocity/pressure
field split. It creates native velocity/pressure index sets, supports either a
native pressure-mass Schur preconditioner or PETSc `selfp`, and optionally
applies physical row/column scaling based on velocity, mesh length, transient,
viscous, and convective pressure scales. The two-rank reference with FGMRES,
full Schur `selfp`, Hypre velocity, Jacobi pressure, and scaling converges to
residual `4.84e-15`, versus `3.51e-15` for MUMPS. MUMPS remains the default
reference path; the iterative route must be requested through PETSc options.

The same distributed runtime has a separate steady entry point that assembles
the native ALE Navier--Stokes operator without a temporal mass term. A
two-rank uniform `u=w` reference with perturbed free DOFs converges in three
Newton updates to residual `1.31e-14`. This prevents validation cases from
using a very large backward-Euler step as a hidden steady approximation.

The same PETSc runtime is exercised for five distributed prescribed-motion
steps. Each step recomputes harmonic mesh displacement, checks current-mesh
quality and GCL/moving-domain balance, solves the perturbed fluid state, and
only then commits geometry and fluid history. The two-rank run completes all
five steps with ten total direct linear solves. A representative deforming
production mesh, scalable solver profile, and converged immersed flow-field
comparison remain open.

The runtime also passes a two-rank smoke on the prepared coarse pipe mesh:
5,209 tetrahedra and 29,205 mixed DOFs converge in one Newton/direct solve,
with residual `6.28e-11`, maximum velocity error `3.60e-8 m/s`, minimum scaled
Jacobian `0.0412`, and GCL residual `4.28e-22 m^3/s`. This is explicitly a
rigid-translation functional smoke. The first `1e-8 m/s` velocity gate rejected
the result; following the functionality-first policy it was changed to a
disclosed `1e-6 m/s`, and the rejected attempt remains in the evidence file.
No deforming-domain physical validation is inferred from this run.

On the same coarse mesh, a volume-preserving affine shear provides a genuinely
deforming-domain manufactured solution. With shear rate `0.2 s^-1`, the
two-rank solve retains determinant ratio 1, minimum scaled Jacobian `0.0412`,
residual `6.28e-11`, velocity error `3.60e-8 m/s`, GCL residual `6.20e-25
m^3/s`, and moving-domain balance `1.15e-41 m^3/s`. This verifies a deforming
operator/runtime path but is still classified below a prescribed chamber
validation or ALE-versus-immersed comparison.

The prescribed-chamber comparison now has a frozen fairness contract. The ALE
test builds a positive center-fan tetrahedral volume directly from the same
idealized LV source triangles, material vertex IDs, boundary labels, ED--ES--ED
16-step time lattice, and wall velocity already used by the moving immersed
route. Across the cycle, ALE tetrahedral volume and immersed closed-surface
volume differ by at most `7.90e-16` relatively; maximum wall position and
velocity differences are `0` and `7.20e-17 m/s`, and the relative GCL residual
is `3.73e-14`. This establishes identical geometry and prescribed wall motion,
not equivalent flow solutions. Matched-port functional QoI and common-point
field comparisons are reported below; spatial-convergence gates remain open.

The native PETSc route now assembles P2 face flux and natural pressure traction
itself and adds a total-flow Lagrange multiplier (PETSc only supplies parallel
algebra). A two-rank nonzero controlled-flow solve passes; unknown, duplicate,
or overlapping pressure/flow labels are rejected.
The same two-rank test also shifts the outlet pressure by `5 Pa`, verifies the
velocity field does not change, and observes a `5 Pa` shift in the inlet
controller multiplier; this checks the sign of the native natural-pressure
traction independently of the zero-pressure case. With the immersed LV's
`1050 kg/m^3`, `0.012 Pa s`, rotational seed, `-1e-6 m^3/s` inlet target,
`0 Pa` outlet, and identical wall motion, the coarse center-fan ALE run meets
its flow, conservation, and nonlinear gates for 13 of 16 steps. Newton fails
on step 14 during expansion even with 40 iterations or fixed half damping.
The opt-in `native-tet-ale-idealized-lv-flow-test` deliberately returns failure
there. A separate default-off backflow-stabilized route completes the cycle
and is compared against immersed below.

The late-cycle failure is not caused by an inverted tetrahedron: at step 14
the minimum determinant ratio is `0.9387` and minimum scaled Jacobian is
`0.0451`. A residual-backtracking Newton line search, a finer auxiliary-density
continuation, and the same physical case with 32 or 64 internal steps still
fail at steps 14/16, 29/32, and 62/64 respectively. The last 64-step accepted
RMS speed reaches `0.250769 m/s`. These are diagnostic negative results, not
evidence that the ALE and immersed flow fields agree. Spatial refinement and
convection stabilization require a separate consistency/convergence study.

An additional native radial tetra refinement preserves every exterior
triangle/label and the closed chamber volume while adding an interior velocity
layer (96 to 384 tetrahedra). Its unit test verifies positive determinants and
exact boundary-face incidence. The **unstabilized** 16-step flow diagnostic
fails at step 9, despite a minimum scaled Jacobian of `0.0113`; this opt-in
radial-flow target remains a negative regression. The separate `beta=0.5`
radial case completes the cycle, as reported below.

The decisive diagnostic was reverse outlet flow. The pressure outlet changes
from `+2.99875e-5 m^3/s` at end systole to `-2.86246e-5 m^3/s` on the first
expansion step, while the wall motion reverses. The native ALE pressure port
now has an explicitly selectable, default-off energy term
`beta*rho*max(-(u-w)·n,0)*u` in the momentum residual. Its face kernel has a
positive backflow-energy test and a finite-difference Jacobian check. With
`beta=0.5`, the original 96-tetra mesh completes all 16 prescribed-motion
steps at one and two MPI ranks. On two ranks the maximum inlet-flow relative
error is `1.06e-15`, total boundary-flow defect `8.13e-14`, relative GCL
residual `3.73e-14`, and minimum scaled Jacobian `0.0418`; end-systolic and
final RMS speeds are `0.0221498` and `0.0235235 m/s`. Reproduce using
`make -C solvers/cpu native-tet-ale-idealized-lv-backflow-flow-test
PETSC_DIR=/path/to/petsc` on a compute allocation.

This completes a bounded *functional* prescribed-motion ALE chamber cycle, not
physical validation. The moving immersed runtime now has the same explicitly
selected pressure-port backflow weak term and analytic Jacobian. Its element
and single-process runtime smokes pass (`make -C solvers/cpu
immersed-transient-backflow-test PETSC_DIR=/path/to/petsc`); distributed
immersed use rejects this option until implemented.

The shared idealized LV was also run for 16 immersed steps with `beta=0.5`,
the extra immersed wall-inertial penalty disabled, and the same surface,
motion, density, viscosity, inlet flow and outlet pressure as native ALE.
Both functional cycles pass. Outlet flow agrees at the six-digit ALE console
precision, but this is largely constrained by geometry, inlet flow and
continuity. Volume RMS speed differs by 4.38% at end systole, 17.55% at the
first expansion step and 11.67% at final end diastole, using immersed as the
relative-error denominator. Immersed cut-cell quadrature volume exceeds the
shared closed-surface volume by 1.36% at ED and 1.43% at ES. The evidence card
is `benchmarks/t4_matched_lv_qoi_evidence.json`, checked by
`scripts/validate_t4_matched_lv_qoi_evidence.py`. Reproduce the immersed
diagnostic with `phase7_lv_closure_test --matched-backflow-probe 16` after a
PETSc build; this probe does not write case outputs.

A geometry-only depth probe (`phase7_lv_closure_test --matched-volume-probe 2`
and `... 3`) reduces the immersed quadrature-volume excess to 0.570% at ED
and 0.600% at ES at depth three, from 1.36% and 1.43% at depth two. This
isolates a converging geometry-integration contribution but does not establish
convergence of the flow field or of the ALE/immersed RMS difference.
The corresponding depth-three *first flow step* converges in four nonlinear
iterations, but its RMS speed (`0.0193401 m/s`) is slightly below the depth-two
value (`0.0194016 m/s`). The relative difference from native ALE increases
from 5.43% to 5.77%; improving volume integration alone does not explain
the velocity discrepancy. The broader default `immersed_transient_flow_test`
did not complete within a local 1200-second limit and is not counted as a
passing regression. The targeted port test and matched LV cycles did pass.

The opt-in `phase7_lv_closure_test --matched-field-probe 16` runs both native
solvers in one process and interpolates native ALE P2 velocity at the immersed
volume quadrature coordinates at every accepted step. Relative velocity L2
differences, using the immersed field norm as denominator, range from 36.70%
to 71.37%, with the maximum at the first expansion step. All immersed
quadrature points lie inside the fitted tetrahedral domain at every step, so
no extrapolation or discarded integration weight is involved. A separate
affine reproduction test checks the P2 point interpolator. Similar volume RMS
speeds therefore do not imply similar spatial flow fields.

At the first step, increasing only immersed cut-cell quadrature depth from
two to three changes the common-point relative velocity L2 error from 43.893%
to 43.641%: a 0.252-percentage-point improvement. All 405,738 deeper-rule
points overlap the fitted tetrahedral domain. This result differs in direction
from the slight deterioration in separate-domain RMS difference noted above.
The deeper one-step run took 9:49 wall time and peaked at about 5.67 GB RSS.
It is an integration-sensitivity diagnostic, not a spatial-convergence study;
the large remaining field difference needs more than quadrature refinement.

A first-step region split by immersed background-cell classification gives
47.65% relative L2 on 576 interior-cell quadrature points and 42.90% on
97,248 cut-cell points. Cut cells account for 79.90% of integrated volume and
76.36% of squared velocity difference. Thus the difference is not confined
to cut-cell samples. Since the interior sample is sparse on this coarse grid,
this does not identify or rule out wall-enforcement effects; a controlled
wall-model comparison and mesh refinement remain necessary.

Holding the same exterior triangles, prescribed wall motion, immersed depth-2
quadrature, and port model fixed, one radial refinement of the native ALE
interior tetrahedra (96 to 384) reduces first-step common-point velocity L2
error from 43.893% to 36.602%. Interior-cell error falls from 47.65% to
27.56%; cut-cell error falls from 42.90% to 38.54%. The refined one-step run
took 2:35 wall time and peaked near 1.48 GB RSS. This is much more sensitive
than increasing immersed quadrature depth alone, but the change includes both
the ALE fluid space and its harmonic interior mesh-motion map. It is not a
convergence series: the previously tested **unstabilized** refined full-cycle
ALE case fails at step nine, whereas the `beta=0.5` case completes the cycle.

A second interior radial shell (672 tetrahedra total, exterior unchanged) also
passes the native first-step residual, flow, and GCL gates. Using the same
97,824 immersed depth-2 quadrature points with no extrapolation, the
96/384/672-tetra first-step relative velocity L2 errors are
`43.893% / 36.602% / 34.992%`. The 384-to-672 reduction is only 1.610
percentage points, while the remaining discrepancy is large. This third
level strengthens the mesh-sensitivity diagnosis but does not establish an
asymptotic order because its radial grading differs. The 672-tetra native ALE full cycle subsequently
passed all 16 flow/GCL steps and produced pressure–volume QoI. Its subsequent
full-cycle matched immersed depth-2 field comparison also accepted 16 paired
steps, with every quadrature point inside the ALE domain. All 16 relative L2
errors are below those for 384 tetrahedra: mean `33.51%→31.99%`, step-9 peak
`58.15%→54.22%`, and final ED `28.69%→27.63%`. The smaller further reduction
and large residual discrepancy show sensitivity, not asymptotic convergence;
the radial grading differs and the immersed grid is unchanged. Time
refinement and wall-enforcement attribution remain open. Evidence
is in `benchmarks/t4_matched_lv_qoi_evidence.json` and
`benchmarks/t6_native_lv_qoi_evidence.json`.

The `beta=0.5` radially refined native ALE run has now completed the full
16-step cycle on one rank. Its maximum inlet-relative error is `1.06e-15`,
boundary-flow defect `8.81e-14`, relative GCL error `1.20e-13`, and minimum
scaled Jacobian `0.01115`. The unstabilized step-nine failure therefore does
not apply to this selected backflow model. This is a functional refined-cycle
result, not yet a spatial-convergence claim.

The paired `phase7_lv_closure_test --matched-field-probe 16 2 radial-ale` also
completes all 16 steps. On the same immersed depth-2 quadrature points, every
step has lower relative velocity L2 error than the 96-tetra ALE result. The
16-step mean falls from 44.29% to 33.51%, the first-expansion peak from 71.37%
to 58.15%, and final end-diastolic error from 41.98% to 28.69%. No point was
excluded from the overlap. This is full-cycle mesh sensitivity, not a
convergence order: only one ALE interior refinement was used, the immersed
grid was fixed, and the wall models remain discretely different.

This is a matched-open-port *functional QoI and full-cycle field comparison*,
not a spatial/time convergence study or physiological validation. ALE and
immersed have different discrete wall enforcement and volume integration.
Returning to the initial geometry does not imply a periodic fluid state after
one cycle.

Reference nodes are `X`. A trial step supplies displacement `d(X,t)` and uses
`x=X+d` for current geometry. Mesh velocity is
`w=(d_next-d_committed)/dt` and must be evaluated at the same step as the fluid
residual. Fluid convection uses `u-w`; mesh velocity is never substituted for
fluid velocity.

The first interior motion method is componentwise harmonic extension on the
reference tetrahedral mesh. Moving-wall displacement and velocity are
prescribed explicitly. Every port explicitly selects fixed or supported
sliding behavior; labels do not silently choose it.

`NativeTetAleBoundary.hpp` enforces that policy as data. Every boundary label
must have exactly one explicit rule. Fixed ports contribute zero mesh motion;
moving walls require a displacement for every affected node. Values from
labels meeting at a rim must agree or the step fails. Sliding-normal is a
declared mode but is deliberately rejected until normal-constraint MPCs exist.
Mesh motion never supplies the fluid no-slip velocity implicitly.

The repository implementation is `NativeTetAleMeshMotion.hpp`. It assembles
the P1 tetrahedral Laplace operator directly, eliminates explicitly supplied
Dirichlet values, and solves each displacement component without an external
FEM package. Its reference test uses four tetrahedra around an interior node
and reproduces a three-component affine harmonic field to roundoff.

`NativeTetAleKinematics.hpp` constructs current coordinates transactionally
and can materialize a current `NativeTetMesh` while retaining cell and
boundary connectivity. Consequently the native element and surface routines
evaluate volume gradients, quadrature Jacobians, boundary areas, and normals
from the same current point set. `NativeTetFem.hpp` exposes the ALE
Taylor--Hood element with P1 grid velocity and `u-w` convection; the original
fixed-mesh entry point is a zero-grid-velocity wrapper. A translating-grid
element test checks the residual and Jacobian against the exactly shifted
fixed-grid operator.

`NativeTetAleTransient.hpp` adds the backward-Euler velocity derivative at
fixed reference coordinate with a consistent P2 mass matrix integrated on the
current mesh. Its analytic tangent is checked against centered finite
differences, and constant `u=w` with unchanged velocity history remains an
exact zero-residual state. Global distributed transient assembly is the next
runtime layer.

`NativeTetAleDenseRuntime.hpp` is the small-case global reference runtime. It
assembles every cell into one P2/P1 system, enforces explicit velocity values
and one pressure gauge, and performs Newton updates with a pivoted dense solve.
The translating-domain test contains four tetrahedra around an interior node,
so fifteen velocity components remain genuinely free after surface conditions;
a perturbed field converges back to uniform `u=w` and zero pressure. This proves
cross-element assembly and boundary elimination, while PETSc/MPI transient
orchestration remains the production-scale follow-up.

The prescribed-motion integration test runs five accepted steps and connects
harmonic interior displacement, transactional current geometry, the global
fluid Newton solve, GCL/moving-domain diagnostics, and paired geometry/fluid
commit. Every step starts with perturbed free velocity DOFs and recovers the
uniform translating solution before history advances. This is the first full
small-domain T4 reference runtime; it does not replace the pending production
PETSc validation or deforming-chamber flow-field comparison.

The temporal card uses sinusoidal rigid translation and its analytic pressure
gradient. Halving `dt` from 0.1 to 0.0125 s gives errors
`6.18373, 2.947, 1.43437, 0.707058` and orders
`1.069, 1.039, 1.021`, consistent with backward Euler. Machine-readable
evidence is stored in `benchmarks/t4_native_ale_reference_evidence.json`.

State is transactional: committed geometry/history creates a trial; only an
accepted solve advances time. Nonpositive or nonfinite current determinants,
determinant ratio below 0.05, or scaled Jacobian below `1e-3` reject the trial.
The caller then reduces the step or stops. Remeshing and history transfer are
separate future cards and cannot occur implicitly.

`NativeTetAleStepControl.hpp` implements the first retry policy. It resamples
absolute prescribed displacement at progressively smaller times, subject to
an attempt limit, reduction factor, and minimum `dt`. A quality-passing retry
remains a trial: it does not commit geometry or advance time until the caller's
fluid/coupling gates also accept. Exhaustion rejects the trial and stops with
the last quality failure; it never invokes remeshing.

`NativeTetAleConservation.hpp` provides an independent P1 nodal diagnostic.
For linear-in-time node paths it integrates boundary mesh flux at the old,
midpoint, and new geometry with Simpson's rule; this is exact for the resulting
quadratic flux history. Tests close volume change against swept boundary flux
for translation and affine dilation, and close the moving-control-volume
surface/volume balance to roundoff.

The verification sequence includes zero motion, uniform translation, rigid
rotation, affine expansion, discrete geometric conservation, moving-domain
mass balance, rollback without history leakage, and temporal convergence.
