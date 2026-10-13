# T2 fixed-wall flow validation contract

The first T2 card is frozen in
`benchmarks/t2_fixed_flow_contract.json`. Its physics and QoI schema remain
backend-neutral, but the implementation decision is deliberately native:
the same circular-pipe geometry, SI material values, boundary signs, analytic
reference, mesh sequence, QoI tolerances, failure classes, result fields, and
resource budget apply to the native FEM implementation.

The steady solution is
`u_x(r) = 2 Q/(pi R^2) (1-r^2/R^2)` and
`Delta p = 8 mu L Q/(pi R^4)`.

The inlet receives that analytic profile, the wall is no-slip, and the outlet
uses zero gauge pressure with the backend's consistent natural traction. Flow
is outward-positive, so the inlet value is `-Q` and the outlet value is `+Q`.
This unidirectional solution has zero convective acceleration and is valid for
the declared steady incompressible Navier--Stokes equation as well as Stokes.

Three successively refined meshes are mandatory. Conservation and backflow
gates apply to every level; the absolute velocity, pressure-drop, and WSS gates
apply to the finest level, while all three positive errors determine the
observed convergence order. Geometry approximation error
must be reported separately from the volume velocity norm and may not be hidden
inside a combined error. A backend must name its velocity and pressure spaces
and either use a stable pair or document and test its equal-order stabilization.
Timing begins only after the correctness gates pass; assembly and solve time
remain separate.

For this verification card, outlet backflow beyond the frozen area-fraction
threshold is a failed case, not something silently stabilized away. A later
production boundary may add an explicitly configured energy-stable backflow
treatment, but that is a distinct validation card.

Run `make t2-contract-audit` to recalculate every analytic value and reject
changes to signs, labels, refinement structure, required evidence, or the
fail-closed backend decision. Passing this audit proves the test definition,
not a solver result.

The native solver must write the three-level result schema enforced by
`scripts/validate_t2_fixed_flow_result.py`. That validator independently
recomputes mass imbalance, pressure error, and observed velocity convergence
order; binds results to the exact contract bytes; requires one unchanged
backend/discrete formulation across all levels; checks provenance hashes,
convergence reasons, geometry/label evidence, QoI gates, and separate assembly,
solve, and memory measurements. Its healthy data are synthetic test fixtures
only and are never recorded as physical solver evidence.

A Git commit alone is insufficient provenance for a dirty worktree. Each level
therefore records `source_tree_state` and a SHA-256 over the exact source-file
manifest in addition to the full Git object ID.

There are deliberately two execution classes. `functional` uses a looser
nonlinear stopping tolerance to exercise mesh loading, native assembly,
boundary conditions, Newton/KSP, QoIs, and provenance on local hardware. It
writes `result_classification=functional_smoke` and
`validation_gates_enforced=false`. `validation` uses strict tolerances and
per-level physical gates. Only three strict level results collected as
`result_classification=physical_validation` are accepted by the validator;
functional output cannot be promoted merely because its QoIs happen to fit a
threshold.

`scripts/generate_circular_pipe_surface.py` creates the closed labelled VTP
input without storing generated case data in Git. Wall, inlet, and outlet use
labels 0, 1, and 2. Each mesh level must refine both the axial/circumferential
surface and the tetrahedral volume; the polygonal-cylinder volume error and the
circle sagitta are analytic geometry-error evidence. The regression sends this
surface through the same C++ manifold/orientation/self-intersection preflight
used by the immersed and FEM routes.

After building `surface_fem_preflight`, prepare the complete untracked mesh
series with:

```bash
python3 scripts/prepare_t2_poiseuille_meshes.py /tmp/t2-poiseuille-meshes
```

`mesh-series.json` records both the exact contract SHA-256 and a canonical hash
of only the mesh-affecting geometry/label/target-size subset. A later result
schema or solver-policy edit therefore cannot masquerade as a geometry change,
while a real mesh-contract change invalidates reuse. The index also records
surface/mesh hashes, segment counts,
tetrahedron quality, boundary preservation, circle sagitta, and polygonal-volume
error. The circumferential edge target is one half of the volume target so the
polygonal circle does not dominate the PDE series. Its `kind` explicitly
prevents treating mesh preparation as solver
validation evidence.

The local WSL gate on 2026-09-18 used Gmsh 4.8.4, contract SHA-256
recorded by each generated `mesh-series.json`,
and produced:

| target size (m) | circumference | tetrahedra | min scaled Jacobian | max circle distance (m) | relative volume error |
|---:|---:|---:|---:|---:|---:|
| 0.0025 | 26 | 5,209 | 0.0411828 | 3.64556e-5 | 9.70496e-3 |
| 0.00125 | 51 | 31,721 | 0.00719376 | 9.48336e-6 | 2.52778e-3 |
| 0.000625 | 101 | 211,177 | 0.00227314 | 2.41859e-6 | 6.44884e-4 |

All levels retained exactly labels 0/1/2, matched the tetrahedral boundary,
had positive determinants, and passed the frozen `1e-3` scaled-Jacobian floor.
The declining minimum quality is caused by the current deterministic polar cap
triangulation and is recorded rather than hidden; replace that topology before
production if it impairs the solver convergence gate.

The initial result schema incorrectly applied final accuracy to every coarse
convergence level. It was corrected to the standard finest-level absolute gate
without loosening any tolerance. No external FEM solve is accepted as native
solver evidence; the full three-level series must come from the repository's
own discretization and assembly implementation.

## Native FEM implementation status

The selected implementation is `native_cpp_petsc_tetrahedral_fem`. The project
owns the ASCII Gmsh 4.1 reader and boundary validation, P2 edge/P1 vertex DOF
topology, Taylor--Hood basis, affine tetrahedron mapping, degree-five-capable
quadrature, steady Navier--Stokes residual and consistent Jacobian, strong
velocity conditions, Newton assembly, and QoI integration. PETSc supplies
distributed matrices/vectors, Krylov/factor solvers, and MPI only. The element
Jacobian is checked against centered finite differences by
`native_tet_fem_test`.

The 2026-09-18 native coarse run contained 1,461 vertices, 7,787 unique edges,
5,209 tetrahedra, and 29,205 mixed DOFs. Its Newton residual decreased from
`2.04345e-7` to `1.53285e-9` and `7.20237e-14` in two updates. The native QoI
integrator reported inlet/outlet flows `-9.95091e-7/+9.95091e-7 m^3/s`, pressure
drop `1.4466 Pa`, velocity relative L2 `1.00774%`, WSS relative L2 `3.50116%`,
mass imbalance `6.35e-16`, divergence closure `4.21e-16`, and zero outlet
backflow. The MUMPS direct path took 221.6 s and is evidence only for the coarse
correctness gate. Initial fully iterative Schur/Hypre experiments did not
converge within their iteration limit. The functional path below uses a Schur
field split with a direct velocity subsolve. Its PETSc field index sets,
boundary rows, and pressure-preconditioner ownership follow the distributed
matrix layout; each result records both requested and observed MPI rank counts
and the maximum per-rank RSS. This is a working distributed path, but the
direct block factorization is not claimed to be a scalable solver.

The local functional pipeline can be exercised after mesh preparation with:

```bash
make t2-native-functional-coarse \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real \
  T2_MESH_SERIES=/tmp/t2-poiseuille-meshes/mesh-series.json \
  T2_FUNCTIONAL_RESULT=/tmp/native-t2-functional-level-1.json \
  T2_MPI_RANKS=2
```

Any level can instead be selected with `make t2-native-functional-level`,
`T2_FUNCTIONAL_LEVEL=1|2|3`, and `T2_MPI_RANKS=N`. The Python runner launches
MPI without a shell, verifies that the solver-reported rank count matches the
request, and only promotes the result after successful exit and provenance
checks. A timeout reports an explicit failure and leaves no output at the
requested result path.

The checked-in functional observations are in
`benchmarks/t2_native_functional_evidence.json` and are fail-closed by the T2
audit:

| level | tets | mixed DOFs | ranks | assembly (s) | linear solve (s) | total (s) | velocity L2 | WSS L2 | mass imbalance |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 5,209 | 29,205 | 2 | 6.31 | 2.95 | 9.26 | 1.0075% | 3.5010% | 2.85e-8 |
| 2 | 31,721 | 160,607 | 4 | 148.54 | 36.41 | 184.96 | 0.4153% | 1.9712% | 2.01e-8 |

Both completed one Newton update. The velocity and WSS errors decrease on the
medium mesh, and its pressure drop is `1.42565 Pa`. These observations remain
`functional_smoke`: their loose solve policy does not enforce strict physical
gates, and the port imbalance remains above the frozen `1e-8` validation gate.

The 211,177-tet fine level was attempted with both eight physical-core ranks
and 16 Open MPI hardware-thread ranks; both reached the 600 s local wall-time
limit. The runner promoted no result and no solver
process remained. A single Hypre AMG V-cycle failed to converge in 600 outer
iterations on the coarse level; inner GMRES+AMG converged but took 128.6 s,
versus the roughly 9--10 s retained direct-block coarse observations. Therefore the fine level and
the strict three-level validation remain open pending a stronger scalable
block preconditioner; the timeout is not presented as partial success.

## Native Y-bifurcation functional case

The first native bifurcation smoke uses the labelled Y fixture, the exact
boundary-preserving Gmsh adapter, and the same project-owned P2/P1 element
assembly. Its mesh has 501 vertices, 1,332 tetrahedra, 8,886 mixed DOFs,
minimum scaled Jacobian `0.00631`, and relative volume error `3.35e-16`.
Label 0 is wall, label 1 inlet, and labels 2/3 are the symmetric outlets.

A two-rank, five-step low-Reynolds run reports inlet flow
`-3.67717e-7 m^3/s`, outlet flows `1.83392e-7` and `1.84325e-7 m^3/s`, upper
flow fraction `0.498732`, and relative mass imbalance `1.44e-16`. Average port
pressures are `0.64171`, `0.004363`, and `0.004894 Pa`. The first formulation
also deleted one pressure/continuity row and produced `1.54e-3` mass error; it
was rejected. With natural outlets, all continuity rows are retained and no
extra pressure gauge is imposed.

This result is explicitly `functional_smoke_not_mesh_convergence`. It proves
native labelled bifurcation flow, pressure reporting, symmetric flow division,
and roundoff mass closure. Multiple mesh/time levels plus traction and WSS
convergence remain required before the T2 bifurcation checklist item can be
closed. Machine-readable evidence is in
`benchmarks/t2_native_bifurcation_evidence.json`.

Current-source caveat (2026-09-20): the Y runtime regression now imposes
`0 Pa` on both outlet labels for the solved vessel→tissue numerical smoke.
The archived bifurcation QoI in that evidence card came from the earlier
natural-outlet branch. They remain historical measurements, **not** a fresh
current-source convergence study or directly comparable outlet-pressure QoI.
The evidence validator checks this classification and the current source
hash; a new matched refinement study is needed before claiming current-source
spatial/temporal convergence.

The same canonical surface was then held fixed while the volume target size
was reduced from `0.004` to `0.002` and `0.0015 m`, producing 1,332, 3,599,
and 8,994 tetrahedra for the original five-step study. Its medium-to-fine
upper-flow-fraction change is `7.71e-5`, below its `1e-3` functional gate, but
inlet pressure, mean/RMS WSS, and integrated wall-traction norm change by
`9.2%`, `11.9/12.1%`, and `12.2%`.

A second fixed-time (`1 x 0.05 s`) spatial study adds a fourth exact-boundary
volume with target `0.001 m`: 6,793 vertices, 39,157 tetrahedra, and 166,405
mixed DOFs. The native P2 field is integrated over every cell using the same
degree-five tetrahedron quadrature, rather than treating outlet flow as a
surrogate for velocity. Fine-to-finest relative changes are `0.0365%` for
flow fraction, `0.993%` for inlet pressure, `0.00894%` for volume-mean speed,
`0.860%` for volume-RMS speed, and `1.73%` for kinetic energy. Those integral
flow, pressure, and velocity QoIs pass their functional spatial gates.
Mean/RMS WSS still change by `22.2/21.6%`, and the wall integrated-traction
norm changes by `72.0%`; therefore traction/WSS remain explicitly
nonconverged and the T2 bifurcation checklist item stays open. The fourth run
completed on four local MPI ranks with residual `1.63e-10` and roundoff mass
imbalance. These failures are retained rather than hidden by a looser solve
tolerance.

The integrated wall-traction vector above is cancellation-sensitive on a
curved, nearly symmetric wall. It remains recorded, but is not an appropriate
standalone relative convergence scalar. The native QoI implementation now
also reports signed/RMS normal traction, mean/RMS total traction, integrated
tangential traction, and raw mean/RMS WSS. Uniform pressure on a tetrahedron
verifies the traction sign, magnitude, zero tangential component, and optional
surface-region filtering analytically.

`scripts/prepare_t2_bifurcation_meshes.py` provides a reproducible second
series that refines both the triangulated OCC surface and tetrahedral volume.
Its first three completed levels contain 661, 3,599, and 36,517 tetrahedra;
surface triangle counts increase from 490 to 924 and 3,144. Medium-to-fine
changes are `0.0868%` flow fraction, `1.53%` inlet pressure, `0.360%` volume
RMS speed, `3.51%` kinetic energy, and `1.12%` wall RMS total traction. A fixed
smooth trunk-wall diagnostic (`label=0`, `0.01 <= x <= 0.025 m`) gives `1.36%`
RMS total-traction change. These integral quantities pass their functional
5% gates (flow fraction uses 0.1%).

Raw WSS does not yet pass: global RMS WSS changes `24.0%`, while the same
smooth trunk-wall region changes `22.3%`. This demonstrates that neither the
Boolean-union junction singularity nor traction-vector cancellation alone
explains the remaining near-wall gradient error. A fourth coupled mesh with
5,124 surface triangles, 92,179 tetrahedra, and minimum scaled Jacobian
`0.00842` passed geometry gates but timed out after 300 s on eight local MPI
ranks; no result was promoted. The remaining work is a scalable linear solver
and/or a verified near-wall gradient-recovery route, not a looser nonlinear
tolerance.

A native patch-recovery experiment volume-weights the P2 velocity gradient at
every shared vertex/edge node over incident tetrahedra, then interpolates the
recovered tensor to boundary quadrature points. It exactly reproduces an
affine vector field. On this Y series, however, medium-to-fine recovered RMS
WSS changes are `26.2%` globally and `24.7%` on the fixed smooth trunk region,
worse than the raw-gradient changes. Recovered WSS is therefore retained as a
negative result and is not promoted to the acceptance metric.

The project now has a distinct native steady PETSc solve entry point; it uses
the same autonomous P2/P1 topology, quadrature, residual/Jacobian, distributed
assembly, constraints, and PETSc algebra, but omits the backward-Euler mass
term rather than approximating steady flow with an artificially large time
step. On the fixed 924-triangle Y surface, 1,332/3,599/8,994-tetrahedron
steady solves all close mass to roundoff. Medium-to-fine changes are `0.0443%`
flow fraction, `0.225%` outlet-referenced pressure drop, `0.0410%` volume RMS
speed, `0.0819%` kinetic energy, `0.724%` wall RMS total traction, and `0.126%`
wall RMS WSS. The fixed smooth trunk region also passes traction and WSS gates.
This completes the idealized bifurcation numerical-convergence item; it is not
a physiological validation claim. The unresolved startup-boundary-layer study
remains separately classified as nonconverged.

The ALE transient runtime now exposes a native PETSc field-split path rather
than requiring monolithic MUMPS. On the 3,599-tetrahedron level, scaled FGMRES
with full Schur `selfp`, Hypre velocity, and Jacobi pressure completes two
Newton updates in 493 total linear iterations. Its mass imbalance is
`6.79e-11`, and flow fraction, pressure, RMS velocity, RMS total traction, and
RMS WSS differ from the direct reference by at most `1.24e-10` relative. The
native pressure-mass Schur attempt reached `DIVERGED_ITS` at 400 iterations.
Both scaled and unscaled `selfp` attempts on 36,517 tetrahedra exceeded the
300 s local limit, so this is a verified medium functional solver path, not a
scalability claim.

Fixed-domain temporal verification is separate from the ALE card. The native
P2/P1 backward-Euler operator is exercised on a stationary four-tetrahedron
domain with the analytic spatially uniform field
`u_x(t)=0.01 sin(2t) m/s` and `dp/dx=-rho du_x/dt`. For time steps
`0.1, 0.05, 0.025, 0.0125 s`, pressure-gradient errors are
`1.33723, 0.693547, 0.352798, 0.177878 Pa/m`; observed orders are
`0.947, 0.975, 0.988`, passing the frozen `0.85–1.15` backward-Euler gate.
`benchmarks/t2_native_temporal_evidence.json` binds these values to the exact
native element/runtime/test sources. This completes the temporal half of the
Poiseuille/manufactured checklist item, but not the still-missing three-level
strict Poiseuille spatial-validation half.

Fixed-domain spatial verification is now independent as well. The native
momentum element accepts an affine body acceleration (constant plus `x/y/z`
coefficients for each component), with finite-value validation and a zero
default that preserves earlier cases. On a unit cube, the manufactured steady
solution `u_x=1e-4 y^3 m/s`, `u_y=u_z=0`, `p=0` is driven by
`b_x=-6 mu 1e-4 y/rho`. Exact velocity is prescribed on the whole boundary;
all interior P2 velocity and P1 pressure unknowns are solved by the autonomous
dense reference runtime. With 2/3/4 subdivisions per axis, relative velocity
L2 errors are `1.13759e-2`, `3.37066e-3`, and `1.42203e-3`; the nonuniform
refinement ratios give orders `2.99998` and `2.99991`. This verifies the
expected P2 L2 spatial behavior of the solved operator. It does not replace the
large frozen Poiseuille fine level, whose strict three-level gate remains open.
The source-bound record is
`benchmarks/t2_native_spatial_evidence.json`.
