# Native moving-tetra species reference path

`NativeTetMovingSpeciesTransport.hpp` adds a project-owned P1 concentration
transport step on the native tetrahedral ALE mesh. It accepts the native P2
Taylor–Hood flow state (or P2 velocity nodes directly) and the same P1 mesh
velocity used by ALE flow. Concentration is in `mol/m³`, diffusion in `m²/s`,
source in `mol/m³/s`, and outward boundary transfer in `mol/s`.

The backward-Euler weak form uses **separate previous and current mass
matrices**, not a fixed-domain mass matrix:

`(M_current c_current − M_previous c_previous)/dt + advection(u−w)
+ diffusion = source`.

Advection is integrated conservatively: the volume term uses `u−w`, outward
faces use the solved interior concentration, and inward faces require an
explicit label-to-inflow-concentration value. Missing inward data fail closed.
Mesh velocity must agree with `(x_current−x_previous)/dt` at each vertex;
inconsistent geometry/velocity input is rejected before assembly.
With the natural zero diffusive flux used here, the reported budget is
`(inventory_current−inventory_previous)/dt + outward_advective_flux−source`.

`make t7-native-moving-species-test` runs two-cell tests for stationary
concentration, rigid mesh translation with `u=w`, 20% mesh expansion with
unchanged inventory but correctly diluted concentration, through-flow,
inflow loading, a translating-mesh `u=2`, `w=1` case that matches a stationary
`u=1` inflow case, diffusion, uniform source, native P2-state adaptation, and
invalid-input rejection. The expansion case conserves `2/3 mol` while its
concentration changes from `2` to `5/3 mol/m³`; this would fail if the old
mass matrix were silently evaluated on the current mesh. The build is
PETSc-free, and the FEM basis, quadrature, assembly, and scalar dense solver
are native code, not DOLFIN.

`NativeTetMovingSpeciesPetscRuntime.hpp` uses the **same sparse native
assembly**, partitioning tetrahedra and boundary faces by MPI rank. PETSc is
used for KSP/matrix/vector algebra only. On this WSL workstation, 1-, 2-, and
4-rank runs matched the dense small-case concentration and mass budget. A
160-vertex test (beyond the dense limit) passed uniform-state, nonzero-source,
and balance checks; a two-step translated-then-expanded mesh preserved
inventory. Build with
`make t7-native-moving-species-petsc-build PETSC_DIR=/path/to/petsc` and run
`make t7-native-moving-species-petsc-test PETSC_DIR=/path/to/petsc
T7_SPECIES_MPI_RANKS=2` on an allocated compute resource.
An independent scalar spatial-convergence test uses an exact unit cube,
zero flow, homogeneous diffusive Neumann boundaries, and initial field
`c₀=2+cos(πx) mol/m³`. At `D=0.1 m²/s`, `dt=0.1 s`, the continuous-in-space
backward-Euler reference is
`c₁=2+cos(πx)/(1+Dπ²dt)`, so the measured error isolates P1 spatial
discretization from time integration. On 2/4/8 subdivisions per axis,
Galerkin relative L2 errors are `0.221088/0.0596235/0.0152552`, orders
`1.891/1.967`; opt-in monotone graph diffusion gives
`0.393917/0.124696/0.0345526`, orders `1.659/1.852`. Both preserve the
`2 mol` inventory and have step balance defect below `1e-9 mol/s` at
1/2/4 MPI ranks. The monotone option is more diffusive at these resolutions;
this no-flow smooth-solution result does not validate high-Péclet accuracy.
Run `make t7-native-moving-species-spatial-convergence-test
PETSC_DIR=/path/to/petsc T7_SPECIES_MPI_RANKS=2`.
The same 2/4/8 subdivision test also covers a zero-diffusivity (`Pe=∞`),
constant `uₓ=1 m/s` inflow front, starting from `c₀=0` with inlet
`c=1 mol/m³`. At `dt=0.2 s`, the continuous-in-space backward-Euler
reference is `c₁(x)=exp(−x/(uₓdt))`; side faces carry zero normal flux and
the outlet is natural. Galerkin relative L2 errors are
`0.190895/0.0587482/0.0148265` (orders `1.700/1.986`), but the coarsest
mesh undershoots to `−0.01330 mol/m³`. The monotone mode stays nonnegative
at all three levels and is conservative; errors
`0.428964/0.281438/0.172991` decrease with orders `0.608/0.702`.
This quantifies the substantial low-order diffusion cost rather than
claiming high-order accuracy or positivity for all geometry and source data.
For the same front at every level, a rigid `+0.02 m` ALE translation uses
`wₓ=0.1 m/s` and `uₓ=1.1 m/s`, preserving the fixed-case relative velocity
`u−w=1 m/s`. Galerkin and monotone moving solutions both match their fixed
counterparts: across 1/2/4-rank tests the largest nodal concentration
difference was `1.34e-15 mol/m³`, inventory difference
`4.45e-16 mol`, and outward-flux difference `2.23e-16 mol/s` (gate `1e-10`
in each quantity). This tests ALE translation invariance of the native
transport operator, not deforming-mesh accuracy or geometric conservation
under arbitrary motion.
The spatial regression also stretches the cube by 20% in `x` over `0.2 s`
with prescribed `u=w=x m/s` on the previous mesh. Zero relative advection
and zero source require the uniform `2 mol/m³` initial field to become
`5/3 mol/m³` while retaining `2 mol` inventory and zero boundary species
flux. Galerkin and monotone modes on 2/4/8 subdivisions pass at 1/2/4 MPI
ranks; the largest nodal concentration error is `8.91e-12 mol/m³` (gate
`1e-10`). This isolates the changing-geometry mass matrix and ALE species
accounting with a prescribed velocity field; `u=w=x` is not claimed as an
incompressible Navier–Stokes solution.
The same native face quadrature now reports outward-positive advective species
flux by boundary label, including negative prescribed inflow. A labelled
two-tetra test gives `-2 mol/s` on the inlet and `+2 mol/s` across the other
faces; dense and MPI/PETSc paths agree, and labelled flux sums to the global
budget. A translating `u=2`, `w=1 m/s` case matches the fixed `u=1` case
label by label. This is the port-level accounting used by the graph adapter below.

`NativeTetMovingSpeciesPorts.hpp` binds one logical graph species to unique
tetra boundary labels using the existing `CouplingPort`, `PortBoundaryData`,
`PortState`, and `SpeciesStepAccounting` types. It requires inlet concentration
only for inward flow, rejects concentration imposed on outward ports, maps
native outward-positive relative `m³/s` and species `mol/s` to observations,
and converts labelled flux into accepted-step port amounts in mol. A boundary
label with simultaneous material inflow and outflow is rejected: a single
graph-port donor concentration cannot faithfully represent that mixture.
The contract test verifies two-port exchange and inventory balance, a
successful reversed-flow exchange, missing/wrong-direction inputs,
duplicate/missing labels, unported species transfer, and
mixed-direction rejection. The staged runtime below uses this contract.
The test also constructs a real `SimulationGraph` pressure/flow/species edge,
checks its acyclic component plan and conservative edge residual, and exercises
reversed flow. At zero flow it publishes the current-mesh area-mean boundary
concentration so a graph executor retaining its prior donor still has a valid
concentration. Invalid negative nodal concentration is rejected before a port
observation or accepted-state commit; high-Péclet positivity is not yet solved.

`NativeTetAleGraphPorts.hpp` supplies the matching native hydraulic port
contract for the staged domain runtime. It translates an outward-positive
graph *relative* flow target to the native ALE solver's *absolute* velocity
controller by adding the mesh's labelled normal flow, and observes relative
`m³/s`, area `m²`, and area-mean P1 pressure `Pa` from a native P2/P1 state.
An analytic one-tetra test verifies fixed and translating meshes, pressure and
flow input selection, missing/wrong-time rejection, and `u=2,w=1` matching
fixed `u=1`. A separate PETSc test uses the adapter's boundary conditions in
the actual native nonlinear ALE Navier–Stokes solver: fixed and translating
four-tetra cases converge, and measured relative flow matches its graph
target (`0.1` fixed, `0.075 m³/s` moving, with `0.1 m³/s` absolute controller
in both). The solved P2 velocity then drives the native P1 species PETSc
solver; hydraulic and scalar port relative flows agree, the species flux is
`2 mol/m³` times the signed port flow, and the moving-domain inventory remains
`1/3 mol`. Run `make t7-native-ale-graph-petsc-test PETSC_DIR=/path/to/petsc
T7_SPECIES_MPI_RANKS=2`. This test validates the solved native
flow-to-species port contract independently of the graph executor.

`NativeTetAleFlowTransportDomainAdapter.hpp` now implements the existing
`CoupledDomainRuntime` and `StagedFlowTransportDomainRuntime` interfaces with
native ALE/PETSc P2/P1 hydraulics and native P1 transport. It supports
prescribed trial mesh motion, pressure/flow graph inputs, optional no-slip
wall labels, transport concentration routing, port observations, amount
accounting, hydraulic/transport rollback, abort, and accepted-step commit.
The standalone staged test retries both hydraulic and transport trials during
a fixed first step, rejects and aborts a
moving second step with a missing inlet concentration without changing the
committed scalar checkpoint, then retries and accepts the moving step. It
passes on 1, 2, and 4 MPI ranks via `make t7-native-staged-ale-test
PETSC_DIR=/path/to/petsc T7_SPECIES_MPI_RANKS=2`. This verifies the runtime's
local lifecycle. A separate graph test wires this runtime to two real
`OneDFlowTransportDomainAdapter` instances through
`SpeciesPressureFlowComponentExecutor`: configured-open-loop 1D → native
tetra ALE → coupled-root 1D, with two pressure/flow and tracer edges. The
explicit test accepts a first fixed step, then deliberately tries a moving
second step. The latter produces local inflow on boundary label 2 despite
its net outflow, so the scalar solve rejects the missing backflow
concentration. The executor aborts the whole graph without changing the
committed native flow/species, 1D flow steps, or donor ownership; a fixed-mesh
retry of that same second step succeeds. A separate rank-0-only precommit
failure is propagated by the graph's collective execution policy to every
rank; it also leaves committed state unchanged before the successful retry.
Both accepted steps check edge and domain amounts, global conservation, and
accepted native state. The same test executable also constructs a 24-tetra
channel with separate inlet, outlet, and no-slip wall labels. Its first fixed
step and second rigidly translating ALE step both commit through the complete
1D → native tetra → 1D graph on 1, 2, and 4 MPI ranks. A roundoff-scale
normal flux on a no-slip wall is treated as numerical zero using the global
velocity scale and machine precision; the distinct, finite local backflow in
the four-tetra negative case remains rejected. Run
`make t7-native-ale-species-graph-1d-test PETSC_DIR=/path/to/petsc
T7_SPECIES_MPI_RANKS=2`. The target also runs fixed-point pressure/flow
iteration on that moving channel. The 1D configured-open-loop adapter now
rearms its sampled inlet after hydraulic rollback, so each iteration has the
same inlet contract; its unit regression and the full graph pass. With
pressure relative tolerance `1e-4`, flow relative tolerance `1e-6`, and
relaxation `0.5`, the moving second step converges in 15 hydraulic iterations
on 1, 2, and 4 ranks. This proves functional strong-coupling convergence at
the stated loose tolerance, not strict physical accuracy, resolved
mixed-direction boundary exchange, or arbitrary failures inside a PETSc
collective. A separate 1D source → native tetra ALE → terminal RCR species
graph now passes explicit and fixed-point coupling on the same moving
24-tetra channel at 1, 2, and 4 ranks. Its second fixed-point step takes 15
iterations in both ordinary and distal-reverse RCR loads; the reverse load
requires an explicit distal donor concentration. Rank-0-only native FEM/RCR
precommit rollback, retry, branch amount, RCR distal amount, and global
species balance are checked by
`make t7-native-ale-species-graph-rcr-test PETSC_DIR=/path/to/petsc
T7_SPECIES_MPI_RANKS=2`. This remains a low-complexity functional example,
not physiological validation.

`NativeTetMovingSpeciesCheckpoint.hpp` adds an in-memory accepted-step
checkpoint for the same native scalar path. It binds the reference mesh,
boundary labels, and diffusivity to a model SHA-256 identity, stores the
current moving coordinates, concentration, accepted-step count, and time,
and verifies payload checksum and state identity before restore. The commit
occurs only after the PETSc solve and mass-budget gate succeed; a rejected
trial leaves the previous state unchanged. The two-tetra restart test resumes
after translation and matches an uninterrupted expansion step, including
`2/3 mol` inventory. It also rejects truncated/corrupted payloads and wrong
models. Run `make t7-native-moving-species-checkpoint-test
PETSC_DIR=/path/to/petsc T7_SPECIES_MPI_RANKS=2` on an allocated resource.
This standalone codec is not itself a graph-level or distributed file
checkpoint; a fixed-test five-shard graph bundle is described below.

The native scalar assembler now also offers an explicit `monotone=true`
low-order option. For each owned tetrahedron and boundary face it adds the
smallest symmetric graph-Laplacian diffusion needed to eliminate positive
off-diagonal entries in that contribution. The row and column sums of this
correction are zero, so the constant-state and inventory identities remain
unchanged. In a two-tetra, zero-diffusivity inlet-front test, the original
Galerkin step has a `-0.109 mol/m³` nodal undershoot; the low-order option is
nonnegative (to `1e-12 mol/m³`) and conservative. A rigidly translated mesh
with the same relative velocity reproduces the fixed-mesh concentration and
budget. Dense and PETSc solutions agree at 1/2/4 ranks. This is deliberately
diffusive. The staged native FEM runtime now selects it explicitly and binds
that choice into its checkpoint model identity; a default-mode checkpoint
cannot restore into a monotone-mode runtime. A labelled fTetWild pipe with
zero initial FEM concentration and positive 0D inlet donor passes the
source→native FEM→RCR graph with the option enabled, while the same
unstabilized front is rejected for negative port concentration. This is a
bounded functionality test, not a proof of positivity for arbitrary moving
geometries, signed body sources, or underresolved graph cases.
The generated fTetWild pipe front also passes a second, rigidly translated
ALE step with mesh displacement `-1e-5 m`, nonnegative FEM concentration,
conservative edge/global budgets and 1/2/4-rank same-mesh field agreement.
The same fTetWild graph driver accepts optional
`--species-output-dir=/new/output/path`. It publishes immutable rank-owned
linear-tetra VTU snapshots after each accepted step with the native P1
`concentration_mol_m3`, `reference_position_m`, and `displacement_m`; VTK
`Points` are current coordinates. The local regression verifies all tetrahedra,
all vertex concentrations against the accepted native state, and the second
step's `-1e-5 m` translation at 1/2/4 ranks. This is a reusable native field
exporter wired into a fixed pipe functional driver, not yet a configurable
species production CLI or a generic organ workflow.
The 24-tetra channel additionally exercises a moving second step from a
zero-concentration front, rank-local precommit rollback, and five-shard
cross-process continuation at 1/2/4 ranks. The monotone option remains
explicit in the staged FEM model identity; this is a functional ALE/graph
restart test, not a spatial-convergence study.

This is still not production ALE species support. The dense reference is
limited to 128 vertices; the sparse PETSc path does not materialize an `N²`
matrix, but currently replicates geometry, previous concentrations, and the
solved field on every rank. The optional low-order mode does not yet provide
a general high-Péclet positivity guarantee; there is no reaction or wall-exchange term, and only a fixed-test graph/file bundle rather than a generic production checkpoint coordinator. The channel graph test solves native flow on a
prescribed rigidly moving mesh. None is a
coupled LV species simulation. Those gaps keep T7's
parent ALE-species item open.
