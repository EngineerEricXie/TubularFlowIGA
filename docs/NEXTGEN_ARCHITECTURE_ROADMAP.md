# TubularFlowIGA Next-Generation Architecture and Implementation Roadmap

**Repository:** `EngineerEricXie/TubularFlowIGA`
**Primary audience:** Codex coding agent and human reviewers
**Plan date:** 2026-08-31
**Purpose:** Evolve TubularFlowIGA from a tube-specialized IGA flow code into a multiscale, geometry-agnostic cardiovascular / biomedical simulation framework supporting 0D/1D/3D coupling, arbitrary anatomical surfaces, moving domains, transport, and eventually fluid-structure interaction (FSI).

---

# 0. Instructions to the Codex Agent

This document is an **engineering roadmap**, not permission for a big-bang rewrite.

The highest priority is to preserve existing validated behavior while progressively introducing reusable abstractions.

## Non-negotiable rules

1. **Do not rewrite the current solver stack from scratch.**
   - Preserve the validated C++/PETSc CPU path.
   - Preserve the existing 1D solver.
   - Preserve the CUDA backend unless a change is explicitly needed for a phase.
   - Preserve existing source-only examples and validation behavior.

2. **Do not start with heart FSI.**
   - First build the coupling architecture using existing 1D and body-fitted 3D solvers.
   - Then generalize geometry through immersed IGA.
   - Then add moving boundaries.
   - Only then add two-way FSI.

3. **Each implementation step must be independently testable.**
   - Prefer small PR-sized changes.
   - Add unit tests and regression tests before moving to the next phase.
   - No phase is complete only because it compiles.

4. **Maintain backward compatibility.**
   - Existing 1D schema-v3 cases must continue to work.
   - Existing 3D schema-v4 cases must continue to work.
   - Existing `.ntiga` databases should remain readable unless a new version is explicitly introduced.
   - Old standalone executables remain usable even after reusable runtime libraries are introduced.

5. **Use SI consistently at coupling interfaces.**
   - Pressure: Pa
   - Flow: m^3/s
   - Area: m^2
   - Velocity: m/s
   - Time: s
   - Species flux units must be explicit and consistent with the species definition.

6. **Do not silently hide conservation errors.**
   - Coupling residuals, mass imbalance, species imbalance, and failed nonlinear solves are first-class diagnostics.
   - Fail loudly when configured convergence criteria are not met.

7. **Keep geometry, discretization, physics, and coupling separate.**
   - A Navier–Stokes kernel must not know whether the geometry came from SWC sweeping, a Cartesian background spline, or a future external mesher.
   - A coupling algorithm must not special-case "1D-to-3D" when a generic port/subsystem interface is sufficient.

8. **Use the existing repository as the source of truth.**
   - Before modifying any component, inspect the current implementation and tests.
   - This plan describes intended architecture; if repository details have changed, adapt the implementation while preserving the design principles.

---

# 1. Executive Summary

TubularFlowIGA currently has a strong but specialized architecture:

```text
SWC / radius-annotated centerline
        ↓
smoothing + resampling
        ↓
swept hexahedral control mesh
        ↓
spline / Bézier extraction
        ↓
.ntiga database
        ↓
MPI/PETSc or CUDA
        ↓
3D Navier–Stokes / transport
```

It also contains a separate native 1D vascular solver with compliant flow, pressure/flow network formulations, junction models, terminal models, transport, physiology, and VCA coupling.

The core limitation is **not IGA itself**. The major limitation is that the current 3D geometry pipeline assumes a tube-like rooted network and generates a body-fitted swept control mesh. That representation becomes unnatural for chambers, ventricles, valves, highly irregular anatomy, and large moving interfaces.

The target architecture is therefore:

```text
                   MULTISCALE SIMULATION GRAPH
                              │
        ┌─────────────────────┼─────────────────────┐
        │                     │                     │
       0D                    1D                    3D
 circulation/circuits   vascular networks     resolved ROIs
        │                     │                     │
        └────────────── Coupling Ports ─────────────┘
                              │
                    P / Q / species flux
                              │
                   common time/coupling driver
```

For 3D geometry:

```text
                     Geometry Frontends
                           │
             ┌─────────────┴─────────────┐
             │                           │
       Tubular / centerline         Surface anatomy
             │                           │
       swept body-fitted            immersed geometry
          spline mesh             on background B-splines
             │                           │
             └───────────┬───────────────┘
                         ↓
                 common physics kernels
```

The framework should ultimately support:

- Standalone 1D vascular simulation.
- Standalone 3D IGA simulation.
- 0D/1D/3D heterogeneous coupling.
- Arbitrary sequences such as `1D → 3D → 1D → 3D → 0D`.
- Multiple independent 3D "islands" embedded in a large 1D network.
- Conservative flow and pressure coupling.
- Conservative species transport coupling.
- Arbitrary anatomical surfaces using immersed / unfitted IGA.
- Prescribed moving-wall cardiac flow.
- Deformable vessels and valve FSI.
- Optional advanced divergence-conforming flow spaces for valve/no-leak applications.
- Existing tubular workflows as a fast specialized geometry backend rather than the only supported geometry path.

---

# 2. Verified Current Repository Capabilities and Constraints

The following observations should guide refactoring.

## 2.1 Existing 3D path

Current documentation describes the 3D code as a research toolkit specialized for tube-like networks.

The present 3D path includes:

- SWC or radius-annotated line-OBJ geometry.
- C++ smoothing and resampling.
- Swept hexahedral control-mesh generation.
- Spline and Bézier extraction.
- Partition-aware `.ntiga` database.
- CPU MPI/PETSc flow and transport.
- Single-GPU CUDA backend.
- Steady and transient incompressible Navier–Stokes.
- Stabilized VMS formulation.
- Generic multispecies transport.
- 3D VCA coupling.
- Boundary flow, pressure, and species-flux measurement.

The tubular mesher currently assumes a connected rooted tree and currently supports one child or binary branching at nonterminal locations. This assumption must **remain valid for the tubular backend**, but it must stop being a global framework assumption.

## 2.2 Existing spline representation is more general than the tubular mesher

The spline data structures already contain concepts for:

- boundary entities,
- extraordinary entities,
- T-junctions / T-edges,
- unstructured spline construction,
- variable element connectivity,
- Bézier extraction.

Therefore:

> Do not equate "current tubular mesher" with "the maximum capability of the spline / solver layer."

However, validate general unstructured control-mesh behavior before depending on it for production general-anatomy meshing.

## 2.3 Existing `.ntiga` / 3D element representation is a useful abstraction boundary

The current 3D element representation already separates:

- global basis connectivity,
- extraction coefficients,
- Bézier geometry points,
- owner/partition information,
- boundary labels.

Element connectivity is variable-size.

The currently hard-coded geometric/discretization assumption is primarily:

- cubic tensor-product Bernstein basis,
- 4 x 4 x 4 = 64 Bézier points,
- full-element tensor-product quadrature.

This makes **quadrature abstraction** a key step for immersed IGA.

## 2.4 Existing 3D runtime already contains coupling-relevant infrastructure

`TransientFlowRuntime` already contains concepts such as:

- `FlowPortMeasurements`,
- outlet flow measurement,
- area-averaged pressure measurement,
- species-flux aggregation,
- pressure-traction/outlet model handling,
- iterative outlet coupling,
- PETSc state vectors,
- transient state management,
- convergence and mass-balance diagnostics.

This code should be reused and generalized, not duplicated.

## 2.5 Existing 1D solver is already substantial

The native 1D solver supports:

- rooted vascular networks,
- steady rigid Poiseuille flow,
- compliant A/Q formulations,
- explicit Rusanov finite-volume flow,
- multiple implicit PETSc formulations,
- compliant wall laws,
- junction flow conservation,
- static or total pressure balancing,
- junction loss models,
- pressure, resistance, and RCR terminals,
- transport of arbitrary species,
- physiology helpers,
- checkpoint/restart,
- VCA coupling.

The 1D documentation explicitly states that 1D and 3D are not yet directly coupled. The next-generation architecture should close this gap before attempting heart FSI.

---

# 3. Product / Scientific Vision

The final project should no longer be defined as:

> centerline → swept mesh → IGA flow solver

It should be defined as:

> **a multiscale spline-based cardiovascular / biomedical simulation engine in which different regions of one physical system can use different spatial fidelities and geometry representations.**

Four defining properties:

## 3.1 Multiscale

```text
0D + 1D + 3D
```

Use each dimensional model where appropriate.

## 3.2 Multi-fidelity

```text
large cheap network
+
small high-fidelity local regions
```

Example:

```text
systemic 1D network
    │
    ├── 3D aortic root / valve
    ├── 3D aneurysm
    ├── 3D stenosis
    └── 3D organ ROI
```

## 3.3 Multiphysics

Long-term support:

```text
blood flow
+ species transport
+ reduced physiology
+ structural mechanics
+ FSI
+ external 0D circulation / devices
```

## 3.4 Geometry-agnostic core

Support both:

```text
centerline/SWC
→ tubular body-fitted spline
```

and

```text
STL/VTP/anatomical surface
→ immersed background B-spline
```

without changing the physics layer.

---

# 4. Target User-Visible Capabilities

The completed architecture should support the following classes of simulation.

| Problem | Target support |
|---|---|
| Straight vessel | Yes |
| Vascular bifurcation | Yes |
| Large vascular tree | Yes |
| Neuron material transport | Preserve |
| 1D compliant arterial network | Yes |
| 1D multispecies transport | Yes |
| 1D → 3D → 1D vessel | Yes |
| Multiple 3D regions inside one 1D network | Yes |
| 0D → 1D → 3D circulation | Yes |
| Patient-specific irregular vessel | Yes |
| Aneurysm CFD | Yes |
| Aortic root CFD | Yes |
| Rigid LV / LA chamber CFD | Yes |
| Prescribed moving LV / LA CFD | Yes |
| 1D–3D species transport | Yes |
| Compliant artery FSI | Yes |
| Aneurysm FSI | Yes |
| Aortic/mitral valve FSI | Long-term |
| Valve contact / near-contact flow | Long-term |
| Full myocardial FSI | Long-term |
| Full electromechanics + blood flow | Future research, not near-term scope |

---

# 5. High-Level Software Architecture

The target architecture should separate five layers.

```text
┌─────────────────────────────────────────────────────────────┐
│                    Simulation Configuration                 │
│         domains + geometry + physics + couplings           │
└─────────────────────────────┬───────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│                  Multiscale Simulation Graph                │
│       nodes = subsystems, edges = coupling relations       │
└─────────────────────────────┬───────────────────────────────┘
                              │
             ┌────────────────┼────────────────┐
             ▼                ▼                ▼
        0D subsystem      1D subsystem      3D subsystem
             │                │                │
             └────────── Coupling Ports ───────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│                  Coupling / Time Integrator                 │
│ fixed point / Aitken / subcycling / commit / rollback     │
└─────────────────────────────┬───────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│              Geometry + Discretization Backends            │
│ tubular spline | immersed spline | future other backend    │
└─────────────────────────────┬───────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│                       Physics Kernels                       │
│ Navier–Stokes | transport | structure | constitutive laws │
└─────────────────────────────────────────────────────────────┘
```

---

# 6. Core Abstraction: CoupledSubsystem

The most important near-term refactor is to make the existing standalone solvers reusable inside a coupled simulation.

A conceptual interface:

```cpp
class CoupledSubsystem {
public:
    virtual ~CoupledSubsystem() = default;

    virtual void Initialize() = 0;

    // Begin a physical macro-step and preserve the committed state.
    virtual void BeginStep(double time, double dt) = 0;

    // Apply current trial coupling boundary data.
    virtual void SetPortInput(
        const std::string& port_id,
        const PortBoundaryData& input) = 0;

    // Solve from the last committed state to the trial state at time + dt.
    virtual void SolveTrial() = 0;

    // Return interface quantities generated by the trial solve.
    virtual PortState GetPortState(
        const std::string& port_id) const = 0;

    // Throw away the current trial result and return to the state at BeginStep().
    virtual void RollbackTrial() = 0;

    // Close an opened physical step without accepting it.  Unlike rollback,
    // this returns to the idle committed phase and is idempotent there.
    virtual void AbortStep() = 0;

    // Validate that the trial is commit-ready without publishing it.
    virtual void PrepareCommitStep() = 0;

    // Publish a prepared commit.  This operation must not throw.
    virtual void FinalizeCommitStep() noexcept = 0;

    // Compatibility shorthand for prepare followed by finalize.
    virtual void CommitStep() = 0;
};
```

Possible implementations:

```cpp
class OneDFlowDomain;
class ThreeDBodyFittedFlowDomain;
class ThreeDImmersedFlowDomain;
class ZeroDCircuitDomain;
class StructuralDomain; // later
```

## 6.1 Why rollback / commit is mandatory

For strongly partitioned coupling, one physical step may require multiple trial solves.

Example:

```text
state at t_n
   │
   ├── trial iteration 0 → reject
   ├── trial iteration 1 → reject
   ├── trial iteration 2 → converged
   │
   └── commit state at t_(n+1)
```

A subsystem must not advance its internal physical history multiple times during coupling iterations.
For a graph of subsystems, committing sequentially is also insufficient: a
late validation failure would otherwise leave earlier domains committed and
later domains uncommitted. The graph executor must prepare every domain before
it invokes any nonthrowing finalizer. On failure it must attempt to abort every
opened domain, preserving the primary error and reporting cleanup failures.

The implementation should distinguish:

- committed previous state,
- trial current state,
- trial auxiliary states,
- persistent terminal/capacitor state,
- transient transport state.

Do not use disk checkpoint/restart as the normal in-memory coupling rollback mechanism. Reuse its state-serialization concepts where useful, but coupling rollback should be efficient and in-memory.

---

# 7. Core Abstraction: CouplingPort and PortState

Every dimension should expose the same conceptual interface.

```cpp
enum class PortQuantity {
    FlowRate,
    MeanPressure,
    MeanNormalTraction,
    TotalPressure,
    Area,
    SpeciesConcentration,
    SpeciesFlux
};

struct PortState {
    double time = 0.0;

    double area = 0.0;
    double flow_rate = 0.0;

    double mean_pressure = 0.0;
    double mean_normal_traction = 0.0;
    double total_pressure = 0.0;

    std::map<std::string, double> concentration;
    std::map<std::string, double> species_flux;
};

struct CouplingPort {
    std::string id;
    std::string subsystem_id;

    // Dimension-specific selector:
    // 1D: network terminal/node/segment
    // 3D: boundary label / immersed interface id
    // 0D: circuit node/branch
    PortLocator locator;

    PortOrientation orientation;

    std::set<PortQuantity> provides;
    std::set<PortQuantity> requires;
};
```

## 7.1 Sign convention

Define one convention globally and test it.

Recommended:

- Positive `flow_rate` means flow **out of the subsystem through the port outward normal / oriented terminal direction**.
- Coupling edges apply the corresponding sign transformation so mass transferred out of subsystem A enters subsystem B.
- Never rely on implicit label ordering to determine sign.

The coupling manifest must record port orientations.

---

# 8. 1D–3D Coupling Formulation

This is the first major new scientific capability and should be implemented before immersed anatomy.

## 8.1 Flow conservation

For a 3D interface surface \(\Gamma\):

\[
Q_{3D} =
\int_{\Gamma} \mathbf{u}\cdot\mathbf{n}\,dA.
\]

The 1D model naturally carries \(Q_{1D}\).

At a conservative interface:

\[
Q_{1D} = Q_{3D}
\]

after sign convention/orientation is accounted for.

## 8.2 Pressure / traction transfer

The simplest initial pressure measurement is:

\[
\bar p_{3D}
=
\frac{1}{A}
\int_\Gamma p\,dA.
\]

For applying 1D pressure back to 3D, use a natural pressure traction:

\[
\mathbf t = -P_{1D}\mathbf n.
\]

This is compatible with the existing 3D pressure-traction mechanism and avoids imposing pressure by deleting a continuity row.

Later implementations may support:

- mean normal traction,
- total pressure,
- energy-consistent coupling conditions.

Do not block the initial coupling architecture on the most sophisticated interface formulation.

## 8.3 Coupling Mode A: 3D flow → 1D pressure → 3D traction

This should be the first/default coupling mode.

```text
3D solve
   │
   └── measure Q_3D
              │
              ▼
           1D solve
              │
              └── return P_1D
                        │
                        ▼
              3D pressure traction
```

Use cases:

- 3D aneurysm connected to distal 1D tree.
- 3D bifurcation with 1D upstream/downstream vessels.
- 3D heart/aortic root connected to systemic 1D circulation.

## 8.4 Coupling Mode B: 3D pressure → 1D flow → 3D velocity profile

Support later.

```text
3D mean pressure
      ↓
   1D solve
      ↓
     Q_1D
      ↓
3D inlet velocity reconstruction
```

The reconstructed profile must satisfy:

\[
\int_\Gamma u_n\,dA = Q_{1D}.
\]

Profile options may eventually include:

1. uniform,
2. parabolic,
3. Womersley,
4. profile recycling,
5. data-driven or reduced-order profiles.

Do not implement all profile models in the first PR.

---

# 9. Strong Partitioned Coupling

The first prototype may use explicit staggered coupling, but the production path should support strongly coupled partitioned iterations.

For a macro time step \(t_n \rightarrow t_{n+1}\):

```text
initial port guess
       │
       ▼
  solve subsystem A
       │
  obtain port output
       │
       ▼
  solve subsystem B
       │
  obtain updated port output
       │
       ▼
 compute interface residual
       │
    converged?
     /     \
   no       yes
   │         │
relax       commit
   │
   └────────→ repeat
```

## 9.1 Initial convergence metrics

Make thresholds configurable.

Examples:

\[
r_Q =
\frac{|Q_A + Q_B|}
{\max(Q_{\mathrm{ref}}, |Q_A|, |Q_B|)}
\]

and

\[
r_P =
\frac{|P^{k+1}-P^k|}
{\max(P_{\mathrm{ref}}, |P^{k+1}|)}.
\]

Do not use raw dimensional differences alone.

Recommended first default target:

```text
coupling_relative_tolerance = 1e-4
```

This value must remain configurable and should be revisited during verification.

## 9.2 Relaxation

Implement in this order:

1. fixed under-relaxation,
2. Aitken Δ² dynamic relaxation.

Conceptually:

\[
x^{k+1}
=
x^k + \omega_k
(\tilde x^{k+1}-x^k).
\]

Aitken relaxation will later also be reusable for FSI.

## 9.3 Failure behavior

If maximum coupling iterations are reached:

- report the interface,
- report residual history,
- report last P/Q values,
- exit nonzero unless explicitly configured for diagnostic continuation.

Never silently accept an unconverged strongly-coupled step.

---

# 10. Multirate Time Integration and Subcycling

Different subsystems should not be forced to use the same internal time step.

Define a coupling macro-step:

\[
\Delta T
\]

while individual subsystems may use internal substeps.

Example:

```text
3D: one step of ΔT = 1e-3 s

1D:
  substep 1
  substep 2
  ...
  substep 10
```

The existing explicit 1D compliant solver already has CFL-driven internal substeps, so the coupling architecture should expose a macro-step without breaking internal stability logic.

The coupling driver should eventually support:

```cpp
AdvanceTrial(time, macro_dt, port_boundary_history);
```

rather than assuming one internal solver step per coupling step.

Future work may require interpolation of coupling data over a macro-step. The first implementation may hold trial boundary values constant over the subcycles if documented and verified.

---

# 11. Multidomain Simulation Graph

Do not design the framework around exactly one 1D and one 3D solver.

Model the complete simulation as a graph.

```cpp
struct SimulationGraph {
    std::map<std::string, std::unique_ptr<CoupledSubsystem>> domains;
    std::vector<CouplingEdge> couplings;
};
```

Examples:

```text
1D → 3D → 1D
```

```text
1D → 3D → 1D → 3D → 1D
```

```text
              ┌→ 3D aneurysm → 1D
1D network ───┼→ 3D renal ROI → 1D
              └→ 1D branch
```

This graph abstraction is essential for future organ/system digital twins.

---

# 12. Configuration Schema Evolution

Do not break schema v3/v4.

Introduce a new schema only when the first multidomain case is ready.

Tentative next-generation concept:

```json
{
  "schema_version": 5,

  "time": {
    "dt": 0.001,
    "steps": 1000
  },

  "domains": [
    {
      "id": "upstream",
      "dimension": "1d",
      "kind": "network_flow",
      "case": "domains/upstream"
    },
    {
      "id": "aneurysm",
      "dimension": "3d",
      "kind": "body_fitted_iga_flow",
      "database": "aneurysm.ntiga",
      "case": "domains/aneurysm"
    },
    {
      "id": "downstream",
      "dimension": "1d",
      "kind": "network_flow",
      "case": "domains/downstream"
    }
  ],

  "couplings": [
    {
      "id": "upstream_to_aneurysm",
      "a": {
        "domain": "upstream",
        "port": "outlet"
      },
      "b": {
        "domain": "aneurysm",
        "port": "inlet"
      },
      "mode": "pressure_flow",
      "algorithm": {
        "kind": "strong_partitioned",
        "relaxation": "aitken",
        "relative_tolerance": 1e-4,
        "maximum_iterations": 30
      }
    }
  ]
}
```

The exact schema should be refined during implementation.

## Backward compatibility requirement

Standalone schemas must continue to dispatch directly to existing solver behavior.

Do not force old examples into graph format.

---

# 13. 1D–3D Transport Coupling

After flow coupling is verified, extend the same port abstraction to scalar transport.

The 3D conservative species flux may be represented as:

\[
\Phi_c^{3D}
=
\int_\Gamma
\left(
c\mathbf u - D\nabla c
\right)\cdot\mathbf n\,dA.
\]

The 1D solver already uses conservative \(A C\)-type state variables.

The interface should conserve transported quantity:

\[
\Phi_c^{1D} + \Phi_c^{3D} = 0
\]

under the selected port orientation.

Support:

- concentration transfer for inflow boundaries,
- conservative advective/diffusive species flux,
- arbitrary configured species names,
- multiple species on one interface,
- diagnostic mass balance per species.

Do not hard-code oxygen as the only coupled species.

---

# 14. 0D Coupling

Existing R/RC/RCR and VCA mechanisms already demonstrate reduced-order external coupling.

The long-term design should expose 0D circuits as normal graph nodes:

```text
3D ── port ── 0D
1D ── port ── 0D
0D ── port ── 1D ── port ── 3D
```

Potential 0D components:

- Windkessel,
- heart chamber elastance,
- systemic / pulmonary circulation,
- device circuits,
- reservoirs,
- oxygenator/dialyzer-style devices.

Avoid creating a separate coupling framework for every reduced-order model.

---

# 15. Geometry Refactor: Stop Making Sweeping a Global Assumption

The current tubular path remains valuable and should be preserved as a high-efficiency backend.

Target geometry architecture:

```text
GeometrySource
     │
     ├── CenterlineGeometry
     │       ↓
     │   TubularSweepMesher
     │       ↓
     │   body-fitted spline
     │
     └── SurfaceGeometry
             ↓
        ImmersedGeometry
             ↓
       background spline
```

Potential future geometry sources:

- SWC,
- radius-annotated OBJ line network,
- STL,
- VTP/VTK surface,
- segmentation-derived surface,
- moving surface sequence,
- displacement field on reference surface.

Physics code must not parse SWC.

---

# 16. Quadrature Abstraction: Critical Step for Immersed IGA

Current 3D element assembly assumes full tensor-product quadrature.

Refactor the element kernel conceptually from:

```cpp
BuildNavierStokesElement(element, state, parameters);
```

where quadrature is hard-coded internally, to:

```cpp
BuildNavierStokesElement(
    element,
    state,
    parameters,
    volume_quadrature);
```

Define:

```cpp
struct VolumeQuadraturePoint {
    std::array<double, 3> parametric;
    double weight;
};

struct SurfaceQuadraturePoint {
    std::array<double, 3> parametric;
    std::array<double, 3> physical;
    std::array<double, 3> normal;
    double weight;
    int boundary_id;
};
```

Providers:

```text
FullCellQuadrature
CutCellVolumeQuadrature
ImmersedSurfaceQuadrature
```

## Acceptance requirement

Existing body-fitted cases must produce numerically equivalent results after replacing hard-coded loops with `FullCellQuadrature`.

This refactor should happen **before** implementing cut-cell geometry.

---

# 17. Cartesian Background B-Spline Discretization

For immersed anatomy, introduce a structured Cartesian cubic B-spline background.

Initial scope:

- uniform Cartesian grid,
- cubic tensor-product B-splines,
- regular basis connectivity,
- analytic Bézier extraction,
- existing PETSc assembly model where possible.

Do not initially require:

- T-spline local refinement,
- hierarchical splines,
- arbitrary-order splines.

The first goal is a robust reference implementation.

Concept:

```text
bounding box
    ↓
Cartesian cells
    ↓
cubic B-spline connectivity
    ↓
Bezier extraction
    ↓
existing element basis evaluation
```

---

# 18. Surface Geometry Engine

Initial arbitrary-anatomy input should be a clean closed triangulated surface.

Support:

```text
STL / VTP
```

Required infrastructure:

- surface validation,
- orientation / outward normal consistency,
- AABB/BVH acceleration structure,
- point-in-domain classification,
- background-cell / triangle intersection,
- cut-cell identification,
- boundary ID support.

Cell classes:

```cpp
enum class CellClassification {
    Outside,
    Inside,
    Cut
};
```

Do not immediately attempt automatic patient-specific segmentation. The solver should consume an already segmented/cleaned surface.

---

# 19. Cut-Cell Volume Integration

Start with the most robust implementation, not the most optimized.

Recommended first method:

```text
cut background cell
      ↓
recursive octree subdivision
      ↓
classify subcells
      ↓
inside → regular Gauss quadrature
outside → discard
cut → subdivide again
```

Later optimize with:

- direct polyhedral clipping/tessellation,
- moment fitting,
- adaptive quadrature caching.

Validation quantities:

- immersed volume,
- immersed surface area,
- integration of constants,
- integration of linear/polynomial fields,
- convergence under background refinement,
- conservation in CFD cases.

---

# 20. Immersed Boundary Conditions

For arbitrary non-body-fitted walls, strong nodal Dirichlet constraints are not generally available.

Implement weak velocity boundary enforcement using Nitsche-type terms.

Required conceptual features:

- prescribed no-slip wall,
- prescribed moving-wall velocity,
- consistent weak terms,
- penalty/stabilization parameter with physically/numerically motivated scaling,
- immersed surface quadrature.

Do not treat an immersed wall as a fake body-fitted element face.

---

# 21. Cut-Cell Stabilization

Small physical intersections can create ill-conditioning.

Add cut-cell stabilization after the first correct immersed integration prototype.

Candidate path:

- ghost penalty / face-based stabilization,
- support-based stabilization appropriate to B-splines.

Requirements:

- conditioning diagnostics,
- robustness versus arbitrarily small cut fractions,
- no loss of expected convergence in regular cases.

Do not compensate for cut-cell instability merely by increasing Nitsche penalty to extreme values.

---

# 22. Prescribed Moving-Domain Cardiac Flow

Before two-way FSI, support moving anatomical surfaces with prescribed motion.

Input concepts:

```text
surface sequence Γ(t)
```

or

```text
reference surface + displacement field d(X,t)
```

Applications:

- LV filling/ejection,
- LA flow,
- aortic root motion,
- prescribed valve motion,
- patient-specific cine MRI / 4D CT wall trajectories.

Immersed background mesh remains fixed.

At every time step:

1. update surface position,
2. update BVH / acceleration data,
3. reclassify affected background cells,
4. update cut-cell volume quadrature,
5. update immersed surface quadrature,
6. impose wall velocity,
7. solve transient Navier–Stokes,
8. output flow / derived quantities.

Potential outputs:

- velocity,
- pressure,
- vortical structures,
- kinetic energy,
- energy loss,
- wall shear metrics,
- residence time,
- washout,
- transported species / particles.

This is the recommended first "heart" application.

---

# 23. Fluid–Structure Interaction

Only begin after:

- 1D–3D coupling is verified,
- immersed static flow is verified,
- prescribed moving boundaries are verified.

## 23.1 Structural subsystem options

### Vessel wall

Possible shell or solid mechanics.

Applications:

- compliant arteries,
- aneurysm deformation,
- pressure-induced wall motion.

### Valve leaflet

IGA shell is attractive because high-continuity splines are well suited to thin-shell formulations.

Required future capabilities:

- nonlinear shell mechanics,
- immersed fluid-shell coupling,
- opening/closing,
- leaflet contact,
- near-contact flow.

### Myocardium

Full myocardium requires 3D nonlinear anisotropic solid mechanics and possibly active contraction.

Do not approximate the long-term myocardium problem as only a thin shell without explicit scientific justification.

## 23.2 FSI coupling architecture

Reuse `CoupledSubsystem`, `PortState` concepts where appropriate, but FSI interface data is field-valued rather than only scalar P/Q ports.

Need a separate distributed interface coupling representation for:

```text
fluid traction → structure
structure displacement/velocity → fluid boundary
```

Reuse:

- rollback/commit,
- strong partitioned iteration,
- Aitken relaxation,
- time-step synchronization.

---

# 24. Divergence-Conforming Flow Space: Advanced Valve Phase

Current stabilized equal-order VMS flow should remain the primary development path.

For thin immersed valves, numerical leakage across a closed leaflet may motivate a divergence-conforming spline formulation.

Treat this as an advanced later phase.

Possible benefits:

- stronger discrete mass conservation,
- reduced leakage through closed valves,
- improved incompressibility behavior.

Do not block general immersed CFD on this advanced formulation.

---

# 25. Database / File-Format Strategy

## 25.1 Near-term 1D–3D coupling

Do not change `.ntiga` merely to add subsystem coupling.

Use existing databases and runtime port APIs.

## 25.2 Immersed IGA

A future database version may need to represent:

- background cell state,
- cut-cell metadata,
- cached volume quadrature,
- immersed surface quadrature,
- boundary/interface IDs.

Consider either:

### Option A: `.ntiga` v6

Extend the database while retaining v3–v5 readers.

### Option B: topology database + separate cut-geometry cache

Example:

```text
case.ntiga
case.cutcache
```

This may be preferable because moving geometry invalidates cut quadrature much more often than the spline topology.

Do not choose the format prematurely. Prototype the in-memory data model first.

---

# 26. Validation Strategy

Every phase needs four layers of testing.

## 26.1 Unit tests

Examples:

- sign conventions,
- port orientation,
- Aitken update,
- coupling residual,
- state rollback,
- scalar flux balance,
- point-in-surface,
- cell classification,
- quadrature polynomial integration.

## 26.2 Regression tests

Existing examples must continue to match accepted references within current tolerances.

Important existing paths:

- straight 3D vessel,
- Y-bifurcation,
- transient flow cases,
- 1D compliant bifurcation,
- multispecies transport,
- VCA coupling.

## 26.3 Numerical verification benchmarks

Use problems with known/reference behavior.

## 26.4 Application validation

Only after numerical verification.

Do not use a realistic-looking patient geometry as proof that the method is correct.

---

# 27. Required 1D–3D Benchmarks

## Benchmark A: straight vessel replacement test

Construct physically equivalent systems.

### Case A1

```text
all 1D
```

### Case A2

```text
1D → 3D straight section → 1D
```

Compare:

- inlet/outlet Q(t),
- pressure waveforms,
- mean pressure drop,
- phase shift,
- mass conservation,
- pulse transit behavior,
- coupling residual history.

Purpose:

> Verify that inserting a 3D region does not create unphysical global behavior.

## Benchmark B: 3D bifurcation embedded in 1D

```text
                 upstream 1D
                      │
                 ┌────▼────┐
                 │ 3D Y    │
                 │ branch  │
                 └──┬───┬──┘
                    │   │
                   1D   1D
```

Verify:

- branch flow split,
- pressure loss,
- global mass conservation,
- waveform continuity,
- sensitivity to coupling tolerance.

## Benchmark C: aneurysm ROI

```text
upstream 1D
     │
     ▼
┌───────────────┐
│  3D aneurysm  │
└──────┬────────┘
       │
 downstream 1D
```

Demonstrate:

Global:
- pressure,
- flow,
- pulse transmission.

Local:
- vortex structure,
- WSS,
- recirculation,
- residence-time-related quantities.

This should be the first scientifically representative multiscale showcase.

---

# 28. Required Immersed Benchmarks

Before heart anatomy:

1. planar/cubic embedded domain integration,
2. immersed straight channel/tube,
3. curved static vessel,
4. sphere/cavity-like domain,
5. irregular static chamber,
6. immersed aneurysm geometry.

Verify convergence with background refinement.

Only then attempt ventricular geometry.

---

# 29. Required Moving-Boundary Benchmarks

Before patient-specific heart:

1. oscillating wall/channel,
2. prescribed expanding/contracting cavity,
3. pulsatile chamber with known volume change,
4. mass balance versus boundary displacement,
5. moving curved vessel.

Then:

```text
prescribed LV wall motion
```

---

# 30. Implementation Roadmap

The phases below are ordered by architectural dependency and scientific risk.

---

## Phase 0 — Preserve baseline and create runtime boundaries

### Goal

Make existing solvers reusable without changing numerical results.

### PR 0.1 — Baseline audit and validation manifest

Tasks:

- enumerate existing validated examples,
- record expected test commands,
- record current database/schema compatibility,
- document current 1D/3D/VCA runtime state,
- add missing lightweight regression harness if needed.

Acceptance:

- clean clone can run the agreed baseline validation set.

### PR 0.2 — Common coupling data types

Add:

```text
PortState
CouplingPort
PortBoundaryData
PortOrientation
CouplingResidual
```

Do not yet modify solver behavior.

Acceptance:

- unit tests for signs/serialization/config validation.

### PR 0.3 — Generalize current 3D port measurement

Refactor existing `FlowPortMeasurements` functionality behind a reusable port API.

Preserve VCA behavior.

Acceptance:

- VCA tests unchanged,
- 3D flow/pressure/species measurement parity.

### PR 0.4 — 3D trial/rollback/commit runtime API

Wrap/refactor `TransientFlowRuntime`.

Need:

```text
BeginStep
SetPortInput
SolveTrial
RollbackTrial
CommitStep
GetPortState
```

Acceptance:

- standalone flow trajectory unchanged,
- repeated rollback/resolve gives identical state,
- no hidden time advancement.

### PR 0.5 — Extract reusable 1D runtime

Refactor the 1D executable so solver logic can run as a library/runtime object.

CLI remains a thin driver.

Need:

- initialization,
- port boundary override,
- trial advance,
- rollback,
- commit,
- port measurement.

Acceptance:

- all existing 1D examples unchanged,
- checkpoint/restart remains valid,
- no duplicated 1D solver implementation.

---

## Phase 1 — First 1D–3D coupling

### Goal

Run one body-fitted 3D region inside a 1D system.

### PR 1.1 — Explicit staggered straight-vessel prototype

Implement:

```text
1D → 3D → 1D
```

with simplest stable coupling.

Purpose:

- validate port mapping,
- validate signs,
- validate state management,
- validate units.

Not production-complete.

### PR 1.2 — Strong partitioned coupling

Add:

- iteration loop,
- convergence metrics,
- fixed relaxation,
- failure diagnostics.

### PR 1.3 — Aitken relaxation

Add dynamic relaxation and residual histories.

### PR 1.4 — Time subcycling

Permit 1D internal subcycling under one 3D macro-step.

### Phase 1 exit criteria

Straight-vessel benchmark passes.

Required diagnostics:

- P/Q time series,
- interface residual,
- mass imbalance,
- iteration count,
- relaxation factor.

---

## Phase 2 — Multidomain graph

### Goal

Remove the assumption of exactly two coupled solvers.

### PR 2.1 — SimulationGraph

Implement domain registry and coupling edges.

### PR 2.2 — Schema-v5 prototype

Add `domains` and `couplings`.

Keep v3/v4 backward compatibility.

### PR 2.3 — Multiple coupled interfaces

Support:

```text
1D → 3D → 1D
```

with independent upstream/downstream interfaces.

### PR 2.4 — 3D bifurcation benchmark

Support one 3D junction attached to multiple 1D branches.

### PR 2.5 — Multiple 3D islands

Support at least two separate 3D regions in one global simulation graph.

### Phase 2 exit criteria

A single case can contain arbitrary 1D and body-fitted 3D domain nodes connected by explicit coupling edges.

---

## Phase 3 — 1D–3D species transport

### Goal

Conserve transported quantities across dimensional interfaces.

### PR 3.1 — General species port data

Add concentration and flux to graph coupling.

### PR 3.2 — Conservative interface transfer

Implement 1D↔3D species exchange.

### PR 3.3 — Multispecies test

At least two independent scalar species across a 1D→3D→1D case.

### Exit criteria

Per-species global balance is reported and verified.

---

## Phase 4 — Quadrature abstraction

### Goal

Remove full-cell quadrature from the Navier–Stokes/transport kernel assumptions.

### PR 4.1 — `QuadratureRule` API

Refactor full tensor quadrature through generic providers.

### PR 4.2 — Navier–Stokes parity

Existing full-cell flow results remain equivalent.

### PR 4.3 — Transport parity

Existing transport results remain equivalent.

### Exit criteria

Physics kernels can integrate over arbitrary provided quadrature points without knowing how they were generated.

---

## Phase 5 — Background B-spline + immersed static geometry

### Goal

Support arbitrary closed surfaces without body-fitted volumetric meshing.

### PR 5.1 — Cartesian cubic B-spline background

Generate regular background spline topology/extraction.

### PR 5.2 — Surface reader and validator

Initial STL/VTP support.

### PR 5.3 — BVH and cell classification

Inside / outside / cut.

### PR 5.4 — Adaptive cut-cell volume quadrature

Octree reference implementation.

### PR 5.5 — Immersed surface quadrature

Generate boundary quadrature with physical points/normals.

### PR 5.6 — Nitsche wall boundary

Static no-slip immersed wall.

### PR 5.7 — Cut-cell stabilization

Ghost penalty or selected consistent stabilization.

### Phase 5 exit criteria

Static arbitrary surface → background cubic B-spline → converged incompressible flow.

---

## Phase 6 — 1D + immersed 3D integration

### Goal

Replace the body-fitted 3D island with an immersed 3D island without changing the coupling layer.

Required architecture proof:

```text
CouplingGraph
      │
      ├── OneDFlowDomain
      │
      └── ThreeDImmersedFlowDomain
```

No coupling algorithm should need to know that the 3D discretization changed.

### Benchmark

```text
1D → immersed 3D aneurysm → 1D
```

---

## Phase 7 — Prescribed moving anatomy

### Goal

Support heart-like time-dependent domains without two-way FSI.

Add:

- moving surface representation,
- geometry interpolation in time,
- moving cut-cell update,
- prescribed wall velocity,
- moving-interface diagnostics.

### First cardiac target

```text
patient-specific or idealized LV
with prescribed wall motion
```

Outputs:

- 3D velocity/pressure,
- vortex metrics,
- washout/residence-time-related metrics,
- mass conservation.

---

## Phase 8 — Structural mechanics and FSI

### Goal

Introduce two-way fluid-structure coupling.

Recommended order:

1. simple deformable wall benchmark,
2. compliant tube FSI,
3. aneurysm wall FSI,
4. thin-shell valve benchmark,
5. valve opening/closing,
6. leaflet contact,
7. patient-specific valve.

Reuse strong coupling/Aitken infrastructure from 1D–3D.

---

## Phase 9 — Whole-circulation multiscale cardiovascular model

Target:

```text
                 0D heart/circulation
                         │
                    3D heart ROI
                         │
                   3D aortic root
                         │
                  systemic 1D tree
                  /      |       \
             3D ROI   1D branch  3D ROI
                 \        |        /
                      terminal 0D
```

This is the long-term digital-twin architecture.

---

# 31. Suggested Repository Organization

Do not reorganize everything immediately. Migrate incrementally.

Target conceptual organization:

```text
core/
  coupling/
    PortState.hpp
    CouplingPort.hpp
    CouplingEdge.hpp
    CouplingAlgorithm.hpp
    AitkenRelaxation.hpp
    SimulationGraph.hpp
  time/
    CoupledTimeIntegrator.hpp

geometry/
  tubular/
    ...
  surface/
    SurfaceMesh.hpp
    BVH.hpp
    CellClassifier.hpp
  immersed/
    CutCell.hpp
    CutQuadrature.hpp

discretization/
  spline/
    ...
  background_bspline/
    ...

physics/
  flow/
    NavierStokesElement.hpp
  transport/
    ...
  structure/
    ... future

domains/
  one_d/
    OneDFlowDomain.hpp
  three_d/
    ThreeDFlowDomain.hpp
  immersed_three_d/
    ImmersedFlowDomain.hpp
  zero_d/
    ZeroDCircuitDomain.hpp

solvers/
  one_d/
    existing CLI
  cpu/
    existing CLI
  cuda/
    existing CLI
```

This is a direction, not a required immediate directory rename.

Avoid mass file movement before functionality is stable because it makes review and regression diagnosis harder.

---

# 32. Numerical / Software Diagnostics Required in Output

Every coupled run should eventually write a machine-readable coupling manifest.

Suggested fields:

```json
{
  "step": 120,
  "time": 0.120,
  "coupling_iterations": 5,
  "interfaces": [
    {
      "id": "aorta_to_roi",
      "flow_a": 1.23e-5,
      "flow_b": -1.23e-5,
      "pressure_a": 12000.0,
      "pressure_b": 11998.0,
      "relative_flow_residual": 2.0e-6,
      "relative_pressure_update": 8.0e-6,
      "relaxation": 0.61
    }
  ]
}
```

Also track:

- cumulative mass imbalance,
- cumulative species imbalance,
- maximum coupling iterations,
- mean coupling iterations,
- failed trial solves,
- nonlinear iterations,
- KSP iterations,
- wall-clock time per subsystem.

This is essential for later performance research and V&V.

---

# 33. Performance Strategy

Correctness first, then performance.

## 33.1 Near term

- Reuse PETSc objects across trial solves where mathematically valid.
- Avoid rereading geometry/database during each coupling iteration.
- Keep rollback states in memory.
- Reuse 1D network topology.
- Cache port mappings.

## 33.2 Immersed geometry

Cache:

- background spline topology,
- BVH,
- static cut-cell classification,
- static quadrature.

For moving surfaces:

- only update cells affected by motion where practical,
- cache reusable topology,
- later investigate incremental BVH updates.

## 33.3 Parallelism

Initial coupling driver may execute subsystems sequentially.

Later possibilities:

- MPI communicator splitting,
- concurrent independent 3D islands,
- concurrent graph branches,
- GPU 3D + CPU 1D execution.

Do not introduce distributed concurrency before functional coupling is verified.

---

# 34. CUDA Strategy

Do not block the architecture on CUDA parity.

Priority order:

1. CPU reference implementation.
2. Numerical verification.
3. CUDA lowering only after the API is stable.

The coupling abstraction itself should be backend-independent.

A future 3D domain may choose:

```text
cpu_petsc
cuda
```

behind the same subsystem/port API.

Immersed CUDA work should be a separate optimization phase.

---

# 35. What Not to Do

The Codex agent should explicitly avoid the following.

## Do not try to make sweeping handle a heart

Do not spend major effort forcing a centerline parameterization onto:

- ventricles,
- atria,
- valve structures,
- general chambers.

Sweeping remains a specialized tubular backend.

## Do not replace IGA with generic tetrahedral FEM merely because anatomy is complex

The architecture goal is to remove body-fitted spline-meshing dependence through immersed IGA, not discard the spline solver.

A future FEM backend may be possible but is not the current objective.

## Do not implement FSI before the coupling infrastructure

FSI requires the same state-management / strong-coupling ideas at higher complexity.

## Do not combine 1D–3D coupling and immersed geometry in the first implementation

First verify:

```text
1D → existing body-fitted 3D → 1D
```

Then replace the 3D geometry backend.

## Do not change all schemas at once

Introduce graph schema only when a runnable multidomain case exists.

## Do not hard-code cardiovascular species

Transport remains generic.

## Do not assume interface pressure equivalence without documenting the chosen definition

Mean static pressure, normal traction, and total pressure are not interchangeable in every situation.

## Do not silently remesh / alter input anatomy

Geometry transformations must be explicit and recorded.

---

# 36. Recommended First Codex Work Package

**Do only this package first. Do not start Phase 1 implementation until this is reviewed.**

## Objective

Prepare TubularFlowIGA for direct 1D–3D coupling with **no intended change to existing numerical results**.

## Tasks

1. Audit the current 1D runtime architecture.
2. Audit `TransientFlowRuntime` and the VCA port path.
3. Identify all mutable state that must participate in trial/rollback/commit.
4. Propose exact C++ interfaces for:
   - `PortState`
   - `CouplingPort`
   - `PortBoundaryData`
   - `CoupledSubsystem`
5. Identify the smallest refactor that exposes the 3D runtime as a reusable subsystem.
6. Identify the smallest refactor that exposes the 1D runtime as a reusable subsystem.
7. Define port sign conventions.
8. Define SI interface units.
9. Define the straight-vessel 1D→3D→1D benchmark configuration.
10. Add architecture documentation and tests for the data types.
11. Preserve every existing example and executable.
12. Do **not** implement immersed geometry, new meshing, or FSI.

## Required deliverable before coding the full coupling

Create/update an architecture document containing:

- current state ownership,
- proposed state ownership,
- subsystem lifecycle,
- coupling iteration sequence,
- rollback semantics,
- port sign convention,
- backward compatibility risks,
- exact files expected to change.

Then implement only the minimal Phase 0 refactor.

---

# 37. Suggested PR Sequence Summary

```text
PR 0.1  Baseline validation manifest
PR 0.2  Common port/coupling data types
PR 0.3  Generalize 3D port measurement
PR 0.4  3D trial / rollback / commit
PR 0.5  Reusable 1D runtime

PR 1.1  Explicit 1D→3D→1D straight vessel
PR 1.2  Strong fixed-point coupling
PR 1.3  Aitken relaxation
PR 1.4  Multirate/subcycling

PR 2.1  SimulationGraph
PR 2.2  Multidomain schema
PR 2.3  Multi-interface coupling
PR 2.4  3D bifurcation + multiple 1D branches
PR 2.5  Multiple 3D islands

PR 3.x  Conservative 1D–3D species coupling

PR 4.x  Generic quadrature abstraction

PR 5.x  Background B-spline + surface + cut cells
        + Nitsche + cut-cell stabilization

PR 6.x  1D + immersed 3D

PR 7.x  Prescribed moving anatomy

PR 8.x  Structural solver + FSI

PR 9.x  Whole-circulation multiscale cases
```

---

# 38. Definition of the Major Milestones

## Milestone M1 — Multiscale core

```text
1D → body-fitted 3D → 1D
```

- pulsatile,
- conservative,
- strongly coupled,
- rollback-safe,
- validated.

This is the **first major milestone**.

## Milestone M2 — Arbitrary heterogeneous graph

```text
multiple 1D / 3D / 0D domains
```

with reusable coupling ports.

## Milestone M3 — General anatomy

```text
STL/VTP
→ immersed cubic B-spline
→ static 3D CFD
```

## Milestone M4 — Multiscale general anatomy

```text
1D → immersed 3D anatomy → 1D
```

## Milestone M5 — Beating anatomy

```text
1D / 0D
↔ prescribed moving 3D heart
↔ 1D systemic network
```

## Milestone M6 — FSI

```text
blood ↔ deformable wall / valve
```

## Milestone M7 — Whole cardiovascular multiphysics platform

```text
0D + 1D + 3D
+ transport
+ moving anatomy
+ FSI
```

---

# 39. Research / Publication Opportunities Created by the Architecture

The software roadmap should not be driven only by feature count. It creates separable research questions.

Potential research layers:

1. **Robust multiscale 1D–3D coupling for spline-based vascular CFD**
   - interface conditions,
   - strong coupling,
   - wave reflection,
   - multirate integration.

2. **Immersed spline CFD for arbitrary biomedical anatomy**
   - cut-cell treatment,
   - boundary accuracy,
   - stabilization,
   - efficiency.

3. **Multiscale immersed cardiovascular simulation**
   - 1D systemic network + local 3D anatomy.

4. **Moving-domain cardiac CFD**
   - patient-specific prescribed wall motion,
   - washout,
   - vortical dynamics.

5. **Immersed IGA FSI**
   - compliant wall,
   - valves,
   - contact.

6. **Automated fidelity selection / simulation agents**
   - choosing which network regions deserve 3D,
   - adaptive 1D↔3D model allocation,
   - cost-aware V&V.

The implementation should preserve diagnostics and provenance so these scientific questions can be studied later.

---

# 40. Final Target Architecture

```text
                           INPUT
                             │
            ┌────────────────┼────────────────┐
            │                │                │
       Centerline         Surface       Moving Surface
       SWC / OBJ         STL / VTP        Γ(t)
            │                │                │
            ↓                └───────┬────────┘
      Tubular Sweep                  ↓
            │                  Immersed Geometry
            ↓                        │
   Body-fitted spline       Background B-spline
            │                        │
            └───────────┬────────────┘
                        │
                        ▼
                  3D Flow Domain
               Navier–Stokes / VMS
                        │
            ┌───────────┼─────────────┐
            │           │             │
        Transport     moving wall     FSI
            │                         │
            └───────────┬─────────────┘
                        │
                     Ports
                        │
        ┌───────────────┼────────────────┐
        │               │                │
      0D domain       1D domain       other 3D
        │               │                │
        └──────────── SimulationGraph ────┘
                        │
                        ▼
              Coupled Time Integrator
          strong iteration / Aitken / subcycling
                        │
                        ▼
               Multiscale Simulation
```

The key architectural transition is:

> **Do not keep asking how to generate a more complicated swept volumetric mesh. Make the solver framework stop requiring one geometry/discretization strategy for every anatomical region.**

And the key implementation transition is:

> **Before adding immersed anatomy or FSI, convert the existing 1D and 3D solvers into reusable rollback-safe subsystems connected by conservative ports.**

That sequence minimizes technical risk and allows each scientific capability to be verified independently.

---

# 41. Codex Agent Operating Strategy

For this repository, model choice should depend on the task.

## Default implementation work

Use **GPT-5.6 Terra with High reasoning effort** for:

- normal refactoring,
- adding well-specified interfaces,
- tests,
- config/schema plumbing,
- moving code into reusable runtime components,
- routine PETSc/C++ integration,
- documentation,
- straightforward benchmark setup.

This should be the default because the project contains many long but not individually frontier-hard implementation tasks.

## Architecture / numerically sensitive work

Use **GPT-5.6 Sol with High reasoning effort** for:

- the initial subsystem/coupling architecture,
- state ownership and rollback semantics,
- pressure/flow interface formulation,
- strong partitioned coupling,
- Aitken implementation review,
- quadrature abstraction,
- Nitsche formulation,
- cut-cell stabilization,
- moving-boundary formulation,
- FSI coupling,
- difficult nonlinear/PETSc debugging,
- major code reviews before merging a milestone.

## Max reasoning

Use **Sol Max selectively**, not as the default.

Good uses:

- a coupling algorithm fails in a subtle way,
- conservation is violated despite apparently correct code,
- numerical stability depends on interacting discretizations,
- a large refactor crosses 1D/3D/runtime/config/PETSc ownership boundaries,
- immersed integration or FSI results are difficult to diagnose,
- reviewing a milestone before committing to the next architecture.

Avoid Max for:

- adding boilerplate tests,
- renaming/refactoring obvious types,
- formatting/config updates,
- simple CLI plumbing,
- repetitive documentation.

## Ultra

Do not use Ultra as the default development mode.

Reserve it for a rare architecture/research review where parallel subagents can independently inspect:

- numerical formulation,
- software architecture,
- tests/V&V,
- performance.

For ordinary implementation it is likely excessive.

## Recommended cost/performance policy

```text
~70–80% of implementation:
    Terra High

~15–25% architecture/review/hard debugging:
    Sol High

~5% or less:
    Sol Max

Ultra:
    exceptional milestone reviews only
```

If only one model can be used for an entire long run, use:

```text
Terra High
```

for a tightly scoped PR with explicit acceptance criteria.

If giving one agent broad autonomy over a whole phase with architectural freedom, use:

```text
Sol High
```

and still restrict it to one phase at a time.

Never ask an autonomous agent to implement the entire roadmap in one run.
