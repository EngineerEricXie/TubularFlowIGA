# Native 3D VCA species-coupling implementation plan

## Objective

Complete native CPU 3D `vca_closed_loop` support for transient
Navier--Stokes plus transported species. Each physical step must use one
consistent set of in-memory states:

```text
reservoir arterial state at t
  -> flow inlet profile and scalar inlet Dirichlet values
  -> Navier--Stokes step
  -> transport step using the new flow velocity
  -> MPI-reduced outlet flow, pressure, and species fluxes
  -> conservative reservoir/device advance to t + dt
```

The existing native 1D `iga_1d` implementation is the reference for port
contracts, explicit-staggered ordering, reservoir behaviour, and manifest
semantics. The local reference implementation is
`../HexSim/Hex_Sim_FEniCS` branch `feature/vca-vascular-integration`.

## Current state

Already implemented:

- Shared SI VCA contracts, external-circuit model, replay reader, oxygen
  capacity conversion, and 1D closed loop.
- Explicit 3D port declaration:
  `coupling.three_d_ports.inlet_label` and `outlet_labels`.
- CPU transient Navier--Stokes flow-only VCA bridge. It scales the inlet
  `initial_velocityfield.txt` profile, aggregates outlet flow and
  area-averaged pressure across MPI ranks, advances a species-free reservoir,
  and writes `coupling_manifest.json`.
- 3D transport contract checks and dynamic VCA inlet-species materialization
  in `include/ThreeDVcaCoupling.hpp`.

Not yet implemented:

- In-process 3D scalar transport step owned by the same VCA runner as flow.
- Outlet species flux integration and reservoir species advance from 3D data.
- VCA checkpoint/restart, replay/open-loop execution in 3D, CUDA support, and
  a committed packed 3D VCA smoke case.

The current flow-only path deliberately rejects reservoir species, oxygenator,
dialyzer, and infusion rates. Do not remove that protection until the species
loop below is complete.

## Required refactor

### 1. Extract a reusable transient flow runtime

Move the step-owned portions of
`solvers/cpu/src/iga_navier_stokes.cpp` into a new CPU header/source pair,
for example `solvers/cpu/include/TransientFlowRuntime.hpp` and
`solvers/cpu/src/TransientFlowRuntime.cpp`.

The runtime should own:

- flow `OwnedRowAssembler`, PETSc state/previous/update/RHS/Jacobian, KSP,
  ghost scatter, and outlet-model state;
- `Advance(step_configuration, physical_time)`, which performs exactly one
  nonlinear flow step and returns a distributed flow state;
- a port-measurement method returning MPI-global outlet flow and
  area-averaged pressure for configured labels;
- checkpoint import/export hooks, but do not implement VCA restart until all
  state types below are serialized.

Keep CLI parsing, output scheduling, and final summaries in
`iga_navier_stokes.cpp`. Existing non-VCA and R/RC/RCR behaviours must retain
their current code path and numerical results.

### 2. Extract a reusable transient transport runtime

Move the core of `solvers/cpu/src/iga_solve.cpp` into a corresponding
`TransientTransportRuntime`.

The runtime should own:

- transport assembler, `left`/`previous` matrices, forcing/current/next/RHS
  vectors, KSP, field order, and fixed Dirichlet rows;
- `Advance(step_configuration, velocity_at_new_time)`, which reassembles
  advection-dependent operators, materializes scalar boundaries, and solves
  exactly one backward-Euler step;
- a method to expose a ghosted or gathered scalar state for port integration;
- total species mass and source-integral diagnostics.

Initial implementation constraints are acceptable and should be explicit:

- exactly one `linear_transport` system;
- its `time.dt` and `time.steps` equal the transient flow system;
- no independent velocity snapshot source in VCA mode; velocity must come
  directly from the current in-memory flow state;
- all reservoir species are transported fields and have inlet Dirichlet
  conditions at the declared VCA inlet.

Avoid a per-step text/VTK velocity round-trip. It is non-conservative, slow,
and breaks MPI ownership. A temporary all-rank PETSc gather is acceptable only
as a clearly marked correctness-first implementation; replace it with a
shared ghost-vector evaluation before performance validation.

### 3. Add boundary species-flux integration

Extend `solvers/cpu/include/BoundaryFlow.hpp` with a quadrature helper for
the outward-normal advective species flux:

```text
integral_boundary concentration * (velocity dot n) dA
```

It must use the same face orientation and quadrature as
`IntegrateBoundaryFlow`. For every configured VCA outlet, integrate each
transported reservoir species on locally owned faces, then use `MPI_Allreduce`.

Create `VascularOutletState` entries with:

- positive outward flow convention;
- area-averaged pressure from the flow state;
- `species_flux[field]` from the boundary integral;
- flux-weighted concentration only when `abs(flow)` exceeds
  `coupling.flow_epsilon_m3_s`.

Use `AggregateVascularOutlets` without changing its signed conservation rules.

### 4. Add the coupled runner

Either add `iga_vca_3d` or make `iga_navier_stokes` dispatch to a dedicated
`RunThreeDVcaClosedLoop(...)` function. A dedicated function is preferred so
the ordinary flow solver remains simple.

Per step, it must:

1. Call `VcaExternalCircuit::InletState(t)`.
2. Materialize normal temporal functions at the new flow time.
3. Call `ApplyThreeDVascularInlet` to set the pump-flow velocity scale.
4. Call `ApplyThreeDVascularSpeciesInlet` to set scalar inlet values.
5. Advance flow.
6. Advance transport with that flow velocity.
7. Measure every outlet and construct `VascularStepResult`, including total
   species mass, source integrals, and balance residuals.
8. Aggregate outlets and call `VcaExternalCircuit::Advance`.
9. Append all three records to `CouplingHistoryWriter`.

Maintain the established explicit lag: the vascular step consumes arterial
state at `t`, and the reservoir consumes the resulting venous return after
that step.

## Files expected to change

| File | Change |
| --- | --- |
| `solvers/cpu/src/iga_navier_stokes.cpp` | Retain CLI; dispatch VCA work to the coupled runner. |
| `solvers/cpu/src/iga_solve.cpp` | Retain CLI; use extracted transport runtime. |
| `solvers/cpu/include/TransientFlowRuntime.hpp` | New reusable one-step flow runtime. |
| `solvers/cpu/include/TransientTransportRuntime.hpp` | New reusable one-step scalar runtime. |
| `solvers/cpu/include/BoundaryFlow.hpp` | Add species-flux quadrature helper. |
| `include/ThreeDVcaCoupling.hpp` | Keep port/config checks; add result construction and diagnostics helpers as needed. |
| `include/CouplingHistory.hpp` | Extract the generic manifest writer from the 1D-only header before the 3D runner uses it. |
| `solvers/one_d/include/OneDCoupling.hpp` | Include the generic history header after extraction; preserve 1D output format. |
| `solvers/cpu/tests/test_boundary_flow.cpp` | Unit-test signed scalar boundary flux and flux-weighted concentration. |
| `solvers/cpu/tests/test_simulation_config.cpp` | Test VCA transport-system and inlet-boundary validation. |
| `solvers/cpu/tests/test_vca_3d_runtime.cpp` | New small runtime test using a packed one-element/two-outlet fixture. |
| `solvers/cpu/Makefile` | Build new runtime test and any new executable/source. |
| `docs/PDE_CONFIGURATION.md` | Replace the flow-only limitation with supported 3D species configuration and remaining limitations. |
| `examples/vascular_flow/...` | Add a small two-outlet VCA case, source geometry, config, and run instructions; never commit generated `.ntiga` or result files. |

## Checkpoint design (after the coupled runtime works)

Extend the flow checkpoint metadata or introduce a VCA checkpoint envelope that
contains:

- flow PETSc state and prior flow state;
- every transport PETSc state;
- R/RC/RCR outlet state;
- reservoir volume, species concentrations, temperature, hematocrit, and the
  latest arterial state needed by devices;
- completed step, physical time, time step, schema/config fingerprint, port
  labels, field order, and transport-system name.

Reject a restart when any of these identities differ. Do not silently restart
only flow while resetting reservoir or species state.

## Validation order and acceptance criteria

1. Unit tests: signed scalar flux, outlet aggregation, dynamic inlet species,
   and failure cases for missing ports/fields/boundaries.
2. One-element/two-outlet transient smoke case: outlet flows sum to the pump
   flow within the incompressibility tolerance; reservoir volume remains
   bounded.
3. Passive tracer case: with no source and equal initial/inlet concentration,
   combined vascular plus reservoir tracer mass is constant to the selected
   tolerance.
4. Oxygenator case: reservoir/device oxygen source is reported separately and
   matches vascular-plus-reservoir mass change.
5. MPI: repeat the smoke case on one and two ranks; compare port histories,
   final reservoir state, and field norms.
6. Restart: compare uninterrupted and checkpoint/resumed coupled runs.
7. Regression: run `make -C solvers/cpu test`, PETSc flow/transport tests,
   `make -C solvers/one_d core-test`, and `git diff --check`.

Do not claim full 3D VCA support until items 1--5 pass. Do not enable CUDA VCA
until it has the same explicit port, species-flux, reservoir, and restart
semantics.
