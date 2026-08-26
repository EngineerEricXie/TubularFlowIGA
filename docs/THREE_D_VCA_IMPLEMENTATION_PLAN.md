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

Implemented in `feature/native-vca-vascular-integration`:

- Shared SI VCA contracts, external-circuit model, replay reader, oxygen
  capacity conversion, and 1D closed loop.
- Explicit 3D port declaration:
  `coupling.three_d_ports.inlet_label` and `outlet_labels`.
- CPU backward-Euler VCA loop in `iga_navier_stokes`: it applies the arterial
  inlet state, solves flow, solves one in-memory transport system using the
  current flow velocity, reduces outlet flow/pressure/species fluxes across
  MPI, advances the external circuit, and writes `coupling_manifest.json`.
- `TransientTransportRuntime`, including gathered scalar state, total-mass and
  volume-source diagnostics, and PETSc binary state import/export.
- Signed outward species-flux quadrature in `BoundaryFlow.hpp`, dynamic scalar
  inlet materialization, and port-result construction in
  `ThreeDVcaCoupling.hpp`.
- A coupled checkpoint envelope: the existing flow checkpoint retains flow and
  outlet state; `<prefix>.vca.json` retains step/time/dt, transport field
  order, and reservoir state; `<prefix>.vca_transport.state` retains the
  transport PETSc vector.
- Source-only bifurcation VCA case in `examples/vascular_flow/vca_bifurcation`.
  Generated meshes, databases, partitions, results, and checkpoints are
  intentionally not committed.
- Unit coverage for signed fluxes, VCA configuration checks, checkpoint
  metadata, and a one-element/single-rank PETSc transport runtime including
  state restore.
- A one-element/two-port end-to-end smoke test that runs flow, transport,
  external-circuit advance, manifest generation, checkpoint, and restart. Its
  uninterrupted and resumed single-rank flow/transport PETSc state files
  compare byte-for-byte; the one- and two-rank velocity/pressure fields and
  reservoir state agree to `1e-10` relative tolerance.

Not yet demonstrated or intentionally unsupported:

- A completed coupled flow-and-transport smoke run through the source
  bifurcation case on two ranks. The locally prepared case is large enough
  that it must be run on an allocated CPU resource, not a shared login node.
- One-versus-two-rank agreement for the physical source case, passive-tracer conservation, oxygenator
  source accounting, and uninterrupted-versus-restart equivalence for the
  coupled 3D runner.
- 3D VCA replay/open-loop modes, CUDA VCA, more than one transport system,
  and a transport source other than the current in-memory flow state.
- Flow-only VCA checkpoint/restart. It remains rejected because there is no
  scalar state to make the checkpoint contract uniform.

The species-free flow bridge remains available, but reservoir species require
the transport path above. Do not extend the accepted circuit/device models
without adding the corresponding conservation test.

## Required refactor

### 1. Extract a reusable transient flow runtime — implemented

The step-owned portions of `solvers/cpu/src/iga_navier_stokes.cpp` now live in
the header-only `solvers/cpu/include/TransientFlowRuntime.hpp`, matching the
existing transport runtime's deployment style. It owns PETSc state and ghost
scatter, Newton/KSP work, R/RC/RCR fixed-point state, VCA port measurements,
velocity gather, and final norms. The CLI retains option parsing, VCA circuit
and transport ordering, output, and checkpoint metadata.

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

### 2. Extract a reusable transient transport runtime — implemented

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

### 3. Add boundary species-flux integration — implemented

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

### 4. Add the coupled runner — implemented in `iga_navier_stokes`

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

## Files changed and remaining modifications

| File | Status and follow-up |
| --- | --- |
| `solvers/cpu/src/iga_navier_stokes.cpp` | **Changed.** Retains CLI, VCA circuit/transport ordering, output, and checkpoint metadata; flow step work delegates to `TransientFlowRuntime`. |
| `solvers/cpu/include/TransientFlowRuntime.hpp` | **Changed.** Owns transient Navier--Stokes state, ghost scatter, Newton/KSP, outlet fixed-point coupling, port measurements, velocity gather, and summary norms. |
| `solvers/cpu/include/TransientTransportRuntime.hpp` | **Changed.** One-step in-memory scalar solver. MPI total-mass and source diagnostics integrate only element owners. Next: replace its correctness-first all-rank velocity/state gather with shared ghost evaluation. |
| `solvers/cpu/include/BoundaryFlow.hpp` | **Changed.** Signed scalar flux and scalar/area integrals. |
| `include/ThreeDVcaCoupling.hpp` | **Changed.** Port and transport validation plus dynamic inlet/result helpers. |
| `include/CouplingHistory.hpp` | **Changed.** Generic manifest writer with the established explicit-staggered schema, units, port validity, and JSON escaping. |
| `solvers/one_d/include/OneDCoupling.hpp` | **Changed.** Uses the generic history writer; 1D-specific inlet/mass/port calculations remain local. |
| `include/VcaCheckpoint.hpp` | **Changed.** VCA metadata parser/validator. Next: add configuration/port fingerprint fields before treating restart as production-ready. |
| `solvers/cpu/tests/test_boundary_flow.cpp` | **Changed.** Covers positive and negative signed species flux. |
| `solvers/cpu/tests/test_simulation_config.cpp` | **Changed.** Covers ports, transport-system validation, and inlet-species update. |
| `solvers/cpu/tests/test_vca_checkpoint.cpp` | **Changed.** Metadata round-trip and mismatch rejection. |
| `solvers/cpu/tests/test_vca_3d_runtime.cpp` | **Changed.** One-element/single-rank PETSc transport, mass, and state restore. |
| `solvers/cpu/tests/test_vca_3d_smoke.cpp` | **Changed.** One-element/two-port flow+transport+reservoir smoke test; it compares uninterrupted/restarted state and one-/two-rank flow fields plus reservoir state. |
| `solvers/cpu/Makefile` | **Changed.** Adds `vca_checkpoint_test` and PETSc `vca_3d_runtime_test`; tracks new VCA headers. |
| `docs/PDE_CONFIGURATION.md` | **Changed.** Records supported CPU 3D VCA configuration and limitations. |
| `examples/vascular_flow/vca_bifurcation/` | **Changed.** Source-only two-outlet case and run instructions; do not commit generated data. |

## Checkpoint design — initial implementation complete

The current VCA checkpoint envelope contains:

- flow PETSc state and reconstructed prior flow state;
- the one supported transport PETSc state;
- R/RC/RCR outlet state in the flow checkpoint;
- reservoir volume, species concentrations, temperature, hematocrit, and the
  latest arterial state needed by devices;
- completed step, physical time, time step, schema version, field order, and
  transport-state filename.

The reader rejects mismatched step/time/dt/field order and never resets the
reservoir or species state silently. Before production use, add a configuration
fingerprint, port labels, transport-system name, and device-model identity to
the metadata validation.

## Handoff checklist and validation order

1. **Done locally:** unit tests for signed scalar flux, outlet aggregation,
   dynamic inlet species, VCA checkpoint metadata, one-element transport, and
   a one-element/two-port flow+transport+checkpoint/restart smoke run on one
   and two ranks. Run
   `make -C solvers/cpu test` and, with PETSc,
   `make -C solvers/cpu PETSC_DIR=... petsc-test`.
2. **Next required run:** execute the source bifurcation case on one and two
   allocated CPU ranks. Check port histories, field norms, and finite reservoir
   and balance values. Do not commit generated case data.
3. Passive tracer case: with no source and equal initial/inlet concentration,
   combined vascular plus reservoir tracer mass is constant to the selected
   tolerance.
4. Oxygenator case: reservoir/device oxygen source is reported separately and
   matches vascular-plus-reservoir mass change.
5. MPI: repeat the physical smoke case on one and two ranks; compare port histories,
   final reservoir state, and field norms.
6. Restart: repeat the one-rank restart equivalence test on two ranks.
7. Regression: run `make -C solvers/cpu test`, PETSc flow/transport tests,
   `make -C solvers/one_d core-test`, and `git diff --check`.

Do not claim full 3D VCA support until items 1--5 pass. Do not enable CUDA VCA
until it has the same explicit port, species-flux, reservoir, and restart
semantics.
