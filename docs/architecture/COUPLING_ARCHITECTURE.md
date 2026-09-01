# Coupling Architecture

Status: Phase 0 design baseline. No solver behavior is changed by this document.

This document fixes the runtime and interface contracts that precede direct
1D--3D coupling. The implementation remains incremental: existing standalone
executables continue to own configuration and output while reusable subsystem
objects acquire explicit trial, rollback, and commit semantics.

## Current runtime boundaries

### Three-dimensional flow

`solvers/cpu/src/iga_navier_stokes.cpp` owns configuration materialization,
physical time and step counting, VCA circuit state, transport orchestration,
checkpoints, and output. `TransientFlowRuntime` owns the PETSc flow objects,
resolved boundary data, outlet-model state, assembly, nonlinear solves, and
boundary measurements.

The current `Advance` method is not a trial solve. For `step > 0` it copies
`state_` to `previous_`, then mutates `state_`, resolved boundary data,
pressure tractions, outlet flow/pressure, capacitor pressure, and accumulated
linear-iteration diagnostics. Calling it repeatedly during one physical step
would therefore advance history incorrectly.

`MeasurePorts` already performs distributed outward-normal flow integration,
area-weighted pressure integration, and advective species-flux integration.
It is coupled to `ThreeDVascularPortDefinition`, exposes only VCA outlet
labels, and discards the integrated area after forming mean pressure. The
integration implementation is the basis of the generic 3D port adapter; it
must not be duplicated.

### One-dimensional flow and transport

`solvers/one_d/src/iga_1d.cpp` is both CLI and runtime orchestrator. It owns the
network, flow state, transport state, VCA/replay providers, dynamic-radius
updates, physical time and steps, checkpoints, histories, and output.
Numerical routines are already reusable header functions, but their mutable
objects are not collected behind a lifecycle boundary.

The native 1D inlet is flow-controlled. Outlet leaves use pressure,
resistance, or RCR closures. Direct 1D--3D coupling will eventually require
explicit per-port boundary overrides, but Phase 0 must not introduce a new
case schema or silently reinterpret existing boundaries.

### Existing VCA path

VCA is an explicit vascular-to-circuit bridge, not direct 1D--3D coupling.
The common implementation should reuse its SI validation and measurement
logic while keeping the generic port layer independent of reservoir-specific
types and the fixed inlet/outlet layout.

## State ownership and rollback inventory

The committed snapshot is the complete physical state at `t_n`. Solver
matrices, Krylov work vectors, ghost scatters, topology, and immutable
configuration are rebuildable or reusable scratch and are not physical state.

| Owner | Physical or trial state | Required handling |
|---|---|---|
| 3D flow runtime | PETSc `state_` | Snapshot at `BeginStep`; restore on rollback; accept once on commit. |
| 3D flow runtime | backward-Euler `previous_` | Set from the committed `state_` once per macro-step; never from a rejected trial. |
| 3D flow runtime | `outlet_models_`, including flow, pressure, and capacitor pressure | Snapshot and restore with the flow vector. |
| 3D flow runtime | `boundaries_` and `pressure_tractions_` | Treat materialized trial inputs as trial state; restore the committed values. |
| 3D flow runtime | linear/nonlinear diagnostic counters | Trial diagnostics remain attributable to the attempt; physical cumulative counters update only on commit. |
| 3D transport runtime | PETSc `current_` and logical `steps_` | Add to the coupled snapshot when transport enters the lifecycle; do not advance on rejected flow trials. |
| 3D CLI/VCA | `VcaExternalCircuit` state, last arterial species, and previous species mass | Advance exactly once after the coupled step converges. |
| 1D flow runtime | area, flow, pressure, nodal pressure, segment flow, inlet flow, outlet/RCR states | Copy as one committed `OneDFlowState` snapshot. |
| 1D flow runtime | physical time, completed step, internal substep count | Advance in the trial copy; publish only on commit. |
| 1D network | dynamic `radius0`, `area0`, and resistance from vasodilation | Snapshot when dynamic geometry is enabled. |
| 1D transport | every species concentration and mutable inlet value/waveform | Snapshot with the trial state. |
| 1D configuration | hematocrit/hemoglobin fields currently changed by `ApplyOneDCoupledInlet` | Move future trial boundary values into runtime-owned data or include the mutated fields in the snapshot. |
| CLI output | writers, coupling history, prior mass, checkpoints | Side effects occur after `CommitStep`, never during `SolveTrial`. |

Disk checkpoint/restart remains a persistence interface. It is not the
in-memory rollback mechanism.

## Common port contract

The first common types will be dependency-free C++17 values under `include/`.
They will not depend on PETSc, MPI, JSON parsing, VCA, or either solver.

Conceptually:

```cpp
enum class PortQuantity {
	Area,
	FlowRate,
	MeanPressure,
	MeanNormalTraction,
	TotalPressure,
	SpeciesConcentration,
	SpeciesFlux
};

struct PortOrientation {
	// Converts the dimension-specific positive direction to positive outward.
	int native_to_outward_sign = 1; // must be +1 or -1
};

struct PortState {
	double time_s = 0.0;
	double area_m2 = 0.0;
	double outward_flow_m3_s = 0.0;
	double mean_pressure_pa = 0.0;
	double mean_normal_traction_pa = 0.0;
	double total_pressure_pa = 0.0;
	std::map<std::string, double> concentration;
	std::map<std::string, double> outward_species_flux;
};

struct PortBoundaryData {
	double time_s = 0.0;
	std::optional<double> outward_flow_m3_s;
	std::optional<double> mean_pressure_pa;
	std::optional<double> mean_normal_traction_pa;
	std::map<std::string, double> concentration;
	std::map<std::string, double> inward_species_flux;
};

struct CouplingPort {
	std::string id;
	std::string subsystem_id;
	std::string locator_kind;
	std::string locator;
	PortOrientation orientation;
	std::set<PortQuantity> provides;
	std::set<PortQuantity> requires;
};
```

Field presence in `PortBoundaryData` is explicit. A numeric zero is valid data
and must not mean "not provided." Species concentration units are defined by
each configured species; species flux uses concentration unit times m^3/s.
The first implementation validates finite values, positive area when area is
provided, a sign of exactly `+1` or `-1`, unique non-empty identifiers, and no
quantity appearing in both `provides` and `requires` for the same direction.

## SI units and sign convention

All coupling interfaces use:

| Quantity | Unit |
|---|---|
| time | s |
| area | m^2 |
| flow rate | m^3/s |
| velocity | m/s |
| pressure and normal traction | Pa |
| species concentration | the species' explicitly configured concentration unit |
| species flux | concentration unit m^3/s |

`PortState::outward_flow_m3_s > 0` means material leaves the subsystem through
the port. In 3D this is `integral(u dot n dA)` with the geometric outward
normal. In a rooted 1D vessel, root-to-leaf flow is outward at a distal port
but inward at the root port. `PortOrientation::native_to_outward_sign`
performs only this local conversion.

A conservative edge joins ports A and B using

```text
Q_A,out + Q_B,out = 0
Phi_A,out(species) + Phi_B,out(species) = 0.
```

The edge applies the sign reversal explicitly. Boundary labels, graph order,
and inlet/outlet names never imply a sign.

## Subsystem lifecycle

The reusable interface is deliberately small:

```cpp
class CoupledSubsystem {
public:
	virtual ~CoupledSubsystem() = default;
	virtual void Initialize() = 0;
	virtual void BeginStep(double time_s, double dt_s) = 0;
	virtual void SetPortInput(const std::string&, const PortBoundaryData&) = 0;
	virtual void SolveTrial() = 0;
	virtual PortState GetPortState(const std::string&) const = 0;
	virtual void RollbackTrial() = 0;
	virtual void CommitStep() = 0;
};
```

The legal sequence is:

```text
committed -> BeginStep -> trial-ready
trial-ready -> SetPortInput -> SolveTrial -> trial-solved
trial-solved -> RollbackTrial -> trial-ready
trial-solved -> CommitStep -> committed
```

Invalid transitions fail clearly. `SolveTrial` always starts from the same
committed `t_n` snapshot, including internal 1D subcycling and terminal-model
state. `CommitStep` succeeds exactly once. A failed solve leaves a rollback
path and cannot emit output or advance external circuit state.

For fixed-point coupling, the driver performs:

```text
BeginStep(A, B)
  -> apply trial interface data
  -> SolveTrial(A), measure A
  -> SolveTrial(B), measure B
  -> compute dimensionless residuals and relax
  -> RollbackTrial(A, B) and repeat, or CommitStep(A, B)
```

## Smallest implementation boundaries

The 3D refactor will first be flow-only. A wrapper or additions to
`TransientFlowRuntime` will make snapshot ownership explicit and remove the
`step > 0` history copy from trial execution. Generic measurements will accept
an explicit set of boundary-label locators and return area as well as flow and
pressure. The VCA adapter will call that implementation so its output remains
unchanged. Transport rollback follows only when flow lifecycle parity passes.

The 1D refactor will collect existing orchestration into a runtime object
without duplicating numerical routines. The CLI becomes a thin driver only
after parity tests cover rigid, compliant, transport, checkpoint/restart, and
VCA paths. Dynamic network and configuration mutations are included in the
snapshot or moved into runtime-owned trial data.

## Straight-vessel benchmark definition

Phase 1 will compare two physically equivalent cases:

```text
reference: upstream 1D segment -> replacement 1D segment -> downstream 1D segment
coupled:   upstream 1D domain  -> body-fitted 3D straight tube -> downstream 1D domain
```

Both cases use identical SI length, radius, density, viscosity, wall law,
terminal model, inlet waveform, macro-step interval, and total axial length.
The 3D interface planes use their geometric outward normals. The benchmark
records P/Q waveforms at both interfaces, pressure drop, phase and pulse
transit, per-step mass imbalance, coupling residual histories, iteration
counts, and tolerance sensitivity. Phase 1 cannot start until Phase 0 proves
rollback determinism and standalone/VCA parity.

## Compatibility matrix

| Interface | Phase 0 requirement |
|---|---|
| 1D schema v3 | Parse and execute unchanged. |
| 3D schema v2/v3/v4 accepted by current readers | Parse and execute unchanged. |
| `.ntiga` databases | No format or reader change. |
| standalone `iga_1d` | Same CLI, numerical path, checkpoints, and output. |
| standalone `iga_navier_stokes` | Same CLI, numerical path, checkpoints, and output. |
| VCA replay/closed loop | Same signs, SI values, histories, and reservoir advancement. |
| CUDA | No Phase 0 behavior change; common interfaces remain backend-neutral. |

Schema v5 is deferred until a runnable multidomain graph exists. No file
format version is introduced for runtime-only coupling data.

## Validation manifest

Fast dependency-free baseline:

```bash
make cpu-test
make -C solvers/one_d core-test
```

PETSc baseline when the local PETSc environment is configured:

```bash
make cpu-petsc PETSC_DIR=/path/to/petsc PETSC_ARCH=arch
make -C solvers/cpu petsc-test PETSC_DIR=/path/to/petsc PETSC_ARCH=arch
make one-d-test PETSC_DIR=/path/to/petsc PETSC_ARCH=arch
```

Configuration checks after building `iga_1d`:

```bash
./solvers/one_d/iga_1d examples/one_d/rigid_straight --check
./solvers/one_d/iga_1d examples/one_d/compliant_bifurcation --check
./solvers/one_d/iga_1d examples/one_d/multispecies_physiology --check
./solvers/one_d/iga_1d examples/one_d/vca_pfc_closed_loop --check
```

The two-rank 3D VCA regression follows
`examples/vascular_flow/vca_bifurcation/README.md` on an allocated compute
resource. Large numerical cases do not run on a shared cluster login node.

## Planned Phase 0 touch list

PR 0.2 is limited to common value types and dependency-free tests:

- `include/CouplingPort.hpp` (new)
- `solvers/cpu/tests/test_coupling_port.cpp` (new shared test location)
- `solvers/cpu/Makefile`

Later 3D lifecycle/measurement work is expected to touch:

- `solvers/cpu/include/TransientFlowRuntime.hpp`
- `solvers/cpu/src/iga_navier_stokes.cpp`
- `include/ThreeDVcaCoupling.hpp`
- focused CPU/PETSc tests and Makefile dependencies

Later 1D runtime extraction is expected to touch:

- `solvers/one_d/include/OneDRuntime.hpp` (new)
- `solvers/one_d/src/iga_1d.cpp`
- `solvers/one_d/include/OneDCoupling.hpp`
- focused 1D tests and Makefile dependencies

Changes outside these boundaries require a documented reason and renewed
compatibility review.
