# Coupled checkpoint state contract

Bundle v1 stores accepted state for native 0D, 1D, and body-fitted 3D graph
runtimes, complete graph history, and restart metadata for a new MPI job.
Moving-flow and bounded moving-FSI library runtimes also support redistribution
by stable global ID. User-facing commands are documented in
[Native graph checkpoint and restart](../COUPLED_RESTART.md); the wire format is
documented in [Checkpoint bundle v1](COUPLED_CHECKPOINT_BUNDLE.md).

## 1. Checkpoint boundary

A checkpoint may be captured only after every domain has completed an accepted
macro-step. For accepted step count `N >= 1`, the following must have finished:

1. all trials and conservation checks;
2. the before-commit callback;
3. prepare and finalize for every domain;
4. executor donor publication;
5. accepted-result and next-pressure-guess updates; and
6. accepted clock advancement.

This boundary occurs at the end of the native graph runner loop, before the
next `BeginStep`. Capture is forbidden inside Newton or FSI iterations, between
hydraulic and transport phases, between prepare and finalize, or from a
before-commit callback. No unaccepted trial state is stored.

Every participating rank agrees on the domain catalog, epoch, accepted clock,
and idle lifecycle before capture. A scheduler signal sets a pending flag;
normal control flow services the flag at the next checkpoint boundary. Signal
handlers never call PETSc, MPI, or file I/O. If the job ends before that
boundary, restart uses the previous complete generation.

## 2. Time, identity, and representation

| Field | Contract |
|---|---|
| `accepted_macro_steps` | Integer `N`; conversions are range- and overflow-checked |
| `accepted_time_s` | The accumulated binary64 `EndTime()` value, never reconstructed as `N*dt` |
| `macro_dt_s` and configured end step | Must match the case; bundle v1 does not silently change the schedule |
| Next graph context | `step_index=N`, `start_time_s=accepted_time_s`; the first graph step uses index 0 |
| 0D clock | `CommittedStepIndex=N-1`, `CommittedStepCount=N`, and the accepted time |
| Immersed transient clock | Backend `index=N`, accepted time, and no active trial |
| 1D clock | Per-domain `completed_step`, `physical_time`, and `internal_substeps`, validated against configured subcycling |
| 3D transport clock | The runtime's actual `Steps()` value |
| Epoch | One unique generation shared by metadata, shards, and history, with the previous complete generation recorded |
| Model identity | Domain/edge/port IDs, directions, field order, units, parameters, coupling method/tolerances, and routing controls |
| Distribution identity | Stable communicator membership, rank map, owned IDs, layout, and partition identity; never an MPI handle |
| Execution identity | Schema/runtime versions, scalar/index widths, byte order, PETSc and numerical options, and source revision |

External inputs include `.ntiga` files and partitions, 1D networks, surfaces,
FSI reference/patch maps, waveforms, replay data, and velocity histories. Their
content SHA-256 and byte length establish identity; paths are only locators.
Effective numerical PETSc options are compatibility inputs. Output locations
and checkpoint frequency are operational metadata and do not affect numerical
identity.

Floating-point state must be finite and round-trip as binary64. Optional state
preserves presence, so an absent pressure is not converted to zero. Maps use
stable IDs; sequences declare their order and size. A loader checks bounds,
duplicates, missing IDs, enum values, finite values, and allocation limits
before allocating large arrays. Checksums do not replace physical and identity
validation.

## 3. Cross-step coupling state

| State | Save or reconstruct |
|---|---|
| Next pressure guess | **Save** `pressure[edge_id]` after it has been updated from the accepted measured pressure |
| Pressure-flow Aitken memory | **Reconstruct** because `PressureFlowComponentExecutor::Advance` initializes residual and relaxation state for each macro-step |
| Species hydraulic Aitken memory | **Reconstruct** for the same step-local reason |
| Species donor hysteresis | **Save** `CommittedDonorOwnership()` for every `(edge_id, species_id)` pair |
| Species routes and transport order | **Reconstruct** from the next hydraulic result plus restored donor ownership |
| Strong-FSI Aitken memory | **Reconstruct**; the next predictor uses committed displacement, velocity, `dt`, clamps, and reference normals |
| Accepted coupling history | **Save** as an immutable versioned prefix |

The accepted history includes iteration counts, residuals, applied/measured
pressures, relaxation, ports, species amounts and donors, 0D accounting, and
FSI accepted identities. A resumed output contains exactly one record for each
accepted step: the stored prefix followed by new records. Wall time, MPI time,
and memory use are measurements of each job and are not presented as a
continuous-run history.

## 4. Runtime payloads

A graph is checkpointable only when every required runtime provides a complete
capture and restore path. Missing runtime support rejects checkpoint creation;
the writer never publishes a partial graph as complete.

### 4.1 0D domains and outlets

[`ZeroDFlowDomain`](../../include/ZeroDFlowDomain.hpp) stores:

- `stored_pressure_pa`;
- committed index, count, and time;
- optional committed `PortState`; and
- `ZeroDFlowStepAccounting`, including initial/final storage, source, distal
  sink, integrated graph-port amount, and residual.

The model, role, ports, and model digest must match. Base, input, trial, and
staged states restore as inactive.

[`OutletCheckpoint`](../../include/OutletCheckpoint.hpp) stores the outlet
label/kind, flow, pressure, and capacitor pressure. Restore also validates all
resistance, capacitance, reference/distal pressure, and boundary mappings. RCR
capacitor pressure is dynamic state and cannot be replaced by a graph pressure
guess.

### 4.2 1D hydraulic and transport state

[`OneDRuntime`](../../solvers/one_d/include/OneDRuntime.hpp) is checkpointable
in its Ready phase after trial frames and rollback snapshots have been cleared.
The payload includes:

- complete `OneDFlowState`: cell area/flow/pressure, node pressure, segment
  flow, dynamic outlet values, inlet flow, completed step, time, and internal
  substeps;
- each `OneDSpeciesState`: concentration, root/outlet fluxes, validity flags,
  initial mass, boundary/source amounts, and accounting validity;
- dynamic inlet values and waveform state;
- complete `LastInlet()`, including optional flow/pressure, species,
  temperature, hematocrit, and metadata;
- accepted vasodilation radius and mutable physiology/perfusate scalars.

Network/node/segment/cell order is part of the input identity. Resistance and
reference area derived from an accepted dynamic radius are reconstructed with
the normal formula and validated; the baseline radius remains immutable.
Hydraulic frames, trial inputs, outlet overrides, and rollback caches are
recreated for the next step.

Standalone `OneDCheckpoint v2` does not contain the complete graph payload, so
it is not simply wrapped as a graph checkpoint. Graph providers use typed
accepted capture and fresh-candidate restore with a verified identity covering
the complete configuration, selected system, network, and external inputs.

### 4.3 Body-fitted 3D flow and transport

[`TransientFlowRuntime`](../../solvers/cpu/include/TransientFlowRuntime.hpp)
stores all owned `(u,v,w,p)` rows, dynamic outlets, accepted boundary controls
and pressure tractions, accumulated linear iterations, and the graph clock.
Backward Euler needs one accepted field: `BeginStep` recreates the committed
and previous snapshots from it. Steady graphs store the field for the next
Newton warm start.

[`TransientTransportRuntime`](../../solvers/cpu/include/TransientTransportRuntime.hpp)
stores the current vector for every ordered field and the true accepted step
count. Its `previous_` member is an assembled matrix, not a second
concentration-history vector. Matrices, KSPs, scatters, forcing, RHS, and trial
vectors are rebuilt.

Metadata and replicated boundaries are stored once per domain. Owned fields
are sharded by bundle world rank with bounded stream buffers. Fresh runtime
shape and ownership determine decoder expectations. A candidate is published
only after every payload, checksum, consumer finalization, runtime validation,
and communicator agreement succeeds.

Native body-fitted graph restart currently retains the same rank count because
`.ntiga` partitions are rank-specific.

### 4.4 VCA and external circuits

[`VcaExternalCircuit`](../../include/VascularCoupling.hpp) requires reservoir
volume, all species concentrations, temperature, hematocrit, and
`last_arterial_species_`. Pump, oxygenator, dialyzer, infusion, and oxygen-state
definitions belong to configuration identity.

The current dedicated VCA checkpoint formats do not provide a complete atomic
graph provider. Closed-loop 1D restart remains outside bundle-v1 graph support.
Future replay/open-loop providers must resume from accepted time using a
checksummed time series rather than a process-local file cursor.

### 4.5 Fixed immersed flow

`ImmersedDistributedNewtonRuntime` requires its owned accepted vector,
layout/partition, controller/gauge identity, port controls, and logical commit
count. Newton candidates, updates, prepared vectors, ghosts, matrices, and KSP
objects are rebuilt.

`ImmersedTransientDistributedRuntime` additionally stores committed time and
index. The next trial recreates history and freezes force/motion input for the
new interval. Callback state must be reproducible from verified model inputs;
arbitrary mutable process-local callback captures are not restartable.

The graph adapter's accepted controls, publications, and clock are part of the
payload. Trial inputs and rollback copies are not. Serial immersed runtimes
retain their own lifecycle state and cannot be represented by omitting fields
from the distributed schema.

### 4.6 Moving geometry

`MovingImmersedTransientFlowRuntime` stores the complete accepted
`ImmersedGlobalFlowState`: time, index, geometry identity, stable node IDs,
four field coefficients, port IDs/multipliers, gauge state, port controls, and
all committed conservation transition primitives.

Material-surface state includes reference vertices, current source vertices
and velocities, ordered triangle IDs/vertices/labels, canonical-corner
provenance, evaluation time, step interval, and material/topology/content/epoch
identities. The original reference and prescribed-motion model remain
configuration inputs.

Restore rebuilds cut geometry, quadrature, ghost, extension, layout, and
integration catalogs, then recomputes geometry identity. Publication identity
also includes predecessor geometry and old/new cell classifications. The
loader validates and reconstructs that predecessor relationship; it cannot
substitute a new genesis publication or accept an arbitrary hash override.

Moving fields can be redistributed by stable node ID when the target rank
count changes. The target receives a new partition identity after complete
source coverage is verified.

### 4.7 Structures and FSI

`PretensionedMembrane` stores scalar displacement and velocity indexed by
owned global surface-node IDs. Reference geometry, normals/areas, clamped mask,
material constants, layout, partition, and model identity are verified inputs.
Matrices, factorization, and trial handles are rebuilt.

The FSI pair additionally stores committed material kinematics, traction and
diagnostics, composition/patch identities, and the complete surface stamp.
The stamp includes time, step, accepted coupling iteration, reference/layout/
partition identities, and producer-state identities. Restore returns the trial
lifecycle to Idle and restores only the committed output.

Moving flow, rank-owned traction slices, and owner membrane state are published
under one manifest. The pair becomes visible only after every candidate and
clock validates. Surface ownership may be redistributed by stable node ID; the
restored target receives a new partition provenance. Native graph CLI wiring
for moving FSI is separate from this library-level bundle support.

## 5. Restore sequence

1. Read a complete manifest and validate schema, epoch, file sizes/checksums,
   and compatibility without changing a live runtime or output directory.
2. Construct models, target ownership, geometry, and PETSc objects for the
   compatible domain membership. Verify complete stable-ID coverage with no
   duplicate ownership or zero-filled gaps.
3. Load every domain state, clock, publication, pressure guess, donor map, and
   accepted-history prefix into an isolated candidate. Recompute identities and
   validate finite values, positivity, clocks/counters, donor coverage, FSI
   pairing, and the final history record.
4. Publish the complete candidate only after communicator-wide success.
   Otherwise destroy it and leave the bundle and output unchanged.
5. Resume with human-readable step `N+1`, producing no duplicate or missing
   accepted records.

Native 0D/1D/body-fitted graphs currently restore with the same rank
membership. Moving flow and bounded moving-FSI support explicit global-ID
redistribution. Unsupported runtime versions, backends, or rank transitions
are rejected before new shards are written.

Legacy standalone checkpoint formats keep their original identity and are not
relabeled as atomic graph bundles. CUDA raw-state checkpoints are not PETSc
vectors and cannot be loaded by a CPU graph provider. Visualization snapshots
are output products, not restart state.

## 6. Validation requirements

| Test | Required result |
|---|---|
| Invalid lifecycle capture | Trial, prepare, partial-finalize, staged-species, and active-FSI states reject capture without publishing a manifest |
| Full-history restart | A new MPI job resumes from nonzero `N`; accepted records, clocks, IDs, dynamic 0D/outlet state, 1D counters, and 3D fields/species match an uninterrupted run within their defined numerical tolerances |
| Coupling memory | Explicit, fixed-point, and Aitken cases distinguish the next pressure guess; near-zero species flow preserves donor hysteresis and donor reversal remains reproducible |
| Dynamic history | Tests include nonzero transient fields, RCR capacitor state, species accounting, reservoir composition, FSI displacement/velocity, and a real moving-cell transition |
| Publication failure | Termination before commit, during shard writes, or before manifest publication leaves only the previous generation loadable; truncated, corrupt, mixed-epoch, or incompatible bundles are rejected |
| Publication identity | Moving predecessor provenance, conservation transitions, committed FSI stamp, and next-step predictor reproduce without substituting a genesis identity |
| Scale and resources | Exercise 1/2/4 ranks and split communicators; record solve, I/O, synchronization, communication, and per-rank peak memory separately; use a scheduler for large or cross-node runs |

These gates define the validation scope of each checkpoint provider. The
corresponding regression target and recorded output establish whether a route
meets the contract.
