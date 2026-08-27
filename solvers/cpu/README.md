# TubularFlowIGA CPU Backend

A project-owned C++17 isogeometric-analysis pipeline for stabilized steady and
backward-Euler Navier-Stokes flow plus configurable transient multi-field transport. PETSc provides distributed
sparse matrices and Krylov solvers; Bezier extraction, basis derivatives,
quadrature, VMS/SUPG weak forms, boundary conditions, nonlinear iteration, and
time integration are implemented in this repository.

This backend replaces the legacy solver. The matching single-GPU implementation is in [../cuda](../cuda).

## Capabilities

| Executable | Purpose |
| --- | --- |
| `iga_pack` | Validate sparse cache or legacy text and create indexed `.ntiga` |
| `iga_inspect` | Print database metadata and partition statistics |
| `iga_case_check` | Validate case inputs and resolved boundary roles without PETSc |
| `iga_mesh_check` | Evaluate every element at 4x4x4 quadrature points and reject non-positive Jacobians |
| `iga_config_check` | Strictly validate `simulation_config.json` without PETSc |
| `iga_flow_validate` | Integrate boundary flow, check divergence-theorem closure, mass balance, and cycle repeatability |
| `iga_womersley_reference` | Generate a time-aligned analytical straight-tube Womersley velocity series |
| `iga_solve` | Solve configured one-to-many-field transport/operator systems |
| `iga_assembly_smoke` | Exercise distributed owned-row sparse assembly |
| `iga_navier_stokes` | Solve the four-field stabilized steady or transient velocity-pressure system |
| `iga_transport` | Transitional old-input CLI, lowered through the generic assembler |

`iga_solve` writes fields in configured order plus a `.fields` name file.
`iga_navier_stokes` reads flow physics and boundary values from that same config.
`iga_solve` can also read a configured time-resolved velocity manifest and
linearly interpolate flow snapshots onto the transport time grid.

## Why this version is faster and smaller

The legacy MPI programs made every rank scan large ASCII extraction files,
stored dense extraction data, preallocated 1,000 matrix entries per row, and
relied on PETSc off-process insertion. This version changes the data and
assembly path:

- Pack the sparse preprocessing cache, mesh, and partition data into a binary
  database with direct element offsets and per-rank touching-element indices.
- Store only nonzero Bezier extraction rows and seek directly to records needed
  by each rank.
- Assemble only owned matrix rows and derive exact diagonal/off-diagonal
  preallocation from symbolic node adjacency.
- Turn unexpected sparse allocation into an error with
  `MAT_NEW_NONZERO_ALLOCATION_ERR` instead of silently growing matrices.
- Keep element matrices in temporary contiguous buffers and release them after
  insertion.
- Assemble time-independent transport operators once, reuse KSP and
  preconditioners, and reuse the current state as the next initial guess.
- Validate input lengths, partition consistency, Jacobian signs, allocation,
  and `KSPConvergedReason` collectively before accepting results.
- Run the full Navier-Stokes nonlinear loop instead of the legacy hard-coded
  single update.

See [ARCHITECTURE.md](ARCHITECTURE.md) for implementation details.

## Measured improvement over the legacy solver

All measurements below were collected on PSC Bridges-2 with optimized PETSc and
OpenMPI 4.0.5. The cylinder contains 4,221 nodes and 3,600 elements.

| Cylinder measurement | Legacy release | TubularFlowIGA CPU |
| --- | ---: | ---: |
| Comparable first Newton update | 199 s | about 12 s |
| First-update speedup | 1.0x | about 16.6x |
| Peak memory | 6.72 GiB | 1.51 GiB |
| Memory reduction | baseline | about 77.5% |
| Complete converged solve | not implemented | 94 s, 7 completed updates |

The legacy timing is one hard-coded update and is not a converged nonlinear
reference. The CPU solver completes all seven updates in less than half the
legacy single-update time while using about 4.45x less peak memory.

Two cylinder transport steps reproduce the legacy fields to relative L2 errors
of `4.18e-7` for `N0` and `3.13e-7` for `Nplus`. The converged default
block-Jacobi/ILU Navier-Stokes run at the time of this benchmark used 94 s and
1.51 GiB. A different measured Schur field-split alternative used 180 s and
1.67 GiB. These historical measurements do not characterize the current
size-aware GAMG Schur default described below; rebenchmark it before making a
new Bridges-2 performance claim.

### Larger validated CPU case

`NMO_66748_subtree` contains 57,456 nodes and 50,940 elements. Its 16-way
database is 507 MiB and all quadrature samples are positive
(`min(detJ)=8.26e-9`).

| Stage, 16 MPI ranks | Result |
| --- | ---: |
| Navier-Stokes convergence | 2 updates, 1,929 Krylov iterations |
| Navier-Stokes solver time | 401 s |
| Navier-Stokes peak RSS | 10.0 GiB |
| Two-step transport assembly | 32.5 s |
| Two-step transport solves | 5.74 s, 372 iterations |
| Transport total wall time | 45 s |
| Transport peak RSS | 6.34 GiB |

`NMO_54499_new` was not accepted for this CPU validation: 6 of 31,500 elements
contained 43 non-positive quadrature samples with minimum determinant
`-5.67e-6`. The legacy solver silently assembled signed volumes; this solver
fails early. Detailed numerical checks are recorded in
[VALIDATION.md](VALIDATION.md).

## Repository layout

- `include/`: database, case input, element kernels, and owned-row assembly.
- `src/`: packer, inspector, mesh checker, smoke test, and production solvers.
- `slurm/`: Bridges-2 PETSc build workflow.
- `ARCHITECTURE.md`: data layout and numerical design.
- `VALIDATION.md`: measured regression and large-case evidence.

Generated executables, `.ntiga` databases, simulation results, and Slurm output
are intentionally ignored.

## Requirements

- A C++17 compiler.
- MPI; Bridges-2 validation used OpenMPI 4.0.5.
- PETSc with C++ and MPI support; the validated build uses optimized
  `arch-linux-c-opt`.
- METIS/`mpmetis` to create a partition matching the MPI rank count.
- Spline outputs: `bzmeshinfo.txt` and `spline_cache.igacache`.
- Optional legacy outputs: `cmat.txt` and `bzpt.txt`.

## Build

From the repository root:

~~~bash
make cpu
make cpu-petsc PETSC_DIR=/path/to/petsc PETSC_ARCH=your-petsc-arch
~~~

The first command builds `iga_pack`, `iga_inspect`, `iga_case_check`,
`iga_config_check`, `iga_flow_validate`, and `iga_womersley_reference` with a
standard C++17 compiler.
PETSc targets require an MPI-aware PETSc configuration.
`PETSC_ARCH` may be omitted for an installed PETSc layout.

## Prepare a case

Keep large cases outside this repository:

~~~bash
export IGA_CPU_ROOT="$(pwd)/solvers/cpu"
export CASE_DIR=/path/to/case
export DATABASE=/path/to/work/case-8.ntiga
export RANKS=8

mpmetis "$CASE_DIR/bzmeshinfo.txt" "$RANKS"
"$IGA_CPU_ROOT/iga_pack" "$CASE_DIR" "$RANKS" "$DATABASE"
"$IGA_CPU_ROOT/iga_inspect" "$DATABASE"
~~~

The partition suffix and `mpiexec -np` count must match the rank count passed
to `iga_pack`. The packer prefers the sparse cache. Add `--legacy-text` after
`$DATABASE` to force the legacy reader for comparison. New databases include
per-element boundary-face labels; version 3 databases remain readable but do
not provide surface metadata.

## Run

Use an allocated compute resource on shared clusters:

~~~bash
mpiexec -np "$RANKS" "$IGA_CPU_ROOT/iga_mesh_check" "$DATABASE"
mpiexec -np "$RANKS" "$IGA_CPU_ROOT/iga_navier_stokes" \
  "$DATABASE" "$CASE_DIR" --max-newton 12 --output velocity.txt
mpiexec -np "$RANKS" "$IGA_CPU_ROOT/iga_solve" \
  "$DATABASE" "$CASE_DIR" --system neuron_transport --output neuron.txt
~~~

Use the Navier--Stokes command with a prepared `vascular_flow/*` case and the
transport command with a prepared `neuron_transport/*` case. The
application-specific source cases and preparation command are documented in
the [examples catalog](../../examples/README.md). Most small configurations
contain one matching equation system; `vascular_flow/multispecies_pulse`
intentionally runs flow first and then six-species transport from its generated
velocity snapshot manifest.

The canonical [`simulation_config.json`](../../docs/PDE_CONFIGURATION.md) names
fields, systems, operators, viscosity, time integration, and per-field boundary
conditions. A velocity Dirichlet condition can use a three-component value or
`initial_velocityfield.txt` profile plus scale. The old two-file input remains
accepted only as a transition adapter. Navier-Stokes uses nonlinear relative
tolerance `1e-5`, absolute residual RMS tolerance `1e-10` per equation, and
boundary mass-imbalance tolerance `1e-3`; configured transport uses relative
tolerance `1e-8`. Override the flow criteria with `--nonlinear-rtol R`,
`--nonlinear-atol A`, and `--mass-rtol R`. The solver accepts a Newton state
only when the nonlinear and mass criteria both pass and prints L2/RMS,
continuity, and boundary-flow diagnostics for every iteration.

The default Navier--Stokes preconditioner is size-aware. Cases below 1,000
control points retain block Jacobi, which is efficient and deterministic for
small validation problems. Larger mixed velocity-pressure systems use a full
Schur field split with the stabilized pressure block (`A11`) and PETSc GAMG for
both subproblems. Standard PETSc options still take precedence; for example,
`-pc_type bjacobi` explicitly selects the legacy preconditioner.

Backward-Euler flow supports time-indexed text output and a PETSc checkpoint
state with validated JSON metadata:

~~~bash
mpiexec -np "$RANKS" "$IGA_CPU_ROOT/iga_navier_stokes" \
  "$DATABASE" "$CASE_DIR" --output velocity.txt --output-every 10 \
  --checkpoint flow-checkpoint --checkpoint-every 10

mpiexec -np "$RANKS" "$IGA_CPU_ROOT/iga_navier_stokes" \
  "$DATABASE" "$CASE_DIR" --output resumed.txt \
  --restart flow-checkpoint
~~~

Alongside text output, both production solvers write a final ParaView XML
`.vtu`. With `--output-every N`, they also write step-indexed `.vtu` files and
a `.pvd` time collection. Flow arrays are `velocity` and `pressure`; transport
arrays include every configured species and requested derived physiology
field. Transient collections begin with the initialized `step000000` state, so
animations preserve the configured initial values before the first solve.

The metadata validates node count, completed step, physical time, `dt`, density,
and viscosity before restart. A final checkpoint is written even when the final
step is not an interval boundary. Legacy positional `MAX_NEWTON` and `OUTPUT`
remain accepted; the default maximum is 12 nonlinear iterations.
For a reproducible interruption test, add `--stop-after-step N`; the checkpoint
retains the original configured final step, so the same case can resume without
editing its time configuration.

Configured transport uses equivalent flags and validates ordered field names,
system name, and velocity-source name in addition to time metadata. A
snapshot-series restart therefore resumes interpolation at the completed
physical time:

~~~bash
mpiexec -np "$RANKS" "$IGA_CPU_ROOT/iga_solve" \
  "$DATABASE" "$CASE_DIR" --system neuron_transport \
  --output neuron.txt --output-every 10 --stop-after-step 10 \
  --checkpoint neuron-checkpoint --checkpoint-every 10

mpiexec -np "$RANKS" "$IGA_CPU_ROOT/iga_solve" \
	"$DATABASE" "$CASE_DIR" --system neuron_transport \
	--output neuron-resumed.txt --restart neuron-checkpoint
~~~

Add `--memory-report transport-memory.jsonl` to collect one JSON object per
stage after database loading, required-element expansion, matrix
preallocation, first operator assembly, first KSP setup, and the completed
solve loop. The report
contains current and peak RSS, PETSc process memory, and PETSc allocator usage
for every rank plus maximum and aggregate values. Allocation tracking is
enabled before `PetscInitialize` only when this option is present, so normal
runs do not pay the PETSc malloc-debug overhead.

Configured transport derives separate left/previous field-coupling patterns
from the compiled terms and Robin/Dirichlet requirements. PETSc scalar AIJ
preallocation and reusable element `nen x nen` blocks include only those active
equation/trial pairs; `MAT_NEW_NONZERO_ALLOCATION_ERR` remains enabled.

Validate one velocity snapshot, or a complete manifest with an optional cardiac
period for cycle-to-cycle comparison:

~~~bash
"$IGA_CPU_ROOT/iga_flow_validate" "$DATABASE" velocity.txt
"$IGA_CPU_ROOT/iga_flow_validate" "$DATABASE" \
  --compare velocity-reference.txt velocity-current.txt
"$IGA_CPU_ROOT/iga_flow_validate" "$DATABASE" \
  --manifest "$CASE_DIR" velocity_series.csv 0.8
"$IGA_CPU_ROOT/iga_flow_validate" "$DATABASE" \
  --compare-manifests cpu-output velocity_series.csv \
  "$CASE_DIR" velocity_series.csv
"$IGA_CPU_ROOT/iga_flow_validate" "$DATABASE" \
  --womersley examples/validation/womersley/womersley_reference.json \
  "$CASE_DIR" velocity_series.csv
~~~

The reported relative mass imbalance is
`2*abs(sum(Q))/sum(abs(Q))`, using outward flow on every packed boundary label.
The divergence-theorem error is also normalized by the larger of total
absolute boundary flow and absolute volume divergence, making its gate portable
across geometric and unit scales.
Manifest mode also reports the maximum velocity relative L2 between snapshots
separated by the supplied period. The utility accepts either CPU or CUDA
three-column velocity output. `--compare-manifests` requires every reference
time to appear in the numerical manifest (which may contain additional startup
snapshots) and reports the per-snapshot and maximum velocity relative L2,
suitable for CPU/CUDA coefficient parity or an independently projected
coefficient reference. Raw analytical point samples are not IGA control
coefficients. `--womersley` instead evaluates the numerical spline and the
analytical field at volume quadrature points and integrates their physical
relative L2 norm, which is the Womersley acceptance metric.
The generator configuration and pressure-gradient convention are documented in
[`examples/validation/womersley`](../../examples/validation/womersley/README.md).

For PSC module, interactive-node, and Slurm examples, see the
[Bridges-2 guide](../../docs/BRIDGES2.md).

## Reproducibility notes

- Record compiler, MPI, PETSc commit/configuration, rank count, node type, and
  PETSc options with every benchmark.
- Compare numerical fields as well as wall time; changing partition or
  preconditioner can alter reduction order.
- Run `iga_mesh_check` before timing assembly or solves.
- Benchmark from allocated compute nodes and report both stage timers and peak RSS.
- The large validation databases and outputs are not versioned; use matching
  case snapshots when reproducing the published numbers.

## Current limitations

- The CPU solver depends on PETSc for distributed sparse linear algebra.
- Scaling has been validated to 16 MPI ranks, not exhaustively across nodes.
- The binary database currently stores partition-specific touching-element
  indices and should be repacked when the rank count changes.
- CPU Navier-Stokes supports steady and backward-Euler rigid-wall flow,
  per-step velocity Dirichlet and pressure-traction waveforms,
  checkpoint/restart, and time-indexed output.
- Configured scalar flux/Robin surface assembly is available in `iga_solve`
  with a version 4 `.ntiga` database.
- Natural pressure traction and R/RC/RCR outlets are available on CPU. The
  outlet models apply `-p n` without replacing continuity rows; compliant walls
  are not implemented.
- The stabilized formulations follow the legacy model; this repository does not
  claim a new physical model or discretization order.
