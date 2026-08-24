# TubularFlowIGA CPU Backend

A project-owned C++17 isogeometric-analysis pipeline for stabilized steady
Navier-Stokes flow and configurable transient multi-field transport. PETSc provides distributed
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
| `iga_solve` | Solve configured one-to-many-field transport/operator systems |
| `iga_assembly_smoke` | Exercise distributed owned-row sparse assembly |
| `iga_navier_stokes` | Solve the four-field stabilized steady velocity-pressure system |
| `iga_transport` | Transitional old-input CLI, lowered through the generic assembler |

`iga_solve` writes fields in configured order plus a `.fields` name file.
`iga_navier_stokes` reads flow physics and boundary values from that same config.

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
block-Jacobi/ILU Navier-Stokes run used 94 s and 1.51 GiB. A measured Schur
field-split alternative used 180 s and 1.67 GiB, so block-Jacobi/ILU remains
the default.

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

The first command builds `iga_pack`, `iga_inspect`, and `iga_case_check` with a
standard C++17 compiler. PETSc targets require an MPI-aware PETSc configuration.
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
  "$DATABASE" "$CASE_DIR" 8 velocity.txt
mpiexec -np "$RANKS" "$IGA_CPU_ROOT/iga_solve" \
  "$DATABASE" "$CASE_DIR" tracer_transport tracer.txt velocity.txt
~~~

The canonical [`simulation_config.json`](../../docs/PDE_CONFIGURATION.md) names
fields, systems, operators, viscosity, time integration, and per-field boundary
conditions. A velocity Dirichlet condition can use a three-component value or
`initial_velocityfield.txt` profile plus scale. The old two-file input remains
accepted only as a transition adapter. Navier-Stokes uses nonlinear relative
tolerance `1e-5`; configured transport uses relative tolerance `1e-8`.

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
- Navier-Stokes is steady and uses fixed-in-time velocity/pressure boundaries;
  configured `time` currently applies only to transport.
- Flux/Robin surface assembly, pulsatile inflow, physiological outlet models,
  and compliant walls are not implemented.
- The stabilized formulations follow the legacy model; this repository does not
  claim a new physical model or discretization order.
