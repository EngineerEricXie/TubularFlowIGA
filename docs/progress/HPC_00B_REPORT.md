# HPC-00B: phase timing for the selected baselines

Status: **complete for HPC-00B**. The repeated scaling matrix (HPC-00C),
remaining numerical-budget acceptance (HPC-00D), and HPC-01–09 remain open.
Session: 2026-09-08 UTC, using the workstation/toolchain in
[HPC-00A](HPC_00A_REPORT.md). This is a modified worktree, not clean HEAD.

## Accounting contract and coverage

The selected body-fitted flow, transport, depth-2 immersed aneurysm and
compliant-channel FSI now produce per-rank application phase records with
`IGA_PROFILE=1`. The supported CUDA flow/transport paths produce the same
records plus their requested project-buffer allocation high-water mark.
The [benchmark guide](../HPC_BENCHMARKS.md#rank-local-application-phases)
defines each interval and its limitations.

`PhaseProfile.hpp` records monotonic main-thread wall time. An enclosing phase
receives inclusive time and exclusive time after subtracting nested scopes.
Repeated/reentrant phases work the same way; inclusive geometry totals can
therefore exceed exclusive geometry time. Exclusive phase sums plus unscoped
time partition application elapsed time. MPI initialization/finalization and
process launch lie outside application elapsed but inside wrapper process wall.

Coverage includes configuration/field loading, geometry/catalog/layout work,
element/surface/ghost/traction/membrane assembly, linear solver setup and solve,
strong FSI macro-step coupling, explicitly wrapped distributed assembly/ghost
exchange/output gathers, and field/checkpoint output. In-memory serial
immersed/FSI fixtures have no field-output or inter-rank communication phase;
zero calls there are N/A, not evidence of a scalable distributed path.
Console diagnostics and remaining control work are reported as unscoped time.

`PetscPhaseProfile.hpp` calls `KSPSetUp` and `KSPSetUpOnBlocks` before
`KSPSolve`. The latter is the
[PETSc-supported profiling boundary](https://petsc.org/release/manualpages/KSP/KSPSetUpOnBlocks/)
for deferred block preconditioner initialization. Return codes propagate to
the existing solver error path. The solver strategy and tolerances are
unchanged. Legacy console solve intervals still include explicit setup;
the new exclusive phase ledger separates it. Library-internal communication
stays inside the library call's phase; the communication bucket measures
wrapped application collectives, including local work/waits in those calls.
It is not an isolated measurement of network latency.

CUDA assembly scopes finish at existing device synchronization. Profiled
GMRES adds synchronization before entry and after block-inverse construction
to avoid timing launch latency as completed GPU work. Setup is nested inside
linear solve; use exclusive time. Those additional synchronizations are off
when profiling is disabled. Transfer work stays in the enclosing phase.
Membrane linear time includes its small dense factorization.

## Executed evidence

Every row below has successful process records, per-rank logs and a validated
`hpc_profile_summary.py` result under `outputs/hpc00/initial/`. Times are
maximum **exclusive seconds** over ranks. These runs overlapped other work and
are instrumentation/correctness evidence, not speedup measurements.

| Evidence directory | Ranks | Geometry | Assembly | Setup | Linear solve |
|---|---:|---:|---:|---:|---:|
| `profile-flow1` | 1 | 0.111725 | 48.193339 | 3.658016 | 0.019454 |
| `profile-flow2` | 2 | 0.084911 | 44.675699 | 2.681911 | 0.015334 |
| `profile-transport-block2` | 2 | 0.054945 | 2.434215 | 0.009386 | 0.046720 |
| `profile-immersed` | 1 | 2.042037 | 83.773833 | 0.238187 | 0.003501 |
| `profile-fsi` | 1 | 11.223351 | 809.932922 | 0.963057 | 0.020465 |
| `profile-cuda-flow` | 1 | 0.130899 | 6.047296 | 0.002253 | 38.563803 |
| `profile-cuda-transport` | 1 | 0.100038 | 0.089259 | 0.000357 | 1.073160 |

FSI coupling inclusive time was `820.841022 s`; exclusive coupling control
was `0.913322 s`, with nested geometry/assembly/solve removed. Application
elapsed was `823.181449 s` and unscoped time `0.128332 s`. The 28 assembly
calls include fluid, traction and membrane work, and the 12 linear calls
include eight fluid solves and four membrane solves.

The two-rank transport PETSc log separately records `PCSetUpOnBlocks` at
approximately `0.009261 s`, consistent with the explicit setup bucket
(`0.009386 s`, including parent setup and wrapper overhead). That event is
nested with other PETSc events and must not be added to the phase total.

Every rank's peak RSS is retained, not only rank zero. For example,
`profile-flow2` recorded `[195379200, 201969664] bytes` and
`profile-transport-block2` recorded `[80642048, 71835648] bytes`.
Summaries retain min/mean/max, per-rank values, and max/mean imbalance.
The sum of individual RSS peaks is explicitly labeled as such, not a measured
simultaneous aggregate peak. CUDA reports project-buffer requested peaks of
`52,207,740 bytes` for flow and `24,744,428 bytes` for transport; library/driver
allocations and rounding are excluded. Whole-device `gpu_used_gib` stays separate.

## Numerical and failure-path checks

- Body-fitted CPU one/two-rank velocity and pressure fields are byte-identical
  to their corresponding pre-instrumentation direct references.
- Two-rank transport fields are byte-identical to the earlier reference.
- Profiled CUDA flow passes CPU comparison: velocity `4.80008e-8`, pressure
  below the unchanged `1e-5` relative L2 limit. CUDA transport also passes
  the unchanged `1e-5` gate. Individual comparisons include file hashes.
- The complete immersed FD/Newton/conservation gate passes with the same
  printed numerical values as its preceding run (including FD defect
  `1.41286e-11` and normalized open-port balance `2.08177e-15`).
- The complete FSI native gate passes: four strong iterations, eight fluid
  Newton/linear iterations, final RMS `2.71104e-8 m`, center displacement
  `3.80669e-5 m`, normalized moving mass/wall leakage `0.0155613`, continuity
  `1.82225e-14`. Force and moment projection assertions remain active.
- Deterministic C++ phase tests and all 17 Python collection/comparison tests
  pass, including nested/reentrant accounting, exception unwinding, disabled
  profiling, missing/failed ranks, tampered logs and invalid time partitions.
- Existing membrane convergence, membrane-runtime, strong-FSI control and
  fluid-traction tests pass. The traction test required MPI local-socket
  access after a sandbox initialization failure; this was not a numerical failure.

Reproduction uses the commands in each `rank-N/run.json`, with the existing
baseline inputs and fresh output directories. Compile with the system PETSc
prefix from HPC-00A and CUDA `CUDA_ARCHS=89`; enable profiling before the
launcher. The full FSI wrapper used a 1,500-second timeout and the immersed
wrapper a 300-second timeout. Numerical tolerances were not loosened.

These measurements identify fluid assembly as the dominant cost in the
selected CPU immersed/FSI cases and GPU linear solve as the dominant CUDA
flow cost. Choosing an optimization still requires the isolated repetitions
and fixed-resource comparisons of HPC-00C; this report makes no speedup claim.
