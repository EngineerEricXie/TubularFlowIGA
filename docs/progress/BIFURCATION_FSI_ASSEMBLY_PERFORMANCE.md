# Bifurcation FSI assembly performance

User authorized efficiency improvements, pausing the original job only if needed.
The heartbeat job 45971412 retains its archived binary on r325. Benchmarks and
regressions use separate Slurm allocations and local scratch. No physical or
convergence settings are relaxed by this optimization.

## Candidate

The current assembly parallelizes volume kernels but executes conservative trace
and Nitsche wall kernels in the main-thread consumption callback. The candidate
computes volume, trace and wall systems together in bounded worker batches.
PETSc insertion, wall-diagnostic accumulation, port assembly and ghost assembly
retain their original main-thread order. No MPI/PETSc calls enter workers.

Each worker can retain two additional element systems (trace and wall), bounded
by the existing batch capacity, default at most eight cells. Geometry, quadrature,
forcing and time integration are unchanged.

## Validation and evidence

- 45981392: candidate regression, 4 allocated CPU cores on r388. All six cases
  (expanded/compact × 1/2/4 threads) passed exact residual/Jacobian action and
  retry/abort checks. Added nonfinite material-wall fault verifies that a wall
  kernel exception is propagated and a subsequent assembly matches exactly.
- 45981464: Y-vessel assembly comparison pending; baseline/candidate at 8 threads,
  candidate at 1 thread, three states each. Same 16×16×5 grid, depth-1 cut
  quadrature with depth-3 empty-cell rescue, density 1060, viscosity .0035,
  dt .05, wall inertial gamma 1. Reports geometry/runtime/begin/assembly costs,
  residual-plus-Jacobian-action SHA256 and PETSc matrix reallocations.
- 45981241 was rejected before any assembly due to the benchmark's incorrect
  fixed-geometry target time. v2 sets an initial state at -.05 and advances to
  the fixed geometry's time 0, with the same dt .05. Geometry took 13.5675 s and
  runtime construction 2.30989 s; 4672 algebraic rows. These are setup timings,
  not a full moving FSI profile.
- 45981204 was cancelled after source staging into `.nsvms_diagnostics` failed
  with a transport-endpoint error. Benchmark staging now uses `benchmarks/`;
  all computation uses Slurm node-local scratch.

Archives and stdout/stderr are under `benchmarks/heartbeat-{performance,regression}-*`.
They are ignored by Git. Candidate has been applied to the worktree after the six passing regressions;
no running or queued heartbeat binary was replaced. Measured speedup and a
moving-wall integration check remain pending before performance closure.

## Follow-up execution

Shared-filesystem header reads stalled compilers with nearly zero CPU time in
45981464 and 45981570; these tests were cancelled without using their wall times
as performance measurements. 45981681 was cancelled before testing to correct
the benchmark's nonnegative time contract. Final benchmark uses geometry evaluated
at .05, zero initial state at 0, and a .05 step. All spatial/physical settings
remain equal between variants. Job **45981715** embeds its test payload in the
Slurm script, requests r253 (where the first build succeeded), and sequences the
assembly comparison and both versions of the existing compliant-channel FSI test.
Its output archive will be `benchmarks/heartbeat-performance-45981715.tar.gz`.

The local candidate header is byte-identical to the header tested by regression
45981392. The production regression test now includes that same wall-failure
retry check. No measured speedup is claimed yet. The previous one-step geometry
setup timings do not establish the full FSI bottleneck distribution.

## First measured assembly result (45981715)

Same Y-vessel and 8 threads: baseline median 28.543797581 s versus candidate
15.982731132 s, a 1.7859148943× assembly speedup (44.0% less time). Three samples
per variant; residual and Jacobian-action hashes match exactly across baseline,
candidate and candidate-serial. PETSc matrix reallocations remain zero.
Candidate-serial median is 58.190114988 s. This is assembly-only, not a claim of
1.79× full-run speedup. The original simple 27-cell complete FSI step took
473.56193825 s; the optimized counterpart is still running.

Pending old-binary retry 45980201 is held while the full FSI comparison finishes,
to avoid starting a lengthy rerun just before a validated optimized binary is
available. Original heartbeat 45971412 is not stopped by this hold.

## Completed comparison (45981715)

Slurm completed with exit 0. See `benchmarks/heartbeat_assembly_comparison.json`.
The 27-cell complete coupled step passed both variants, with identical reported
high-precision coupling residual histories and final diagnostics. Each variant
ran once on the same 8-core allocation: baseline 473.56193825 s, candidate
445.392055263 s, 1.06324738× speedup, 5.9485% less elapsed time. This single-pair
observation is not a statistical guarantee. Assembly consumed 436.9328 s versus
420.0329 s in this volume-heavy depth-3 case. The Y-vessel depth-1 assembly result
remains 1.7859× from three samples. Do not transfer that factor to full heartbeat
wall time. Both tests retain their original numerical criteria.

Original heartbeat 45971412 reached its 4h limit and timed out. The held old
retry will be replaced with the tested parallel assembly version; its physical
case, twenty frames and visualization-only conservation policy remain unchanged.

Optimized full heartbeat job **45982619** submitted with 8 CPU cores and 12h.
It builds from the original 45971412 source snapshot, replacing only the
ImmersedTransientFlowRuntime.hpp archived by passing comparison job 45981715.
Old held retry 45980201 was cancelled after the new submission succeeded.
Outputs: `benchmarks/heartbeat-optimized-45982619/results/`; wrapper:
`benchmarks/heartbeat-optimized-launch.sbatch`. Full heartbeat speedup and
expansion/contraction verification are still pending.

## Completion evidence and next optimization plan (2026-09-14 PSC time)

Job 45982619 completed all 20 solver steps in 13603.010422705 s. Its final
Slurm FAILED status came from native VTK validation of a missing
NumberOfTuples attribute on display_displacement_scale, after the solve and
heartbeat motion checks had succeeded. Export repair job 46006012 completed
with exit 0, preserving numerical values; native VTK validation and the
ParaView archive are now available under heartbeat-optimized-45982619.
The case remains visualization_only, with 17 conservation-failing steps.

Exclusive full-run costs: assembly 65.43%, geometry 16.61%, solver setup
10.80%, other coupling work 6.34%, linear solve 0.17%. These supersede the
small-fixture profile as the guide to heartbeat optimization. No paired
full-cycle old/new speedup has been established.

The next work is specified in
[the PSC time optimization plan](../plans/BIFURCATION_FSI_PSC_TIME_OPTIMIZATION_PLAN.md)
and [the gpt-sol handoff](../plans/BIFURCATION_FSI_GPT_SOL_HANDOFF.md).
Those new tasks are planned, not implemented by this documentation update.
