# HPC-00B/D: numerical references and phase instrumentation

Status: **partial; both checklist items and the complete HPC goal remain open**.
Later update: HPC-00B has subsequently passed its selected-baseline acceptance;
see [HPC-00B report](HPC_00B_REPORT.md). HPC-00D has also subsequently passed
the selected baseline gates, independent species budgets, and complete
immersed/FSI state comparisons; see [HPC-00D report](HPC_00D_REPORT.md).
The observations and remaining-work list below preserve the preceding collection state.
Session: 2026-09-07 America/New_York / 2026-09-08 UTC. Toolchain and source
HEAD are recorded in [HPC-00A](HPC_00A_REPORT.md); the worktree is modified.
Generated evidence below lives in ignored `outputs/hpc00/initial/`.

## Numerical reference results

All comparisons retain the predeclared catalog limits: CPU rank coefficient
relative L2 `1e-6`, CPU/CUDA `1e-5`, and zero-reference absolute L2 `1e-12`.
No pressure offset was removed. All field files contain 1,005 nodes.

| Comparison | Velocity relative L2 | Pressure relative L2 | Result |
|---|---:|---:|---|
| Default CPU 1 vs 2 ranks | 9.63026e-7 | 1.20546e-6 | Pressure failed |
| Direct-reference CPU 1 vs 2 ranks | 4.29358e-15 | 1.34034e-15 | Passed |
| CPU vs CUDA before wall-trace correction | 0.317297 | 0.101225 | Failed |
| CPU vs CUDA after correction, steady | 4.80008e-8 | 2.74249e-8 | Passed |
| CPU vs CUDA, waveform transient step 1 | 7.62746e-9 | 3.24176e-9 | Passed |
| CPU vs CUDA, waveform transient step 2 | 6.57355e-9 | 4.40615e-9 | Passed |

The small CPU flow reference uses
`PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps'`
on both rank counts and `--nonlinear-rtol 1e-8 --nonlinear-atol 1e-12
--mass-rtol 1e-6`. This is a correctness reference using the installed MUMPS
backend, not a scalable solver selection or a production-default change.
The tightened iterative serial attempt timed out at 240 s; its logs remain in
`flow-strict1`. Default-solver failures are not overwritten by the accepted runs.

The packed spline geometry gate also passed on two ranks:
`iga_mesh_check flow-baseline/straight_tube-2.ntiga` reported 720 elements,
minimum sampled `detJ=9.01032e-5`, zero bad elements and zero bad samples.
This is distinct from the control-mesh quality check in HPC-00A.

The independent `iga_flow_validate DATABASE VELOCITY` tool reports, rather
than enforces, its conservation metrics. The corrected steady CUDA field has
zero wall flux, relative mass imbalance `3.43214e-8` and relative
divergence-theorem error `1.67600e-8`. CPU and CUDA transient accepted steps
have maximum mass imbalance `4.54177e-8` and maximum relative
divergence-theorem error `2.17603e-8`. Solver convergence additionally enforces
the requested mass threshold `1e-6`. Initial step-zero boundary initialization
is compared as an initial condition, not asserted to be a converged flow.

Transport CPU 1/2 ranks passed at `1.81155e-7`; CPU/CUDA passed at
`4.06867e-6`. The immersed depth-2 native FD/conservation gate passed as
recorded in HPC-00A. The complete serial compliant-channel FSI native gate
also passed: four strong iterations, final weighted displacement RMS
`2.71104e-8 m` below `4.80669e-8 m`, center displacement `3.80669e-5 m`,
eight Newton and eight linear iterations. Normalized moving mass and wall
leakage were both `0.0155613`, and discrete continuity was `1.82225e-14`.
The native force/moment projection assertions were executed, including their
`1e-11 N` / `1e-11 N m` component limits. Complete log:
`fsi-output.log`; process wall `858.66 s`, peak RSS `74,432,512 bytes`.

These runs validate correctness and instrumentation. Some overlapped other
simulations or compilation; none establish isolated scaling or speedup.

## CUDA boundary correction and memory scope

CPU flow constrains all spline bases with nonzero trace on a wall face.
CUDA previously constrained only the control-point label mask. That omitted
wall-trace bases at the inlet/wall intersection and produced a different
discrete problem even when both nonlinear solvers converged.

CUDA now reads the same labeled hex connectivity and calls the existing
`WallTraceBasis(database, mesh)` helper. It applies zero velocity to those
bases at initialization and after each configured boundary/waveform refresh.
This preserves the CPU discretization, pressure gauge, and file formats.
The steady and two-step transient comparisons above cover both paths.

Reproduction after preparing the HPC-00A flow case and building CPU/CUDA:

```bash
mkdir /tmp/hpc-cuda-walltrace
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
LD_LIBRARY_PATH=/home/tsungyeh/anaconda3/envs/tubularflow-cuda/targets/x86_64-linux/lib \
python3 scripts/hpc_rank_run.py --output-dir /tmp/hpc-cuda-walltrace \
  --expected-ranks 1 --timeout 240 -- ./solvers/cuda/iga_cuda navier-stokes \
  /tmp/tubularflow-hpc-flow-baseline/straight_tube-2.ntiga \
  /tmp/tubularflow-hpc-flow-baseline --output /tmp/hpc-cuda-walltrace/velocity.txt \
  --max-newton 30 --nonlinear-rtol 1e-8 --nonlinear-atol 1e-12 --mass-rtol 1e-6
```

For the transient regression, copy `controlmesh.vtk`, `initial_velocityfield.txt`
and `geometry_transform.json` into a fresh case directory. Copy the baseline
configuration and set `time={"dt":0.1,"steps":2}` and flow
`time_integration="backward_euler"`. Add this temporal function and set the
inlet velocity condition's `waveform` to `hpc_inlet`:

```json
{"name":"hpc_inlet","kind":"sinusoid","units":"dimensionless",
 "mean":1.0,"amplitude":0.25,"period":0.4,"phase":0.0}
```

Use the same nonlinear tolerances, add `--output-every 1
--visualization-format vtu`, and run the CPU reference with the one-rank packed
database and direct-solver environment above. CUDA can read the two-rank
database into its single GPU. Compare initial, step-1, step-2, and final
velocity and pressure files independently. Exact executed argv and modified
configuration are retained in `transient-cpu/rank-0/run.json`,
`transient-cuda/rank-0/run.json`, and `flow-transient/simulation_config.json`.
`transient-cuda/comparisons.json` contains all eight field comparisons.

`DeviceBuffer` now tracks requested current/peak allocation bytes, with
ownership transferred on moves and removed after successful release. Byte
count overflow is rejected before allocation. The actual GPU allocation test
passed allocation, move construction/assignment, zero-size, overflow, and
final-release checks (peak 108 bytes; final live zero):

```bash
conda run -n tubularflow-cuda make -C solvers/cuda device_allocation_test CUDA_ARCHS=89
LD_LIBRARY_PATH=/home/tsungyeh/anaconda3/envs/tubularflow-cuda/targets/x86_64-linux/lib \
  ./solvers/cuda/device_allocation_test
```

Measured project peaks: steady flow `52,207,740 bytes`; transport
`24,744,428 bytes`; both ended with zero live project buffer bytes. These
exclude CUDA/cuBLAS internal allocation and allocator rounding. The existing
`gpu_used_gib` remains a separate whole-device usage sample.

## Rank-local phase accounting

`PhaseProfile.hpp` and opt-in `IGA_PROFILE=1` instrumentation now cover the
configured CPU transport path. The [benchmark guide](../HPC_BENCHMARKS.md)
defines phase boundaries, inclusive/exclusive accounting, and unscoped work.
`hpc_profile_summary.py` checks per-rank log hashes, rank completeness,
successful status, finite timing, and the exclusive-time partition. It reports
every rank, extrema, mean, imbalance and process RSS without adding nested time.

Validation executed:

```bash
make -C solvers/cpu phase-profile-test iga_solve \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
python3 -m unittest discover -s scripts/tests -p 'test_hpc*.py' -v
```

The deterministic C++ clock test covers nesting, reentrant categories, explicit
stop, exception unwinding, and disabled profiling. All 17 Python tests passed,
including missing/failed rank, tampered log, nonfinite timing, and overlap
rejection. Actual MPI transport runs with profiling enabled on 1 and 2 ranks,
and disabled on 1 rank, passed. Their fields are byte-identical to each
corresponding pre-instrumentation reference. The disabled run emits no profile.
Raw results and summaries are in `profile-transport1`, `profile-transport2`,
and `profile-disabled`.

The two-rank instrumented run measured assembly exclusive times
`[1.99808, 1.95790] s` and explicitly wrapped collective times
`[0.00352810, 0.0437022] s`. This illustrates observable rank imbalance; it
is not a performance comparison against the concurrent one-rank run.

## Remaining acceptance work

- Complete phase coverage for CPU flow, immersed/moving flow, FSI coupling,
  and remaining output/control work. Separate library PC setup from solve
  using verified PETSc event scope; explicit `KSPSetUp` alone can miss deferred
  setup in `KSPSolve`.
- Complete the repeated isolated resource matrix and machine-readable
  repetition aggregation. No speedup or memory-regression conclusion yet.
- Complete an explicit transport mass/species-budget gate and accepted
  reference provenance for every selected fixture. A matching field norm and
  converged KSP alone do not demonstrate species conservation.
- Preserve all failed observations; keep HPC-00B/C/D and HPC-01–09 open until
  their own full implementation and acceptance criteria are satisfied.
