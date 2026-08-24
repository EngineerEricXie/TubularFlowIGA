# Configurable neuron-transport and cardiovascular-flow roadmap

This roadmap records the intended direction of TubularFlowIGA as of 2026-08-24.
The target is one configuration model for tubular and branching simulations
that can describe neuron transport or cardiovascular flow without embedding
application-specific field names in numerical kernels.

## Current baseline

| Capability | Status | Current boundary |
|---|---|---|
| Configured multi-field transport | Working on CPU and CUDA | Scalar systems; CUDA supports 1–8 fields |
| Legacy neuron transport expressed as generic operators | Working on CPU and CUDA | Uses a prescribed velocity field |
| Configured cardiovascular flow | Working on CPU and CUDA | Stabilized, rigid-wall, steady Navier–Stokes |
| Shared `simulation_config.json` | Working | Transport and flow still use separate CLI entry points |
| Spatial inlet velocity profile | Working | `initial_velocityfield.txt` times a constant scale |
| Pulsatile inlet waveform | Working for CPU and CUDA configured transport | Cardiovascular flow remains steady |
| Transient Navier–Stokes | Not implemented | No velocity time derivative or previous flow state |
| Flux and Robin surface terms | Working on CPU and CUDA | V100 parity relative L2 `6.03e-8` |
| Physiological outlet models | Not implemented | No resistance or Windkessel model |
| Compliant wall / FSI | Not implemented | Current geometry is rigid |

Existing compute-node smoke evidence covers configured one-field and two-field
transport, the legacy adapter, steady Navier–Stokes, and CPU/CUDA comparisons.
Those results establish the baseline but do not validate pulsatile flow.

## Milestone 0: stabilize the current rewrite

1. Repair and exercise all Makefile targets from a clean build.
2. Reconcile the local branch with the benchmark-documentation commits on
   `origin/main`, then commit the config-driven rewrite in reviewable units.
3. Keep `simulation_config.json` as the canonical schema and retain the v1 text
   inputs only as a documented transition adapter.
4. Run parser, generic-element, preprocessing, packing, two-rank CPU smoke, and
   single-GPU smoke gates. Record commands, job IDs, norms, and field errors.
5. Reject warnings introduced by the rewrite and verify cache and legacy-text
   packing paths.

Exit gate: a fresh clone can build and run the public configured transport case
on CPU, and the same case produces an agreed CPU/CUDA relative L2 result.

## Milestone 1: complete configuration and boundary execution

1. Assemble scalar `flux` and `robin` surface integrals on CPU. **Complete.**
2. Add the same surface operators to CUDA and compare against CPU. **Complete.**
3. Ensure every field/boundary combination is executed or rejected; a
   schema-valid condition must not be silently ignored. **Complete for the
   advertised schema v2 boundary kinds; waveform execution remains an explicit
   Milestone 2 rejection.**
4. Decide whether to add a unified CLI dispatcher. The schema supports multiple
   systems, but `iga_solve` currently dispatches only linear transport while
   Navier–Stokes has a separate executable or command.

Exit gate: every boundary type advertised by schema version 2 has a tested
numerical implementation, or the parser explicitly rejects it.

## Milestone 2: time-dependent inlet waveforms

Extend the schema with named temporal functions instead of a special pulse
branch. Candidate waveform sources are:

- constant;
- sinusoid with mean, amplitude, period, and phase;
- periodic tabulated CSV data with declared interpolation;
- Fourier coefficients for measured or fitted cardiac cycles.

The dependency-free schema, strict periodic CSV reader, and evaluators for all
four sources are complete. CPU and CUDA configured transport materialize
waveform-backed Dirichlet data at every step; their constant-factor regression
has relative L2 `3.52e-8`.
Steady Navier–Stokes intentionally rejects waveform references until transient
CPU and CUDA time stepping is implemented.

A velocity inlet combines a spatial profile and temporal amplitude:

```text
u_inlet(x, t) = spatial_profile(x) * waveform(t)
```

Pressure boundaries should use the same temporal-function interface. Parsing,
interpolation, periodic wrapping, units, and out-of-range behavior need
dependency-free unit tests before solver integration.

Exit gate: constant waveforms reproduce the current fixed-boundary result, and
tabulated/Fourier waveforms pass value, periodicity, and restart-time tests.

## Milestone 3: transient cardiovascular flow

1. Add density and a velocity time-derivative term to the Navier–Stokes weak
   form; begin with backward Euler and preserve a path to BDF2.
2. Add previous-state storage, physical-time stepping, nonlinear checks per
   step, restart/checkpoint metadata, and time-indexed output.
3. Re-evaluate stabilization parameters to include the temporal scale.
4. Apply temporal inlet and pressure values at each physical time step.
5. Implement CPU first, then CUDA after the CPU formulation is stable.

Validation includes the constant/steady limit, temporal refinement, mass
balance, and a straight-tube pulsatile benchmark such as Womersley flow.

Exit gate: CPU and CUDA complete multiple cardiac cycles with bounded mass
imbalance, repeatable cycle-to-cycle behavior, and documented relative L2
agreement.

## Milestone 4: physiological outlets and coupled transport

1. Implement pressure waveforms and resistance/RC/RCR (Windkessel) outlets.
2. Couple time-resolved velocity snapshots to configured transport without
   manually selecting one steady velocity file.
3. Define interpolation when flow and transport time steps differ.
4. Preserve the existing offline one-way workflow for reproducibility.

Exit gate: a branching case balances inlet/outlet flow, advances configured
transport with time-dependent velocity, and restarts consistently.

## Scope boundary: pulsatile inflow versus pulse-wave propagation

A time-varying inlet in a rigid tube produces pulsatile flow. Physical arterial
pulse-wave speed and reflections also depend on wall compliance. Claiming
pulse-wave propagation requires a compliant-wall or FSI milestone, additional
material parameters, and separate validation. Until then, documentation should
use “pulsatile inflow” or “pulsatile rigid-wall flow.”

## Bridges-2 validation policy

When Codex or a developer is already inside a Slurm compute allocation, use it
for lightweight serial builds, unit tests, dependency checks, and the smallest
smoke case after loading required modules. Do not request a nested interactive
allocation.

Use `sbatch` for large geometry, long simulations, benchmarks, GPU work
without an existing GPU allocation, or several independent jobs that should run
in parallel. Never run simulations or expensive builds on login nodes. Record
job ID, hardware, module versions, rank/thread count, assembly and solve times,
memory, convergence reason, norms, and CPU/CUDA relative L2 for numerical gates.
