# Validation Notes

All simulations below ran on Bridges-2 CPU nodes with the optimized PETSc build and OpenMPI 4.0.5. Generated databases and results are kept in `.nsvms_diagnostics/` and are intentionally ignored by Git.

## Cylinder regression

- Exact owned-row preallocation produced 17,022,960 used/allocated nonzeros with zero allocation growth.
- Two transport steps matched the legacy fields to relative L2 errors of `4.18e-7` (`N0`) and `3.13e-7` (`Nplus`). Compacting the block-diagonal previous-state matrix produced identical output.
- Steady Navier–Stokes reached relative residual `5.75e-6` in seven completed updates. The default block-Jacobi/ILU run took 94 seconds and 1.51 GiB; the Schur field split took 180 seconds and 1.67 GiB. Their converged velocity fields agree to relative L2 `9.13e-10`.
- The first default v2 update took about 12 seconds. The release legacy solver used 199 seconds and 6.72 GiB for its single hard-coded update; that legacy output is not a converged nonlinear reference.
- Feeding the converged velocity directly to two transport steps completed successfully. The result differed from transport using the prescribed initial velocity by relative L2 `2.40e-3`, as expected.

## Larger-case selection

`NMO_54499_new` was rejected: 6 of 31,500 elements have non-positive quadrature-point Jacobians (43 bad samples; minimum `-5.67e-6`; first element 8184). The legacy solvers silently assemble these signed volumes.

`NMO_66748_subtree` was selected instead. It has 57,456 nodes and 50,940 elements—about 14 times the cylinder node count. All quadrature samples passed, with minimum determinant `8.26e-9`. Its generated 16-way database is 507 MiB and has 1.41x touching-element duplication.

The 16-rank Navier–Stokes run converged after two updates: residual `7.51` to `3.28e-5`, 1,929 total Krylov iterations, 401 solver seconds, and 10.0 GiB peak Slurm RSS. The generated 57,456-row velocity/pressure fields contain no nonfinite values. All 11,436 wall velocities, 161 inlet velocities, and 805 pressure-outlet values satisfy their constraints exactly.

Using that velocity directly, two transport steps took 32.5 seconds to assemble and 5.74 seconds to solve (372 iterations); total job wall time was 45 seconds and peak RSS was 6.34 GiB. Reusing the current state as the next KSP initial guess reduced solve time from 9.64 seconds, while changing the result only `4.54e-7` in relative L2. Both concentration fields are finite and positive: `N0` range `[9.62e-6, 2]`, `Nplus` range `[9.95e-7, 1.033]`. All inlet values match the prescribed `N0=2`, `Nplus=1` exactly.

## Transient restart regression

A local two-rank regression on 2026-08-24 used a 2,211-node, 1,800-element
straight tube, backward Euler with `dt=0.1`, two physical steps, density `1`,
and dynamic viscosity `1`. Step 1 converged at residual L2 `1.11153e-8` and
step 2 at `6.87843e-8`. The uninterrupted final norms were velocity L2
`0.248348` and pressure L2 `44.9562`.

The first run wrote a checkpoint after step 1. A second run loaded that state,
resumed at physical time `0.1`, and completed step 2. `cmp` reported identical
uninterrupted and restarted velocity and pressure files. Both time-indexed
velocity outputs contained exactly 2,211 rows. This regression validates state
transfer and restart time; it is not a Womersley, temporal-refinement, or
compute-node performance result.

## Time-resolved transport coupling smoke

The same local two-rank case used the two time-indexed flow outputs as a named
`snapshot_series` velocity source for configured tracer transport. Both the
exact snapshot-time run and a run whose first transport step lay halfway
between two flow snapshots completed two steps with 44 total Krylov iterations
and final L2 `14.1464`. The interpolated and exact-time fields differed by L2
`5.76381e-5`, confirming that the interpolated field reached assembly rather
than being silently replaced by a steady velocity. Operator assembly took
3.64 seconds locally. This is an execution smoke, not a coupled-physics
accuracy validation.

## Resistance outlet smoke

The local 2,211-node straight tube replaced its zero-pressure outlet with
`resistance=0.001` and `reference_pressure=0`. Surface quadrature measured
outward flow `0.00370264`; the coupled pressure converged to `3.70264e-6`,
matching `p=R*Q`. The relaxed fixed-point loop reached its `1.01e-6` pressure
tolerance in three outlet iterations. Pressure was applied as the natural
traction `-p n`, with zero pressure-constrained nodes. The final relative mass
imbalance was `4.88703e-7`, and the relative divergence-theorem error was
`2.48410e-7`. This checks execution, sign convention, coupling, and continuity
retention; physiological parameter calibration remains a separate task.

## Mass-balance diagnostic

`iga_flow_validate` independently integrated all 760 packed boundary faces of
the local straight-tube snapshot. For the final two-step field it measured wall,
inlet, and outlet outward flows `1.54465e-4`, `-3.82790e-3`, and
`4.58733e-3`. The net outward flow was `9.13893e-4`, giving relative mass
imbalance `0.213285` under `2*abs(sum(Q))/sum(abs(Q))`.

The volume integral of `div(u)` was `9.13891e-4`; its difference from the
surface integral was `2.06468e-9`. This verifies the postprocessor geometry,
orientation, and quadrature closure and shows that the imbalance belongs to the
coarse startup solution. It is not a passing branching or multi-cycle
mass-balance result.

## Temporal-refinement diagnostic

The 2,211-node straight tube reached physical time `0.2` with three
backward-Euler grids: `dt=0.1` (2 steps), `0.05` (4 steps), and `0.025`
(8 steps). All nonlinear and linear solves converged. At matching time `0.1`,
successive velocity relative L2 differences were `2.69120e-2` and
`1.14446e-2`, giving observed self-convergence order `1.2336`. At time `0.2`,
the differences were `4.06530e-3` and `8.03039e-3`, so that endpoint was not a
monotone self-convergence sequence. The result supports the expected first-order
behavior at the earlier sample but does not pass the planned Womersley or
multi-cycle gate.

The final relative mass imbalance decreased from `0.213285` to `0.205853` and
`0.186729` as `dt` was halved. Divergence-theorem errors remained near `2e-9`,
confirming that these values are solution diagnostics rather than surface
postprocessing error.

## Womersley and multi-cycle acceptance

A local two-rank straight-tube run on 2026-08-24 used 2,211 nodes, 1,800
elements, backward Euler with `dt=0.1`, two periods of `0.4`, dynamic viscosity
`1`, density `1e-4`, and radius `0.5` in normalized database coordinates. The
sinusoidal pressure-gradient/inlet amplitude was `0.2` about mean `1`; a
zero-valued natural pressure traction closed the outlet. All eight physical
steps and their linear/nonlinear solves converged.

Across the run, `iga_flow_validate --manifest` reported maximum relative mass
imbalance `4.6845851769939491e-7`, maximum relative divergence-theorem error
`2.3807004951043619e-7`, and maximum cycle-to-cycle velocity relative L2
`6.7578940198139713e-6`. On second-cycle times `0.5`, `0.6`, `0.7`, and `0.8`,
`--womersley` integrated physical-volume relative L2 values between
`0.0464696362` and `0.0464865370`; the maximum `0.04648653701365682` passed the
`0.05` gate. Analytical point samples were not treated as spline coefficients:
both numerical and analytical fields were evaluated at 4x4x4 element
quadrature points.

## Branching snapshot transport and restart

A version 4 repack of the 14,565-node, 12,780-element bifurcation contains
3,380 boundary faces. Its wall-supported velocity field measured inlet flow
`1.232562e-2`, outlet flows `4.181540e-3` and `8.142626e-3`, relative mass
imbalance `1.18025e-4`, and divergence-theorem error `1.67875e-10`.

Configured tracer transport then interpolated between two named branch velocity
snapshots for two steps on four MPI ranks. It completed 44 Krylov iterations
with final L2 `24.9893`; all time-indexed fields contained 14,565 rows. A
separate one-step run wrote a transport checkpoint. Resuming step 2 from that
checkpoint produced a final tracer field and `.fields` metadata that were
byte-identical to the uninterrupted run.
