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
