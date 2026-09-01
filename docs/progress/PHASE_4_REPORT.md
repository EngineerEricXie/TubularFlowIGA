# Phase 4 report

Status: **in progress**.

## PR 4.1: immutable CPU quadrature rules

The CPU backend now has an element-local quadrature API with providers whose
rules remain valid by immutable reference for the provider lifetime.

- Volume points use `(u,v,w)` in `[0,1]^3` and a dimensionless reference-cube
  weight `du dv dw`. Future volume kernels multiply the raw geometric
  determinant `det(dx/dxi)` exactly once.
- Surface points store mapped metres, a unit geometric outward `normal`, a
  nonnegative boundary label, and an already-physical area weight `dA` in
  square metres. Future surface kernels must not multiply a Jacobian again.
- The concrete volume rule preserves the existing `qz/qy/qx` 4x4x4 order.
  The body-fitted surface rule preserves existing face then `qi/qj` ordering
  and the six reference-face conventions.
- `BasisValues::determinant` remains the historical `raw_detJ/8` assembly
  measure. `BasisValues::raw_determinant` exposes `raw_detJ` separately; no
  solver loop was migrated in this slice.

Rule validation rejects invalid coordinates, weights, surface data, mapping or
normal inconsistencies, non-positive volume Jacobians, and an empty active
volume rule. Generic surface validation permits immersed interior-parametric
points because an element map alone cannot establish their outward orientation.
The body-fitted provider additionally verifies reference-face membership and
the mapped outward `normal` convention.

## PR 4.2: flow-kernel migration to explicit rules

Navier--Stokes assembly now has a primary overload that consumes a
`VolumeQuadratureRule`.  It evaluates each supplied point once and uses
`point.weight * basis.raw_determinant` as the volume measure.  The established
full-cell entry point remains a compatibility wrapper around the immutable
4x4x4 provider; stabilization and VMS expressions are otherwise unchanged.

The transient flow runtime now passes explicit full-cell volume rules into its
element assembly and body-fitted surface rules into flow, pressure/area, and
pressure-traction diagnostics.  Their primary surface paths consume the
already-physical point `weight`, the supplied unit `normal`, and `boundary_id`
directly.  A repeated boundary label is integrated once over all matching
points rather than once per face.  The legacy face-based interfaces remain for
existing callers.  Transport/species boundary integration has intentionally
not migrated yet.

Focused checks compare explicit Navier--Stokes assembly against a test-only
copy of the pre-refactor 4x4x4 tensor loop (including its stabilization
calculation) for nonconstant steady and transient states on unit and curved
elements.  Boundary flow, scalar/area, traction, and divergence compare the
explicit rules with preserved legacy face paths.  An off-centre single volume
point analytically isolates the temporal velocity residual
`-rho*N_a*(u-u_previous)/dt*weight*raw_detJ`; an immersed synthetic surface
point independently checks direct use of supplied physical weight and normal.
Surface-point permutation leaves integrals unchanged.
Evidence: `make -C solvers/cpu navier_stokes_test boundary_flow_test
PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real`, both focused
binaries, `make -C solvers/cpu test`, and `vca_3d_runtime_test` all passed
warning-clean on the local PETSc 3.15 build.
