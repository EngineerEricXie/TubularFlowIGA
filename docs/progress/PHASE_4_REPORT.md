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
