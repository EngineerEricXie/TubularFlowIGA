# Phase 4 report

Status: **complete**. The CPU Navier--Stokes and transport kernels, their
surface terms and diagnostics, and transport conservation accounting now
consume explicit immutable quadrature rules. Existing body-fitted behavior is
preserved by the full-cell and body-fitted providers.

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

## PR 4.3: transport-kernel migration to explicit rules

The legacy and generic CPU transport volume assemblers now take an explicit
`VolumeQuadratureRule` primary path and apply `weight * raw_detJ` once per
supplied point.  Generic Flux and Robin assembly, along with species flux and
advective--diffusive species diagnostics, have explicit surface-rule paths
that use the point label, physical area weight, and supplied unit normal
directly.  Legacy entry points remain available; runtime transport assembly
and flow-runtime species diagnostics use the explicit paths.  SUPG remains the
existing background-metric formula in this slice.

Focused transport and boundary tests cover full-rule compatibility, curved
geometry, non-tensor volume rules, immersed surface points, and permutation
invariance.

Transport mass and source accounting now uses the same runtime-owned,
per-element volume-rule catalog as assembly.  The default catalog owns one
validated full-cell rule per local element; an injectable catalog must cover
the local element ids exactly and is validated once when the runtime is
constructed.  Field mass and physical-volume helpers use
`point.weight * raw_detJ` once.  Focused affine/skewed off-centre checks verify
constant field mass and source volume against analytic raw-determinant values.

## Phase closure

The phase-closing Sol review accepted the reference/physical weight
conventions, outward-normal and traction/flux signs, unchanged VMS/SUPG
formulas, arbitrary-point kernel paths, legacy compatibility wrappers, and the
shared assembly/accounting rule ownership. The review also confirmed that the
remaining hard-coded tensor loops are limited to geometry preflight,
reference-only diagnostics, compatibility wrappers, visualization, and the
unchanged CUDA backend; they are not CPU physics-kernel integration paths.

Final validation on the local PETSc 3.15 installation:

```text
make -C solvers/cpu test
make -C solvers/cpu petsc PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
make -C solvers/cpu petsc-test PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
make -C solvers/coupling multidomain-test PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
```

All builds were warning-clean and all dependency-free, PETSc, restart, and
one/two-rank MPI tests passed. The multidomain two-species closure retained
maximum integrated residuals of `2.1684043449710089e-19` per edge,
`1.2878566756861082e-14` per domain, and `1.2434504652065331e-14` globally.

Phase 4 therefore meets its exit criterion: CPU physics assembly and
conservation diagnostics integrate arbitrary supplied quadrature points
without knowing whether they originated from a tensor full cell or a future
cut/immersed geometry provider.
