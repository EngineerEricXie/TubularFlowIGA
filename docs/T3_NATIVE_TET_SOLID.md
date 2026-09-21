# Native tetrahedral hyperelastic solid

T3-B now has a repository-owned element vertical slice for total-Lagrangian
three-dimensional solid mechanics. The frozen scope is recorded in
[`benchmarks/t3_native_tet_solid_contract.json`](../benchmarks/t3_native_tet_solid_contract.json).
No DOLFIN, external FEM discretization, or external assembly framework is
used. PETSc is reserved for a later distributed global runtime.

The initial element is a four-node constant-strain tetrahedron with three
displacement unknowns per reference node and a compressible isotropic
neo-Hookean material:

```text
W(F) = mu/2 (tr(F^T F)-3) - mu ln J + lambda/2 (ln J)^2
P(F) = mu (F-F^-T) + lambda ln J F^-T
```

[`NativeTetHyperelasticSolid.hpp`](../solvers/cpu/include/NativeTetHyperelasticSolid.hpp)
computes deformation gradient, positive deformation Jacobian, strain energy,
first Piola stress, internal residual, analytic consistent tangent, and
consistent mass. It rejects invalid material/state data and inverted current
configurations.

[`NativeTetSolidStaticSolver.hpp`](../solvers/cpu/include/NativeTetSolidStaticSolver.hpp)
adds repository-owned global assembly and a bounded dense reference Newton
runtime. It accepts explicit nodal dead loads and displacement-DOF constraints,
uses a pivoted solve plus backtracking, and returns displacements, constraint
reactions, free residual, strain energy/load potential, iteration count, and
minimum deformation Jacobian. A manufactured homogeneous finite-deformation
case converges from zero in four Newton updates with a `2.83e-14 N` free
residual and closes total applied force with support reaction.

[`NativeTetSolidPrestress.hpp`](../solvers/cpu/include/NativeTetSolidPrestress.hpp)
adds an explicit inverse-elastostatics initializer for loaded image geometry.
It requires the material, dead nodal loads, and zero-displacement supports;
none are inferred from an organ name or silently set to zero. Starting from
the loaded tetrahedral coordinates, it repeatedly runs the full native
nonlinear forward equilibrium solve and relaxes the unloaded reference
coordinates until the forward configuration reconstructs the image. It
returns the estimated unloaded mesh, loaded displacement, equilibrating nodal
internal force, and first-Piola prestress for every tetrahedron. The verified
single-tetra finite-deformation case converges in 16 inverse iterations with
`9.27e-11 m` loaded-geometry error and `9.21e-11 m` reference-coordinate error.

Near incompressibility uses a separate formulation rather than relaxing the
compressible P1 guard. `NativeTetMixedHyperelasticSolid.hpp` implements a
total-Lagrangian equal-order P1 displacement/P1 pressure tetrahedron with an
isochoric neo-Hookean response, pressure-volume constraint, finite bulk
compliance, and element-length-scaled pressure-gradient stabilization. Its
analytic consistent tangent has relative centered-finite-difference error
`6.03e-10`. `NativeTetMixedSolidStaticSolver.hpp` provides the corresponding
bounded global Newton reference runtime.

The locking card bends a `4 x 1 x 1` block at `nu=0.4999`. Across six
tetrahedral refinement levels, mean tip displacement progresses from
`-5.62069e-4` to `-2.2098e-3 m`, with a `4.76%` change on the last refinement.
On the third mesh, changing `nu` from `0.49` through `0.499` to `0.4999`
changes displacement by at most `0.166%` instead of collapsing the response;
maximum `|J-1|` at `nu=0.4999` is `1.20e-4`. This is a numerical locking
verification, not tissue calibration or a claim that equal-order P1/P1 is the
only suitable mixed formulation.

Run the local gate with:

```bash
make t3-solid-contract-audit
```

The test covers zero state, rigid translation, finite rigid rotation,
homogeneous finite deformation, global internal-force and moment balance,
analytic-tangent finite differences, tangent symmetry, consistent total mass,
inversion rejection, inverse reference reconstruction, and prestress recovery.

This completes the functional small dense tetrahedral subset, but is not a
production-scalable T3-B backend. The displacement-only P1 tetrahedron remains
limited to Poisson ratio 0.45; near-incompressible work must explicitly select
the verified mixed formulation.
The current input geometry must be an explicitly declared unloaded reference,
or must pass through the bounded inverse initializer with explicit loads,
material, and supports. Unknown patient-specific loads or material parameters
are not inferred. Hexahedra, fibres,
active contraction, volume IGA, global Newton/PETSc assembly, boundary loads,
distributed/PETSc assembly, surface traction integration, spatial convergence,
and the other listed element/material options remain explicit later work.
