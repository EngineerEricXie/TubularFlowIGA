# T3 native IGA Kirchhoff--Love shell

The first structural backend is repository-owned C++ isogeometric
Kirchhoff--Love shell code. PETSc may provide distributed algebra, but no
external package may provide shell basis evaluation, kinematics, integration,
residuals, or tangents. External shell/solid solvers are cross-validation
references only.

The frozen machine-readable card is
`benchmarks/t3_native_iga_shell_contract.json`. Run `make t3-contract-audit`
before changing the implementation or its tests.

## Initial scope

The unknowns are three translational components at tensor-product B-spline or
NURBS control points. Rotations follow from midsurface derivatives. The first
formulation is total-Lagrangian Kirchhoff--Love with membrane and bending
resultants, geometric nonlinearity, consistent tangent, and isotropic
St. Venant--Kirchhoff material. Transverse shear is neglected, so the initial
route is recommended only when thickness divided by characteristic in-plane
length is at most 0.05.

Polynomial degrees must be at least two and the analysis surface must be at
least C1 between element boundaries inside a patch. Version one accepts one
patch and rejects non-C1 seams. A later interface card may add a verified
penalty, Nitsche, bending-strip, or mortar coupling; silently treating a C0
seam as a Kirchhoff--Love connection is forbidden.

Inputs use metres, pascals, and kilograms per cubic metre. The midsurface,
positive weights, thickness, material, loads, constraints, and reference
configuration are explicit. An organ name cannot choose material or thickness,
and image geometry is not assumed stress-free unless its reference-state
policy says so.

## Verification sequence

Implementation proceeds through basis/derivative identities, rigid-motion
invariance, constant membrane strain, constant change of curvature, centered
finite-difference tangent agreement, then a three-level manufactured spatial
convergence series. A rendered deformation is demonstration evidence only.
The old P1 normal pretension membrane remains a separate reduced model and is
not renamed or counted as this shell.

The card freezes numerical tolerances, failure classes, output fields, and a
4 GiB/300 s local small-case budget before element implementation begins.

## Implemented native foundation

`NativeIgaShellSurface.hpp` now owns open tensor-product B-spline/NURBS basis
evaluation through second derivatives, rational quotient derivatives,
midsurface tangents, second derivatives, normal, and surface Jacobian. It
rejects degree below two, malformed/open-knot violations, internal continuity
below C1, inconsistent control grids, nonpositive weights, parameters outside
the patch, and singular metrics. Tests cover partition of unity, derivative
identities, an exact plane, finite-difference first/second derivatives, and a
rational quarter cylinder.

`NativeIgaKirchhoffLove.hpp` computes reference/current metrics and curvature,
total-Lagrangian Green membrane strain, curvature change, plane-stress
St. Venant--Kirchhoff membrane resultants and bending moments, strain energy,
and surface mass. Its analytic checks cover rigid motion, finite membrane
stretch, and pure bending.

`NativeIgaShellPatchSystem.hpp` integrates every nonzero knot span with
four-point tensor Gauss quadrature and returns patch energy, internal residual,
consistent tangent, and consistent mass. The initial residual/tangent are
obtained by centered differentiation of the repository-owned energy and are a
small-case reference implementation, not the intended production-performance
path. Tests enforce rigid-translation energy invariance, tangent symmetry and
directional finite-difference agreement, and exact total mass. Analytic or
automatic-differentiation assembly, loads/constraints, Newton orchestration,
multi-patch coupling, and convergence cards remain open.

`NativeIgaShellStaticSolver.hpp` supplies the first end-to-end autonomous
small-patch runtime. It integrates constant dead surface traction into
control-point forces, accepts explicit prescribed DOFs, reduces the full
three-component tangent to free DOFs, uses pivoted dense elimination, and runs
Newton with potential-energy line search. The regression solves a constrained
nonlinear extension in two updates, reaches a `6.8e-6 N` free residual, and
balances the applied 20 N with support reactions. This dense runtime is a
verification vehicle; larger patches will use the same native assembly with
PETSc algebra after the analytic/AD tangent replaces numerical differentiation.

The curvature approximation regression uses a smooth manufactured
`sin(pi u) sin(pi v)` field, cubic C2 splines, Greville collocation, and
2/4/8 element levels. Its relative H2-like curvature errors are
`7.8557e-2`, `1.7476e-2`, and `4.1455e-3`, giving observed orders 2.168 and
2.076. This proves the bending-field approximation path, but it is not yet the
full loaded-shell PDE convergence card. Rigid rotations are checked at 0.43
and 2.40 radians with zero membrane and bending energy.

The loaded simply-supported unit plate uses the analytic solution
`w=sin(pi x) sin(pi y)` and load `q=4 pi^4 D w`. Native cubic IGA assembly on
2/4/8 element levels gives displacement L2 errors `8.433e-3`, `7.010e-4`, and
`3.657e-5` (orders 3.59 and 4.26), and H2-like curvature errors `6.598e-2`,
`1.734e-2`, and `4.142e-3` (orders 1.93 and 2.07). A deformed pure-bending
state also retains its membrane and bending energies after a superposed
2.10-radian rotation and translation, checking large-rotation objectivity.

`NativeIgaShellPatchInterface.hpp` implements the initial compatible-interface
diagnostic for a left `u=max` edge joined to a right `u=min` edge with the same
seam parameter direction. It samples both reference and current configurations
and checks position, along/across-seam tangents, normals (the thickness
direction), membrane resultants, and bending moments. A two-patch stretched and
bent surface passes with action--reaction jumps below `1e-8`; injected gaps and
C0 kinks fail closed. This proves interface compatibility diagnostics, not a
general nonmatching or weak multi-patch coupling. The runtime remains
single-patch for nonlinear orchestration, but the compatible C1 assembly
foundation is implemented: explicit C0 displacement and span/degree-scaled C1
derivative MPCs, an elimination transformation, and `T^T r` / `T^T K T`
reduction. Arbitrary independent vectors satisfy every interface constraint;
the virtual-work identity and reduced-tangent symmetry are regression-tested.
Nonmatching or weak multi-patch coupling remains a later card.
