# Native matching ALE–solid interface

The first T5 interface slice is frozen in
[`benchmarks/t5_native_matching_interface_contract.json`](../benchmarks/t5_native_matching_interface_contract.json)
and implemented by
[`NativeTetMatchingFsiInterface.hpp`](../solvers/cpu/include/NativeTetMatchingFsiInterface.hpp).
It is repository-owned and uses no external coupling or FEM framework.

“Matching” has a strict meaning here: fluid ALE vertices and solid reference
vertices use the same stable node IDs, and the selected interface triangles
have identical connectivity. Solid displacement creates the current boundary
positions. Solid velocity becomes ALE boundary velocity and, under no slip,
fluid wall velocity. The native fluid P2 edge-node value is the average of its
two P1 solid endpoint velocities. Interior ALE velocity is not assigned by
this transfer; the native harmonic mesh-motion stage extends boundary motion.

The dynamic input is explicitly `traction_on_structure_pa`. A caller starting
from `sigma * n_fluid` must negate it. Constant traction on a current triangle
is integrated into three equal P1 nodal forces. The transfer independently
reports total force, moment about the global origin, current area, nodal power,
and triangle-quadrature power.

The same header also evaluates the native P2/P1 fluid state on each current
interface triangle, area-averages `sigma*n_fluid`, and performs that sign flip.
A pressure-only regression on the `z=0` face confirms that `p=2 Pa` produces
`-2 Pa` in the structural z direction.

`NativeTetSolidFsiRuntime.hpp` connects the matching publication contract to
the repository-owned hyperelastic tetrahedral solver. It embeds the same
`FsiTrialLifecycle` used by the existing partitioned strong-coupling
coordinator; it does not introduce a second state machine. An exact stamped
`SurfaceTraction.consistent_nodal_force_n` is scattered directly to the native
solid external-force vector. A successful trial publishes interface
displacement and backward-Euler velocity. Rejecting an iteration, or aborting
after prepare, leaves the committed numerical state and committed publication
unchanged. Finalization prevalidates both owners before no-throw swaps.

Run:

```bash
make t5-matching-interface-audit
```

The transfer test exercises finite in-plane deformation and affine velocity.
Total force and moment match direct integration and the discrete power mismatch
is `1.67e-16 W`. Missing/extraneous/duplicate traction, missing labels,
nonfinite state, degenerate current faces, and missing P2 topology fail closed.
The runtime test additionally rejects stale iteration stamps and verifies
trial rejection, prepared-step abort, successful commit, and committed output
identity. This first native solid adapter is deliberately single-partition.
The same test instantiates the existing generic `StrongFluidStructureCoupling`
and `DynamicWeightedAitkenRelaxation` with a constant-load test fluid and the
real native solid adapter. It converges in two iterations after one rejected
iterate and unit relaxation. An injected failure after both sides prepare but
before either finalizes leaves the solid committed identity unchanged. This
validates coordinator wiring and rollback; it is not a fluid-physics case.

`NativeTetAleFsiRuntime.hpp` supplies the corresponding small-case fluid-side
adapter. Its zero-state closed-domain invariant test executes the complete
repository-owned path: matching kinematics, harmonic volume-mesh motion,
native P2/P1 transient ALE flow, current-configuration Cauchy traction,
consistent solid nodal force, native hyperelastic solve, and paired commit.
The ALE kinematics owner now has a fallible prepare followed by a no-throw
finalize, so an injected failure after both domains prepare preserves fluid
time/state and solid state. The adapter also exposed and fixed a multi-tetra
bug: native edge DOF order follows first encounter and is not globally sorted,
so matching lookup must not use binary search.

The invariant smoke begins from rest and must remain at rest. It establishes
operator wiring and transaction semantics only. This adapter uses the dense
reference algebra and one partition; production PETSc/distributed integration
remains separate.

The coordinator no longer has to accept a step from displacement alone. Five
new finite-or-disabled controls independently gate area-weighted traction,
consistent nodal force, fluid nonlinear residual, solid free-DOF residual,
and interface power defect in addition to the existing displacement gate.
Their default is infinity so existing validated immersed cases retain their
prior policy. The native invariant smoke explicitly enables finite (loose,
function-first) values and requires two publications before traction/force
fixed-point convergence can pass. Added-mass-sensitive nontrivial cases are
still required before these tolerances can be called physically calibrated.

A second end-to-end test is nonzero. A small internal P2 velocity pulse decays
in the closed native ALE domain, producing computed wall Cauchy traction that
loads a supported native hyperelastic solid. It converges in two coupling
iterations with committed solid displacement norm `3.01552e-7 m`. The accepted
displacement/traction/fluid/solid/power diagnostics are respectively
`8.33843e-9 m`, `0.016315 Pa`, `1.97091e-9`, `8.71882e-11 N`, and
`3.30872e-24 W`. This proves a nonzero two-way operator path, but its stiff
solid and deliberately loose function-first tolerances do not make it a
compliant-channel, elastic-tube, added-mass, or physiological validation.

The first open-boundary functional case is an idealized `20 x 4 x 4 mm`
channel. It prescribes a cross-sectional parabolic P2 inlet velocity, treats
the outlet as a natural velocity boundary, fixes every non-interface mesh
boundary, and couples the matching top wall to supported hyperelastic tetrahedra.
The fluid adapter requires an explicit value for every essential non-interface
P2 node, rejects absent natural labels and interface/boundary conflicts, and
does not silently turn an unspecified boundary into a wall. The first `0.01 s`
step converges in two coupling iterations; four further steps advance to
`0.05 s`. Each accepted step passes displacement, traction, nodal-force,
fluid, solid, and interface-power gates. The first-step maximum committed wall displacement
is `2.89519e-9 m`; accepted displacement, traction, nodal-force, fluid, solid,
and power diagnostics are `2.89216e-10 m`, `3.26841e-4 Pa`, `3.29236e-9 N`,
`2.0769e-9`, `1.55224e-11 N`, and `8.07794e-28 W`. This closes an autonomous
end-to-end functionality gap, but the coarse mesh, supported
surrogate wall, and loose thresholds make it neither a physical compliant-wall
benchmark nor added-mass, spatial, or temporal validation.

Distributed ownership, a production PETSc fluid adapter, traction varying
within a triangle, nonmatching projection/mortar, the NURBS-shell-to-ALE
parametric map, elastic-tube/bifurcation cases, and physically calibrated
coupled convergence remain separate work.

`NativeTetFsiCheckpoint.hpp` adds a paired accepted-step checkpoint for this
single-partition vertical slice. One bounded deterministic payload stores the
fluid algebraic state, committed ALE and interface displacements, solid
displacement and velocity, and both model/state identities. Capture while a
macro step is active is rejected, and a SHA-256 checksum detects payload
corruption before state construction. The compliant-channel test restores both
runtimes and advances steps two through five along restored and uninterrupted
paths; their fluid and solid committed identities are exactly equal after
every accepted step. This does
not provide distributed shard ownership or repartition restart.
