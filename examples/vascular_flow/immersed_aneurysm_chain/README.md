# Immersed aneurysm chain

This deterministic steady, quasi-static Phase 6 closure case couples a rigid
Poiseuille-compatible 1D source and sink through an immersed low-Reynolds
number surface-of-revolution bulge.  The planar octagonal inlet and outlet
caps are labels 1 and 2; all lateral facets are wall label 0.  It is not a
transient or backward-Euler-equivalence example.

Run the focused depth-2 Jacobian/conservation regression from
`solvers/cpu` with the system PETSc configuration:

```bash
cd solvers/cpu
make phase6-aneurysm-depth2-regression PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
```

This target makes an isolated depth-2 copy of `immersed/` and runs the
aneurysm Jacobian plus conservation regression. It is not a full coupled-chain
closure target. The checked-in depth-3 fixture is the completed Phase 6
closure benchmark; its composed closure target evaluates three quasi-static
3D load samples for each coupling mode, rather than transient time
integration.
