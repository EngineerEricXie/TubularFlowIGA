# Native 3D VCA bifurcation smoke case

This source-only case exercises the CPU 3D explicit-staggered VCA path with a
flow-controlled pump, one reservoir species, and two venous outlets. Generate
the mesh/database on an allocated resource; generated files are intentionally
not committed.

```bash
RANKS=2 ./scripts/prepare_example.sh vascular_flow/vca_bifurcation /tmp/vca-bifurcation
mpiexec -np 2 ./solvers/cpu/iga_navier_stokes \
  /tmp/vca-bifurcation/vca_bifurcation-2.ntiga /tmp/vca-bifurcation \
  --output /tmp/vca-bifurcation/vca-flow.txt
```

Inspect `/tmp/vca-bifurcation/coupling_manifest.json`. Its inlet and venous
histories should contain four records, both outlet labels, oxygen species flux,
and reservoir volume changes close to zero for an incompressible closed loop.
