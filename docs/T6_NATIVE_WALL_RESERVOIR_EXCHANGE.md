# Native vessel–0D wall reservoir exchange

`NativeTetWallReservoirExchange.hpp` couples labelled, zero-relative-flow
walls of the project-owned moving P1 tetra species FEM to fixed-volume,
well-mixed 0D species stores. The standalone CLI and reusable calculation API
both accept multiple distinct wall labels, each with its own volume and
committed amount. PETSc supplies sparse algebra and MPI only; the
element and face weak forms belong to this repository. fTetWild or Gmsh may
create the labelled tetra mesh, not the discretization or solver.

The vessel's outward species flux on the selected label is
`F = ∫wall h(c_blood−C_reservoir) dA` in mol/s. `h` has units m/s; both
concentrations are mol/m³. The reservoir holds amount `M` in mol and fixed
volume `V` in m³. Backward Euler solves both at the same time:

```text
vessel inventory change / dt + advective outflow + reaction sink + F = body source
(M_new − M_old) / dt = F
C_reservoir,new = M_new / V
```

Positive `F` transfers species from vessel to reservoir; negative `F` is the
reverse direction. There is **no hydraulic volume flow** across this wall,
and this fixed-volume store is distinct from the variable-volume 0D
pressure/flow species reservoir. The boundary label must exist and have zero
relative normal velocity; missing labels, invalid storage/transfer parameters,
and flow through the exchange wall fail closed.

The FEM matrix does not depend on the external concentration. The bounded
reference coupling uses native FEM/PETSc solves at external concentrations 0
and 1 to obtain the exact linear wall-flux response, solves the scalar 0D
backward-Euler balance algebraically, then solves the FEM with that accepted
concentration. It checks the vessel, reservoir and combined species balances
before returning a new state; the caller's committed state is never mutated.
This is a three-solve functional coupling, not a performance-optimized
production algorithm. For multiple regions, the same calculation obtains
the full cross-label flux-response matrix with one zero-concentration solve
and one unit perturbation per label, then solves the small 0D backward-Euler
Schur complement. Each region and the combined vessel–regions budget are
checked before returning. `monotone=true` is used in the bounded tests because
ordinary Galerkin may locally overshoot on coarse, low-diffusion meshes.

Local regression:

```bash
make t6-native-wall-reservoir-exchange-test PETSC_DIR=/path/to/petsc
make t6-native-wall-reservoir-exchange-test PETSC_DIR=/path/to/petsc \
  FTETWILD_BIN=/path/to/FloatTetwild_bin
```

The 1/2/4-rank two-tetra checks cover forward and reverse exchange, two
regions with opposite-direction transfer and cross-label response, rigid
ALE translation, body source plus first-order decay, positivity for the
selected monotone case, rejected invalid inputs and paired mass balance. The
optional fTetWild run generates a newly labelled 0.03 m × 0.005 m pipe,
checks two consecutive steps, reverse transfer and sensitivity to a doubled
`h`. Tetra counts vary between meshing runs; they are not frozen assertions.

This is a single-tracer computational component. Its one-label form is
available through the standalone [species CLI](T9_NATIVE_TET_SPECIES_CLI.md)
as `finite_wall_reservoir`; the multi-label form is available as
`finite_wall_reservoirs_by_label`. Both publish a separate checksummed tissue
sidecar at each accepted step and support paired vessel/tissue restart. The older
`wall_exchange_by_label` option still prescribes an external concentration
and does **not** update finite storage. 3D Darcy perfusion,
tissue mechanics, fluid/graph coupling and liver physiology remain outside
this component.
