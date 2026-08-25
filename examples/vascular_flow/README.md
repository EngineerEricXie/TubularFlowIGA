# Vascular-flow examples

These small cases solve steady, rigid-wall, three-dimensional incompressible
Navier--Stokes flow. They are intended for first runs, boundary-condition
checks, and CPU/CUDA correctness comparisons rather than physiological
calibration.

## Cases

| Case | Geometry | Boundary labels |
|---|---|---|
| `straight_tube` | Short constant-radius vessel | wall 0, inlet 1, outlet 2 |
| `bent_tube` | Planar quarter bend | wall 0, inlet 1, outlet 2 |
| `y_bifurcation` | One inlet splitting into two arms | wall 0, inlet 1, outlets 2 and 3 |

Every configuration contains only the `blood_flow` system with velocity and
pressure fields. Walls are no-slip, the generated inlet profile is multiplied
by scale `0.01`, and outlets use zero natural pressure traction. The viscosity
and density values are numerical example parameters, not patient-specific
blood calibration.

Spline extraction normalizes the domain coordinates before packing. Treat the
committed values as a numerical smoke configuration; a physiological study
must define a consistent length, time, velocity, pressure, density, and
viscosity scaling.

## Prepare and run on CPU

From the repository root:

```bash
VASCULAR_WORK="$(mktemp -d /tmp/tubularflowiga-vascular.XXXXXX)"
RANKS=2 ./scripts/prepare_example.sh \
  vascular_flow/straight_tube "$VASCULAR_WORK"
VASCULAR_DB="$VASCULAR_WORK/straight_tube-2.ntiga"

mpiexec -np 2 ./solvers/cpu/iga_mesh_check "$VASCULAR_DB"
mpiexec -np 2 ./solvers/cpu/iga_navier_stokes \
  "$VASCULAR_DB" "$VASCULAR_WORK" \
  --output "$VASCULAR_WORK/velocity-cpu.txt"

./solvers/cpu/iga_flow_validate \
  "$VASCULAR_DB" "$VASCULAR_WORK/velocity-cpu.txt"
```

Build the MPI/PETSc solver first if `solvers/cpu/iga_navier_stokes` is absent:

```bash
make cpu-petsc PETSC_DIR="$PETSC_DIR" PETSC_ARCH="${PETSC_ARCH:-}"
```

The requested output contains three velocity columns per node. Pressure is in
the neighboring `velocity-cpu.txt.pressure` file. `iga_flow_validate` reports
flow on every boundary, net outward flow, relative mass imbalance, volume
divergence, and divergence-theorem error.

## Run on CUDA

```bash
./solvers/cuda/iga_cuda mesh-check "$VASCULAR_DB"
./solvers/cuda/iga_cuda navier-stokes \
  "$VASCULAR_DB" "$VASCULAR_WORK" \
  --output "$VASCULAR_WORK/velocity-cuda.txt"
```

Build with `make cuda CUDA_ARCHS=YOUR_GPU_ARCH` and verify
`./solvers/cuda/iga_cuda device-info` before running. CUDA additionally writes
`velocity-cuda.txt.vtk` for direct visualization.

## Customize

- Change centerline coordinates and radii in `skeleton_initial.swc`.
- Change segment length, smoothing, or bifurcation refinement in
  `mesh_parameter.txt`.
- Change viscosity, density, time integration, inlet profile scale, temporal
  functions, or outlet models in `simulation_config.json`.

These geometries are rigid. A time-dependent inlet produces pulsatile
rigid-wall flow, not compliant-wall pulse-wave propagation. See the
[boundary-condition guide](../../docs/BOUNDARY_CONDITIONS.md) and
[PDE configuration guide](../../docs/PDE_CONFIGURATION.md) before changing the
physics.

Measured geometry, flow, mass-balance, and CPU/CUDA results are recorded in
[`VALIDATION.md`](VALIDATION.md). Its tracer table is retained as historical
evidence from before application-specific example configurations were split.
