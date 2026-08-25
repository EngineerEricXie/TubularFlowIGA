# Simple-geometry validation

All three committed cases were prepared and solved on 2026-08-24. CPU runs
used two MPI ranks with PETSc 3.15.5. CUDA runs used one NVIDIA GeForce RTX
4080 SUPER (SM 89) with CUDA 12.6.3. Results are correctness evidence for the
example inputs, not cross-hardware performance benchmarks.

## Geometry and packing

| Case | Nodes | Elements | Control-mesh minimum scaled J | Packed minimum detJ |
|---|---:|---:|---:|---:|
| `straight_tube` | 1,005 | 720 | 0.763084 | `9.00989e-5` |
| `bent_tube` | 2,010 | 1,620 | 0.750036 | `5.92202e-5` |
| `y_bifurcation` | 7,329 | 6,300 | 0.305198 | `1.52566e-5` |

Every sampled Jacobian was positive. The schema-v2 case checker found the
expected wall/inlet/outlet labels and successfully resolved both configured
equation systems without any legacy input file.

## Historical two-step tracer transport

These transport results were recorded before the public cases were separated
into application-specific configurations. They remain useful numerical
evidence for the same geometries, but tracer transport is no longer part of
the vascular-flow example configuration.

| Case | CPU iterations | CUDA iterations | Final L2 | CPU/CUDA relative L2 |
|---|---:|---:|---:|---:|
| `straight_tube` | 35 | 105 | 13.7781 | `3.75856e-8` |
| `bent_tube` | 36 | 128 | 13.8932 | `5.10308e-8` |
| `y_bifurcation` | 42 | 136 | 13.7895 | `5.97178e-8` |

All solves reported zero singular diagonal blocks on CUDA.

## Steady Navier–Stokes

The inlet profile scale was `0.01`, viscosity and density were `1`, walls were
no-slip, and terminal pressures were natural zero-pressure tractions.

| Case | CPU linear iterations | CUDA linear iterations | Velocity relative L2 | Pressure relative L2 | Relative mass imbalance |
|---|---:|---:|---:|---:|---:|
| `straight_tube` | 404 | 1,378 | `1.99420e-10` | `1.78270e-10` | `2.83042e-7` |
| `bent_tube` | 592 | 3,823 | `8.99887e-11` | `7.78402e-11` | `8.39831e-8` |
| `y_bifurcation` | 1,509 | 3,440 | `6.52563e-10` | `6.88461e-10` | `8.36495e-8` |

The symmetric Y case produced outward terminal flows `1.868446e-3` and
`1.868436e-3`, confirming that both outlet labels execute independently and
balance the inlet.
