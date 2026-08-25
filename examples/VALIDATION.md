# Public-example validation

These checks were rerun on 2026-08-25 after separating the public neuron and
vascular configurations. The local x86_64 workstation used an Intel Core
i9-14900KF, GCC 11.4.0, OpenMPI 4.1.2, PETSc 3.15.5, and two MPI ranks.
ParaView 5.13.1 rendered the committed gallery images in offscreen batch mode.

## Preparation gates

All three cases below completed mesh generation, Bezier extraction, two-way
METIS partitioning, `.ntiga` packing, schema validation, and packed boundary
resolution. Every sampled control-mesh Jacobian was positive.

| Case | Nodes | Elements | Minimum scaled J | Packed boundary labels |
|---|---:|---:|---:|---|
| `vascular_flow/straight_tube` | 1,005 | 720 | `0.763084` | wall 0, inlet 1, outlet 2 |
| `vascular_flow/y_bifurcation` | 7,329 | 6,300 | `0.305198` | wall 0, inlet 1, outlets 2 and 3 |
| `neuron_transport/branched_neurite` | 7,329 | 6,300 | `0.305198` | wall 0, inlet 1, outlets 2 and 3 |

## Steady vascular flow

Both rigid-wall Navier--Stokes cases converged after six Newton updates.

| Case | Linear iterations | Velocity L2 | Pressure L2 | Relative mass imbalance | Relative divergence error |
|---|---:|---:|---:|---:|---:|
| `straight_tube` | 404 | `0.166093` | `11.8216` | `2.83042e-7` | `1.45381e-7` |
| `y_bifurcation` | 1,509 | `0.332098` | `58.8361` | `8.36495e-8` | `4.24139e-8` |

The Y case produced outward flows `1.868446e-3` and `1.868436e-3` at the two
terminal labels. The final net outward flow was `3.20202e-10`.

## Branched neuron transport

The configured `N0`/`Nplus` system completed two physical steps in 88 Krylov
iterations. Assembly took `12.0796 s`, the time-loop solves took `0.4497 s`,
and the final coefficient L2 norm was `29.6961`. The output contained exactly
7,329 node rows and the `.fields` file recorded `N0` followed by `Nplus`.

The coefficient ranges were:

- `N0`: `[-9.21959e-3, 1]`;
- `Nplus`: `[-2.01918e-1, 2]`.

This intentionally short case is an execution and coupling smoke test, not a
positivity-preserving or biologically calibrated transport benchmark. The
reported undershoot must not be interpreted as a physical concentration.

## ParaView rendering

The gallery PNGs were produced with `scripts/render_example.py` and
`pvbatch --force-offscreen-rendering`. The renderer deep-copies the legacy VTK
mesh, validates result row counts and field ordering, attaches CPU coefficients
as point arrays, and colors a center-plane slice.

These images are piecewise-linear control-mesh visualizations. Quantitative
flow values above come from `iga_flow_validate` using the packed IGA basis and
quadrature rather than measurements from the rendered pixels.

CUDA was not rerun for this release check because the current execution
environment did not expose an NVIDIA device. Existing scheduled CPU/CUDA parity
evidence remains in the solver validation documents.
