# Public-example validation

These checks were rerun on 2026-08-25 after separating the public neuron and
vascular configurations. The local x86_64 workstation used an Intel Core
i9-14900KF, GCC 11.4.0, OpenMPI 4.1.2, and PETSc 3.15.5. Earlier gallery
checks used two MPI ranks; the pulse multispecies case used eight.
The updated branched-neurite animation case also used eight ranks.
ParaView 5.13.1 rendered the committed gallery images in offscreen batch mode.

## Preparation gates

All five small showcase cases below completed mesh generation, Bezier extraction, METIS
partitioning, `.ntiga` packing, schema validation, and packed boundary
resolution. Every sampled control-mesh Jacobian was positive.

| Case | Nodes | Elements | Minimum scaled J | Packed boundary labels |
|---|---:|---:|---:|---|
| `vascular_flow/straight_tube` | 1,005 | 720 | `0.763084` | wall 0, inlet 1, outlet 2 |
| `vascular_flow/y_bifurcation` | 7,329 | 6,300 | `0.305198` | wall 0, inlet 1, outlets 2 and 3 |
| `vascular_flow/iga_wordmark` | 24,924 | 22,140 | `0.149585` | wall 0, inlet 1, outlet 2 |
| `vascular_flow/multispecies_pulse` | 2,010 | 1,620 | `0.750036` | wall 0, inlet 1, outlet 2 |
| `neuron_transport/branched_neurite` | 7,329 | 6,300 | `0.305198` | wall 0, inlet 1, outlets 2 and 3 |

The large `neuron_transport/nmo_06840_bifurcation` regression was validated
separately on 2026-08-26. It produces 29,238 nodes and 25,920 elements; its
12-rank packed geometry check reported `minimum_detJ=4.5706e-11`, zero bad
elements, zero bad quadrature samples, and wall/inlet/outlet labels 0/1/2/3.
See its [case README](neuron_transport/nmo_06840_bifurcation/README.md) for the
resource warning, exact commands, provenance, and transport evidence.

## Steady vascular flow

The straight and Y cases converged after six Newton updates. The wordmark is a
Stokes-limit showcase with density `1e-5`, viscosity `1`, and inlet-profile
scale `0.005`; it converged after five updates. It is a numerical visualization
case, not a physiological blood-flow parameter set.

| Case | Linear iterations | Velocity L2 | Pressure L2 | Relative mass imbalance | Relative divergence error |
|---|---:|---:|---:|---:|---:|
| `straight_tube` | 404 | `0.166093` | `11.8216` | `2.83042e-7` | `1.45381e-7` |
| `y_bifurcation` | 1,509 | `0.332098` | `58.8361` | `8.36495e-8` | `4.24139e-8` |
| `iga_wordmark` | 1,304 | `0.428108` | `2439.12` | `9.79398e-8` | `9.20554e-8` |

The Y case produced outward flows `1.868446e-3` and `1.868436e-3` at the two
terminal labels. The final net outward flow was `3.20202e-10`.

The wordmark run took `755.45 s` on two MPI ranks and used PETSc FGMRES with
restart 200 plus overlapping additive Schwarz (`PCASM`, overlap 2, local ILU
level 1). Its final nonlinear residual was `1.78285e-7`, below the
`3.7574e-7` tolerance, and its net outward flow was `1.98302e-11`.

## Branched neuron transport

The configured `N0`/`Nplus` showcase completed 40 physical steps on eight MPI
ranks in 1,215 Krylov iterations. Assembly took `8.57693 s`, the time-loop
solves took `2.51929 s`, and the final coefficient L2 norm was `66.0504`. The
output contained exactly 7,329 node rows and the `.fields` file recorded `N0`
followed by `Nplus`.

At `t=4`, `N0` ranged from `2.6032e-4` to `1`, and `Nplus` from `9.1568e-4`
to `2.36306`. The worst early transient undershoots were `-1.87144e-2` for
`N0` and `-9.48942e-2` for `Nplus`. This is an execution, coupling, and
visualization case, not a positivity-preserving or biologically calibrated
transport benchmark; the undershoots must not be interpreted as physical
concentrations.

## Large NMO neuron transport

The two-step `neuron_transport/nmo_06840_bifurcation` acceptance run used 12
OpenMPI ranks. Assembly took `26.6397 s`, the time-loop solve took `109.36 s`,
and 10,520 Krylov iterations produced final coefficient L2 `330.035`. The
output contained exactly 29,238 finite rows with field order `N0`, `Nplus`.
Final ranges were `0.0403333` to `1.0001424` for `N0` and `1.7740167` to
`2.0012732` for `Nplus`. This short two-step horizon validates execution and
output, not biological steady state. CUDA runtime validation remains pending
for this large public case.

## Pulse multispecies showcase

The curved `multispecies_pulse` case was packed for eight MPI ranks and ran
eight backward-Euler flow steps followed by eight six-species transport steps.
It is a transit-matched, low-density visualization case (`dt=1`, pulse period
`8`, inlet-profile scale `2`, `density=1e-8`, and `viscosity=1`), not a
dimensional blood calibration. All flow steps converged. At the last step the
nodal speed had mean `0.860707` and maximum `2.22747`, so the eight-second
animation spans an advective distance comparable to the curved vessel length.

The PVD collections contain the actual initial state at `t=0` and all eight
solved states. At `t=8`, glucose ranged from `4.57843` to `6.02892`, oxygen
from `0.0552356` to `0.14`, and lactate from `0.492938` to `0.95707`.
Relative to useful visualization thresholds, 1,880 of 2,010 nodes had glucose
above `5.1`, 742 had oxygen above `0.122`, and 1,968 had lactate below `0.9`.
Thus the inlet signal reaches the vessel interior instead of remaining on the
inlet face. Increasing the showcase diffusivities reduced the earlier
high-Peclet overshoot to less than 0.5% for glucose and 1.5% for lactate.

The CPU and SM 89 CUDA transport paths were compared after one step on an
NVIDIA GeForce RTX 4080 SUPER; the combined six-field relative L2 difference
was `1.11722e-7`. Both backends wrote readable VTU/PVD output containing the
six solved species plus `pO2`, `pCO2`, `pH`, `SaO2`, `total_oxygen`, and
`hematocrit`.

## ParaView rendering

The gallery and README hero PNGs were produced with
`scripts/render_example.py` and `pvbatch --force-offscreen-rendering`. The
renderer deep-copies the legacy VTK mesh, validates result row counts and field
ordering, attaches CPU coefficients as point arrays, and colors a center-plane
or explicitly positioned slice. The wordmark hero uses `--slice-z 0` to show
the letter strokes while leaving its out-of-plane connector arcs as a
translucent context surface.

The two multispecies GIFs were rendered from the native PVD collections with
`scripts/render_multiphysics_gif.py`. Their transfer functions use each
field's range over the entire time series, so a fixed color has the same value
in every frame and transport fronts cannot be hidden by per-frame rescaling.
The branched-neurite Nplus animation uses the same fixed-range rule through
`scripts/render_transport_gif.py`.

These images are piecewise-linear control-mesh visualizations. Quantitative
flow values above come from `iga_flow_validate` using the packed IGA basis and
quadrature rather than measurements from the rendered pixels.

CUDA was rerun only for the focused one-step pulse showcase above. The older
gallery cases were not rerun on CUDA; their scheduled parity evidence remains
in the solver validation documents.
