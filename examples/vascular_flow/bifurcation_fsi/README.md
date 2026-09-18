# Idealized bifurcation FSI showcase

This example is being validated. It couples a three-dimensional Y-shaped blood
lumen to the existing small-displacement pre-tensioned membrane model. The inlet
and both outlet rims are clamped; the remaining wall participates in two-way
traction/kinematics coupling with dynamic Aitken relaxation. The structural
model is a scalar normal-displacement membrane, not a calibrated hyperelastic
arterial wall. The driver uses one MPI rank with optional OpenMP assembly.

The parent diameter is 4 mm and the daughter diameters are approximately 3.3 mm.
The inlet and outlets lie at x = -9 and +9 mm. A smooth union connects the parent
and two oblique daughter cylinders. The generated triangulation is closed and
oriented, with wall label 7 and inlet/lower/upper port labels 1/2/3.

Default fluid parameters are density 1060 kg/m³ and viscosity 0.0035 Pa·s.
A constant 20 Pa inlet-to-outlet pressure difference drives startup from rest.
This is a transient startup demonstration, not a cardiac pulse. Membrane
parameters are areal mass 0.3 kg/m², damping 30 kg/(m²·s), pretension 10 N/m,
and foundation 2e7 N/m³; these are illustrative, not patient-calibrated.

The transient wall condition uses the existing Nitsche impedance stabilization:
viscous gamma 2 and inertial gamma 1, adding `16 * rho * h_n / dt` to the wall
penalty coefficient. `-demo_wall_inertial_gamma 0` reproduces the viscous-only
configuration that exceeded the leakage gate at step 3. The stabilized eight-step
run is being validated. This numerical stabilization retains the prescribed
material wall velocity and physical fluid/membrane parameters; see
[the operator definition](../../../docs/architecture/MOVING_DOMAIN_ARCHITECTURE.md#transient-immersed-wall-impedance).

Generate inputs in a separate workspace using Python with NumPy and scikit-image:

```bash
python3 examples/vascular_flow/bifurcation_fsi/generate_surface.py /shared/case/surface.txt
make -C solvers/cpu bifurcation_fsi PETSC_DIR="$PETSC_DIR" PETSC_ARCH="$PETSC_ARCH"
# Run on an allocated compute node, with thread counts matching the allocation.
OMP_NUM_THREADS=8 IGA_ASSEMBLY_THREADS=8 \
  mpiexec -np 1 --map-by slot:PE=8 --bind-to core \
  solvers/cpu/bifurcation_fsi /shared/case/surface.txt /shared/case/run \
  -demo_steps 8 -demo_grid 16 -demo_dt .01 -demo_pressure 20
```

The output directory must be new. On a successful run, `fsi.pvd` references
the accepted fluid and membrane VTU snapshots. Failed runs can retain partial
output; earlier driver builds also wrote snapshots before the conservation
gates and therefore retained rejected steps for diagnosis. File presence alone is not acceptance; require
a complete `run.json` marked `passed` and the independent checks below. `history.csv` contains displacement, port fluxes and moving-domain
conservation, including raw volumetric defects and their normalization scale;
`coupling.csv` contains the actual coupling residual histories.
Fields are exported from committed solver states. Fluid VTU data are integration
point samples, so continuous slice/streamline views require interpolation.

Acceptance requires nonzero wall motion below 0.1 mm, strong-coupling convergence,
normalized moving mass defect and wall-relative leakage below 0.03, and discrete
moving-wall continuity defect below 1e-8. Mesh/time convergence for quantitative
physiology is separate from this demonstration.

The mass-defect and wall-relative-leakage diagnostics divide their absolute
volumetric defects by the maximum of the configured reference flow (1e-6 m³/s),
the absolute backward-Euler volume rate, and the absolute total fluid and
material surface fluxes. Their values are therefore not percentages of inlet
flow. The discrete moving-wall continuity diagnostic has its own normalization;
it checks the assembled continuity balance and does not replace geometric
mass/leakage checks or a mesh-convergence study.

On Bridges-2, submit the dedicated batch wrapper from the repository root. It
requests a 24-hour wall-time limit, builds and runs inside the allocation, and
preserves source, binary hashes,
stdout/stderr, hardware, and accounting in a fresh job directory. Explicit
exports are necessary when the login environment sets `SBATCH_EXPORT=NONE`:

```bash
export IGA_REPO_ROOT="$PWD"
export IGA_FSI_SHOWCASE_ROOT=/ocean/projects/PROJECT/USER/TubularFlowIGA-fsi-showcase
# PETSC_DIR and PETSC_ARCH must identify an installation compatible with
# Open MPI 4.0.5 / GCC 10.2.0 on Bridges-2.
sbatch -A PROJECT \
  --export=IGA_REPO_ROOT,IGA_FSI_SHOWCASE_ROOT,PETSC_DIR,PETSC_ARCH,CXX,HDF5_CFLAGS,HDF5_ALL_CFLAGS,HDF5_LIBS,PKG_CONFIG_PATH \
  solvers/cpu/slurm/bifurcation_fsi.sbatch
```

After a successful run, independently verify saved fields and generate the
presentation (also on compute):

```bash
python3 examples/vascular_flow/bifurcation_fsi/verify.py /shared/job/results \
  --output /shared/job/verification.json
# With the VTK Python package available:
python3 examples/vascular_flow/bifurcation_fsi/verify_vtk.py /shared/job/results \
  --output /shared/job/vtk-verification.json
python3 examples/vascular_flow/bifurcation_fsi/render.py /shared/job/results /shared/job/visuals
```

The renderer only accepts a complete `run.json` marked `passed`. Animation
frames correspond to actual accepted time steps; no intermediate solution
states are invented. Center-plane views interpolate the sampled fluid fields.
The wall view labels its deformation magnification explicitly, while its
colour bar always reports the actual displacement in micrometres.

For an in-progress or failed run, a single numerically accepted step can be
inspected with `--preview-step 1`. This checks that step's recorded convergence,
flow signs and conservation gates and produces a PNG prominently labelled as
an incomplete-run preview. It does not generate a showcase animation or mark
the overall run as passed.

At finite coupling tolerance, the raw membrane solution and relaxed fluid
interface can differ slightly. Each retains its own accepted displacement
history: the fluid boundary velocity is the backward-Euler derivative of its
own accepted geometry. The structural state keeps its membrane time integrator
history. This prevents the accepted residual from becoming an inconsistent
fluid boundary velocity on the next step.

## Heartbeat visualization mode

For a functional ParaView demonstration, `bifurcation_heartbeat.sbatch` uses
one 1-second illustrative pressure cycle, 20 actual coupled time steps, a 16×16×5
background grid and cut-quadrature depth 1. This coarse case prioritizes a visible
animation; it is not the formal grid-16 conservation acceptance case.

The inlet pressure is `p(t) = 5 + 15 * (1 - cos(2*pi*t/1s))` Pa, ranging from
5 to 35 Pa. Both outlet pressures remain zero. The small-displacement membrane
and two-way traction/kinematic coupling compute the wall response; displacement
is not prescribed from the pulse. The cycle starts from rest, so it includes
startup and is not a periodic steady-state solution.

```bash
# Set PETSC_DIR/PETSC_ARCH as in the standard example, then submit from repo root.
sbatch -A mch260002p --export=ALL solvers/cpu/slurm/bifurcation_heartbeat.sbatch
```

Outputs are under `.nsvms_diagnostics/heartbeat/job-JOBID/results/`:

- Open **wall-display.pvd** for a wall-only animation with displacement magnified
  100 times. Press Apply and Play. The geometry is already magnified; do not add
  another Warp By Vector. Color by `normal_displacement_m` for signed physical
  displacement (metres). `display_displacement_scale` records the magnification.
- Open **fsi.pvd** for actual-scale fluid and membrane blocks. Fluid data are
  integration-point samples; use Points to see them directly.
- Keep every `step-*` directory together with the PVD files when downloading.
- `history.csv` records the actual pressure and `conservation_passed` at each step.
- `heartbeat-verification.json` checks finite fields, time series, exact display
  scaling and observed increases/decreases of mean signed wall displacement.

`-demo_visualization_only true` records conservation failures and continues.
Its `run.json` status is always `visualization_only`, never `passed`. Coupling
convergence, finite diagnostics and the small-displacement bound remain required;
port-flow direction is reported without imposing the startup-only sign gate.
The existing formal verifier and renderer still reject this status.

Parameters: `-demo_pressure` is the minimum pressure;
`-demo_pulse_amplitude` is the peak-minus-minimum excursion (default 0 preserves
constant forcing); `-demo_period` is seconds; `-demo_display_scale` changes only
exported display geometry; `-demo_cut_depth` defaults to the original depth 3.
