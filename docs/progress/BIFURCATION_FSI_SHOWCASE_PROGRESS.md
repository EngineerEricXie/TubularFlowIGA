# Bifurcating vessel FSI showcase

Status: **in progress**. The deliverable is an eight-step, 3D Y-vessel two-way
FSI run, independently checked fields, ParaView time series, and result-derived
images/animation. A validated first-step preview exists; it is not full-run
acceptance. The case and reproduction commands are in
[the example README](../../examples/vascular_flow/bifurcation_fsi/README.md).

The illustrative model uses a 4 mm parent and two approximately 3.3 mm daughters,
blood-like Newtonian fluid, constant 20 Pa pressure drive, and a small-displacement
normal membrane with clamped port rims. It models startup, not a cardiac cycle or
a patient-calibrated artery. Mesh/time convergence has not been established.

## Execution and provenance

- Base revision: `db8c66a2ccdbe90101a296b27eae6bac97c0d93b`, with the showcase and
  runtime fixes in the working tree. Per-job source copies, patches, binary
  hashes, module lists and linkage records identify the actual executed code.
- PSC account: `mch260002p`; workspace:
  `/ocean/projects/mch260002p/thsieh1/TubularFlowIGA-fsi-showcase`.
- One MPI rank, eight OpenMP threads, 15,200 MB allocation; all builds, simulations
  and substantial postprocessing run inside Slurm compute allocations.
- Job `45934683` ran on `r177`, a dual-socket AMD EPYC 7742 host
  (128 hardware cores; this job reserves eight). Steps 1 and 2 at 10/20 ms
  passed after six coupling rounds each; step 3 failed the mass/leakage gate. Step 2 RMS residual
  is 2.328153e-10 m, maximum membrane displacement 8.312715e-7 m, normalized mass
  defect 0.02770345, wall leakage 0.02770347 and continuity defect 7.080232e-14.
  Independent VTK 9.4.1 checks of both frames passed, including finite fields,
  metadata, speed norms, membrane geometry and inter-step backward-Euler velocity.
  Step 2 maximum fluid speed is 0.05120900 m/s and wall velocity consistency error
  1.439956e-20 m/s. Evidence: `job-45934683/two-frame-verification.json`.
  Full eight-step results remain pending.
- Build/linkage: GCC 10.2.0, Open MPI 4.0.5, PETSc 3.22.1 with real double
  scalars and 32-bit indices; C++17, `-O3`, OpenMP. This driver writes ASCII VTK
  XML and does not link HDF5. Evidence: `job-45934683/build.stdout`, `modules.txt`,
  `linkage.txt`, and the archived PETSc version/configuration headers.

Job `45934683` retained its original 16-hour limit. An attempt to extend it to
24 hours was denied by Slurm (`job-45934683/extend-time.stderr`); the job
was not interrupted by that request. The reproduction wrapper now requests 24 hours, allowing
more margin for the observed approximately 2.6-hour first step. The job later failed a numerical gate after 7:01:15;
it did not hit its wall-time limit.

## Evidence obtained before the replacement run

Job `45923794` on `r248` passed step 1 at 10 ms after six strong-coupling rounds,
then failed when starting step 2. Its `final.scheduler.json` records exit 1.
`results/history.csv` and `results/coupling.csv` contain the accepted first-step
metrics below; no complete-run `run.json` was produced.

| First-step measurement | Value |
| --- | ---: |
| Coupling RMS residual / threshold | 2.450596e-10 / 1.100000e-9 m |
| Maximum membrane displacement | 7.977435e-7 m |
| Maximum sampled fluid speed | 0.02765046 m/s |
| Sampled pressure range | 0.961953–19.437063 Pa |
| Outward inlet flux | -1.120167e-7 m³/s |
| Lower / upper outward outlet flux | 4.973569e-8 / 4.977545e-8 m³/s |
| Normalized moving mass defect | 0.01454020 |
| Normalized wall-relative leakage | 0.01454185 |
| Discrete moving-wall continuity defect | 6.106578e-14 |
| Elapsed time through accepted step | 9394.02 s |

The mass and leakage diagnostics use a normalization scale at least 1e-6 m³/s;
they are not percentages of inlet flow or a claim of mesh-independent accuracy.
The example README states the normalization definition. All original acceptance
thresholds are retained.
VTK 9.4.1 successfully read both datasets and their multiblock container:
723,546 fluid sample points and 806 membrane nodes / 1,582 triangles. Checks
covered finite arrays, time/unit metadata, speed norms, and membrane geometry
and backward-Euler velocity. Evidence: `job-45923794/first-frame-vtk-verification.json`.
The visually inspected `job-45923794/preview/preview-step-1.png` explicitly marks
its incomplete-run scope and 100× displayed deformation. No synthetic time
frames or complete-run animation were generated from that failed run.

## Fixes and validation

The initial slow attempt `45918820` lost PETSc matrix preallocation during an
empty final assembly. Complete stencil preallocation removed matrix growth;
`preallocation/comparison.json` records bit-identical residual/Jacobian actions
against the original implementation and zero matrix reallocations. The serial full transient-flow regression reached its final diagnostic checks
after passing the nonlinear solve/rollback invariants, then failed the obsolete
bitwise flux-sum assertion described below.

The second-step failure in `45923794` exposed distinct accepted histories:
the fluid commits the relaxed interface, while the membrane commits its raw
response. Deriving the next fluid boundary velocity from the structural history
violated the fluid geometry's backward-Euler check. The coordinator now uses the
fluid adapter's own accepted displacement history for predictor and relaxed
velocities, including that history in publication identity. Structural dynamics
and acceptance thresholds are unchanged. A two-step test deliberately retaining
a nonzero accepted residual fails on the original coordinator and passes on the
fix; the complete existing strong-coordinator unit suite also passes. Evidence:
`kinematic-history/comparison.json` and its baseline/candidate logs.

The replacement retains the same physical case and eight-step target, adds
per-iteration residual logging, and archives all modified runtime headers.
Full-run field verification, final rendering, numerical report closure and the
reviewable source commit remain outstanding.

## Third-step leakage rejection and stabilized rerun

Job `45934683` exited 1 at step 3 (30 ms). Coupling converged in four rounds
(RMS 4.010545e-10 m), but normalized moving mass defect 0.03761045 and wall
leakage 0.03761046 exceeded the unchanged 0.03 limits. Discrete continuity was
7.049006e-14. Only steps 1 and 2 passed. Evidence: `job-45934683/final.scheduler.json`,
`solver.stdout`, `solver.stderr`, and `results/{history.csv,coupling.csv,failure.txt}`.

The failing configuration used the default viscous-only wall penalty (gamma 2,
inertial gamma 0). The next complete run enables the already implemented transient
Nitsche impedance with inertial gamma 1 to test whether stronger enforcement of
the material-wall velocity controls startup leakage. Geometry, grid, density,
viscosity, membrane model, drive, eight time steps and acceptance thresholds stay
fixed. This is a candidate correction; success remains unproven until the rerun
passes. Raw conservation rates and their normalization scale are now exported
and independently reconciled by the field verifier.

The stabilized full eight-step job is `45961543`, submitted with one MPI rank,
eight CPU cores and a 24-hour limit. Submission/accounting evidence is
`submission-inertial-v1.{stdout,stderr,scheduler.json}` in the showcase workspace.
It started on `r328` (dual-socket AMD EPYC 7742), built successfully and began
step 1 with viscous/inertial gammas 2/1. The surface is byte-identical to the
failed run; source and five relevant header copies matched the worktree at startup. Evidence:
`job-45961543/startup-provenance.json`, `build.stdout`, `solver.stdout` and
`running.scheduler.json`. Step 1 has now passed after six coupling rounds;
maximum wall displacement is 7.676181e-7 m, normalized mass defect 0.005880372,
wall leakage 0.005882342 and continuity defect 7.747208e-14. Seven steps remain.
Later heartbeat visualization edits in the shared worktree are separate from
this archived constant-pressure executable; see `worktree-divergence.json`.

## Regression reduction-order assertion

The serial regression in step `45918565.0` exited 134 after 12:18:15 at test
line 716: it required a sum of per-label flux totals to match the independently
accumulated total bit-for-bit. An isolated replay of the exact fixture and field
state on both the original and preallocated runtime produced identical values:
label sum -0.00022377905964951722, direct total -0.00022377905964951673, difference
-4.878910e-19. Thus the discrepancy predates preallocation. Evidence:
`preallocation/{transient-regression.stderr,transient-regression.scheduler.json,flux-reduction-comparison.json}`.

The assertion now uses the floating-point accumulation error bound from the
actual summand count and absolute flux sum. Physical conservation gates, finite
difference thresholds and bitwise rollback invariants remain unchanged. The
complete existing OpenMP regression target passed with eight threads in step
`45918565.17` (exit 0:0, elapsed 02:07:58). Its finite-difference maximum error
was 1.74673e-11, zero-block maximum 0, linear relative residual 1.31502e-15,
open balance 1.66994e-7 and wall leakage 3.39035e-6. Logs and accounting:
`preallocation/transient-regression-openmp{,-build}.{stdout,stderr}` and
`preallocation/transient-regression-openmp.scheduler.json`.
