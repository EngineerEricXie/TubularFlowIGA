# CUDA Validation

## Schema-v3 integration regression

After adding the native 1D subsystem and the schema-v3 dimension dispatcher,
the unchanged 3D CUDA backend was rebuilt and rerun on 2026-08-25. The WSL host
used an RTX 4080 SUPER (SM 89, 15.99 GiB), the `tubularflow-cuda` Conda
environment, and `nvcc` 12.6.85. A clean SM 89 build completed without
warnings. `device-info` and CUDA runtime execution required the environment's
`targets/x86_64-linux/lib` directory in `LD_LIBRARY_PATH`.

The source-only public straight-vessel and straight-neurite examples were
freshly meshed and packed into separate 1,005-node, 720-element databases. Both
device `mesh-check` runs reported `min(detJ)=9.0098909917245762e-5` and zero bad
samples. The steady Navier--Stokes solve converged at Newton check 6 with no
singular diagonal blocks. Relative to the two-rank PETSc result:

```text
velocity_relative_l2=1.9942021597698123e-10
pressure_relative_l2=1.7826968415764227e-10
relative_mass_imbalance=2.830318713527568e-7
relative_divergence_theorem_error=1.4538127906482656e-7
```

The two-step configured `N0`/`Nplus` transport solve completed in 318 CUDA
Krylov iterations with no singular diagonal blocks. Its combined CPU/CUDA
relative L2 difference was `3.9995721083720065e-6`, below the `1e-5` gate, and
node and field ordering matched. Restarting from the step-1 CUDA checkpoint
changed the final field by `3.2091905463522199e-13` relative L2, below the
`1e-12` restart gate.

## Transient and coupled gates

The backward-Euler flow, time-indexed output, raw checkpoint/restart,
R/RC/RCR outlet coupling, and snapshot-series transport paths passed on
2026-08-24. The WSL host used one NVIDIA GeForce RTX 4080 SUPER (SM 89,
15.99 GiB) and CUDA Toolkit 12.6.3 installed in the `tubularflow-cuda` Conda
environment. `nvcc` reported release 12.6, V12.6.85. The CUDA target built for
SM 89 without warnings; `device-info` reported the expected `1/64` FP64 ratio.

The 2,211-node, 1,800-element mesh check evaluated every volume quadrature
sample with minimum `detJ=8.1088652284857333e-5` and zero invalid samples.
The eight-step straight-tube case covered two periods (`dt=0.1`, period `0.4`):

```text
preprocess_s=0.20281
assembly_s=55.4066
linear_s=181.265
total_linear_iterations=35990
state_l2=277.886
velocity_l2=1.53637
pressure_l2=277.882
gpu_used_gib=1.5083
maximum_relative_mass_imbalance=4.6845861945309674e-07
maximum_cycle_velocity_relative_l2=6.7578834255151922e-06
maximum_womersley_volume_relative_l2=0.046486537001491426
cpu_cuda_velocity_relative_l2=4.2926346639797907e-11
cpu_cuda_pressure_relative_l2=2.3669895670887319e-11
```

Restarting the run from its step-4 raw checkpoint changed final velocity by
`1.1524231282768408e-14` and pressure by `1.0965863268123281e-15` relative L2.
The result is below the `1e-12` restart gate. This consumer GPU is correctness
evidence, not a performance baseline for V100/H100-class FP64 hardware.

The R outlet converged to `Q=3.70264e-3` and `p=3.70264e-6`; its relative mass
imbalance was `4.88703e-7`. CPU/CUDA velocity and pressure differences were
`3.11089e-14` and `7.82448e-15`. RC and RCR completed the same fixed-point path.
The RCR checkpoint stored capacitor pressure `2.9386022758583951e-6` at step 1,
restored it at `t=0.1`, and advanced step 2 to flow `3.71854e-3`, applied
pressure `3.71825e-6`, and capacitor pressure `2.97454e-6`.

A two-step 2,211-node snapshot-series transport run matched CPU to
`8.170084106094566e-12` and matched the uninterrupted CUDA result after a
step-1 restart to `8.6877684819497343e-20`. The full 14,565-node,
12,780-element bifurcation then completed the same two-snapshot gate in 217
Krylov iterations with final L2 `24.9893`. Assembly and solve took `1.87246`
and `0.280954` seconds, peak device use was `1.64697 GiB`, CPU/CUDA relative L2
was `3.7467328704604612e-12`, and restart relative L2 was
`6.6421720189409289e-20`.

The dependency-free components used by CUDA also have unit coverage for
waveform/manifest parsing, checkpoint metadata, R/RC/RCR updates, boundary
surface-flow quadrature, and natural pressure-traction quadrature under
`make cpu-test`.

`slurm/validate_transient_v100.sbatch` automates the same build, mesh, restart,
mass, transport, optional CPU/CUDA field, multi-cycle, and Womersley checks for
a scheduled cluster GPU.

The script enforces default final-flow limits of `1e-2` relative mass
imbalance and `1e-6` relative divergence-theorem error. With
`IGA_CARDIAC_PERIOD`, `IGA_FLOW_DT`, and `IGA_FLOW_STEPS`, it also constructs a
manifest for all generated snapshots and gates maximum cycle-to-cycle velocity
relative L2 at `1e-3`. Each tolerance has a named environment override in the
[Bridges-2 guide](../../docs/BRIDGES2.md).
`IGA_WOMERSLEY_CONFIG` adds a physical-volume Womersley relative L2 gate
(default `5e-2`) by evaluating the numerical spline and analytical solution at
element quadrature points. The case-directory/manifest variables are retained
only for an independently projected coefficient reference; raw analytical
point samples must not be treated as IGA coefficients.

## Configured Surface-Term Parity

Job `44364589` ran the public 4,221-node, 3,600-element smoke database on one
Bridges-2 Tesla V100-SXM2-32GB with CUDA 12.4. The tracer wall used
`D grad(c) dot n = 0.5(2-c)`; both backends ran two time steps from the same
version 4 `.ntiga` database and configuration.

```text
CPU final L2:  31.3588
CUDA final L2: 31.3588
relative L2:   6.0297227463691488e-08
CUDA assembly: 0.120425 s
CUDA solve:    0.0719532 s
CUDA iterations: 99
singular diagonal blocks: 0
```

The configured surface-term parity gate was `1e-5`; the run passed. A version
3 database with the same surface configuration is rejected because it has no
packed element-face labels.

## Configured Waveform Parity

Job `44366806` ran a two-step constant-waveform tracer case on one
Tesla V100-SXM2-32GB. A factor-two inlet waveform produced exactly twice the
CPU baseline field; CUDA matched the CPU waveform result:

```text
CPU final L2:  28.91
CUDA final L2: 28.91
relative L2:   3.5232387432606833e-08
CUDA assembly: 0.124109 s
CUDA solve:    0.0787549 s
CUDA iterations: 109
singular diagonal blocks: 0
```

This job validates configured transport boundary updates. The later SM 89
acceptance above separately validates backward-Euler Navier-Stokes boundary
updates, multi-cycle integration, and restart.

## Cylinder Baseline

Final gate job `42863716` ran on one Bridges-2 V100-32GB with CUDA 12.4. The test used
`.nsvms_diagnostics/cylinder-8.ntiga` (4,221 nodes, 3,600 elements, maximum 64
basis functions per element).

Geometry preflight evaluated all 230,400 quadrature samples:

```text
minimum detJ: 4.5049436545644524e-05
non-positive samples: 0
```

The complete Navier–Stokes solve converged at Newton iteration 7. Relative to
the CPU v2 8-rank result:

```text
velocity relative L2 error: 8.241423351250334e-10
pressure relative L2 error: 5.199945255187425e-10
transport relative L2 error: 5.616752505130285e-06
transport tolerance: 1e-8 (matching CPU v2)
all VTK field gates: passed
```

Both implementations report state L2 `591.909` and the same Newton residual
sequence to the displayed precision.

## Performance and Solver Tuning

The initial restart-50 V100 result established an assembly speedup but required 12,417 block-Jacobi GMRES iterations. Job `42862078` validated restart 200 on one L40S-48GB:

| Stage | CPU v2, 8 ranks | CUDA, 1 L40S | Ratio |
| --- | ---: | ---: | ---: |
| NS assembly, 8 passes | about 87.7 s | 7.74 s | 11.3× faster |
| NS linear solves | about 7.63 s | 8.57 s | 0.89× |
| Complete NS numerical work | 94.33 s | 16.32 s | 5.78× faster |
| Transport assembly | 4.17 s | 0.170 s | 24.5× faster |

Restart 200 reduced total NS iterations to 3,486 without changing the nonlinear solution. The optimized L40S velocity and pressure relative L2 errors were `8.24e-10` and `5.20e-10`. A polynomial-Jacobi correction was tested and rejected because it did not converge for this saddle-point system.

For one transport step, the final GPU solve has relative L2 difference `5.62e-6` versus CPU v2. It uses the same `1e-8` relative tolerance as CPU v2.

## Reproduction

```bash
export TUBULARFLOWIGA_ROOT=/path/to/TubularFlowIGA
export IGA_CUDA_ROOT="$TUBULARFLOWIGA_ROOT/solvers/cuda"
export IGA_CASE_ROOT=/path/to/case-data
cd "$IGA_CUDA_ROOT"
sbatch --export=ALL,IGA_CASE_ROOT="$IGA_CASE_ROOT" slurm/validate_v100.sbatch
sbatch --export=ALL,IGA_CASE_ROOT="$IGA_CASE_ROOT" slurm/validate_cpu_transport.sbatch
sbatch --export=ALL,IGA_CASE_ROOT="$IGA_CASE_ROOT" slurm/nmo_full_v100.sbatch

export IGA_TRANSIENT_CASE_DIR=/path/to/transient-case
export IGA_TRANSIENT_DATABASE=/path/to/transient-case/case-1.ntiga
sbatch --export=ALL,IGA_TRANSIENT_CASE_DIR="$IGA_TRANSIENT_CASE_DIR",IGA_TRANSIENT_DATABASE="$IGA_TRANSIENT_DATABASE" \
  slurm/validate_transient_v100.sbatch
```


Validation outputs are written under `.nsvms_diagnostics/` and should not be committed. Cylinder regression and the full corrected `NMO_54499_new` coupled solve have passed.

## NMO_54499_new Short Run

Job `42862637` ran the corrected 35,949-node CPU-v2 snapshot on one H100-80GB. All 2,027,520 geometry samples were positive (`min detJ=1.2767547810311726e-07`). The unique block pattern contains 9,680,271 node blocks.

| Large-case stage | CPU v2, 16 ranks | CUDA, 1 H100 | Ratio |
| --- | ---: | ---: | ---: |
| Transport assembly | 21.82 s | 0.411 s | 53.1× faster |
| First transport solve | about 0.80 s average/step | 0.419 s | about 1.9× faster |
| First NS assembly | 45.91 s | 3.76 s | 12.2× faster |
| First NS linear solve | 28.31 s | 13.82 s | 2.05× faster |


## NMO_54499_new Full V100 Run

Jobs `42863526` and `42863716` ran the full corrected 35,949-node
snapshot on one V100-32GB. Times below exclude file preprocessing.

| Stage | CPU v2, 16 ranks | CUDA, 1 V100 | Ratio |
| --- | ---: | ---: | ---: |
| NS assembly | 191.36 s | 34.57 s | 5.54× faster |
| NS linear solves | 110.50 s | 68.28 s | 1.62× faster |
| Complete NS numerical work | 301.86 s | 102.86 s | 2.93× faster |
| Transport assembly | 21.82 s | 1.19 s | 18.3× faster |
| Transport linear solves | 240.00 s | 175.36 s | 1.37× faster |
| Complete transport numerical work | 261.82 s | 176.55 s | 1.48× faster |
| Coupled numerical work | 563.68 s | 279.40 s | 2.02× faster |

GPU wall times, including preprocessing and output, were 110.70 s for NS
and 182.93 s for transport, versus Slurm elapsed times of 308 s and 268 s
for CPU v2. Peak CPU host RSS was 7.05 GiB (NS) and 5.28 GiB (transport).
CUDA peak host/device use was about 1.09/2.69 GiB and 1.09/1.50 GiB,
respectively, reducing peak combined memory by about 46%.

Relative to CPU v2:

```text
velocity relative L2 error: 5.333774343577596e-06
pressure relative L2 error: 6.074109073211555e-06
transport relative L2 error: 7.562349896264679e-06
```

Transport restart 200 and polynomial Jacobi were measured and rejected.
The final restart-50 solve completed 300 steps in 71,493 iterations. Its
block-Jacobi preconditioner is weaker than CPU PETSc's local ILU, but GPU
assembly and sparse operations still make the complete transport solve faster.

All 2,027,520 geometry samples were positive. This case uses only 2.69 GiB
of device memory at peak, so V100 was selected over a longer-wait H100.
The example directory contains an older 35,748-node mesh; reproduce these comparisons with the matching `.nsvms_diagnostics/NMO_54499_new_v2/case-16.ntiga` snapshot.
