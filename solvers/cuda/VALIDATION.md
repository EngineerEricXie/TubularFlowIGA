# CUDA Validation

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
export IGA_CUDA_ROOT=/ocean/projects/mch260002p/thsieh1/NeuronTransportIGA-CUDA
export IGA_CASE_ROOT=/ocean/projects/mch260002p/thsieh1/NeuronTransportIGA
cd "$IGA_CUDA_ROOT"
sbatch --export=ALL,IGA_CASE_ROOT="$IGA_CASE_ROOT" slurm/validate_v100.sbatch
sbatch --export=ALL,IGA_CASE_ROOT="$IGA_CASE_ROOT" slurm/validate_cpu_transport.sbatch
sbatch --export=ALL,IGA_CASE_ROOT="$IGA_CASE_ROOT" slurm/nmo_full_v100.sbatch
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
