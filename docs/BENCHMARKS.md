# Benchmark Results

These measurements were collected on PSC Bridges-2 with FP64 arithmetic.
Packing and text preprocessing are excluded unless stated. Timings are
hardware- and solver-option-specific; reproduce them from compute nodes with
the same case snapshot.

## Legacy versus optimized CPU

Cylinder: 4,221 nodes and 3,600 elements.

| Measurement | Legacy CPU | Optimized CPU |
| --- | ---: | ---: |
| Comparable first nonlinear update | 199 s | about 12 s |
| First-update speedup | 1.0x | about 16.6x |
| Peak memory | 6.72 GiB | 1.51 GiB |
| Memory reduction | baseline | about 77.5% |
| Complete solve | not implemented | 94 s, 7 updates |

The legacy executable performed one hard-coded update and did not establish
nonlinear convergence. The table diagnoses implementation overhead; it must not
be described as a converged end-to-end speedup.

## CPU versus CUDA: cylinder

CPU used 8 MPI ranks; GPU job 42862078 used one L40S-48GB.

| Stage | CPU | CUDA | Speedup |
| --- | ---: | ---: | ---: |
| Navier-Stokes assembly, 8 passes | 87.7 s | 7.74 s | 11.3x |
| Navier-Stokes linear solves | 7.63 s | 8.57 s | 0.89x |
| Complete Navier-Stokes numerical work | 94.33 s | 16.32 s | 5.78x |
| One transport assembly | 4.17 s | 0.170 s | 24.5x |

CPU-to-CUDA relative L2 differences were `8.24e-10` for velocity,
`5.20e-10` for pressure, and `5.62e-6` for one transport step.

## CPU versus CUDA: corrected NMO_54499_new

This snapshot has 35,949 nodes. All 2,027,520 geometry samples were positive;
`min(detJ)=1.2768e-7`. Jobs 42863526 and 42863716 compared 16 CPU ranks with
one V100-32GB.

| Stage | CPU | CUDA | Speedup |
| --- | ---: | ---: | ---: |
| Navier-Stokes assembly | 191.36 s | 34.57 s | 5.54x |
| Navier-Stokes solves | 110.50 s | 68.28 s | 1.62x |
| Navier-Stokes total | 301.86 s | 102.86 s | 2.93x |
| Transport assembly | 21.82 s | 1.19 s | 18.3x |
| Transport solves, 300 steps | 240.00 s | 175.36 s | 1.37x |
| Transport total | 261.82 s | 176.55 s | 1.48x |
| Coupled numerical work | 563.68 s | 279.40 s | 2.02x |

| Peak allocation | CPU RSS | CUDA host | CUDA device | Combined reduction |
| --- | ---: | ---: | ---: | ---: |
| Navier-Stokes | 7.05 GiB | 1.09 GiB | 2.69 GiB | about 46% |
| Transport | 5.28 GiB | 1.09 GiB | 1.50 GiB | about 51% |

Velocity, pressure, and transport relative L2 differences were `5.33e-6`,
`6.07e-6`, and `7.56e-6`.

Assembly benefits most from CUDA. Large-case GMRES remains limited by
block-Jacobi preconditioning, so total transport speedup is lower than assembly
speedup. The CUDA backend is single-GPU; multi-GPU scaling is not implemented.
