# Benchmark Results

## Cross-case end-to-end runtime comparison

The following reported runtime comparison covers three branching-neuron cases. It contrasts the legacy [NeuronTransportIGA](https://github.com/EngineerEricXie/NeuronTransportIGA) CPU workflow, HexSim CPU, the optimized TubularFlowIGA CPU workflow, and the TubularFlowIGA V100 backend. Each value is a case-specific end-to-end runtime in seconds; it is not a hardware-scaling study.

| Case | Nodes / elements | NeuronTransportIGA CPU | HexSim CPU | TubularFlowIGA CPU | TubularFlowIGA V100 |
| --- | ---: | ---: | ---: | ---: | ---: |
| bifurcation | 14,565 / 12,780 | 442.68 s | 71.98 s | 26.75 s | **7.87 s** |
| NMO_54499_shrink | 16,344 / 14,040 | 589.32 s | 69.61 s | 33.03 s | **16.39 s** |
| NMO_66748_subtree | 57,456 / 50,940 | 12,127.80 s | 137.36 s | 98.77 s | **78.52 s** |

Relative to the legacy CPU workflow, the V100 runtimes correspond to **56.2x**, **36.0x**, and **154.5x** speedups for the three cases, respectively. Relative to the optimized TubularFlowIGA CPU workflow, they correspond to **3.40x**, **2.02x**, and **1.26x** speedups. Hardware, rank count, solver options, and case snapshots must be kept fixed when reproducing or extending this comparison.

These measurements were collected on PSC Bridges-2 with FP64 arithmetic.
Packing and text preprocessing are excluded unless stated. Timings are
hardware- and solver-option-specific; reproduce them from compute nodes with
the same case snapshot.

Throughout this document, **legacy** refers to the original workflow and solver
implementation in [NeuronTransportIGA](https://github.com/EngineerEricXie/NeuronTransportIGA).

## C++ control-mesh generation

The dependency-free C++ replacement was compared against the current MATLAB
workflow, not the stale VTK snapshots in some case directories.

| Case | Points | Hexes | C++ wall time | Peak RSS | Minimum scaled J |
| --- | ---: | ---: | ---: | ---: | ---: |
| Cylinder | 4,221 | 3,600 | 0.06 s | 4.5 MiB | 0.763084 |
| Bifurcation | 14,565 | 12,780 | 0.21 s | 6.6 MiB | 0.382533 |
| NMO_54499_new | 35,949 | 31,680 | 0.76 s | 10.8 MiB | 0.0315119 |

For cylinder, bifurcation, and three-bifurcation regressions, connectivity and
labels matched MATLAB exactly; coordinate differences were at most the
`1e-6` VTK text precision. The large NMO mesh had no invalid control elements.
Its extracted IGA geometry also passed all 2,027,520 determinant samples with
`min(detJ)=1.27675e-7`.

The spline extractor was subsequently changed from whole-mesh dense storage to
chunked construction and ordered streaming. On eight OpenMP threads:

| NMO spline extraction | Wall time | Peak RSS | Relative result |
| --- | ---: | ---: | ---: |
| Legacy | 62.50 s | about 2.19 GiB | baseline |
| Chunked OpenMP | 9.34 s | 158.1 MiB | 6.69x faster, 92.9% less memory |

The four legacy output files and packed `.ntiga` database matched byte for
byte. A repeated pre-cleanup optimized run took 7.27 s; the table uses the
conservative final measurement from a different compute node. See the
[mesh validation report](MESH_CPP_VALIDATION.md) for methodology.

The native sparse-cache path removes dense text generation and parsing:

| NMO preprocessing stage | Legacy path | Sparse-cache path | Speedup |
| --- | ---: | ---: | ---: |
| Spline extraction | 62.50 s | 3.74 s | 16.7x |
| Database packing | 30.02 s | 4.47 s | 6.72x |
| Combined | 92.52 s | 8.21 s | 11.3x |

The cache-only extractor used 160.3 MiB peak RSS and the packer used 4.2 MiB.
The 329 MB cache replaces 537 MB of `cmat.txt` and `bzpt.txt`. The resulting
31,680-element `.ntiga` database matched the legacy database byte for byte.

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
