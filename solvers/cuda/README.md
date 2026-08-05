# TubularFlowIGA CUDA Backend

A project-owned FP64 CUDA implementation of the TubularFlowIGA
Navier-Stokes and transport pipeline. CUDA Runtime and cuBLAS provide device
memory and vector primitives; Bezier extraction, quadrature, stabilized weak
forms, sparse block assembly, preconditioning, and restarted GMRES remain
implemented in this repository.

The numerical reference and database packer live in [../cpu](../cpu).
The original implementation is
[NeuronTransportIGA](https://github.com/EngineerEricXie/NeuronTransportIGA).

## Capabilities

One executable exposes four modes:

| Command | Purpose |
| --- | --- |
| `iga_cuda device-info` | Report the selected CUDA device and compiled capabilities |
| `iga_cuda mesh-check DATABASE.ntiga` | Validate all 4x4x4 element quadrature Jacobians |
| `iga_cuda navier-stokes ...` | Run the stabilized four-field steady solve |
| `iga_cuda transport ...` | Run the transient two-field transport solve |

The database can be packed for any CPU MPI rank count. The single-GPU reader
loads each element exactly once and ignores CPU ownership records.

## GPU-specific optimizations

- Flatten sparse Bezier extraction, connectivity, quadrature tables, and
  geometry into contiguous device arrays.
- Precompute geometry transforms at every quadrature point once before assembly.
- Store unique node adjacency as block CSR: 2x2 blocks for transport and 4x4
  blocks for Navier-Stokes.
- Tile element-pair work in groups of 256; evaluate basis data in shared memory,
  accumulate quadrature contributions in registers, and use FP64 atomics only
  for the final cross-element reduction.
- Avoid the repeated multi-gigabyte COO indices and values that a direct PETSc
  GPU assembly port would create.
- Keep state vectors, sparse matrices, and GMRES workspaces resident on device.
- Assemble transport operators once and reuse the block-Jacobi inverse,
  workspace, and state across all time steps.
- Rebuild only the Navier-Stokes Jacobian-dependent data each nonlinear update.
- Use project-owned block-CSR SpMV and block inversion; cuBLAS is limited to
  FP64 dot, norm, copy, scale, and AXPY operations.
- Tune restarted GMRES separately: restart 200 for Navier-Stokes and restart 50
  for transport. Polynomial-Jacobi alternatives were measured and rejected
  after failing or slowing convergence.

See [ARCHITECTURE.md](ARCHITECTURE.md) for kernel and data-layout details.

## Measured CPU-to-GPU improvement

Measurements were collected on PSC Bridges-2. Timings exclude input packing
unless explicitly described. CPU and GPU fields were compared numerically, not
only by runtime.

### Cylinder: 4,221 nodes, 3,600 elements

The CPU baseline used 8 MPI ranks. GPU job `42862078` used one L40S-48GB.

| Stage | CPU, 8 ranks | CUDA, 1 L40S | Speedup |
| --- | ---: | ---: | ---: |
| Navier-Stokes assembly, 8 passes | about 87.7 s | 7.74 s | 11.3x |
| Navier-Stokes linear solves | about 7.63 s | 8.57 s | 0.89x |
| Complete Navier-Stokes numerical work | 94.33 s | 16.32 s | 5.78x |
| One transport assembly | 4.17 s | 0.170 s | 24.5x |

Restart 200 reduced Navier-Stokes GMRES iterations from 12,417 to 3,486 without
changing the nonlinear solution. Relative to the CPU result, velocity and
pressure relative L2 errors were `8.24e-10` and `5.20e-10`. The one-step
transport relative L2 difference was `5.62e-6` at the same `1e-8` solver
tolerance.

For historical context, the legacy CPU solver required 199 s and 6.72 GiB for
one non-converged hard-coded update. That is not an apples-to-apples GPU
baseline; use the converged CPU implementation above for speedup claims.

### Corrected NMO_54499_new: 35,949 nodes

Jobs `42863526` and `42863716` compared 16 CPU ranks with one V100-32GB. All
2,027,520 geometry samples were positive with `min(detJ)=1.2768e-7`.

| Stage | CPU, 16 ranks | CUDA, 1 V100 | Speedup |
| --- | ---: | ---: | ---: |
| Navier-Stokes assembly | 191.36 s | 34.57 s | 5.54x |
| Navier-Stokes linear solves | 110.50 s | 68.28 s | 1.62x |
| Complete Navier-Stokes numerical work | 301.86 s | 102.86 s | 2.93x |
| Transport assembly | 21.82 s | 1.19 s | 18.3x |
| Transport linear solves, 300 steps | 240.00 s | 175.36 s | 1.37x |
| Complete transport numerical work | 261.82 s | 176.55 s | 1.48x |
| Coupled numerical work | 563.68 s | 279.40 s | 2.02x |

| Peak memory | CPU host RSS | CUDA host | CUDA device | Combined reduction |
| --- | ---: | ---: | ---: | ---: |
| Navier-Stokes | 7.05 GiB | 1.09 GiB | 2.69 GiB | about 46% |
| Transport | 5.28 GiB | 1.09 GiB | 1.50 GiB | about 51% |

GPU wall times including preprocessing and output were 110.70 s for
Navier-Stokes and 182.93 s for transport. Relative to CPU, the velocity,
pressure, and transport relative L2 errors were `5.33e-6`, `6.07e-6`, and
`7.56e-6`.

These results use the corrected 35,949-node snapshot, not the older mesh with
negative Jacobians. Full evidence and tuning history are in
[VALIDATION.md](VALIDATION.md).

## Repository layout

- `include/`: CUDA runtime wrappers, block CSR, device mesh, GMRES, and kernels.
- `src/iga_cuda.cu`: CLI, database loading, assembly, solvers, and VTK output.
- `slurm/`: Bridges-2 build, regression, CPU-reference, and large-case jobs.
- `ARCHITECTURE.md`: kernel organization and scaling boundary.
- `VALIDATION.md`: job IDs, accuracy gates, timing, memory, and rejected tuning.

Generated executables, Slurm `.out` files, databases, and simulation outputs are
ignored.

## Requirements and clone layout

- CUDA toolkit 12.4 or compatible.
- A CUDA GPU with FP64 support. The fat binary targets SM 70, 80, 89, and 90.
- A C++17-compatible host compiler.
- cuBLAS.
- The in-tree CPU backend, which provides the shared `IgaDatabase.hpp` format reader.

~~~text
TubularFlowIGA/
|-- solvers/cpu/
|-- solvers/cuda/
`-- meshgeneration/
~~~

## Build

~~~bash
export TUBULARFLOWIGA_ROOT=/ocean/projects/mch260002p/thsieh1/TubularFlowIGA
export IGA_CUDA_ROOT="$TUBULARFLOWIGA_ROOT/solvers/cuda"
module load cuda/12.4.0
make -C "$IGA_CUDA_ROOT"
"$IGA_CUDA_ROOT/iga_cuda" device-info
~~~

Do not run GPU commands on a login node. The build can be submitted with:

~~~bash
cd "$IGA_CUDA_ROOT"
sbatch slurm/build_v100.sbatch
~~~

## Run on a GPU compute node

~~~bash
export IGA_CASE_ROOT=/ocean/projects/mch260002p/thsieh1/NeuronTransportIGA

interact -A mch260002p -p GPU-shared --gres=gpu:v100-32:1 -t 00:30:00
module load cuda/12.4.0

"$IGA_CUDA_ROOT/iga_cuda" mesh-check DATABASE.ntiga
"$IGA_CUDA_ROOT/iga_cuda" navier-stokes DATABASE.ntiga CASE_DIR 8 velocity.txt
"$IGA_CUDA_ROOT/iga_cuda" transport DATABASE.ntiga CASE_DIR 300 concentration.txt velocity.txt
~~~

Each solver writes its text interchange field and `OUTPUT.vtk`. Navier-Stokes
VTK contains velocity and pressure; transport VTK contains `N0` and `Nplus`.

## Reproduce the Bridges-2 validations

Case data and generated databases remain outside Git:

~~~bash
export TUBULARFLOWIGA_ROOT=/ocean/projects/mch260002p/thsieh1/TubularFlowIGA
export IGA_CUDA_ROOT="$TUBULARFLOWIGA_ROOT/solvers/cuda"
export IGA_CASE_ROOT=/ocean/projects/mch260002p/thsieh1/NeuronTransportIGA
cd "$IGA_CUDA_ROOT"

sbatch --export=ALL,IGA_CASE_ROOT="$IGA_CASE_ROOT" slurm/validate_v100.sbatch
sbatch --export=ALL,IGA_CASE_ROOT="$IGA_CASE_ROOT" slurm/validate_cpu_transport.sbatch
sbatch --export=ALL,IGA_CASE_ROOT="$IGA_CASE_ROOT" slurm/nmo_full_v100.sbatch
~~~

V100-32GB is the preferred current production target because all validated
cases fit within 2.69 GiB device memory and it offers strong FP64 throughput.
H100 is supported; L40S is functional but not preferred for FP64 production.

## Reproducibility notes

- Record GPU model, CUDA version, driver, host compiler, and Slurm job ID.
- Run `mesh-check` before timing assembly or Krylov solves.
- Use the same corrected `.ntiga` snapshot and solver tolerances for CPU/GPU
  comparisons.
- Report assembly, linear solve, total wall time, host RSS, and peak device
  allocation separately.
- Compare velocity, pressure, and both transport fields after performance runs.

## Current limitations

- Execution is single process and single GPU; multi-GPU domain decomposition is
  not implemented.
- Block-Jacobi is weaker than CPU PETSc local ILU for the large transport case,
  limiting solve speedup even though assembly is much faster.
- The CUDA backend consumes the in-tree CPU database header; changes to the format must be validated on both backends.
- Validation focuses on FP64 scientific GPUs and the current IGA formulations;
  other discretizations and mixed precision are not implemented.
