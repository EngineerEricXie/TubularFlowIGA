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

One executable exposes five modes:

| Command | Purpose |
| --- | --- |
| `iga_cuda device-info` | Report the selected CUDA device and compiled capabilities |
| `iga_cuda mesh-check DATABASE.ntiga` | Validate all 4x4x4 element quadrature Jacobians |
| `iga_cuda navier-stokes ...` | Run stabilized steady or backward-Euler four-field flow |
| `iga_cuda transport ...` | Run the transient two-field transport solve |
| `iga_cuda solve ...` | Run configured 1–8-field transport through the generic operator graph |

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
- HDF5 development headers and library. Temporal VTKHDF uses serial HDF5 calls
  on the host.
- The in-tree CPU backend, which provides the shared `IgaDatabase.hpp` format reader.

~~~text
TubularFlowIGA/
|-- solvers/cpu/
|-- solvers/cuda/
`-- meshgeneration/
~~~

## Build

From the repository root:

~~~bash
make cuda CUDA_ARCHS="70 80 89 90"
./solvers/cuda/iga_cuda device-info
~~~

Set `CUDA_ARCHS` to the compute capabilities needed at the deployment site.
Run `device-info` and all solver modes only on a CUDA-capable host.

On WSL systems where `nvidia-smi` works but `nvcc` is absent, the Windows
driver is already exposed to WSL but a Linux CUDA toolkit is still required.
One non-root installation used for the SM 89 validation is:

~~~bash
conda create -n tubularflow-cuda -c nvidia cuda-toolkit=12.6
conda run -n tubularflow-cuda make cuda CUDA_ARCHS=89
~~~

Run the executable from the activated environment, or expose its `bin` and
`targets/x86_64-linux/lib` directories through `PATH` and `LD_LIBRARY_PATH`.
Do not install a second Linux NVIDIA driver inside WSL.

~~~bash
conda activate tubularflow-cuda
export LD_LIBRARY_PATH="$CONDA_PREFIX/targets/x86_64-linux/lib:${LD_LIBRARY_PATH:-}"
./solvers/cuda/iga_cuda device-info
~~~

## Run

~~~bash
export IGA_CUDA_ROOT="$(pwd)/solvers/cuda"
export CASE_DIR=/path/to/case
export DATABASE=/path/to/work/case.ntiga

"$IGA_CUDA_ROOT/iga_cuda" mesh-check "$DATABASE"
"$IGA_CUDA_ROOT/iga_cuda" navier-stokes \
  "$DATABASE" "$CASE_DIR" --max-newton 12 --output velocity.txt
"$IGA_CUDA_ROOT/iga_cuda" solve \
  "$DATABASE" "$CASE_DIR" --system neuron_transport --output neuron.txt
~~~

Use the Navier--Stokes command with a prepared `vascular_flow/*` case and the
transport command with a prepared `neuron_transport/*` case. See the
[examples catalog](../../examples/README.md) for the source inputs and
preparation command. Each public configuration intentionally contains only its
matching application system except `vascular_flow/multispecies_pulse`, which
runs flow followed by six-species transport.

Configured CPU and CUDA transport lower the same weak-form terms and write a
neighboring `.fields` file recording output names. Navier-Stokes reads viscosity
and named velocity/pressure boundaries from the same `simulation_config.json`;
open pressure boundaries use the natural traction `-p n`.
Configured CUDA transport also contains the named snapshot-series interpolation
and per-step reassembly path, raw-state checkpoint/restart, and time-indexed
output. Straight-tube and 14,565-node branching parity/restart gates pass.
Both transient solvers accept `--stop-after-step N`, which writes an intentional
mid-run checkpoint while retaining the full configured horizon for restart
validation.
CUDA Navier--Stokes uses the same mesh-independent convergence controls as the
CPU backend: `--nonlinear-rtol R` (default `1e-5`), `--nonlinear-atol A`
(default residual RMS `1e-10` per equation), and `--mass-rtol R` (default
relative boundary-flow imbalance `1e-3`). A Newton state is accepted only when
both the nonlinear residual and mass-balance criteria pass. Iteration output
includes residual RMS, the unmodified continuity residual, net boundary flow,
and relative mass imbalance.
Transient CUDA output defaults to temporal Bézier VTKHDF; steady output
defaults to VTU. CPU and CUDA share the same output implementation. See the
[visualization output guide](../../docs/VISUALIZATION.md).
The fixed `transport` command remains only for old input compatibility. See the
[configurable PDE guide](../../docs/PDE_CONFIGURATION.md).


The scripts in `slurm/` reproduce the published Bridges-2 measurements and
require an external `IGA_CASE_ROOT`. See the
[Bridges-2 guide](../../docs/BRIDGES2.md) for module, allocation, interactive,
and `sbatch` commands.
`slurm/validate_transient_v100.sbatch` is the focused transient/restart gate.

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
- Navier-Stokes includes backward-Euler previous-state storage, temporal
  stabilization, physical-time stepping, per-step boundary waveforms, and
  checkpoint/restart. The multi-cycle Womersley gate passes on SM 89.
- Configured scalar flux/Robin surface assembly matches the CPU weak form and
  requires version 4-or-later `.ntiga` face labels.
- Configured transport supports per-step temporal Dirichlet waveforms.
- Configured transport raw-state checkpoint/restart validates ordered fields,
  system, velocity source, step, time, and `dt`; GPU restart validation passes.
- CUDA raw-state checkpoint/restart, time-indexed flow output, and R/RC/RCR
  outlet coupling pass GPU runtime validation. Outlet flow and natural pressure
  traction use the same host-side surface quadrature as CPU, and outlet values
  retain all continuity rows. Checkpoint metadata also matches CPU. Compliant
  walls are not implemented.
- Block-Jacobi is weaker than CPU PETSc local ILU for the large transport case,
  limiting solve speedup even though assembly is much faster.
- The CUDA backend consumes the in-tree CPU database header; changes to the format must be validated on both backends.
- Validation focuses on FP64 scientific GPUs and the current IGA formulations;
  other discretizations and mixed precision are not implemented.
