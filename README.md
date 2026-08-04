# NeuronTransportIGA CUDA Solver

This directory is the GPU-native successor to `iga_solver_v2`. It keeps Bézier
extraction, quadrature, stabilized weak forms, boundary conditions, time
integration, sparse assembly, and Krylov iteration in project-owned C++/CUDA
code. CUDA Runtime and cuBLAS provide only device memory and vector primitives.

The implementation targets FP64 scientific computing. V100-32GB is the default because it is readily available, has strong FP64 throughput, and fits all current cases. The fat binary also supports A100 (SM 80), L40S (SM 89), and H100 (SM 90); the current allocation does not have the `GPU-dev` QoS required for PSC's A100 node.

Set the independent solver and legacy-case locations once per shell:

```bash
export IGA_CUDA_ROOT=/ocean/projects/mch260002p/thsieh1/iga_solver_cuda
export IGA_CASE_ROOT=/ocean/projects/mch260002p/thsieh1/NeuronTransportIGA
```

## Build

```bash
module load cuda/12.4.0
make -C "$IGA_CUDA_ROOT"
```

## Run on a GPU compute node

```bash
interact -A mch260002p -p GPU-shared --gres=gpu:v100-32:1 -t 00:30:00
module load cuda/12.4.0
"$IGA_CUDA_ROOT/iga_cuda" mesh-check DATABASE.ntiga
"$IGA_CUDA_ROOT/iga_cuda" transport DATABASE.ntiga CASE_DIR 10 output.txt velocity.txt
"$IGA_CUDA_ROOT/iga_cuda" navier-stokes DATABASE.ntiga CASE_DIR 8 velocity.txt
```

When an output path is supplied, each solver writes its text interchange file
and `OUTPUT.vtk`. Navier–Stokes VTK contains velocity/pressure; transport VTK
contains `N0`/`Nplus` and can be opened directly in ParaView.

The packed database may have been created for any CPU MPI rank count: the
single-GPU reader loads every element exactly once and ignores CPU ownership.
Do not run these commands on a login node.

Submit the full validated large case with:

```bash
cd "$IGA_CUDA_ROOT"
sbatch --export=ALL,IGA_CASE_ROOT="$IGA_CASE_ROOT" slurm/nmo_full_v100.sbatch
```

## Production scope

The single-GPU Navier–Stokes and transport workflow is implemented and numerically regressed on cylinder and the full 35,949-node `NMO_54499_new` case. Peak device use is 2.69 GiB, so one V100 is the preferred production configuration. Run Navier–Stokes first and pass its three-column velocity output to transport. Multi-GPU MPI is reserved for future meshes that do not fit or cannot meet runtime targets on one device. See `VALIDATION.md` for measured accuracy and performance.
