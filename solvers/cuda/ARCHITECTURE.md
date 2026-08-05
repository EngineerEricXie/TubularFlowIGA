# TubularFlowIGA CUDA Backend Architecture

## Scope

`iga_cuda` is a project-owned FP64 implementation of the same two-field transport and four-field stabilized steady Navier–Stokes formulations as the CPU backend. The binary `.ntiga` database remains the interchange format. CPU partition ownership is ignored by the single-GPU reader, so an existing database can be reused without repacking.

## Data Layout

The host reader flattens element connectivity and sparse Bézier extraction. The GPU stores geometry points, sparse extraction rows, quadrature reference tables, and fixed geometry transforms. Geometry is evaluated once at all 4×4×4 quadrature points; a non-positive Jacobian aborts before assembly.

The global matrix is node-level block CSR: transport uses 2×2 blocks and Navier–Stokes uses 4×4 blocks. Only unique node adjacency is stored. Element kernels locate blocks by binary search and use FP64 atomics only for the final cross-element reduction.

Element-pair work is tiled in groups of 256. Each tile evaluates extracted basis values and derivatives in shared memory, accumulates all 64 quadrature contributions in registers, then updates global block CSR once. This avoids the multi-gigabyte repeated COO indices and values that a direct PETSc GPU port would require.

## Numerical Solvers

Transport assembles both time-integration matrices once and keeps all state vectors on device. Navier–Stokes rebuilds the VMS Jacobian and residual at each Newton iteration. Strong boundary rows match CPU v2.

The linear solver is restarted left-preconditioned GMRES. cuBLAS supplies only FP64 dot, norm, scale, copy, and AXPY operations. Sparse block matrix-vector products, block inverses, IGA weak forms, and assembly are project-owned CUDA kernels. Transport reuses one device GMRES workspace and its node-block Jacobi inverse across all time steps. Navier–Stokes reuses the workspace, rebuilds the inverse after each Jacobian assembly, and uses restart 200. A polynomial correction was tested and rejected for the saddle-point system.

Transport uses restart 50 and the same `1e-8` relative tolerance as CPU v2. A polynomial-Jacobi correction was also tested on the large transport system and rejected after it failed to converge; the stable default remains node-block Jacobi.

## GPU Strategy

The default fat binary contains SM 70, 80, 89, and 90 code. V100-32GB on `GPU-shared` is the selected production target for current cases: it has strong FP64 throughput, short queue times, and ample memory. SM80 is included for portability, but this allocation currently lacks the `GPU-dev` QoS required for PSC's A100 node. H100 is optional for larger future cases; L40S is supported but not preferred for FP64 workloads.

## Scaling Boundary

The execution model is one process and one GPU. Current cases fit comfortably, avoiding MPI halo traffic and global Krylov reductions. If future meshes exceed one device, row ownership can partition the existing block CSR while preserving element weak forms and local kernels; owned/interior work should overlap ghost exchange, with global GMRES reductions through `MPI_Allreduce`.
