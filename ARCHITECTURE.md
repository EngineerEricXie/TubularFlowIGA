# NeuronTransportIGA-CPU Architecture

## Data path

MATLAB mesh generation writes `controlmesh.vtk`; `spline_src` produces `bzmeshinfo.txt`, `cmat.txt`, and `bzpt.txt`. `iga_pack` validates these files and creates an indexed `.ntiga` database. Each MPI rank seeks only the element records needed for its owned rows instead of scanning every multi-gigabyte ASCII file.

## Distributed assembly

PETSc uses contiguous node/field ownership. The database stores a one-layer touching-element index for each rank, so a rank computes every contribution to its owned rows and never relies on PETSc's off-process stash. A symbolic adjacency pass derives exact diagonal and off-diagonal nonzero counts. `MAT_NEW_NONZERO_ALLOCATION_ERR` turns a sparsity mistake into an immediate error; validated cylinder assembly has zero allocation growth.

This removes the legacy fixed allocation of 1,000 entries per row, repeated full-file parsing on every rank, and dense extraction storage. Element matrices remain temporary contiguous buffers and are released after insertion.

## Numerical methods

The steady incompressible solver uses the original stabilized VMS residual and approximate Jacobian with physical first/second basis derivatives. Its nonlinear loop rebuilds the residual/Jacobian and enforces velocity and outlet-pressure constraints. Block-Jacobi/ILU is the measured default: despite more Krylov iterations, it is faster and uses less memory than Schur field-split setup on these meshes. PETSc options can override it.

Transport retains the two-field reaction–diffusion–advection/SUPG formulation. Its time-independent left and previous-state operators are assembled once. Each time step uses one sparse matrix-vector product and reuses the KSP/preconditioner.

Both solvers perform a collective 4×4×4 Jacobian-sign preflight. Input errors, invalid geometry, allocation errors, and negative `KSPConvergedReason` values terminate with a nonzero status.
