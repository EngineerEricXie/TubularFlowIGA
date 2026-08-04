# NeuronTransportIGA-CPU

This directory contains the project-owned C++ IGA Navier–Stokes and transport solvers. PETSc supplies distributed sparse matrices and Krylov methods; Bézier extraction, basis derivatives, quadrature, stabilized weak forms, boundary conditions, and time integration remain implemented here.

Set the repository and legacy-case locations once per shell:

```bash
export IGA_CPU_ROOT=/ocean/projects/mch260002p/thsieh1/NeuronTransportIGA-CPU
export IGA_CASE_ROOT=/ocean/projects/mch260002p/thsieh1/NeuronTransportIGA
```

## Build

The default PETSc location is the optimized local build in `../petsc/arch-linux-c-opt`:

```bash
module load openmpi/4.0.5-gcc10.2.0
make -C "$IGA_CPU_ROOT"
make -C "$IGA_CPU_ROOT" petsc
```

Override `PETSC_DIR` and `PETSC_ARCH` when using another installation.

## Prepare and check a case

Generate a partition matching the intended MPI size, convert the legacy text extraction files once, then check every element at the 4×4×4 quadrature points:

```bash
mpmetis "$IGA_CASE_ROOT/example/cylinder"/bzmeshinfo.txt 8
"$IGA_CPU_ROOT/iga_pack" "$IGA_CASE_ROOT/example/cylinder" 8 cylinder-8.ntiga
"$IGA_CPU_ROOT/iga_inspect" cylinder-8.ntiga
mpiexec -np 8 "$IGA_CPU_ROOT/iga_mesh_check" cylinder-8.ntiga
```

`iga_pack` rejects truncated or inconsistent `cmat.txt`, `bzpt.txt`, mesh, and partition data. The binary database stores sparse extraction rows, direct element offsets, and per-rank touching-element indices. Solvers also run the geometry check and collectively reject non-positive Jacobians rather than assembling a corrupt system.

## Run the coupled workflow

```bash
mpiexec -np 8 "$IGA_CPU_ROOT/iga_navier_stokes" \
  cylinder-8.ntiga "$IGA_CASE_ROOT/example/cylinder" 8 velocity.txt
mpiexec -np 8 "$IGA_CPU_ROOT/iga_transport" \
  cylinder-8.ntiga "$IGA_CASE_ROOT/example/cylinder" 300 concentration.txt velocity.txt
```

Navier–Stokes writes a three-column velocity file accepted directly by transport and a companion `velocity.txt.pressure`. Mesh label `1` is the velocity/concentration inlet, `0` is a no-slip wall, and every label `>=2` is a zero-pressure outlet. The nonlinear relative tolerance is `1e-5`; reaching `MAX_NEWTON` first returns a failure. Omit the final transport velocity argument to use `CASE_DIR/initial_velocityfield.txt`.

## Bridges-2 CPU execution

Compile on the login node, but run tests and simulations inside `interact` or `sbatch` with allocation `mch260002p`. This OpenMPI build is verified with `mpiexec`; plain `srun` may fail during `MPI_Init_thread`.

```bash
interact -A mch260002p -p RM-shared -t 00:30:00
module load anaconda3
module load openmpi/4.0.5-gcc10.2.0
mpiexec -np 8 "$IGA_CPU_ROOT/iga_mesh_check" cylinder-8.ntiga
```

Use `RM`/`RM-shared` CPU nodes for production and request memory based on a small-case benchmark. Do not launch MPI simulations on a login node.

See `VALIDATION.md` for cylinder numerical regression and the completed 57,456-node coupled benchmark.
