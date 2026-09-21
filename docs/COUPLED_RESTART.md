# Native graph checkpoint and restart

`iga_multidomain_flow` and `iga_1d_3d_bifurcation` can checkpoint completed
macro-steps and resume them in a new MPI job. The native CLI supports 0D, 1D,
and body-fitted 3D flow graphs, plus 1D and body-fitted 3D species graphs.
Restart currently requires the same MPI rank count and communicator membership.

Moving-flow and bounded moving-FSI library runtimes use versioned bundles and
can redistribute fields by stable global IDs. They are not yet wired into the
native graph CLI. Immersed graph, CUDA graph, and body-fitted repartitioned
restart are also outside the CLI path described here.

## Save and resume

Use the same inputs, build, thread configuration, PETSc options, rank count,
and `.ntiga` partitions for both jobs. The Makefile embeds a SHA-256 identity
of the relevant source, including local modifications. Custom builds must set
`IGA_NATIVE_CHECKPOINT_SOURCE_SHA256`; generate it with:

```bash
python3 scripts/native_checkpoint_source_identity.py
```

A binary without that identity can run normally but cannot save or load a
native graph checkpoint.

```bash
make -C solvers/coupling iga_multidomain_flow PETSC_DIR="$PETSC_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps'

# First job: save every three accepted steps and stop after step six.
mpiexec -np 3 solvers/coupling/iga_multidomain_flow \
  --graph-case CASE --output-dir FIRST_OUTPUT \
  --checkpoint-dir CHECKPOINT_ROOT --checkpoint-every 3 --stop-after-step 6

# New job: load the latest complete generation and finish the configured run.
mpiexec -np 3 solvers/coupling/iga_multidomain_flow \
  --graph-case CASE --output-dir RESUMED_OUTPUT \
  --restart-dir CHECKPOINT_ROOT \
  --checkpoint-dir CHECKPOINT_ROOT --checkpoint-every 3
```

`FIRST_OUTPUT` and `RESUMED_OUTPUT` must not exist. The case still defines the
total number of time steps; `--stop-after-step` only shortens the current job
and cannot precede the restored generation. A resumed output contains the
saved accepted-history prefix followed by new records, so it does not depend
on the first output directory. Supplying only `--restart-dir` leaves the source
bundle read-only.

| Option | Behavior |
|---|---|
| `--checkpoint-dir ROOT` | Save generations under a shared POSIX directory, creating missing directories as needed |
| `--checkpoint-every N` | Save every positive `N` accepted steps; default `1` |
| `--restart-dir ROOT` | Load the latest complete compatible generation |
| `--stop-after-step N` | Stop at accepted step `N`; that final step is checkpointed even when it is off-cycle |

Normal completion also writes a final checkpoint. Each generation contains the
complete graph state, the next pressure guess, species-donor hysteresis, and
the full accepted-history prefix. Ranks write their owned 3D fields; the group
root writes small replicated state and history. Step-local Aitken state is
reinitialized by the normal algorithm.

## Scheduler warning signals

When `--checkpoint-dir` is enabled, the process accepts `SIGUSR1` from before
input loading. The handler only sets a flag. After the next successful accepted
macro-step, the full communicator writes a checkpoint and exits normally.
A signal received by any participating rank triggers the same coordinated
action.

Send the warning to solver ranks and allow enough time for one macro-step plus
checkpoint I/O. If the next step fails or the job is forcibly terminated,
restart uses the previous complete generation. Signal forwarding is scheduler
specific; see [HPC deployment](HPC_DEPLOYMENT.md).

## Compatibility and failure handling

Compatibility checks cover:

- graph and domain configuration plus all referenced asset contents;
- source revision, MPI/PETSc version, scalar/index widths, and byte order;
- communicator membership, rank count, compiler, OpenMP, and fast-math identity;
- thread environment, effective PETSc options, Newton controls, and
  conservation tolerances.

Paths are informational; content SHA-256 values establish identity. The first
format version does not migrate state across builds or numerical
configurations.

Writers persist shards and checksums, synchronize them, and publish the
manifest last. Incomplete or corrupt generations remain diagnostic artifacts
but are skipped during loading. Re-saving a step creates a new epoch rather
than overwriting an existing one. A failed compatibility or candidate-state
check produces no new graph output.

The format and durability rules are specified in
[checkpoint bundle v1](architecture/COUPLED_CHECKPOINT_BUNDLE.md). Bundle v1
targets shared POSIX file systems. Accepted history is streamed by record with
a 16 MiB metadata limit per record; the complete prefix and in-memory history
grow with simulation length.
