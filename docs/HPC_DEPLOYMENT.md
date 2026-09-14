# Reproducible HPC deployment and test tiers

This page defines the repository-owned evidence workflow. It does not replace
site policy or turn a workstation run into cross-node acceptance.

## Build compatibility record

Run dependency checks in the same module environment used to build and run the
solver. The optional report records the source state, compiler and MPI wrapper,
PETSc scalar/index and package capabilities, HDF5 configuration, CUDA compiler
and target architectures, Slurm resources, and optional binary linkage:

```bash
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-linux-c-opt
./scripts/check_dependencies.sh cpu --report /path/to/evidence/cpu-build.json

CUDA_ARCHS="70 80 89 90" \
./scripts/check_dependencies.sh cuda --report /path/to/evidence/cuda-build.json
```

On the documented workstation, the CUDA check also finds
`tubularflow-cuda/bin/nvcc` when that Conda environment is not active. Pass the
same path as `NVCC` to a direct Make invocation, or build through `conda run`.
`make hpc-build-manifest HPC_BUILD_MANIFEST=FILE` is the CPU shorthand. Reports
are evidence artifacts and belong outside Git.

Do not mix the MPI implementation behind `mpicxx`, PETSc, parallel HDF5, and
`mpiexec`. A successful dependency report records compatibility inputs; the
numerical tiers still have to pass.

## Test tiers

`hpc_test_tiers.py` gives every command an outer timeout, preserves stdout and
stderr, propagates nonzero exits, and writes `result.json`. A tier that cannot
run writes `status: skipped` and a concrete `skip_reason`; a skip is never
reported as a pass.

| Tier | Scope | Typical command |
|---|---|---|
| `unit` | mesh core, coupling contracts, CUDA host execution contract | `make hpc-test-unit HPC_TEST_OUTPUT=/path/unit` |
| `mpi` | distributed immersed numerical test at 1, 2, and 4 ranks | `make hpc-test-mpi PETSC_DIR="$PETSC_DIR" HPC_TEST_OUTPUT=/path/mpi` |
| `gpu` | CUDA host contract and real `device-info` probe | `make hpc-test-gpu HPC_TEST_OUTPUT=/path/gpu` |
| `scheduled` | two-node MPI execution inside an existing Slurm allocation | `python3 scripts/hpc_test_tiers.py --tier scheduled --output-dir /path/scheduled` |

The GPU tier is a hardware smoke test. Use the CPU/GPU field matrices in
[HPC_BENCHMARKS.md](HPC_BENCHMARKS.md) for numerical backend comparison. Large
checkpoint, FSI, graph, and scaling cases remain scheduler tiers and should not
be added to the fast unit target.

Every output directory must be new. This prevents a rerun from combining stale
logs with current results.

## Cross-node graph job and safe requeue

`solvers/cpu/slurm/multinode_graph.sbatch` requests two RM nodes, maps ranks by
node, binds each rank to the `cpus-per-task` core set, makes OpenMP match that
allocation, forces BLAS libraries to one thread, and records the build and job
metadata. It copies the graph input to node-local `$SLURM_TMPDIR` (falling back
to PSC `$LOCAL`); results and checkpoints remain on a shared filesystem.

Submit from the repository root after building in a compatible compute allocation.
The wrappers use `SLURM_SUBMIT_DIR` because Slurm executes a spool copy; set
`IGA_REPO_ROOT` explicitly when submitting from another directory. They load
Anaconda Python and the same Open MPI module for every attempt.

Build inside the compute allocation:

```bash
make hpc-cross-node-binaries PETSC_DIR="$PETSC_DIR" PETSC_ARCH="$PETSC_ARCH"
```

```bash
export IGA_GRAPH_CASE=/shared/cases/native-graph
export IGA_OUTPUT_ROOT=/shared/results/graph-${USER}
export IGA_CHECKPOINT_ROOT=/shared/checkpoints/graph-${USER}
export IGA_GRAPH_PREPARE=1
export IGA_CHECKPOINT_TEST_DELAY=60
sbatch -A "$PROJECT_ACCOUNT" \
  --export=ALL,IGA_GRAPH_CASE,IGA_OUTPUT_ROOT,IGA_CHECKPOINT_ROOT,IGA_GRAPH_PREPARE,IGA_CHECKPOINT_TEST_DELAY,PETSC_DIR,PETSC_ARCH \
  solvers/cpu/slurm/multinode_graph.sbatch
```

Slurm sends `SIGUSR1` five minutes before the time limit. The batch shell trap
forwards it to `mpiexec`, whose Open MPI 4.0.5
`--mca ess_base_forward_signals SIGUSR1` setting delivers it to solver ranks. The native graph exits after the next accepted macro-step and a
published checkpoint. The script records `checkpointed` and requeues by
default. Immediately before requeue, the wrapper explicitly enables the job
requeue flag because launcher startup on Bridges-2 may reset the submission flag. Set `IGA_REQUEUE_ON_SIGNAL=0` to stop after the safe checkpoint. Each
attempt uses `attempt-$SLURM_RESTART_COUNT`; a resumed attempt adds
`--restart-dir` when a published generation exists.

With `IGA_GRAPH_PREPARE=1`, attempt zero creates a 16,384-element, eight-step
native graph packed for the allocated rank count. `IGA_CHECKPOINT_TEST_DELAY=60`
sends a controlled `SIGUSR1` to the batch shell after one minute, exercising the
same trap used by Slurm's wall-time warning. The first attempt must publish an
accepted-step checkpoint and requeue; the next attempt discovers that generation
and completes from it. Use new shared case, output, and checkpoint paths. Set both
variables to zero when running a separately prepared production graph and relying
only on the scheduler warning.

This job is for shared-mode native 0D/1D/body-fitted-3D graphs. Grouped graph
checkpoint/restart and native moving/FSI graph checkpoint are current
limitations. Validate signal delivery and requeue on the target site because
Slurm and Open MPI policies vary.

The committed Bridges-2 values follow the current
[PSC Bridges-2 user guide](https://www.psc.edu/resources/bridges-2/user-guide/):
RM jobs allocate complete 128-core nodes, while RM-shared cannot span nodes.
The small graph example therefore reserves more cores than its default eight
MPI ranks and should be used for cross-node correctness or replaced with a
larger packed case. Slurm's `B:` signal form targets only the batch shell, which
is why the explicit trap and Open MPI forwarding are both present; see the
[Slurm sbatch signal contract](https://slurm.schedmd.com/sbatch.html) and
[Open MPI launcher options](https://docs.open-mpi.org/en/main/man-openmpi/man1/mpirun.1.html).

## Cross-node moving FSI and pair restart

`solvers/cpu/slurm/cross_node_fsi.sbatch` closes the scheduled acceptance path
for the current bounded moving-FSI library runtime. In one two-node allocation
it runs a one-rank reference, a four-rank strong-coupling writer distributed as
two ranks per node, and a new two-rank read-only process restored from the
four-rank pair bundle. The existing checkers enforce fields, unique owned rows
and surface IDs, force/traction, Aitken history, ports, conservation, iteration
counts, 4-to-2 repartition provenance, and immutable checkpoint bytes.
It first runs the `scheduled` test tier with one MPI rank on each node, so the
same job also provides the required small cross-node test result and skip-free
`result.json`.

```bash
export IGA_FSI_OUTPUT=/shared/results/fsi-cross-node-${USER}
sbatch -A "$PROJECT_ACCOUNT" \
  --export=ALL,IGA_FSI_OUTPUT,PETSC_DIR,PETSC_ARCH \
  solvers/cpu/slurm/cross_node_fsi.sbatch
```

The fixture is intentionally long and is a correctness/restart workload. Its
timings do not replace the body-fitted strong/weak scaling study. The membrane
matrix remains on one bounded owner even though fluid rows, moving geometry,
surface publications, convergence checks, and checkpoint shards span ranks.

## Strong and weak scaling

Prepare a fixed physical case packed for every strong-scaling rank count. Also
prepare one weak case per rank count, increasing global elements so elements per
rank remain within the configured tolerance. `hpc_prepare_scaling_cases.py`
generates C2 duct cases from the repository fixture: the default strong case has
16,384 elements, while every weak case has exactly 256 elements per rank. Their required-element indexes follow the production packer's contiguous
owned-node ranges; they do not replicate every element on every rank. Tiny
fixtures are rejected as evidence by review even if the collector itself can
run them.

```bash
make hpc-prepare-scaling \
  HPC_SCALING_CASES=/shared/cases/scaling \
  HPC_SCALING_RANKS='1 64 128 256'
```

`hpc_cross_node_scaling.py` performs geometry preflight, at least three fresh
solver repetitions, native physical validation, and strong-case velocity and
pressure comparison to one rank. It records per-rank workload, speedup,
efficiency, solver iterations, phase times, maximum rank RSS, summed individual
RSS peaks, communication, and allocation metadata. Weak cases are validated
independently because their physical fields intentionally differ.

The Bridges-2 wrapper stages all cases on each node and runs the collector in an
exclusive allocation:

```bash
export IGA_SCALING_CASE_ROOT=/shared/cases/scaling
export IGA_SCALING_OUTPUT=/shared/results/scaling-${USER}
export IGA_SCALING_RANKS='1 64 128 256'
export IGA_SCALING_PREPARE=1
sbatch -A "$PROJECT_ACCOUNT" \
  --export=ALL,IGA_SCALING_CASE_ROOT,IGA_SCALING_OUTPUT,IGA_SCALING_RANKS,IGA_SCALING_PREPARE,PETSC_DIR,PETSC_ARCH \
  solvers/cpu/slurm/cross_node_scaling.sbatch
```

After all three jobs pass, consolidate their evidence. The command rejects
single-node results, dirty or mismatched revisions, a missing graph requeue,
skipped scheduled tests, FSI runs that do not span two hostnames, incomplete
repetitions, tiny strong/weak cases, missing numerical checks, and absent wall,
iteration, RSS, or communication metrics:

```bash
python3 scripts/hpc_finalize_cross_node.py \
  --graph-output /shared/results/graph-${USER} \
  --fsi-output /shared/results/fsi-cross-node-${USER} \
  --scaling-output /shared/results/scaling-${USER} \
  --output /shared/results/cross-node-acceptance-${USER}.json
```

Only a `status: passed` result from this finalizer closes HPC-05C, HPC-07C/D,
and HPC-09A-D. Preserve the referenced output trees with the acceptance file.

Auto-preparation requires the prebuilt
`solvers/coupling/hpc_duct_solver_fixture`. Set `IGA_SCALING_PREPARE=0` to use
an existing immutable case tree. The default relative patterns are
`strong-{ranks}/root3d`, `weak-{ranks}/root3d`, and `../duct.ntiga`; override
them with the four `IGA_STRONG_*`／`IGA_WEAK_*` variables for production data.
The collector uses distributed FGMRES with rank-local block Jacobi/ILU for this
study, rather than a global direct solve that would obscure MPI scaling.

Choose rank counts and `--ntasks-per-node` for the actual node type. Retain
negative scaling results. Strong efficiency is `T(1)/(p*T(p))`; weak efficiency
is `T(1)/T(p)` at constant elements per rank. Report maximum-rank critical-path phases separately;
phase maxima may occur on different ranks and cannot be summed into wall time.
The sum of individual peak RSS values is also not a simultaneous memory sample.

## Acceptance evidence

For each job retain the source revision and dirty state, module versions,
dependency/build manifest, Slurm job ID and node list, exact launcher and
binding output, staged input hashes, exit status, native convergence and
physical gates, field comparisons, phase/RSS records, checkpoint generation,
and `sacct` or `sstat` accounting. Generated databases, results, JSON evidence,
and scheduler logs remain outside Git; commit only the scripts and a concise
report that points to their location.

## Recorded Bridges-2 acceptance

The [completed acceptance report](progress/HPC_09_BRIDGES2_CROSS_NODE_ACCEPTANCE.md)
records the tested revision, actual jobs, numerical checks, and formal scaling.
The graph used 256 ranks (128/node) and a four-hour allocation per attempt.
Accepted macro-steps exceeded the default five-minute warning margin: request
checkpoints early with `scancel --signal=USR1 --batch JOBID` and leave
enough time for the current macro-step to converge and publish. The recorded
run used early requests on subsequent attempts without reducing numerical or
physical workload requirements.
