# HPC baseline collection

The source catalog is [hpc_baselines.json](../benchmarks/hpc_baselines.json).
It separates four source fixtures, generated inputs, independent numerical
gates, and existing evidence. Its comparison policy is declared before changes
to the numerical kernels. New runs must still pass the fixture's native gates.

The Python tools require Python 3.10 or newer. The per-rank resource wrapper
uses Linux affinity and GNU `/usr/bin/time`; it does not add dependencies to
the native simulation pipeline.

## Record inputs and environment

For a fixed physical core comparison of body-fitted flow, use the optional
OpenMP CLI and prepack the same case for every requested rank count:

```bash
make -C solvers/cpu iga_navier_stokes_openmp iga_mesh_check PETSC_DIR=/path/to/petsc
python3 scripts/hpc_cpu_matrix.py \
  --case-dir /path/to/prepacked-case --database-stem straight_tube --system flow \
  --output-dir /path/to/new-hybrid-results --ranks 1 2 4 --total-cores 4 \
  --repetitions 3 --timeout 600
```

This compares 1 rank × 4 threads, 2 × 2 and 4 × 1 using Open MPI on one Linux
node. It verifies disjoint physical core bindings, actual native worker teams
and batch bounds. SMT siblings do not count as extra physical cores. All modes
use identical physical inputs, solver settings and output scope. The first
process is excluded from repeated statistics; it is not a controlled cold-cache
run. `observed_speedup_vs_pure_mpi` uses the same-core pure MPI configuration as
reference, and no rank-based parallel efficiency is reported for this mode.
Use otherwise idle allocated resources and retain negative performance results.
Omit `--total-cores` for the existing one-thread-per-rank CPU matrix.

The baseline catalog retains its original capability snapshot and policy.

```bash
python3 scripts/hpc_inventory.py \
  --output /tmp/hpc-inventory.json \
  --petsc-prefix "$PETSC_DIR${PETSC_ARCH:+/$PETSC_ARCH}" \
  --ranks 2 --threads 1 \
  --binary solvers/cpu/iga_navier_stokes \
  --artifact /path/to/case/database.ntiga \
  --artifact /path/to/case/controlmesh.vtk \
  -- mpiexec --bind-to core -np 2 ./solvers/cpu/iga_navier_stokes \
  /path/to/case/database.ntiga /path/to/case --output /path/to/run/velocity.txt
```

This records a command without executing it. Pass the actual configured PETSc
prefix containing `include/petscconf.h`; a source installation with a forwarding
Makefile may require its architecture subdirectory. Each binary has its own
hash and `ldd` result, so a system HDF5 probe is not mistaken for the binary's
effective library choice. Repeat `--artifact` for generated inputs and accepted
reference fields. The full working-tree content manifest includes untracked
source additions and missing tracked files, rather than claiming a clean HEAD.

Probe errors retain their error output. A GPU or MPI sandbox failure does not
prove hardware or software is missing. The inventory also checks the documented
`tubularflow-cuda` Conda compiler location when `nvcc` is absent from PATH;
`--nvcc` selects another installation explicitly. Only an allowlist of resource
environment variables is recorded.

## Measure one process per rank

Prepare inputs and compile executables before measuring. Use a new directory
for every repetition. The wrapper verifies launcher rank count, creates one
exclusive directory per rank, preserves stdout/stderr, and writes `run.json`
with wall time, affinity, exit status, and peak RSS. It propagates failures and
timeouts to the launcher.

```bash
mkdir /tmp/hpc-flow-run
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
PETSC_OPTIONS='-log_view :/tmp/hpc-flow-run/petsc.log' \
mpiexec --bind-to core -np 2 \
  python3 scripts/hpc_rank_run.py \
  --output-dir /tmp/hpc-flow-run --expected-ranks 2 --timeout 600 -- \
  ./solvers/cpu/iga_navier_stokes DATABASE CASE_DIR \
  --output /tmp/hpc-flow-run/velocity.txt
```

`PETSC_OPTIONS` supplies profiling independently of application CLI parsing.
The current flow CLI treats single-dash PETSc arguments as legacy positionals;
direct `-log_view` arguments are not supported by that parser yet. Use the
environment interface until the CLI options work is validated.

For a standalone serial fixture, omit `mpiexec` and use `--expected-ranks 1`.
The output parent must already exist. Do not launch MPI inside the measured
command: the wrapper belongs inside the launcher, once per rank.

Accounting rules:

- Wrapper wall time covers executable startup through final output and exit.
  It is not assembly time, and it excludes the outer MPI launch interval.
- PETSc log events can be nested. Report `PCSetUp`, `KSPSolve`, assembly,
  communication, and application timers with their nesting/scope; do not add
  overlapping intervals into an alleged exclusive total.
- Process peak RSS is recorded per rank. Their sum is a sum of individual
  peaks, not a measured simultaneous aggregate peak.
- Rank-zero console timers are not substitutes for all-rank stage timing.
- CUDA allocation and synchronized kernel timings require the CUDA backend's
  own metrics; this wrapper measures host resource use only. The CUDA solver
  now reports `cuda_allocations scope=project_device_buffers` with requested
  peak/live bytes tracked across `DeviceBuffer` allocation, move, and release.
  This excludes allocator rounding and library/driver-owned allocations.
  Existing `gpu_used_gib` is a device-wide usage sample and must be labeled
  separately; it is not the project's allocation high-water mark.
- A zero exit status is named `process_passed`, not `numerically_validated`.
  Require the independent gates and output comparisons before acceptance.

The CPU matrix runner below aggregates repeated body-fitted MPI runs. Record a
fresh-process first run separately from at least
three repetitions, and do not call it a cold filesystem-cache run unless cache
state was controlled. Keep benchmark runs isolated from other simulations and
compilation when making speedup claims.

## Repeat a workstation CPU matrix

Prepare one case directory containing `STEM-1.ntiga`, `STEM-2.ntiga`, etc. Each
database must be packed for its corresponding rank count from the same geometry
and configuration. Compile before running the matrix. The runner uses Open MPI
launcher flags, binds one rank per core, disables oversubscription, and sets
OpenMP and BLAS thread counts to one. Select ranks that fit the allocated cores.

```bash
python3 scripts/hpc_cpu_matrix.py \
  --case-dir /path/to/prepared-transport \
  --database-stem straight_neurite --system neuron_transport \
  --ranks 1 2 4 8 --repetitions 3 --timeout 600 \
  --output-dir /path/to/new-transport-matrix
```

Use `--system flow --database-stem straight_tube` for the flow fixture. The
flow matrix uses the catalog's MUMPS reference and unchanged nonlinear/mass
tolerances; transport uses its existing default solver. Inherited `PETSC_OPTIONS`
are replaced explicitly and recorded. These are benchmark configurations, not
changes to solver defaults. This runner currently supports the configured
`neuron_transport` fixture and body-fitted flow, not arbitrary equation systems.

Each configuration first passes a separate `iga_mesh_check`. The runner then
executes the first solver process for each rank count, followed by at least
three rounds of fresh processes. Configuration order rotates between repeated
rounds; jobs never overlap within this runner. Geometry preflight can warm input
caches, so even the first solver process is **not a cold-cache measurement**.

`matrix.json` records configuration, topology, tolerances, and hashes of inputs,
binaries, and collection scripts. Each observation retains launcher logs,
rank-local logs/resources, profile summaries, and field comparisons against
the fresh one-rank first run. Velocity and pressure are compared separately;
transport node IDs must match. Inputs and binaries must remain unchanged across
measurement. Output directories are exclusive, preventing accidental reuse.

`summary.json` preserves every attempted observation. Failure stops the matrix
and suppresses repeated statistics; missing, duplicate, or failed samples cannot
be silently dropped. Successful summaries report all repeated values, medians,
min/max and population standard deviation, phase maxima, maximum rank RSS, and
the sum of rank RSS peaks. First processes remain in the observations but are
excluded from repeated statistics. Launcher time includes MPI/Python startup
and shutdown; per-rank process and application time are also available.

Speedup and efficiency are descriptive ratios against one rank, not automatic
performance acceptance. Audit numerical/physical gates and measurement isolation
before making claims. Rank-phase maxima can occur on different ranks and must
not be added into a supposed critical path. Summed RSS peaks are not simultaneous
aggregate memory. Unsupported OpenMP/hybrid assembly is N/A; the supported GPU
backend is recorded as a separate pending measurement, not N/A.

## Repeat serial fixtures and a single GPU

`hpc_serial_matrix.py` runs the other catalog paths with a separate first process
and at least three fresh-process repetitions. It pins the host process to the
selected logical CPU and requests one OpenMP/BLAS thread. Compile and prepare
inputs before running; select a CPU within your allocation.

```bash
python3 scripts/hpc_serial_matrix.py \
  --case immersed --case-dir /path/to/isolated-depth2-fixture \
  --cpu 0 --repetitions 3 --timeout 300 --output-dir /path/to/new-immersed-matrix

python3 scripts/hpc_serial_matrix.py \
  --case fsi --cpu 0 --repetitions 3 --timeout 1500 \
  --output-dir /path/to/new-fsi-matrix
```

The immersed fixture runs the full nine-block centered-FD regression and its
Newton, geometry, and conservation gates. `--solve-only` is not this benchmark.
The FSI fixture enforces its compiled force/moment, moving-mass, wall leakage,
continuity, and strong-coupling assertions. The collector requires successful
process/profile status and complete finite diagnostics; printed values are
rounded and do not replace the native full-precision assertions. These two
fixtures remain serial; unsupported MPI/OpenMP modes are N/A.

For CUDA, use an intact accepted CPU matrix for the same body-fitted catalog
case. The collector verifies the recorded input/binary/tool hashes and the
one-rank first field, then compares every GPU field with that reference using
the unchanged `1e-5` relative L2 and `1e-12` zero-reference absolute limits.

```bash
conda run --no-capture-output -n tubularflow-cuda env \
  LD_LIBRARY_PATH=/home/tsungyeh/anaconda3/envs/tubularflow-cuda/targets/x86_64-linux/lib \
  python3 scripts/hpc_serial_matrix.py \
  --case cuda-transport --cpu-matrix /path/to/accepted-cpu-transport-matrix \
  --cpu 0 --gpu 0 --repetitions 3 --timeout 180 \
  --output-dir /path/to/new-cuda-transport-matrix
```

Use `--case cuda-flow` with the accepted flow matrix for GPU flow. It uses the
previously registered 30-Newton maximum and unchanged nonlinear/mass tolerances.
The Conda path above is this workstation's documented installation; adjust the
environment and CUDA library path on another machine. `--gpu` sets
`CUDA_VISIBLE_DEVICES`; a separate `device-info` preflight records the visible
device before timing. A failed device probe remains an error with preserved logs.

Every accepted GPU run must have exactly one allocation report, a positive
requested project-buffer peak, and zero project buffers live at shutdown. The
summary reports repeated allocation peaks separately from host peak RSS. Driver,
library, and whole-device usage are outside that project allocation scope.
Native solver convergence and field comparisons do not replace independent
physical validation.

The serial runner shares the CPU runner's exclusive output directories and
failure-preserving statistics policy. It preserves all observations, excludes
the first process from repeated medians/ranges/standard deviation, and verifies
inputs, reference fields, binaries, and collection scripts again at completion.
There is no speedup acceptance inferred from a successful matrix. Compare timing
scopes and affinity explicitly: CPU MPI launch time and a directly launched
serial process have different startup overheads.

## Rank-local application phases

The CPU `iga_solve` and `iga_navier_stokes`, CUDA simulation commands, and the
selected immersed aneurysm/FSI gate executables accept the opt-in environment
variable `IGA_PROFILE=1`. Use it outside `mpiexec` so every rank inherits the
setting. With the per-rank wrapper above, each successful `stdout.log` ends
with one `hpc_profile` JSON record. Profiling is disabled by default.

```bash
python3 scripts/hpc_profile_summary.py /path/to/completed-run > summary.json
make -C solvers/cpu phase-profile-test
```

The summary verifies complete rank coverage, successful exit, the recorded log
hash, and the timing partition. It retains per-rank values and reports their
minimum, mean, maximum, and maximum/mean imbalance, alongside host peak RSS.
A missing, failed, modified, or unprofiled rank is an error, not a zero sample.

`inclusive_s` includes nested phases; `exclusive_s` subtracts them. Only the
sum of exclusive durations plus `unscoped_s` partitions `elapsed_s`.
Application elapsed time starts after PETSc initialization and ends before
PETSc finalization; wrapper wall time covers the complete child process.
The profiler uses the main thread's monotonic wall clock without inserting
MPI barriers. Worker timers remain disabled so overlapping OpenMP wall times
cannot be summed as process time.

Current configured-transport coverage:

| Phase | Measured scope |
|---|---|
| `input` | Initial configuration/field reads; snapshot loading/interpolation |
| `geometry` | Required-element loading, geometry and surface-boundary preflight |
| `assembly` | Operator integration/insertion and boundary-row work |
| `solver_setup` | Explicit `KSPSetUp` and `KSPSetUpOnBlocks` calls |
| `linear_solve` | `KSPSolve` after explicit setup |
| `communication` | Operator finalization/row collectives and output gather |
| `output` | Field/checkpoint writing, with gather subtracted from exclusive time |

Communication scopes include local work and wait time inside the wrapped
collective operations; they are not a measurement of network transfer alone.
MPI inside PETSc solves remains in the solver phase. PVD/physiology manifest
writing, preallocation, time-step RHS work, startup diagnostics, and other
unwrapped work remain visible as `unscoped_s`. Zero calls denotes an unexecuted
or uninstrumented phase.

Body-fitted flow uses the same accounting for initialization, assembly, ghost
exchange, matrix finalization, convergence diagnostics, solver setup/solve,
and field/checkpoint output. Static immersed assembly includes line-search
assemblies; geometry covers catalog/layout creation. Moving geometry includes
construction of the candidate epoch, history extension, and state mapping.
The FSI `coupling` interval encloses an entire accepted/rejected macro-step;
its exclusive time subtracts nested geometry, fluid/traction/membrane assembly,
and linear solves. Membrane linear time includes its small dense factorization.
Serial immersed paths report no inter-rank communication; console diagnostics
and unwrapped lifecycle work can still appear as unscoped time.

The [PETSc profiling interface](https://petsc.org/release/manualpages/KSP/KSPSetUpOnBlocks/)
explicitly supports moving block-PC setup ahead of `KSPSolve`. Setup includes
the parent and subblock work. Legacy console `linear_s`/aggregate solve timers
retain their original wider intervals; use profile **exclusive** times for
disjoint comparisons. Arbitrary user PC implementations can still perform
work inside application/solve callbacks, so the reported phase describes the
actual wrapped API boundary, not a universal division of every backend's internals.

CUDA assembly finishes at its existing device synchronization. Profiled GMRES
synchronizes before entry and after block-inverse construction, charging setup
as a nested child of linear solve. Those two additional synchronizations are
enabled only with `IGA_PROFILE=1`. CPU/GPU transfer work stays in the enclosing
phase; `communication` does not purport to count all PCIe traffic. The wrapper
and existing allocation counter retain host RSS and requested GPU peak bytes.

The benchmark output records CPU/CUDA wall traces, steady/transient
comparisons, and the one/two-rank instrumentation checks.

## Compare fields

CPU rank comparisons use the catalog's predeclared `1e-6` relative L2 limit;
CPU/CUDA comparisons use `1e-5`. A zero reference uses an explicit absolute
L2 limit of `1e-12`. These supplement all native physical/convergence gates.

```bash
python3 scripts/hpc_compare_fields.py reference-velocity.txt velocity.txt \
  --rtol 1e-6 --zero-atol 1e-12
python3 scripts/hpc_compare_fields.py reference-transport.txt transport.txt \
  --node-ids --rtol 1e-6 --zero-atol 1e-12
```

The comparison streams numeric rows, checks shape and finite values, and checks
strictly increasing matching IDs when `--node-ids` is used. Compare velocity,
pressure, and transport separately with the same ordering, units, and pressure
gauge. These are coefficient norms, not quadrature-weighted physical L2 norms.
Do not silently shift pressure or ignore ID mismatch to make a comparison pass.

## Independent transport budgets

`iga_transport_validate` independently integrates the selected static,
prescribed-velocity backward-Euler transport fields. It checks per-species
volume and surface balances, free-row residuals, and essential values at both
time levels. The surface balance includes the concentration-weighted velocity
divergence required by the native nonconservative advection term.

```bash
make -C solvers/cpu transport-budget-test
./solvers/cpu/iga_transport_validate DATABASE CASE_DIR neuron_transport \
  field.step000000.txt field.step000001.txt --rtol 1e-6 --zero-atol 1e-12
```

The CLI returns 0 for numerical acceptance, 2 for a failed numerical gate, and
1 for invalid input or execution failure. Use
`scripts/hpc_transport_budget_regression.py` to validate complete CPU/GPU
histories against an intact accepted CPU matrix and exercise corrupted inputs.
The selected two-species regression additionally requires zero net linear
transfer. This is an offline small-case checker with full fields in memory;
it does not implement distributed diagnostics for large runs.

The registered normalization, commands, analytic checks, and unsupported
modes are encoded by the transport budget script and catalog.

## Complete immersed and FSI reference states

The native fixtures support optional `--reference-output NEW_DIRECTORY` after
their numerical gates. The [reference collector](../scripts/hpc_reference_states.py)
records native acceptance, complete physical fields, stable IDs, units, and
source/binary/input provenance. It preserves the predeclared fieldwise `1e-6`
relative L2 and `1e-12` zero-reference absolute L2 gates.

```bash
python3 scripts/hpc_reference_states.py collect --case fsi \
  --output-dir NEW_REFERENCE --cpu 0 --timeout 1500
python3 scripts/hpc_reference_states.py collect --case fsi \
  --output-dir NEW_CANDIDATE --cpu 0 --timeout 1500
python3 scripts/hpc_reference_states.py compare NEW_REFERENCE NEW_CANDIDATE
```

For immersed flow, use `--case immersed --case-dir DEPTH2_CASE`. Collection
always runs the complete FD fixture. Comparison requires two accepted
collections with matching fixture definitions and physical inputs; different
binary versions are recorded and may be compared. A failed numerical comparison
returns 2, invalid/incomplete evidence returns 1, and acceptance returns 0.

The reference-state manifests document field coverage and source identity.
These numerical state files are not restart checkpoints. Collection with extra output, particularly concurrent
collection, is separate from an isolated performance matrix.

The immersed transient runtime now has optional OpenMP volume assembly. Build
the FSI fixture with `make -C solvers/cpu compliant_channel_fsi_openmp_test
PETSC_DIR=...`, then collect it with explicit resources:

```bash
python3 scripts/hpc_reference_states.py collect --case fsi \
  --binary solvers/cpu/compliant_channel_fsi_openmp_test \
  --threads 4 --cpus 0 2 4 6 --output-dir NEW_OMP_CANDIDATE
python3 scripts/hpc_reference_states.py compare NEW_REFERENCE NEW_OMP_CANDIDATE
```

Choose available CPU IDs from the actual machine topology. The CPU list
overrides `--cpu`; it must contain unique available IDs and at least one CPU
per thread. The collector sets OpenMP binding, limits BLAS to one thread,
records execution settings separately from physical inputs, and requires
native evidence of the requested assembly team and bounded batches.
`--binary` must implement the same fixture, native gates and complete state
export. It does not enable parallelism in a serial executable. The static
immersed aneurysm fixture has not received this volume-assembly integration.

For an isolated repeated comparison, the OpenMP matrix controller serializes
fresh processes and compares all nine fields against the accepted reference
after every run:

```bash
python3 scripts/hpc_openmp_matrix.py \
  --reference outputs/hpc00/reference-states/fsi-reference \
  --output-dir NEW_ISOLATED_MATRIX --threads 1 4 --cpus 0 2 4 6 \
  --repetitions 3 --timeout 1800
```

Run this after other simulations finish, on otherwise idle allocated CPUs.
The controller cannot exclude unrelated external work. It alternates thread
configurations between repetitions, excludes repetition 0, and preserves each
collection and field comparison. Its `matrix.json` is updated atomically after
each accepted run. Numerical acceptance is separate from speedup eligibility:
the latter requires at least 10% lower median end-to-end launcher time and a
maximum repeated RSS ratio no greater than 1.25. One-rank RSS is also aggregate
RSS here. Negative performance results remain valid observations. A failed or
incomplete run prevents an accepted matrix; fewer iterations or build success
cannot replace complete native and field gates.

## Tests

```bash
python3 -m unittest discover -s scripts/tests -p 'test_hpc*.py' -v
```

The tests cover provenance changes, invalid fixture paths, missing probes,
configured PETSc identification, failed children, timeout, rank mismatch,
non-overwriting result directories, literal argument handling, and field
comparison rejection paths. Actual MPI and numerical runs are separate evidence.
