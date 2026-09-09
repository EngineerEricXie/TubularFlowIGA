# HPC-01A: communicator injection and independent graph execution

Status: **complete for HPC-01A**.
The 1D, body-fitted 3D, flow graph, and species graph split-group tests have
passed, as has the full multidomain CLI compatibility regression. This does not certify
the rest of HPC-01 or distributed immersed/FSI support.

## Implemented contract

The four synchronous entry points in `OneDImplicit.hpp` now accept an optional
`MPI_Comm communicator = PETSC_COMM_WORLD`: pressure network, linearized AQ,
nonlinear AQ, and `AdvanceImplicitOneD` dispatch (including implicit PDE).
Matrix/vector/KSP/SNES creation and the nonlinear MUMPS rank-count check use
the supplied communicator. The nonlinear initial-guess solve receives the
same communicator. `VecScatterCreateToAll` derives its communicator from the
vector, so the existing replication stays inside the supplied group.

These calls **borrow** the communicator for their synchronous duration. They
do not duplicate, retain, or free it. The caller must keep it valid until the
call returns; a runtime callback that captures it must keep its owner alive
for all calls. Standalone CLI call sites explicitly pass their application communicator,
which remains world. The optional default preserves existing external callers.
This changes neither the 1D equations nor their numerical/file formats.

## Current inventory

| Path | Ownership/lifetime and remaining work |
|---|---|
| `OwnedRowAssembler` | Stores a borrowed communicator; caller outlives assembler and PETSc objects; matrix/vector allocation already uses it |
| `TransientFlowRuntime` | Stores a borrowed communicator; domain collectives use it; `COMM_SELF` ghost vectors/index strides are intentionally local |
| `TransientTransportRuntime` | Stores a borrowed communicator; solver, field reductions and checkpoint viewers use it |
| 1D implicit entry points | Newly injected borrowed communicator; four formulations tested below |
| Standalone CPU CLIs | Still choose world as their default application group |
| Immersed runtimes | Sequential PETSc objects use `COMM_SELF`; loaders still reject unsupported MPI sizes |
| FSI | Existing single-partition surface contract remains; no distributed claim |
| Graph factories/executors | `RunMultidomainFlow` borrows its whole-graph communicator, passes it through 1D callbacks and 3D constructors, and destroys adapters before backends; per-domain scheduling remains HPC-08 |

The borrowed lifetime is documented at the 3D public constructors and the new
`include/MultidomainRunner.hpp` API. No duplicate is needed for these synchronous
calls: the caller owns the group, the runner owns runtime lifetimes within the
call, and the tests free groups only after all runtimes have been destroyed.
`DomainRuntimeRegistry`, pressure/species executors, and adapters have no hidden
world collectives. Local ghost vectors and their local stride index sets use
`COMM_SELF` intentionally; scatter contexts derive from the distributed vectors.
The graph runner now uses its supplied group for preflight, geometry and port
checks, construction, reductions, and output status broadcast. The default CLI
selects world and owns PETSc initialization/finalization. PETSc options must be
configured before the embedded call; the runner does not change the options
database. Immersed loaders still require their supplied graph size to be one;
sequential immersed and membrane objects intentionally retain `COMM_SELF`.

## Split-group test and evidence

`one_d_subcommunicator_test` splits four world ranks into groups of sizes 1
and 3. Groups use different viscosity/inlet data and different step counts.
Each rank independently computes a `COMM_SELF` reference, then runs its own
group's distributed solve. All four formulations compare area, cell and nodal
pressure, cell and segment flow against that reference. Different collective
sequences make unintended world participation detectable by the 90-second
timeout; test failures abort the test job, not masquerade as runtime recovery.

```bash
make -C solvers/one_d subcommunicator-test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
```

The target selects the installed direct MUMPS backend for this correctness
test. The test passed on Open MPI 4.1.2 / PETSc 3.15.5. Maximum reported error
was `1.16694e-15` for nodal pressure, below the unchanged `1e-6` relative L2
limit; zero references use `1e-12` absolute L2. The other nonzero maxima were
`9.12576e-16` for cell pressure and `8.81019e-17` for flow. Per-rank commands,
logs and process resources are retained in `outputs/hpc00/initial/subcomm-direct`.

An earlier run using default iterative PCs with `-ksp_rtol 1e-12` failed to
converge in linearized AQ (including the nonlinear initial guess). It exited
through the test's failure path without hanging. The direct reference above
does not establish robustness or scalability of those iterative defaults;
solver configuration/conditioning remains an explicit HPC-04 concern. The
comparison tolerance was not relaxed to accept the direct result.

The existing `one_d_petsc_test` also passed on two world ranks, including its
four-formulation and checkpoint checks. `iga_1d`, that existing test, and the
new test rebuilt with the selected PETSc toolchain.

## Body-fitted split-group acceptance

`body_fitted_subcommunicator_test` splits three ranks into groups of sizes 1
and 2. Each rank computes its own serial reference. The groups use different
inlet/scalar boundary values and one/two transport steps. The test compares
velocity and pressure separately, flow norms and reference inlet flow, the
transport field, integrated mass and sources. Every transport step solves,
rolls back, repeats and commits. Each group writes and restores its own first
step checkpoint; the runtime dies before its borrowed communicator is freed.
The flow runtime also enforces positive Jacobians, nonlinear convergence and
`1e-6` relative boundary mass imbalance.

Inputs are the HPC-00 straight-tube and straight-neurite cases. Their control
meshes both hash to
`b42b159449cfb4b222ff76d7bb08532c4a65aaa6568101c0483870db7b593db9`;
the test uses the flow database with matching transport geometry/labels.
The group of size two reported relative L2 errors:

| Quantity | Error |
|---|---:|
| Velocity | `3.48316e-15` |
| Pressure | `9.94483e-16` |
| Transport | `1.70503e-14` |
| Integrated mass | `1.62206e-14` |

The size-one group matched exactly. Every rank exited zero, with no timeout.
The unchanged comparison limits are `1e-6` relative L2 and `1e-12` absolute L2
for a zero reference. The test uses direct MUMPS as the numerical reference.

The first run exposed a production defect: `TransientTransportRuntime` enabled
nonzero initial guesses after the first step even when options selected
`KSPPREONLY`. PETSc rejects that combination. The runtime now queries the
configured KSP type and enables warm starts only for iterative solvers. No
operator, boundary condition, convergence tolerance or production default
changed. The existing `vca_3d_runtime_test` also passed with its default
iterative configuration, including its multi-step, rollback and restart gates.

## Whole-graph split-group acceptance

`multidomain_subcommunicator_test` links the production graph runner with
`IGA_MULTIDOMAIN_NO_MAIN`. It creates separate group directories and serial
references. Only each group leader runs the serial oracle while the other
members wait in that group, exposing accidental world participation even
before distributed solves begin.

- Flow mode runs schema-v5 explicit coupling in the size-one group and Aitken
  coupling in the size-two group. It compares accepted edge pressure/flow
  columns, edge identities and acceptance iteration counts. Both serial and
  group output must pass the original bifurcation physical/conservation gates.
  Maximum relative L2 was `7.45082e-15`; the size-one comparison was exact.
- Species mode runs schema-v6 flow/transport graphs, with reversed declaration
  ordering in the size-two group. It validates canonical species/edge IDs,
  donor routing, receiver concentration, domain and global integrated budgets,
  external-port signs, units and manifest contents. Nonzero integrated amounts
  and inventories are compared per species and quantity, avoiding mixed-unit
  norms. Maximum relative L2 was `3.35727e-15`; the size-one comparison was exact.

The constant-concentration, source-free species fixture has analytically zero
net outward amounts. Both solves are checked against that analytic zero with
the already declared `1e-12` absolute L2 bound, in addition to the unchanged
native budget gates. An earlier comparator divided cancellation roundoff
(`~1e-19`) by serial cancellation roundoff (`~1e-20`) and reported relative
errors above two; that failed comparison is retained. No computed field or
native tolerance was changed to resolve the invalid normalization. Another
initial graph test incorrectly required the legacy bifurcation manifest from
the multidomain entry; shared physical checks now validate the correct manifest
for each entry, preserving the legacy entry's original marker requirements.

## Reproducible commands and retained evidence

Toolchain: GCC 11.4.0, Open MPI 4.1.2, PETSc 3.15.5 real-double/32-bit indices
with MUMPS, on the local i9-14900KF WSL workstation. Tests use one OpenMP and
BLAS thread and MPI core binding. These concurrent correctness runs are not
isolated timing or scaling evidence. In particular, the graph fixture has one
Bezier element and is not evidence for large-problem load balance.

```bash
make -C solvers/cpu body_fitted_subcommunicator_test vca_3d_runtime_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
make -C solvers/coupling multidomain_subcommunicator_test petsc \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real

# Select a new output parent for every command/repetition.
mkdir /tmp/hpc-body-groups
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps' \
mpiexec --bind-to core -np 3 python3 scripts/hpc_rank_run.py \
  --output-dir /tmp/hpc-body-groups --expected-ranks 3 --timeout 360 -- \
  ./solvers/cpu/body_fitted_subcommunicator_test \
  /tmp/tubularflow-hpc-flow-baseline/straight_tube-1.ntiga \
  /tmp/tubularflow-hpc-flow-baseline/straight_tube-2.ntiga \
  /tmp/tubularflow-hpc-flow-baseline /tmp/tubularflow-hpc-transport-baseline \
  /tmp/hpc-body-groups

mkdir /tmp/hpc-graph-groups
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps' \
mpiexec --bind-to core -np 3 python3 scripts/hpc_rank_run.py \
  --output-dir /tmp/hpc-graph-groups --expected-ranks 3 --timeout 300 -- \
  ./solvers/coupling/multidomain_subcommunicator_test /tmp/hpc-graph-groups
# Repeat with another fresh parent and append "species" for schema-v6.
```

The body test inputs are generated by the HPC-00 preparation/packing workflow;
the graph test generates its own small fixture. Completed per-rank run records,
stdout/stderr, test inputs/outputs, and failed attempts are retained under
`outputs/hpc01/communicator/`. The inventory records the working-tree snapshot
and binary/input hashes; HEAD alone is not the tested source identity.
The body test took about 98.4 seconds per rank including waiting/startup and
used at most 266,858,496 bytes peak RSS per rank. This is test cost, not an
assembly/solve breakdown or a speedup claim. Existing phase measurements are
tracked separately under HPC-00B.

Compatibility checks completed so far:

- Existing bifurcation CLI smoke test: exit zero, including serial explicit,
  Aitken, MPI numerical comparison and failed-run completion-marker rejection.
- Existing explicit CLI smoke test: exit zero, including fixed/Aitken physical
  equivalence and MPI checks.
- Existing VCA runtime test: exit zero with default iterative solver behavior.
- Full multidomain CLI smoke test: the 600-second attempt timed out in
  `immersed_one`, after the serial/two-rank 0D clock cases. Exit 124, input
  fixture and logs are retained in `multidomain-cli-timeout*`. The parent and
  solver were verified absent before starting a new 3600-second attempt with
  line-buffered output. The extended run exited zero and passed all native
  checks. Maximum edge/domain/global species residuals were respectively
  `6.7762635780344027e-20`, `9.556334131140265e-15` and
  `8.890183375706226e-15`. Final output and tested binary hashes are preserved
  in `outputs/hpc01/communicator/multidomain-cli-passed.json` and its log.
  Command: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout 3600
  stdbuf -oL -eL ./multidomain_flow_smoke_test` from `solvers/coupling`.

Subsequent HPC-01B changes rebuilt the body-fitted and graph runtimes. Their
current-build subgroup, numerical, halo and checkpoint checks also passed;
see [HPC-01B](HPC_01B_PROGRESS.md). The full CLI record above identifies the
binary tested for the communicator change rather than attributing that run to
a later rebuild.

## Scope remaining outside HPC-01A

Cross-rank error propagation for arbitrary local failures remains HPC-01C;
test-harness `MPI_Abort` is not its implementation. Options prefixes and scalable
PC selection remain HPC-04. Distributed immersed/moving/FSI computation,
checkpoint publication contracts and per-domain communicator scheduling remain
HPC-03/05/07/08. These subgroup correctness tests do not certify cross-node
operation, which still requires HPC-09 allocation evidence.
