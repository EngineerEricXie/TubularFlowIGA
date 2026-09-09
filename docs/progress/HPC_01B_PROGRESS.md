# HPC-01B: ownership mappings, exact coverage and halo validation

Status: **complete for HPC-01B**. Common mappings, distributed coverage checks,
body-fitted integration and surface publication validation passed their gates.
Existing backend ownership contracts were audited. Actual distributed immersed,
moving and FSI computation remains HPC-03/07; no new distributed solver is claimed.

## Implemented interfaces

- `include/ParallelOwnership.hpp` distinguishes global node IDs from local
  owned offsets. It validates owned ranges, sorted/unique required IDs and
  field-row arithmetic before multiplication or allocation.
- `solvers/cpu/include/ParallelOwnershipValidation.hpp` checks communicator-wide
  node-range coverage, exact entity coverage, and exact `(element, node)`
  incidence coverage. It routes IDs/records to audit shards with `MPI_Alltoallv`
  and sorts only each shard. It does not gather the global ID catalog to every
  rank or accept a count/hash as proof of unique coverage. MPI count/displacement
  limits are checked before exchange. Controlled validation failures reach every
  member of the supplied communicator; this is not the general HPC-01C failure
  protocol and does not claim recovery from process loss or memory exhaustion.
- `SurfaceOwnershipValidation.hpp` adapts the existing surface layout contract.
  It checks interface/domain/subsystem identity, reference geometry, communicator
  partition metadata and exact publication ownership of arbitrary stable global
  node IDs. It preserves local layout/stamp validation and does not implement
  distributed force transfer, quadrature or membrane solves.
- `OwnedRowAssembler::ValidateOwnership()` is an optional collective audit.
  It checks the required-element catalog, database partition ownership,
  minimum-node integration ownership, and availability of every owned-row
  element contribution. File-read errors in the audit are converted to a
  collective validation failure.
- `OwnedRowAssembler::RequiredRows()` constructs checked PETSc rows for the
  flow and transport scatters. Both runtimes use it. In particular, the old
  flow expression `4 * int32_node + field` could overflow a 32-bit intermediate
  even with a 64-bit PETSc build; the new calculation checks and uses 64-bit
  arithmetic before converting to `PetscInt`.

No numerical operator, physical boundary, file format, partitioning algorithm
or convergence tolerance changes. Full ownership audits are explicit validation
operations, not extra communication added to every production time step.
Ordinary runtime construction performs the inexpensive row/required-ID checks.

## Independent element catalogs

The first real-case audit failed because it incorrectly treated database-owned
elements as a subset of the required-element index. The packer writes required
indices from node-row ownership; METIS element ownership is independent.
`TransientFlowRuntime` already calls `LoadOwned` separately for physical
integration. The corrected audit streams the authoritative partition records
from the database, then compares their incidences with the owned rows actually
available in `LoadRequired`. Existing databases are accepted without repacking.

The failed attempt is retained. The architecture document now states this
distinction explicitly. This was a correction to the new audit and its initial
inventory assumption, not a change to the existing packer contract.

## Tests and evidence

Toolchain: local WSL/i9-14900KF, GCC 11.4.0, Open MPI 4.1.2, PETSc 3.15.5,
real double and 32-bit indices. OpenMP and BLAS thread counts are one.
Builds use warnings enabled. The pure arithmetic tests exercise 64-bit row
calculations without allocating large fields; they do not certify an actual
64-bit PETSc installation.

`parallel_ownership_test` passed on three ranks within its 30-second timeout:

1. Global/local mapping inverses, empty owned ranges, invalid offsets,
   duplicate/out-of-order/out-of-range required IDs, and row overflow rejection.
2. Complete global catalogs, an empty catalog, and an empty rank.
3. Missing/duplicate IDs, inconsistent global counts, overlapping ranges, and
   a duplicate/omission pair that preserves both global count and ID sum.
4. Shared-node element incidence coverage, missing remote row contributions,
   duplicate contributions, and wrong node IDs with unchanged counts.
5. Independent groups of sizes one and two, including group-local rejection.
6. A packed algebraic fixture where rank 2 owns the element but no matrix
   rows or required elements. Ranks 0/1 assemble the exact two-by-two matrix
   `[[2,-1],[-1,2]]`. Another fixture removes rank 1's required element while
   both integration-owner coverage checks remain complete; incidence validation
   must reject it collectively. This checks an actual missing assembly
   contribution rather than only synthetic ID vectors.

The existing `vca_3d_runtime_test` passed with the checked row mapping, including
its iterative transport, multi-step, rollback and restart gates.

The extended `body_fitted_subcommunicator_test` additionally audits the actual
720-element database. It overwrites owned flow values with two distinguishable
states, verifies each requested remote/global ID, restores the original state
and verifies that the halo refreshes. Transport required-state buffers are
compared with independently gathered owned values after solve, rollback,
commit and checkpoint read. The same flow/pressure/transport/mass physical
comparisons remain `1e-6` relative L2 or `1e-12` absolute L2 for a zero reference.
The corrected-catalog run passed on all three ranks. The group of size two
reported relative L2 errors of `3.75708e-15` for velocity, `1.1644e-15` for
pressure, `1.70503e-14` for transport and `1.62206e-14` for integrated mass.
The size-one group matched exactly. The earlier audit failure remains retained.

The rebuilt flow and species graph subgroup checks also passed, with maximum
relative L2 errors `7.45082e-15` and `3.35727e-15` respectively. Surface checks
use noncontiguous IDs and reject duplicate publication, wrong partition rank,
changed reference geometry, another interface identity and unknown nodes on a
single rank, with rejection observed by every group member. Existing surface
and FSI runtime contract targets also passed. The final rebuilt VCA runtime
passed its unchanged numerical and lifecycle assertions.

```bash
make -C solvers/cpu parallel-ownership-test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
make -C solvers/cpu body_fitted_subcommunicator_test vca_3d_runtime_test petsc \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
```

The three-rank body test uses the HPC-01A command and the same prepared inputs;
see [that report](HPC_01A_PROGRESS.md). Successful and failed per-rank records,
stdout/stderr, binary/source inventory and build logs are retained under
`outputs/hpc01/ownership/`. These runs overlap compilation or other correctness
tests; their wall/RSS records measure test cost, not isolated scaling.

## Other-backend audit and scope

| Backend | Current ownership and check |
|---|---|
| 1D implicit | Matrix/residual insertion loops use PETSc's owned row range. `OneDGetVectorAll` replicates the result within the supplied group; that copy is not another algebraic owner. Existing four-formulation subgroup tests remain applicable because this code is unchanged. |
| 0D and small 1D graph state | State and graph execution are replicated. Domain accounting combines domains locally; output is published by graph rank zero. No extra sum over identical rank replicas is performed. |
| Immersed static/moving | Current solver objects are sequential on `COMM_SELF`, and graph loading rejects size greater than one. Active-node offsets, physical DOFs, controller multipliers and gauge DOFs are distinct; existing checked row counts and lookup methods enforce their ranges. MPI cell/face ownership and changing-geometry exchanges are required in HPC-03. |
| Surface and membrane | Immutable global material IDs key reference positions; owned IDs alone authorize field publication. Local layouts validate IDs, triangles, positive aligned weights and partition-bound stamps. The new collective adapter validates global publication coverage and identity. Membrane construction still requires one partition; distributed triangle integration, physical area/force/energy checks and exchanges are HPC-07. |
| CUDA | Single process/device; its local arrays do not represent distributed MPI ghosts. No multi-GPU behavior is implied. |

The surface adapter uses the reference catalog already replicated by the
existing layout, contributing it once to the sharded audit. It does not remove
that layout's global-geometry memory cost. It also preserves the current surface
contract's rejection of an empty owned slice; extending surface participation
for ranks without interface nodes belongs to the distributed FSI implementation.
The generic node/incidence validators and body-fitted assembler do support empty
ranks, as tested above.

General cross-rank error propagation outside these validation stages remains
HPC-01C. Oversized audit exchanges currently reject clearly instead of truncating
MPI counts; large-count alternatives and full memory scaling are later work.
The optional audit streams partition records separately and allocates local
incidences/audit shards. It is not an allocation-free production fast path.

The body-fitted ghost mappings are immutable for a runtime's fixed geometry.
`GatherRequiredVelocity` and `GatherRequiredState` synchronously refresh from
the current owned vector on every call. Moving geometry must recreate mappings
and state-transfer rules when geometry/ownership changes; these fixed-geometry
checks do not establish that future behavior.
