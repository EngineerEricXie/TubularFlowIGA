# Parallel ownership and halo contract

HPC-01B contract. Status: **common interfaces and current-backend acceptance complete**.
Row mappings, exact entity/incidence checks, body-fitted integration and surface
publication checks are implemented. Evidence and scope are recorded in the
[progress report](../progress/HPC_01B_PROGRESS.md). This does not certify new
distributed backends.

## Identity and scope

An entity is identified by its domain, entity kind, and stable global ID.
A rank number is relative to the domain communicator, not world. A local
array index, PETSc local row, or position in a ghost buffer is not a global ID.
Each mapping must define its valid range and reject missing entries; an absent
remote value must not silently become zero or refer to the first local entity.

Repartitioning may change owner and local index but must preserve physical
identity. Active-set changes additionally require an explicit state-transfer
rule. Existing `.ntiga` IDs and file layout remain unchanged. New immersed and
surface mappings must state how their stable IDs relate to geometry/history.

## Current body-fitted implementation

| Catalog | Meaning | Current implementation |
|---|---|---|
| Global control-node ID | Database node identity | Connectivity contains `int32_t` node IDs |
| Owned node range | Half-open range of node IDs assigned to this rank | `Database::NodeRange(rank)` |
| Algebraic row | `global_node * fields + field` | `OwnedRowAssembler`, checked against `PetscInt` |
| Local algebraic row | Offset within the owned contiguous PETSc range | Node-range offset followed by field offset |
| Required elements | Elements needed to assemble owned rows; not necessarily all database-owned elements | `Database::LoadRequired(rank)` |
| Database element owner | Unique partition owner from the packed database | `Element::owner`, `Database::LoadOwned(rank)` |
| Minimum-node element owner | Rank owning the smallest connectivity node | `OwnedRowAssembler::OwnsElementByMinimumNode` |
| Required state nodes | Sorted unique nodes referenced by local integration work | Runtime `RequiredNodes()` and ghost-position map |
| Ghost-buffer position | Local position in a gathered required-node vector | Runtime scatter/index sets; not a database node ID |

These are different ownership relations. `TransientFlowRuntime` keeps the
database-owned element catalog for physical boundary integration. Transport
mass and source integrals use the minimum-node predicate on required elements.
Either is a valid unique integration policy when its local catalog contains
every selected entity. Never sum integrals over every required element without
such a predicate: a touching element may appear on multiple ranks.
Do not require the partition-owned catalog to be a subset of required elements:
the packer's rank index is built from node-row ownership, independently of
METIS element ownership. The flow runtime deliberately loads both catalogs.

The collective assembler constructor checks agreement on global node, element
and field counts before PETSc object creation. This is a shape check, not a
content digest or a field-name/order check. The exact ownership/incidence audit
below remains independently useful. Configuration and controlled-failure
coverage are recorded in the [HPC-01C update](../progress/HPC_01C_CONFIGURATION_PROGRESS.md).

An owned-row assembly may evaluate the same element on multiple ranks. Each
rank inserts only its owned test rows, including columns referring to remote
nodes. Those repeated element evaluations are not duplicate algebraic
contributions. Conversely, filtering the assembly to database-owned elements
would omit contributions to rows owned by another rank unless an explicit
off-process assembly protocol replaces the current method.

## Checkable invariants

The following invariants apply to each supported backend, with collectives
restricted to its supplied communicator. Body-fitted checks cover row and
integration catalogs; the surface adapter covers publication ownership.
Distributed quadrature/transfer checks become mandatory when HPC-03/07 add
those currently unsupported modes:

1. Every global algebraic degree of freedom has exactly one owner. Local
   ranges cover the expected global range without overlap or holes, including
   legal empty ranks. Field count, row arithmetic and ID conversions fit their
   declared integer types before allocation or insertion.
2. Global-to-local and local-to-global mappings are inverses on their declared
   catalogs. Required-node IDs are unique and valid. An owned node may also
   appear in the required-state buffer; this does not give it a second owner.
3. For each owned row, every incident element contribution is available and
   inserted exactly once. Symbolic adjacency agrees with the numerical
   insertion pattern, including remote columns.
4. Each global integral selects exactly one contribution per physical entity.
   Check global ID coverage and multiplicity, not only a total count, which can
   hide a duplicate paired with an omission. A chosen integration owner must
   have the entity's geometry and current state.
5. Each requested remote state value has one defined owner and one matching
   exchange entry. Validate request/response identities and counts. Empty
   neighbor lists and empty ranks must participate in required collectives.
6. Cross-partition interface quadrature, ghost penalties and surface triangles
   follow an explicit contribution rule. Summing weights over unique entities
   must recover the global reference area/volume; shared endpoints alone do
   not define ownership of a face or quadrature contribution.

Full ID-coverage checks may be enabled for validation cases rather than every
time step. Their scalable implementation must avoid gathering the entire
catalog onto every rank merely to validate a distributed representation.

## Halo freshness and trial state

Each exchange belongs to a specific field, accepted step/trial, geometry
version and ownership mapping. Refresh the required state after a solve or
boundary update before evaluating dependent integrals. A rollback invalidates
values derived from the rejected trial; either restore the matching halo or
regenerate it from the restored owned state. A geometry/ownership change
invalidates its old scatter and interpolation maps.

Begin/end asynchronous exchanges must complete before buffers, PETSc objects
or communicators are destroyed. OpenMP workers must not share mutable element
scratch or modify PETSc objects concurrently under the existing contract.
Thread-local integration and serial insertion are separate concerns from
MPI entity ownership; their implementation is HPC-02.

## Other backends and outstanding implementation

- Current 1D implicit solvers distribute algebraic objects but gather the
  solution to their own communicator. This replication is intentional for the
  existing small network runtime; it is not distributed ownership of all
  network state. A future distributed catalog must cover segment and junction
  quantities separately.
- Small replicated 0D/1D domains must designate how an aggregate diagnostic is
  contributed once; summing identical replicas would multiply the result by
  communicator size.
- Immersed and moving runtimes are still serial. HPC-03 must define active
  cell/DOF ownership, pressure gauges, port multipliers, ghost-face contributions
  and stable-ID transfer when the active set changes. Removing a size check
  without these mappings is insufficient.
- FSI currently uses one surface partition. HPC-07 must validate owned/ghost
  surface nodes, triangle/area weights, force and moment, discrete work, and
  a single global acceptance decision. Duplicated ghost-node forces must not
  be included twice in global sums or Aitken inner products.
- CUDA remains one process and one device; its local arrays are not MPI halo
  buffers. Multi-GPU ownership is outside this checklist's required scope.

Implemented checks live in `ParallelOwnership.hpp` and
`ParallelOwnershipValidation.hpp`. `OwnedRowAssembler::ValidateOwnership()`
audits both integration catalogs and owned-row incidences, while required-row
construction validates the runtime scatter mapping. Shared elements, an empty
rank and deliberately missing/duplicate IDs have executable tests; real-case
halo/physical acceptance passed as recorded in the report.
`SurfaceOwnershipValidation.hpp` additionally checks partition identity,
reference-geometry agreement and unique surface-node publication. Add the immersed and surface
adapters when their distributed catalogs exist. Serial/MPI physical baselines
remain independent acceptance gates.

Sources: [OwnedRowAssembler](../../solvers/cpu/include/OwnedRowAssembler.hpp),
[Database](../../solvers/cpu/include/IgaDatabase.hpp),
[flow runtime](../../solvers/cpu/include/TransientFlowRuntime.hpp),
[transport runtime](../../solvers/cpu/include/TransientTransportRuntime.hpp),
[surface contract](../../include/DistributedSurfaceInterface.hpp).
