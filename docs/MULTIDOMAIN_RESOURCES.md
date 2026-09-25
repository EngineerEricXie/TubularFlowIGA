# Multidomain resource groups

The multidomain runner keeps the existing shared-communicator mode as its
default. An optional `resources` object assigns disjoint contiguous MPI rank
groups to graph domains:

```json
"resources": {
  "mode": "domain_groups",
  "groups": [
    {"id": "small", "ranks": 1,
     "domains": ["source", "network", "terminal"]},
    {"id": "region_a", "ranks": 8, "domains": ["region_a"]},
    {"id": "region_b", "ranks": 8, "domains": ["region_b"]}
  ]
}
```

Groups consume ranks in manifest order. The example maps `small` to world rank
0, `region_a` to ranks 1–8, and `region_b` to ranks 9–16. A scheduler may place
those world ranks on nodes or NUMA domains through its normal rank-binding
options. The requested rank sum may be smaller than the allocation; remaining
ranks participate only in graph coordination. Every domain must appear exactly
once, group IDs must be unique, and groups cannot overlap.

Only members of a domain's group construct its native PETSc runtime. Every
world rank holds a lightweight proxy with the same graph metadata. The proxy
executes trial, rollback, and commit operations on the owning group and
broadcasts port states and species accounting from the group's first world
rank. A group-local exception is converted to the existing world-level
collective failure protocol before another graph stage begins.

Hydraulic solves use dependency levels from the pressure/flow graph. Domains
in the same level run together when they belong to different groups. Domains
sharing a group run in deterministic graph order. Staged transport rebuilds
the dependency levels after donor direction is resolved, so a flow reversal
changes the transport schedule without changing the coupling scheme.

The completion manifest records `resource_mode` and the reconstruction rule
`contiguous_manifest_group_order`. The graph manifest remains the authoritative
mapping: restart or resubmission must use the same manifest bytes and allocated
rank count. Grouped native checkpoint/restart is currently rejected before
runtime construction; use shared mode when native graph checkpoint/restart is
required. Grouped mode currently supports 0D, 1D, and body-fitted 3D domains.
Immersed 3D domains continue to use shared mode.

Use domain groups when independent, expensive domains can overlap or when
avoiding native state replication reduces memory. Small and sequential graphs
usually gain nothing from the extra broadcasts and synchronization. The
reproducible comparison tool is:

```bash
python3 scripts/hpc_domain_group_scaling.py \
  --output-dir artifacts/benchmarks/hpc08/domain-group-scaling --repeats 3
```

The script uses five ranks: one owner for all small domains and two ranks for
each of two independent 128-element 3D ducts. It compares the same five-rank
shared configuration, applies a `1e-6` numerical gate, and records per-rank
wall time, application phases, affinity, and peak RSS.
