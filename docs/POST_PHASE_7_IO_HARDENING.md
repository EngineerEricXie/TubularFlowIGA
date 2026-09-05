# Post-Phase-7 publication I/O hardening

Phase 7 publishes the basic usable artifact set: one `fields.vtu` and one
`metrics.json` per epoch, plus the PVD collection.  Runtime and console
conservation diagnostics are the authoritative Phase 7 numerical evidence.

The following publication-infrastructure work is explicitly deferred to Phase
9 or a post-Phase-7 hardening effort.  It is not a Phase 7 numerical gate:

- a schema-v4 artifact and full conservation-record serialization;
- logical payload identities and byte-level payload hashes;
- a `commit.json` manifest;
- strict parsers and canonical JSON requirements;
- crash-recovery plus corruption and mixed-epoch validation; and
- extended publication fault injection.

The current schema-v3 retry and recovery behavior is intentionally limited to
the basic snapshot publication contract.  Consumers needing production-grade
artifact integrity or recovery guarantees must not infer them from Phase 7.
