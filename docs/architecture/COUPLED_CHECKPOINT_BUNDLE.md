# Coupled checkpoint bundle v1

HPC-05B implements the shared-storage format, publication, verification, and
recovery discovery in [CoupledCheckpointManifest.hpp](../../include/CoupledCheckpointManifest.hpp)
and [CoupledCheckpointBundle.hpp](../../include/CoupledCheckpointBundle.hpp).
The [state contract](COUPLED_CHECKPOINT_CONTRACT.md) defines the domain payloads
which HPC-05C/D must supply. The production graph does **not** yet expose these
checkpoint/restart APIs through its CLI.

## Format and identity

A caller supplies a unique 64-character lowercase SHA-256-shaped generation ID
for each publication attempt. The root is an existing, durably created POSIX
directory on shared storage. Each generation is a new directory:

```text
checkpoint-root/<generation>/
    <logical-shard-id>.tmp       # incomplete writer, never loaded
    <logical-shard-id>.shard     # immutable envelope and payload
    manifest.tmp                # incomplete publication, never loaded
    manifest                    # sole completion marker
```

Creating an existing generation fails. A failed attempt uses a new generation
on retry; incomplete files are retained for diagnosis. Published generations
are never overwritten, and this layer performs no retention/deletion. There
is no mutable `latest` pointer that can destroy access to the previous version.

The canonical ASCII manifest starts with `IGA_COUPLED_CHECKPOINT 1` and contains:

- generation and optional predecessor IDs;
- case, configuration, and execution compatibility SHA-256 values and rank count;
- accepted macro-step count and exact time/dt binary64 bit patterns;
- sorted, unique logical shard IDs, owner ranks, payload format IDs, payload
  byte counts, and SHA-256 values;
- a checksum of the complete canonical manifest body.

The compatibility digests are supplied by the integration layer and must cover
the state contract's graph/model/input/partition/field/backend/build/options
identities. A format identifier is an opaque payload codec version; this layer
does not reinterpret PETSc or CUDA data. Logical shard IDs use 1–128 portable
ASCII characters (`A-Z a-z 0-9 . _ -`, excluding `.`/`..`); they are not paths.
The provider retains original graph/domain IDs in its metadata and maps them
to these stable file IDs without changing graph interfaces.

Parsing is bounded and flat: manifest ≤4 MiB, at most 4096 shards, token ≤128
bytes, payload ≤1 TiB per shard, and rank count ≤`INT_MAX`. Integer conversion
rejects overflow; canonical reserialization rejects duplicate/unknown fields,
trailing content and alternate numeric spellings. Time/dt use exactly 16 hex
digits representing IEEE binary64; no decimal integer or clock precision is
lost. Nonfinite/nonpositive clocks, empty catalogs, repeated/unsorted IDs,
unknown schema versions and invalid hashes/owners are rejected.

Every shard has a bounded six-line `IGA_CP_SHARD 1` envelope containing its
generation, logical ID, owner, format and payload size. Its checksum covers
both envelope and payload. Copying an otherwise identical shard from another
generation is rejected. Zero-byte payloads are supported for empty ownership;
the envelope and checksum remain mandatory.

## Publication API

1. The graph coordinator reaches the all-domain accepted boundary and agrees
   on the epoch and complete shard catalog. It calls `CreateCoupledCheckpointEpoch`;
   all writers start only after this succeeds. Root creation/durability and
   domain capture/ownership agreement are integration responsibilities.
2. Each owner calls `WriteCoupledCheckpointShard` with a declared size and a
   streaming writer. The helper exclusively creates a temporary file, writes
   the envelope, checks exact payload length, computes a checksum over the
   producer's bytes, fsyncs/closes the file, publishes the immutable shard and
   fsyncs its directory. It returns a writer receipt with ID/owner/format/size/hash.
3. After all owners finish, the coordinator supplies the agreed full catalog
   and all sorted receipts to `PublishCoupledCheckpoint`. It reopens every
   expected shard, verifies envelope/size/hash against the **writer receipt**,
   and fsyncs it. Thus corruption between write and publication is rejected,
   rather than acquiring a new trusted checksum from damaged storage.
4. The coordinator writes/fsyncs/closes `manifest.tmp`, then publishes `manifest`
   last and fsyncs the generation directory. Publication uses same-directory
   `linkat` followed by removal of the temporary link, providing atomic
   no-replace visibility. A conflicting publisher fails without overwriting.

Open descriptors are scoped; writes handle partial progress and EINTR. Readers
require regular files and refuse final-component symlinks/FIFOs/devices.
File and directory sync/write/close/publication errors propagate to the caller.
A failure after linking a final file can leave that file visible; its data was
already synced, but successful durable publication is only reported after
directory fsync returns. Previous generations remain intact in every case.
As with other POSIX crash protocols, power-loss durability requires the
filesystem/storage to honor fsync and atomic link semantics; this batch tests
process termination and injected I/O errors, not physical power loss.

No helper calls MPI. HPC-05C must perform collective catalog/receipt/error
agreement and invoke these local operations in the proper schedule, outside
any lambda that already represents local failure agreement. Missing receipts
or expected files prevent completion. The bundle layer cannot infer missing
domains from a caller that supplies an incomplete catalog: constructing that
catalog from the validated graph is a required provider responsibility.

## Load and recovery

`LoadCoupledCheckpoint(root, generation, compatibility, catalog)` validates the
manifest, compatibility, exact complete catalog and every shard checksum/size.
Explicit selection rejects an incomplete or corrupt generation immediately.
`ReadCoupledCheckpointShard` then streams payload bytes through a consumer,
checking the same envelope/checksum again. The consumer must populate an
unpublished restore candidate because final checksum verification occurs after
streaming; only subsequent all-domain validation may publish accepted state.

`FindLatestCoupledCheckpoint` scans at most 65536 root entries, parses matching
generation manifests, and verifies candidates in descending accepted-step
order. It returns the latest fully verified compatible bundle plus rejection
diagnostics, which callers must report when falling back. Incomplete manifests,
corrupt data, wrong configurations and unreadable generations are excluded.
Two valid candidates at the same greatest step are ambiguous and require
explicit selection. An empty result does not authorize silently starting a
fresh simulation when restart was requested.

Shard processing uses a 64 KiB buffer, without gathering full fields into a
root vector. The coordinator currently verifies all shard bytes serially;
parallel verification/aggregators and retention are later I/O scaling work.
See the [05B report](../progress/HPC_05B_CHECKPOINT_BUNDLE_REPORT.md) for exact
process/fault/streaming evidence and remaining graph integration gates.
