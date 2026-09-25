# Body-fitted runtime cleanup contract

`TransientFlowRuntime` and `TransientTransportRuntime` borrow their MPI
communicator. The caller must keep that communicator and PETSc alive until
cleanup finishes; the runtimes neither own nor free the communicator.

## Normal shutdown

After gathering final fields, writing diagnostics, and saving any checkpoint,
all ranks in the communicator call `Close()` in the same order. Closing is
terminal: only another `Close()` call or destruction is valid afterward. A
retry creates a new runtime.

`Close()` attempts every PETSc destruction operation in a fixed order and
stores the first error and object name without allocating memory. Once all
destruction calls return, `CollectiveLocalStage` coordinates the failing rank
and message. PETSc destruction calls do not run inside a local-only callback.
When several ranks fail, the lowest failing rank supplies the diagnostic and
retains its first failed object.

| Runtime | Destruction order |
|---|---|
| Flow | solver, scatter, destination IS, previous ghost, state ghost, source IS, RHS, update, previous state, committed state, state, Jacobian |
| Transport | solver, scatter, destination IS, state ghost, source IS, RHS, next, committed, current, forcing, previous matrix, left matrix |

Each runtime performs this sequence once. Repeated `Close()` calls coordinate
the stored result; a failed first close does not become successful and does
not retry uncertain destruction calls on only some ranks. The `noexcept`
destructor remains a fallback, but it cannot report success. Embedding callers
therefore call `Close()` explicitly before declaring a job complete.

## Caller and publication order

| Entry point | Final order |
|---|---|
| CPU flow and VCA | Finish fields, indices, and VTKHDF; close optional transport; close flow; publish completion/profile |
| Multidomain runner | Finish accepted-step bookkeeping; close transport and flow in domain-map order; write graph CSV, completion manifest, and stdout |
| Sequential explicit/fixed/Aitken runner | Finish accepted-history bookkeeping; close the 3D runtime; write CSV, completion manifest, and stdout |

Final graph writers use saved accepted histories, configuration, asset
identities, and container counts. They do not query closed PETSc state.
Registry-owned body-fitted adapters borrow these runtimes and do not solve or
abort from their destructors.

An accepted checkpoint or field file may precede a later cleanup failure. It
represents an accepted state, not successful completion of the whole job.
Atomic generations and restart behavior follow the
[coupled checkpoint contract](COUPLED_CHECKPOINT_CONTRACT.md).

## Failure scope

After a construction or solve failure, destruction attempts the same fixed
sequence while preserving the original exit result. Cleanup does not replace
collective ordering inside a runtime: all ranks must complete the coordinated
failure protocol before unwinding.

Regression injection replaces return codes after real destruction calls. It
tests coordinated exit, remaining cleanup, reference counts, and a retry with
a fresh runtime. It does not model a permanent block inside destruction, node
loss, process death, or corrupted PETSc internals. Other CLIs maintain their
own explicit cleanup regressions, including configured and legacy transport,
assembly smoke, and the four implicit 1D formulations.
