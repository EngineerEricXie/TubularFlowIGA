# HPC-00A: baseline provenance and fixture inventory

Status: **complete for inventory tooling and fixture selection**.
HPC-00B/C/D and the remainder of the workstation/HPC goal remain open.
Validation session: 2026-09-07 America/New_York (2026-09-08 UTC).
Source HEAD: `ee8da2ab528849814b6d12c8186476b5f7f59124`, with the current
uncommitted documentation, benchmark tools, and CPU Makefile fix recorded in
the inventory. This is not a claim that HEAD alone reproduces the worktree.

## Delivered artifacts

- [Baseline catalog](../../benchmarks/hpc_baselines.json): four selected source
  fixtures, preparation/run contracts, supported execution modes, generated
  assets, native gate sources, and historical reference evidence.
- [Inventory tool](../../scripts/hpc_inventory.py): tracked and new source
  hashes, fixture hashes, explicit generated artifacts, executable hashes and
  library linkage, selected PETSc configuration, compiler/MPI/HDF5/CUDA/GPU
  probes, requested resources, actual launcher affinity, and resource environment.
- [Per-rank wrapper](../../scripts/hpc_rank_run.py): independent log/resource
  records, launcher rank validation, timeout, nonzero exit propagation, and
  exclusive result directories.
- [Field comparator](../../scripts/hpc_compare_fields.py): streaming coefficient
  L2 comparisons with shape, finite-value, and optional global-node-ID checks.
- [Collection guide](../HPC_BENCHMARKS.md): commands, timer scope, reference
  requirements, and limitations.

Generated evidence is retained locally under `outputs/hpc00/initial/` and is
intentionally ignored by Git. `inventory.json` preserves restricted-sandbox
probe results; `inventory-runtime-access.json` preserves the subsequent probe
with MPI/GPU access. Raw logs and generated cases are evidence, not new
source fixtures or independently verified reference solutions.

The inventory's working-tree content digest at collection was
`f2f497acbf9dd87f57fdeb4dedf97da66614f6dbc3f0bd0c5c65db13e0b20b5b`.
Later report/checklist edits naturally change a newly collected source digest.

| Fixture | Input manifest SHA-256 |
|---|---|
| Body-fitted straight flow | `902a4c00d53d629f11ce166e9015c29445395232971c07e70d33f40b477e3214` |
| Body-fitted straight transport | `7fe0d365d1f39211ef8ac1e39cc4a862f6f6d983aa55cb564a3fe3617b11ddbf` |
| Immersed aneurysm source plus gate | `ca31f9c5c0c8809451c54ae06c5b24f353211d014a24ba00101d66ef6c087178` |
| Compliant-channel FSI fixture plus gate | `8c07a82e971797552c424e49bdad3ee8dfc75c3ea541897c013e690f9d156d9a` |

## Environment

- WSL/Linux x86-64; OS reports Intel Core i9-14900KF, 16 logical CPUs,
  8 exposed cores, and one NUMA node. These are guest-visible resources.
- GNU C++ with the repository's `-O3 -std=c++17 -Wall -Wextra -Wpedantic`
  CPU configuration; spline uses its existing C++11/OpenMP configuration.
- Open MPI 4.1.2 and system PETSc 3.15.5 at
  `/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real`.
  The configured PETSc uses real double scalars and 32-bit indices and includes
  MUMPS and Hypre. Executable linkage is captured separately from pkg-config.
- HDF5 probe reports 1.10.7. The body-fitted visualization binaries also link
  the repository-selected serial HDF5 library; a tool probe is not proof of
  parallel application I/O.
- NVIDIA GeForce RTX 4080 SUPER, compute capability 8.9, driver 591.86,
  reported total memory 16,376 MiB. The Conda CUDA compiler is V12.6.85.
- Initial CPU runs use `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`; MPI
  runs request `--bind-to core`. Actual per-rank affinities are in `run.json`.

MPI initialization and GPU metadata access fail in the restricted sandbox.
The two-rank `iga_1d ... --check` and GPU metadata probe passed with the
required runtime access. This distinguishes sandbox restrictions from absent
hardware or a missing MPI installation.

## Validation and discovered build defect

Commands executed:

```bash
python3 -m unittest discover -s scripts/tests -p 'test_hpc*.py' -v
make cpu-petsc PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
make -C solvers/cpu immersed_aneurysm_jacobian_test compliant_channel_fsi_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
env -u PETSC_DIR -u PETSC_ARCH make cpu
RANKS=2 OMP_NUM_THREADS=2 ./scripts/prepare_example.sh \
  vascular_flow/straight_tube /tmp/tubularflow-hpc-flow-baseline
RANKS=2 OMP_NUM_THREADS=2 ./scripts/prepare_example.sh \
  neuron_transport/straight_neurite /tmp/tubularflow-hpc-transport-baseline
```

All 13 Python tests passed. They exercise corrupted/mismatched inputs and
failure behavior as well as successful collection/comparison. CPU 3D, native
1D, and the two native immersed/FSI gate executables compiled successfully.
Existing spline warnings were observed; no spline numerical code was changed.

The first example-preparation attempt exposed a Makefile default-goal bug:
the earlier shared-header dependency rule selected `navier_stokes_test` when
invoking `make cpu`, which then failed to locate PETSc headers. Setting
`.DEFAULT_GOAL := all` restores the documented PETSc-free CPU tools build.
Both full example preparations subsequently passed. Each has 1,005 nodes,
720 elements, minimum sampled control-mesh determinant `0.000152336`, minimum
scaled Jacobian `0.763084`, no bad elements, and no surface intersections.

## Initial numerical observations, not scaling claims

The per-rank wrapper was exercised with real one- and two-rank flow/transport
runs. They completed successfully and produced the expected 1,005 field rows.
The standalone flow validator measured relative mass imbalance `3.37079e-7`
and relative divergence-theorem error `1.67587e-8` for the default two-rank run.

The initial default-solver velocity comparison passed (`9.63026e-7` relative
L2), but pressure failed the predeclared `1e-6` comparison threshold
(`1.20546e-6`). Both observations are retained. A stricter, identical solver
configuration for both ranks is being evaluated; the comparison threshold has
not been relaxed. Transport passed at `1.81155e-7` relative L2.

The depth-2 immersed Jacobian/conservation gate passed: maximum reported
centered-FD block defect `1.41286e-11`, final residual `2.50713e-17`, inlet
target error `2.71051e-20`, open normalized balance `2.08177e-15`, and wall
leakage `1.573e-5`. Its full FD-and-solve harness reported 95.309 s process
wall time and 61,157,376 bytes peak RSS. Aggregate assembly was 91.9865 s and
aggregate linear solves 0.288852 s; the separately reported line-search
assembly is included in assembly accounting, not an additional disjoint total.

These were initial correctness/instrumentation runs, with other validation
processes active. They are not isolated repetitions or a speedup study.
FSI execution and the tighter flow comparison are ongoing numerical work,
tracked in the next HPC-00 report; their completion is not claimed here.

Subsequent update: FSI and the tighter CPU comparison have now passed, and a
CUDA wall-trace discrepancy was fixed and tested. See
[HPC-00B/D progress](HPC_00BD_PROGRESS.md) for the later evidence and remaining
gates; the observations above retain the original collection context.

## Remaining work

- HPC-00B: add complete all-rank application stage instrumentation; PETSc
  nested events and rank-zero console timings alone do not fulfill this gate.
- HPC-00C: perform isolated repetitions and the supported resource matrix,
  including actual CUDA numerical runs and later cross-node allocation.
- HPC-00D: finish accepting the stronger reference configuration and capture
  all native fixture results. Preserve default-configuration failure evidence.
- Resolve single-dash PETSc CLI argument handling in the later solver-options
  task; profiling currently uses the working `PETSC_OPTIONS` interface.
- No numerical kernel parallelization, distributed immersed/FSI execution,
  coupled restart, or cross-node scalability has been implemented by HPC-00A.
