# PETSc solver options

`iga_multidomain_flow` and `iga_1d_3d_bifurcation` assign an independent PETSc
options prefix to every body-fitted 3D flow/transport runtime and every
implicit 1D runtime. Unprefixed PETSc options form a common baseline; a
domain-specific option takes precedence.

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps'

mpiexec -np 3 solvers/coupling/iga_multidomain_flow \
  --graph-case CASE --output-dir NEW_OUTPUT \
  -domain_island_a_flow_ksp_type fgmres \
  -domain_island_b_flow_ksp_type gmres
```

This example keeps the shared LU/MUMPS preconditioner while selecting FGMRES
for `island_a` and GMRES for `island_b`. It demonstrates option routing, not a
general performance recommendation.

## Prefixes

| Target | Prefix |
|---|---|
| Flow in domain `junction` | `domain_junction_flow_` |
| Transport in domain `junction` | `domain_junction_transport_` |
| Velocity fieldsplit in `island_a` | `domain_island_a_flow_fieldsplit_0_` |
| Pressure fieldsplit in `island_a` | `domain_island_a_flow_fieldsplit_1_` |
| Block-Jacobi sub-KSP in `junction` | `domain_junction_flow_sub_` |

For example, `-domain_junction_transport_ksp_rtol 1e-10` changes only that
transport solve. A child option such as
`-domain_island_a_flow_fieldsplit_0_pc_type gamg` takes effect only when the
parent PC actually creates that fieldsplit. PETSc documents the general
mechanism in
[KSPSetOptionsPrefix](https://petsc.org/release/manualpages/KSP/KSPSetOptionsPrefix/).

Domain IDs containing 1--64 lowercase ASCII letters, digits, or underscores
use the readable form above. Other IDs use
`domainh_<full-domain-id-sha256>_<role>_` to avoid collisions and excessive
prefix length.

Accepted-step output includes a `solver_configuration` JSON record with the
domain, role, prefix, KSP, PC, factor backend, tolerances, iteration limit, and
last KSP iteration count/reason. `last_iterations` belongs to the final KSP
call, not the sum across Newton iterations or rejected trials. A warm-started
solve may report zero iterations. An unselected factor backend is recorded as
`auto`.

## Option precedence and snapshots

Options are applied in this order:

1. runtime defaults;
2. legacy unprefixed PETSc options such as `-ksp_type`, `-pc_type`, and
   `-fieldsplit_0_pc_type` where that runtime supports them;
3. the full domain/role prefix.

Each runtime captures an immutable options snapshot. Values, flags, embedded
spaces, and option-like strings are copied directly rather than reconstructed
from the display form of `PetscOptionsGetAll`. Snapshot content and prefix are
checked for communicator-wide agreement. Applying a domain override does not
mutate the source database, so `-options_left` can still identify options
consumed by child solvers.

PETSc 3.15 child KSP objects may retain a prefix while querying the current
default options database. CoupledFlow therefore installs the private snapshot
during `SetFromOptions`, `SetUp`, and `Solve`, then restores the caller's
database and error handler before coordinating the result. Fieldsplits are
created after assigning the prefix. PETSc calls remain single-threaded within
each MPI process; independent communicators may use separate snapshots.

Embedding callers can pass `solver_options_prefix` to
`TransientFlowRuntime` or `TransientTransportRuntime`. An empty prefix retains
the legacy unprefixed interface. `SolverConfiguration()` is a local query and
must be wrapped by the caller's group error protocol when used collectively.

## Matrix backend options

Factor packages read some settings from the operator matrix rather than the
KSP or PC. An embedding runtime can call
`PetscSolverOptions::Attach(Mat)` before factor setup to bind the same private
database and prefix to the matrix. This path is wired into the shared immersed
MPI Newton core, the serial transient/moving immersed Jacobian, and the scaled
matrix in the distributed extension.

For a distributed moving runtime with prefix `immersed_moving_`:

```text
-immersed_moving_mat_mumps_icntl_14 100
```

This increases the MUMPS factor-workspace allowance relative to its estimate;
it does not change the equation, Jacobian, KSP tolerance, or Newton tolerance.
Use `-immersed_moving_ksp_view` to confirm the effective factor settings.

Matrix-option regressions verify independent snapshots and factor settings on
one, two, and four ranks, including split communicators. A factor failure
diagnostic includes the PC failure reason so structural/numerical zero pivots
can be distinguished from factor-memory exhaustion.

## Immersed, moving, and FSI runtimes

Immersed graph domains use `domain_<id>_flow_`. Their precedence is:

1. runtime defaults;
2. `immersed_static_` or `immersed_transient_` family options;
3. the complete domain prefix.

```bash
mpiexec -np 2 solvers/coupling/iga_multidomain_flow \
  --graph-case CASE --output-dir NEW_OUTPUT \
  -immersed_transient_ksp_type gmres \
  -domain_immersed_flow_ksp_type fgmres
```

Here, the domain named `immersed` uses FGMRES; other transient immersed
domains inherit GMRES. Static-family options use `immersed_static_`.

With no overrides, serial static and distributed immersed runtimes use
GMRES/LU, while serial transient and moving runtimes use FGMRES/LU. Multi-rank
LU retains the existing MUMPS selection. These defaults do not establish that
every PETSc preconditioner is suitable for the saddle-point systems.

Moving geometry reuses one immutable snapshot across committed and trial
epochs. The FSI wrapper derives a prefix from the fluid domain ID unless the
caller supplies one. A shared serial-transient options owner must use the same
communicator and prefix and outlive every attached KSP.

Solver options are part of the transient input identity. The current identities
are `ImmersedTransientInput/v6` and `ImmersedTransientMovingInput/v3`; changing
the effective options changes the trial identity without changing the field
format.

## Body-fitted standalone and VCA CLIs

`iga_navier_stokes` derives a `..._flow_` prefix from the selected
Navier--Stokes system. Legacy boundary configurations without a system name use
`domain_flow_flow_`. VCA transport uses `..._transport_`. `iga_solve` follows
the selected transport system name:

```bash
PETSC_OPTIONS='-ksp_type gmres -pc_type lu -pc_factor_mat_solver_type mumps' \
mpiexec -np 2 solvers/cpu/iga_solve DATABASE.ntiga CASE_DIR \
  --system transport --output result.txt \
  -domain_transport_transport_ksp_type fgmres
```

Place the database and case positional arguments before PETSc options. Both
CLIs accept a single-dash PETSc key followed by its optional value while
retaining the existing double-dash application options.

Each accepted step prints the effective prefix, KSP, PC, factor backend,
iterations, and reason. The configured transport CLI uses a nonzero initial
guess after the first step; its normal GMRES path supports that warm start,
whereas PREONLY does not.

## Implicit 1D and SNES

Graph domains use the same `domain_<id>_flow_` rule. The standalone `iga_1d`
CLI derives its prefix from the selected equation system. For
`blood_flow_1d`:

```bash
mpiexec -np 3 solvers/one_d/iga_1d IMPLICIT_CASE \
  --output-dir NEW_OUTPUT \
  -domain_blood_flow_1d_flow_ksp_type fgmres \
  -domain_blood_flow_1d_flow_snes_type newtontr
```

Only `nonlinear_aq` and `implicit_1d_pde` use SNES. `pressure_network` and
`linearized_aq` use a linear KSP. The linearized initial guess and the KSP
inside SNES share one runtime prefix, so `..._ksp_type` affects both while
`..._snes_type` affects only SNES.

Accepted-step `one_d_solver_configuration` records the effective linear and,
when applicable, nonlinear settings and convergence results. Rigid and
explicit 1D formulations do not create PETSc solvers and do not emit this
record.

Embedding code can create a `OneDPetscSolverContext` and pass it to
`AdvanceImplicitOneD` or the individual solve API. The context must use the
same communicator, remain alive for the complete call, and be constructed
outside local-only preparation stages because construction performs collective
agreement.

## Diagnostics

Useful inspection options include:

```bash
mpiexec -np 3 solvers/coupling/iga_multidomain_flow \
  --graph-case CASE --output-dir NEW_OUTPUT \
  -domain_junction_flow_ksp_view \
  -domain_junction_flow_sub_ksp_converged_reason \
  -domain_source_flow_snes_view
```

`ksp_view` reports the actual KSP, PC, and child solvers. Child convergence
options produce output only when the selected PC creates that child. A viewer
failure is an execution failure and must be detected through the process exit
status.

GAMG smoothers use `..._mg_levels_` prefixes, with a level number added by
PETSc. The coarse solver uses `..._mg_coarse_`. When changing the coarse KSP,
also configure compatible norm behavior; the multilevel regressions use
preconditioned norms for coarse GMRES.

## Legacy transport CLI

`iga_transport` uses `domain_neuron_transport_transport_`, with unprefixed
options as its baseline:

```bash
PETSC_OPTIONS='-domain_neuron_transport_transport_ksp_type fgmres -domain_neuron_transport_transport_ksp_rtol 1e-12' \
mpiexec -np 2 solvers/cpu/iga_transport \
  DATABASE.ntiga CASE_DIR 2 result.txt
```

The default remains GMRES with block Jacobi and a nonzero initial guess after
the first step. The example tolerance is a small-case regression setting, not
a universal tuning recommendation.

Changing solver options does not change `.ntiga`, field-output, or checkpoint
payload formats. It does change source/input compatibility identities, so a
restart must use the same build and effective numerical options.
