# Native tetra ALE flow with 0D source and terminal RCR

`NativeTetAleFlowDomainAdapter.hpp` exposes the project-owned tetrahedral
P2/P1 ALE Navier–Stokes runtime as a hydraulic-only `CoupledDomainRuntime`.
It uses the existing SI, outward-positive graph ports and converts relative
`u−w` flow targets to absolute native flow controls. FEM basis functions,
quadrature, residuals, Jacobians, and nonlinear iterations remain project
code; PETSc provides sparse algebra and MPI.

The integration test connects the existing compliant source reservoir to a
four-tetra native flow domain and then to the existing terminal RCR through
`PressureFlowComponentExecutor`. It runs a fixed first step and a rigidly
translated ALE second step. It checks the two edge flow balances, native
inlet/outlet volume balance, both 0D storage balances, graph commit counts,
and native mesh/time publication. A rank-0-only precommit error must abort all
three domains without altering their accepted states, after which the same
step succeeds. Explicit execution with density `1000 kg/m³` passes on 1, 2,
and 4 MPI ranks.

The native hydraulic adapter now has a bounded, checksummed accepted-step
checkpoint containing its current ALE coordinates, flow state, time, and step
count. Its model identity covers the reference mesh, ports, fluid and solver
parameters, wall labels, and an explicitly supplied motion-model ID. The test
restores all three graph domains into fresh runtimes after step one, advances
the restored graph through step two, and compares its accepted flow and 0D
pressures within `1e-10 × max(1, |reference|)` tolerance, and geometry exactly, with the
uninterrupted graph. It also rejects
corrupted payloads and a different motion-model ID. A dedicated
`NativeTetHydraulicGraphCheckpoint.hpp` bundle now saves the native FEM and
both 0D owners as one atomically published file epoch. All MPI ranks agree on
the accepted epoch and replicated state before rank 0 writes; discovery
checks every shard before any fresh owner is restored. The test deliberately
corrupts a shard, verifies collective rejection without publishing candidate
state, restores the file, and then resumes the second graph step. Callers must
discard the entire fresh candidate graph after any restore error. The
bundle compatibility identity binds the owner models, graph domains, ports,
edges, and pressure/flow iteration controls; the caller also supplies case
and execution identities. A separate two-process test now writes the first
accepted step in one MPI job, starts a new MPI job, restores only from the
published files, and compares the second step's full flow vector, ALE
coordinates, 0D pressures, and outlet flow with an uninterrupted reference
at 1, 2, and 4 ranks. Run it with
`make t7-native-ale-flow-0d-cross-process-test PETSC_DIR=/path/to/petsc
T7_SPECIES_MPI_RANKS=2`. The
existing general graph runner does not yet invoke this dedicated bundle, so
this is not a claim of CLI `--resume` support.

A separate native tetra hydraulic CLI now accepts a versioned case JSON and
labelled Gmsh 4.1 mesh with its own `--checkpoint-dir`/`--restart-dir` path;
see [the CLI contract](T7_NATIVE_TET_HYDRAULIC_CLI.md). It is not the older
IGA multidomain runner.

The same test also exercises fixed-point pressure/flow iteration at density
`1 kg/m³`, pressure relative tolerance `1e-4`, flow relative tolerance
`1e-6`, and relaxation `0.5`. The translated second step converges in 14
hydraulic iterations on 1, 2, and 4 ranks. A separate balanced-impedance
functional case uses density `1000 kg/m³`, source `(C,R,Q) =
(1e-9 m³/Pa, 500 Pa·s/m³, 0.1 m³/s)`, and terminal RCR
`(Rp,Rd,C,pd) = (10 Pa·s/m³, 50 Pa·s/m³, 1e-9 m³/Pa, 0 Pa)`.
Its explicit and fixed-point variants both pass on 1, 2, and 4 ranks; the
translated second fixed-point step converges in 19 iterations with
`0.1 m³/s` throughflow. A preliminary high-density iteration with the
*unbalanced* low-inertia 0D parameters failed in the native Newton solve;
that parameter combination is not claimed as passed. Neither successful
parameter set is physiological calibration. Run all four variants with
`make t7-native-ale-flow-0d-graph-test PETSC_DIR=/path/to/petsc
T7_SPECIES_MPI_RANKS=2` on an allocated compute resource.

The terminal RCR here is hydraulic only. It has pressure storage, distal
outflow, and volume accounting; it has no generic 0D concentration state or
species storage/exchange contract. The separate 1D–tetra–1D graph test is
the current native species path. This four-tetra 0D example does not establish
physiological parameter calibration, strict spatial/temporal accuracy,
broad strong-coupling robustness across impedances, or general-runner file
restart.
