# Generic 0D species reservoir contract

`include/ZeroDSpeciesReservoir.hpp` is a project-owned, dependency-free
well-mixed 0D transport model. It is separate from the hydraulic source and
RCR equations: those supply outward-positive port flows, while this model
tracks one positive compartment volume (`m³`) and nonnegative amount (`mol`)
per named species. It neither infers a pressure/volume law nor converts a
hydraulic RCR capacitor pressure into an anatomical volume.

For each step, the caller supplies one finite flow (`m³/s`) for every declared
port. Negative flow enters the compartment and requires a nonnegative upstream
concentration (`mol/m³`) for every species. Positive flow leaves at the new
well-mixed compartment concentration. Zero/outgoing flow must not carry
unneeded donor data. An optional per-species source rate has units `mol/s`.
The backward-Euler update is

```text
V_new = V_old - dt Σ_ports q_out
M_new = (M_old + dt (source_rate - Σ_inflow q_out c_donor))
        / (1 + dt Σ_outflow q_out / V_new)
```

It rejects missing/extra ports, incomplete or negative inflow concentration,
invalid state/source data, nonfinite results, nonpositive final volume, and
negative trial amount. The result publishes a standard `PortState` for each
port plus explicit volume and per-species amount accounting. The required
conservation equations are `V_new - V_old + dt Σq_out = 0` and
`M_new - M_old - dt source_rate + dt Σ(q_out c_donor) = 0`, with `c_donor`
equal to the upstream concentration on inflow and the new mixed concentration
on outflow. No external framework performs these updates.

Run the local unit contract with `make coupling-zero-d-species-reservoir-test`.
It covers two species, two ports, forward and reversed flow, pure volume
expansion/dilution, zero-flow source, balances, and invalid inputs.

`include/ZeroDSourceReservoirSpecies.hpp` composes the native hydraulic
source-reservoir with that same transport model. The graph port carries the
hydraulic outward flow; the external pump is a separate transport port with
outward flow `-Q_pump`. Positive prescribed pump flow therefore needs an
explicit donor concentration, while graph inflow needs a graph donor. The
composition checks that transport `ΔV` equals hydraulic `C Δp`.
`include/ZeroDSourceReservoirSpeciesDomainRuntime.hpp` exposes this model
through the existing staged graph interface, including separate hydraulic
and transport trials, rollback, abort, and paired prepare/commit. The pump
species amount remains an external boundary term in global accounting.
`make coupling-zero-d-source-reservoir-species-test` covers both graph
directions, reversed pump, volume balance, and missing donors;
`make coupling-zero-d-source-reservoir-species-runtime-test` covers lifecycle
and accounting. `make t7-native-ale-species-graph-zero-d-source-test
PETSC_DIR=/path/to/petsc T7_SPECIES_MPI_RANKS=2` runs the actual
0D source→native tetra ALE FEM→RCR chain in explicit and fixed-point modes. Two steps
(fixed, then translating) passed on 1, 2, and 4 MPI ranks with both edge and
global species balances, pump boundary amount, and a rank-0-only precommit
failure followed by successful retry. The 24-tetra geometry is functional,
not physiological. The fixed-point source uses finite capacitance `1e-4
m³/Pa`, resistance `1e5 Pa·s/m³`, prescribed pump flow `0.1 m³/s`, and
initial stored pressure `1e4 Pa`, putting it near a steady-flow working
point. It converged in 12 outer iterations on the second step with the
original native FEM Newton settings. The earlier low-resistance source
(`0.1 Pa·s/m³`) did not converge in fixed-point mode: pressure feedback
produced a difficult FEM boundary flow. Simply relaxing Newton tolerance
or raising its iteration limit did not resolve that case. The graph
intentionally rejects a direct 0D-to-0D edge.

`include/ZeroDTerminalRcrSpecies.hpp` adds a narrow hydraulic/transport
composition for a terminal RCR. It evaluates the existing native hydraulic
RCR equation, uses its graph-port and distal-branch outward flows as the two
species-reservoir ports, and requires the species volume change to equal the
RCR capacitor's stored-volume change. The absolute mixed volume remains an
explicit transport initial condition; capacitance determines its *change*,
not a physiologically calibrated absolute volume. A reversed distal branch
requires an explicit distal donor concentration. Run
`make coupling-zero-d-terminal-rcr-species-test` for forward/reversed flow,
storage agreement, species balance, and missing-donor rejection.

`include/ZeroDTerminalRcrSpeciesDomainRuntime.hpp` is the first staged graph
adapter for the terminal RCR case. It implements the existing
`CoupledDomainRuntime` and `StagedFlowTransportDomainRuntime` interfaces,
separates hydraulic and transport trials, supports donor lookup before an
outgoing transport solve, and commits pressure/amount/volume only after both
trials are prepared. `make coupling-zero-d-terminal-rcr-species-runtime-test`
checks rollback, abort, retry, reverse distal flow, and no premature commit.
`make coupling-zero-d-terminal-rcr-species-graph-test` runs the production
`SpeciesPressureFlowComponentExecutor` with a **controlled flow source** and
this real RCR adapter. A precommit failure leaves both domains unchanged;
retry and a second, reversed-flow step pass per-edge and global amount gates.

The separate `t7-native-ale-species-graph-rcr-test` target runs a 24-tetra
native P2/P1 ALE flow and P1 species
domain between a configured-open-loop **native 1D source** and this RCR
species runtime. A fixed first step and rigidly translating second step pass
at 1, 2, and 4 MPI ranks in explicit and fixed-point pressure/flow coupling,
including a separate load whose distal RCR branch reverses and receives an
explicit distal donor concentration. Missing that donor is rejected before
commit. The positive and distal-reversed fixed-point cases each take 15 outer
iterations on the second step. A rank-0-only rejected precommit attempt
propagates to all ranks and leaves the native FEM and RCR committed states
and donor ownership unchanged; retry passes. Both graph edges, all domain
amounts, the RCR distal amount, and global balance are checked. This is a
bounded functional chain using project-owned FEM assembly, not a test double
at the 3D or RCR end.

`include/ZeroDFlowSpeciesCheckpoint.hpp` now stores an accepted 0D hydraulic
pressure and well-mixed volume/species amounts together, bound to the
hydraulic model, species catalog/source rates, and configured external donor.
Both source and terminal staged runtimes can capture only after an accepted
commit and restore only into a fresh owner. The bounded binary codec checks
schema, payload digest, and bundle epoch clock. The cross-process
`make coupling-zero-d-flow-species-checkpoint-test` saves source and RCR
shards under one published bundle manifest, loads them in a separate process,
rejects truncation/corruption/wrong domain/model/epoch, then verifies that
the next accepted states match an uninterrupted run. This is the focused
paired-0D checkpoint test; the full three-domain test is described below.

The native tetra staged adapter can now capture and restore its accepted
P2/P1 flow and P1 species checkpoints as a pair. Restore validates model
identity (including a required motion-model ID, species ID, and source rate),
step count, time, current ALE points, field sizes, and positive tetra
geometry before changing either committed field. Its 1/2/4-rank test checks
bitwise-equivalent restored checkpoint bytes and next-step field agreement
within `1e-12`; a 4-rank rerun exposed only `3.55e-15` species roundoff, so
next-step bytes are not claimed identical.

`make t7-native-ale-species-graph-checkpoint-test PETSC_DIR=/path/to/petsc
T7_SPECIES_MPI_RANKS=2` now publishes five accepted-step shards in one
manifest: graph pressure/donor controls, source reservoir, terminal RCR,
native FEM flow, and native FEM species. A separate process validates the
bundle and all five clocks, restores fresh owners, and advances the second,
translating step; 326 state and balance values agree with the uninterrupted
run on 1, 2, and 4 MPI ranks. This verifies a bounded full-graph file restart
for the 24-tetra functional case. Bundle compatibility now hashes the actual
graph topology/ports/species, all three domain model identities, and hydraulic
iteration/routing/amount controls. A changed source model or execution routing
or the graph species unit is rejected before shard loading; each field shard also carries its own model
identity. Epoch naming and the five-shard catalog are still fixture-specific,
so this is not yet a general production checkpoint coordinator for arbitrary
cases or physiological geometry.

For a non-handbuilt mesh, `make t7-native-tet-species-ftetwild-test
PETSC_DIR=/path/to/petsc FTETWILD_BIN=/path/to/FloatTetwild_bin
T7_SPECIES_MPI_RANKS=2` regenerates a labelled 0.03 m × 0.005 m pipe surface,
runs the fTetWild adapter, checks positive/quality-gated tetrahedra and wall,
inlet, outlet labels, then solves one native source→P2/P1 FEM/P1 tracer→RCR
step. In the four ordinary material variants, both edge species amounts must
be nonzero in the physical forward direction, with edge and global balance
residual below `1e-12 mol`. On each
generated mesh the script runs four material settings: baseline, native
diffusivity `1e-5 m²/s` only, body source `0.1 mol/(m³·s)` only, and both.
All use the same nontrivial pump donor (`5 mol/m³`) and `1e-5 m³` source/RCR
mixing compartments. Diffusion alone must change the native nodal-concentration
RMS departure from the initial `2 mol/m³` by more than `1e-8 mol/m³`, while
reporting zero source. In source-only and joint runs, native inventory must
rise by at least one-quarter of the reported positive source amount relative
to the matching no-source run; local runs gave about `2.12e-9 mol` of source
and a matching inventory increase. This isolates functional diffusion and
source effects on identical geometry within each rank run, but does not
establish spatial convergence or calibrate physiological tissue coefficients.
The script also runs a zero-initial-FEM-concentration inlet front. The staged
runtime's opt-in monotone graph-diffusion mode must give a nonnegative native
field, positive native inventory, and global residual below `1e-12 mol`;
the corresponding unstabilized graph must reject its negative port
concentration. The front need not reach the RCR within one step, so its
terminal species edge is allowed to carry zero amount. A checkpoint made in
default mode is rejected by a monotone-mode staged runtime.
The same generated fTetWild mesh is now also advanced through a second
rigidly translated ALE step (`-1e-5 m` in x). The native step count and
actual mesh displacement, both edge/global species residuals, nonnegative
concentration and positive inventory are checked. A fresh 646-tetra
same-mesh 1/2/4-rank run compared all six native concentration fields
(the four ordinary material variants plus fixed and moving fronts); its
maximum relative L2 difference was `1.62e-12`. This adds a non-handbuilt
moving-mesh functional check but does not establish spatial convergence.
On the separate 24-tetra channel, the same opt-in mode now also has a
two-step moving-ALE graph test: a rank-0 precommit failure rolls back all
three domains and donor ownership; a retry accepts the fixed first step;
five checkpoint shards are saved and restored in a new process before the
moving second step. The continued 326-value graph fingerprint matches the
uninterrupted path on 1/2/4 ranks, with nonnegative FEM concentration and
global balance within `1e-8 mol` (2-rank observed minimum
`0.002666 mol/m³`, residual `9.19e-17 mol`). Run
`make t7-native-ale-species-monotone-front-checkpoint-test PETSC_DIR=...
T7_SPECIES_MPI_RANKS=2`. This uses a controlled 24-tetra geometry, not a
physiological high-Péclet validation or a generic production coordinator.
Independent local 1/2/4-rank runs passed. Use
`make t7-native-tet-species-ftetwild-same-mesh-mpi-test PETSC_DIR=...
FTETWILD_BIN=...` for the full-field comparison; its gate is `1e-10` for
ordinary variants and `1e-8` for the near-zero fixed/moving fronts.
fTetWild may produce different valid meshes across repetitions, so independent
runs cannot be used for this field comparison. Neither test is a
mesh-convergence or physiologically calibrated material study. The FEM
Newton tolerance is
`1e-8`: at `1e-3`, the small `~1e-7 m³/s` port flow was inaccurate enough for
the graph's conservative-flow gate to correctly reject the step.

```bash
make t7-native-ale-species-graph-rcr-test PETSC_DIR=/path/to/petsc \
  T7_SPECIES_MPI_RANKS=2
```

Still absent are generalized production graph checkpoint orchestration,
distal-bed physiology, and calibrated absolute mixing volume. The earlier
controlled-source graph test remains a focused RCR lifecycle test, not a
substitute for this native FEM example.
