# PSC FSI time optimization — 20260915T035654Z

## Accepted / verified

P0 passed: 46014291, 225 checks, 1602 s / 1.78 allocated core-hours;
TotalCPU 3910 s, MaxRSS 666544 KiB. Accepted baseline remains P0D.
B0 dirty/untracked source + P0D diagnostics are immutable, see baseline-manifest.json.
Same-binary detail overhead median -0.158% (noise, not speedup). Hot worker work:
volume 760.530, wall 248.229, trace 78.070 s; overlapping, not additive wall.
Caller ghost 39.516 and scatter 27.004 s. W0 expanded/compact/FUNNELED/caller/
SINGLE-negative/worker-wall fault retry passed. Original Y W1 R/J exact archive.

P1-A before-change proof 46021328 passed: 4 identical preassembly/initial input
pairs; unused preassembly 130.313880 s. Removed ONLY adapter moving_.Assemble();
public Assemble and initial SolveTrial assembly preserved, no hidden consumer.
P1ABefore/P1A immutable numerical digests: read their *-source.sha256 files.

P1-A V0/V1/small V2 passed: numerical work 46066555 on r327, 8 CPUs/batch 8;
short reanalysis 46069954 validated saved output, no solver or compiler rerun.
Paired complete adapter processes 46060472 exit 0: zero initial residual,
zero KSP, post-finalize port/traction publication, first assembly throw + exact
same-input retry, original late-publication/prepare/reject/abort checks passed.
Prior test/source/binary hashes and exits reverified in job-46069954/correctness.json.
Three W1 process pairs AB/BA/AB, warmup + 3 samples/state, exact archived R/J;
mallocs 0, correct team/generations. Small complete FSI all native fields and
remaining assembly inputs exactly equal; full calls 20 -> 16. Before 533.813350 s,
after 438.588160 s (~17.8% elapsed reduction, ONE screening pair only).
Formal physics acceptance remains false for the original full cycle.

## P1-A promotion

Full original immersed suite 46070334 PASSED (7234 s, 8 CPUs, 16.076 allocated
core-hours, TotalCPU 42595 s, MaxRSS 673380 KiB); no geometry/quadrature/assertion
changes. Original Y W2 46071831 PASSED (3234 s, 7.187 allocated core-hours,
TotalCPU 9022 s, MaxRSS 500408 KiB). First-two historical history, all native
fields, coupling and remaining assembly inputs EXACT; full calls 70 -> 60.
Before 1682.416 -> after 1541.453 s, 8.379% reduction, ONE screening pair only.
P1-A promoted. Accepted algorithm baseline is P1A; immutable archives retained.

## Live P1-B / terminal P1-A

LATEST: 46096311 PENDING, 8 CPUs/15200 MB/2 h, ONLY missing throw-ready
off/on pair and matched small FSI. Seven completed assertion-bearing pairs
from 46083059 are strictly revalidated by source/binary SHA and selector-only
test delta, never misrepresented as a passed overall process. Prior actual W2
proof PASSED: 20 line-search/next-Newton identical pairs, redundant 315.454 s.
Normal/half/nan-first/failure/nan-all completed pairs: 7 Newton updates exact;
normal full calls 14 -> 8, half/nan-first 15 -> 9, exact operator and retries.
No production changes after immutable P1B. ValidationR1 SHA e43f2fda138e99c0d0107588eeacd2c47725d43afd1919083b83b427d8f49442.
46096441 PENDING, 1 CPU/20 min compile + native-export unit check, NO simulation.
46096555 afterok:46096311, FULL unchanged depth4/depth3 immersed suite, 8 CPU/4 h.
46096556 afterok:46096441, 1 CPU/10 min actual W2 staging-only preflight.
46096574 afterok:46096555:46096556, original Y W2 pair, 8 CPU/2.5 h.
W2 runner SHA d32662dd92575b0c8cdcf1f35db5d907d05aa056a94f3f5b506fc7546dc332fe.
Main calculations strictly sequential; only tiny no-simulation utility parallel.
Sleep THIS conversation 14400 s (4 h), then scheduler-only these FIVE handles.
Terminal-only fresh evidence; failure dependency means inspect prerequisite,
not blind resubmission. P1A still accepted; P1B unpromoted, full plan active.

TERMINAL UPDATE: 46096311 PASSED in 2247 s, 8 CPU, TotalCPU 11356 s,
MaxRSS 836324 KiB. All eight targeted pairs plus matched small FSI exact;
small FSI full calls removed 4, native fields/coupling exact, screen elapsed
370.210 -> 282.967 s (one pair only). 46096441 native build/export PASSED
(137 s, 1 CPU, no simulation); 46096556 W2 preflight PASSED (7 s, no sim).
46096555 FAILED only at original test line 697: obsolete accounting required
attempt_assembly_count >= 2*Newton updates. All prior physical/convergence
checks reached this point; this old pure-work expectation contradicts P1-B.
Updated ONLY its lower bound to initial full + one full candidate per update;
rejections may add calls. No physical/convergence gate or production code changed.
ValidationR2 SHA b5f8449b483ab99bc5234a18f5b2dd4f07c399dd8c1f7f6bd939ad1d2cd0f73e.
46112307 PENDING, complete original immersed suite 8 CPU/4 h.
46113038 afterok:46112307, original Y W2 8 CPU/2.5 h. Old dependency-never
46096574 explicitly cancelled; no lost computation. Conversation sleep 10800 s,
then scheduler-only these two handles, terminal-only evidence. P1B unpromoted.

P1-B PROMOTED: full original suite 46112307 PASSED, elapsed 5730 s,
TotalCPU 33394 s, MaxRSS 667932 KiB. Original Y W2 46113038 PASSED,
elapsed 3294 s, TotalCPU 8455 s, MaxRSS 494620 KiB. Exact original first-two
history, every native field, coupling/conservation/Newton/KSP and reconstructed
full/reuse input sequence. Full assemblies 60 -> 40; one screening pair elapsed
1816.384 -> 1463.325 s (-19.437%, not final repeated claim). P1B accepted.

P2 instrumentation adds default-off result payload bytes and worker completion
spread; unit asserts disabled mode does no profiling work. P2Validation SHA
9352ab9e2aa132e17a33dc3a8b79601685e7c8ed3add5e7e7cf38a67027d1ba0.
MAIN 46126551 submitted, 8 CPUs/15200 MB/3 h: same immutable binary, fixed
8 threads, batch 8/16/32 in forward/reverse process order, exact archived W1
operators, Mat mallocs zero, prepare/consume/spread/payload and GNU MaxRSS.
Select smallest within 2% of fastest only if faster than batch8 in BOTH processes
and RSS<70%; then full W0 serial/OpenMP/FUNNELED/SINGLE-negative/expanded/
compact/ordered scatter/worker throw retry. P2 unpromoted; W2 follows only pass.
Conversation sleep 10800 s, scheduler-only 46126551, terminal-only evidence.

P2 harness correction: 46126551 failed in 9 s BEFORE runner/build/simulation
because common sbatch passes named flags. Added standard --root/--variants/
--case-list parsing and frozen case manifest; numerical candidate unchanged.
P2ValidationR2 SHA e476464325c9ec258c004b2c80b431435c3badb3a4c68aec4c1a6ed9412e4b51.
46138318 RUNNING, same 8 CPU/15200 MB/3 h P2 gate. Next conversation sleep
10800 s; scheduler-only 46138318, terminal-only evidence.

P2 W0/W1 46138318 PASSED in 868 s, TotalCPU 2914 s, MaxRSS 661716 KiB.
Exact archived operators and zero Mat mallocs at all batches. Median W1 assembly:
batch8 16.821 s, batch16 13.453 s, batch32 12.283 s; batch32 faster in BOTH
process orders and selected/fastest. Its MaxRSS 265560 KiB (<70% of 15200 MB),
resident result payload max 50528256 bytes, max worker completion spread 0.4395 s,
prepare 0.0291%, consume 11.805%. All W0 serial/OpenMP/caller/FUNNELED/
SINGLE-negative/expanded/compact/ordered scatter/worker fault retry passed.
P2 remains unpromoted pending W2.

46149676 no-simulation native build PENDING; 46149881 afterok staging preflight;
46149882 afterok original Y matched W2, baseline P1B batch8 vs P2 batch32.
Each sample saves actual batch; exact native/history/coupling and exact complete
assembly/reuse sequence/call count required, plus no >3% timing regression.
P2W2Runner SHA 1fb7037b93258f686972ca07b6802f15eb2f1fd6a6b2232439c226856b5b57e4.
Conversation sleep 10800 s, then scheduler-only these three handles.

P2 PROMOTED: no-simulation build 46149676, preflight 46149881 and original Y
W2 46149882 all PASSED. Exact history/native/coupling and complete full/reuse
sequence/call count; baseline batch8 1326.014 s, batch32 1163.091 s (-12.287%,
one screening pair). P2 accepted configuration is 8 threads, batch32.

P3 implements true ResidualOnly kernels: volume and trace skip dense tangent
allocation/loops; wall skips tangent; ghost uses n×n matrix-free coefficient
action rather than 4n×4n J; ports/gauge never write Mat. Vec, measurements,
controller, diagnostics and accumulation order retained. Independent R/J
generations forbid stale-J Newton solves; line search uses residual-only and
rebuilds full at the next Newton iteration. Full/residual calls recorded.
Preliminary 46159278 compiled then failed a TEST that expected J current after
residual-only convergence; stale J is intentional and generation-guarded.
Test now constructs a fresh full at the same state, requires exact R and fresh J.
Production residual entry is internal; only testing macro exposes it.
P3ValidationR3 SHA 5596d6009419bea4dc820a8894d8f621902b16541e94e7f7412720ae37e59c6c.
46161562 FAILED during Slurm compilation: the production line-search path still
called a removed helper name. This was fixed by calling the private request-based
assembler directly. P3ValidationR4 SHA
cde3c8361ae042a956f549c0912b51b35c9ae4655ee7c6cec42c223be96f1982.
46163904 PASSED targeted eight branch pairs with 8 threads and accepted batch32
in 1179 s. zero, one-step, normal, half, nan-first, failure, nan-all, and
throw-ready all matched bitwise; every check and process exit passed. P3 remains
unpromoted: complete original immersed suite, matched small FSI, and original-Y
W2 gates are still required.

P3VerifierR1 SHA e59f327fb1983b22ee6602e2fe69e79f1714f3f1f48ea940d93fd6061d45fb30
adds state-aligned P2 full+reuse versus P3 residual+rebuild verification without
weakening native/history/coupling/conservation gates. 46168275 runs the complete
original immersed suite. Independent no-simulation native build 46168602 feeds
preflight 46168977. afterok 46168976 is matched small FSI; 46168978 original-Y W2
depends on small FSI and preflight. Initial states: 46168275/46168602 RUNNING,
the rest dependency-pending. P3 remains unpromoted. Conversation sleep 3600 s,
then scheduler-only these five IDs and terminal-only fresh evidence.

Terminal lightweight gates: 46168602 PASSED native exporter plus P2/P3 heartbeat
builds in 107 s, simulation_executed=false. 46168977 PASSED exact source/binary/
surface/history staging preflight in 10 s, simulation_executed=false. 46168275
remains RUNNING (493 s at last scheduler check); 46168976 and 46168978 remain
dependency-pending. Conversation sleep 3600 s, then scheduler-only chain check.

Static queued-input audit at 2026-09-16T19:25:05Z: P3VerifierR1 contains all
three case cards, run_p3_gate.sh, verify_p3.py and P3-aware verify_w2.py; every
archive member hash exactly matches the submitted worktree file. No scheduler
or running log was polled during this audit.

Static P3 contract audit then found R4 public AssembledNegativeResidual() and
AssembledJacobianAction() could expose stale storage after a state generation
change, although Newton itself was guarded. This violates the plan's current-R/J
invariant. R4 46168275 (802 s consumed), pending 46168976 and 46168978 were
cancelled; completed no-sim evidence remains preserved. P3ValidationR5 SHA
f083218030094ec8816bbb1401a606fa92451aff2c3d180ec38f074ee227620e adds
public generation guards and original-suite stale-rejection/fresh-rebuild tests.
P3VerifierR2 SHA 8d1ca7988b3dd8104cb60cbd4c21038a6fec8b392aa84f2999d45f433669c18d.
R5 targeted 46169775 and no-sim build 46169774 RUNNING; full 46169825, small FSI
46169826, preflight 46169827, W2 46169828 dependency-pending. Sleep 1200 s then
scheduler-only these six IDs. P3 remains unpromoted.

R4/R5 source-file hash audit confirms the only numerical/test files changed are
ImmersedTransientFlowRuntime.hpp and test_immersed_transient_flow.cpp. Exact diff
is two accessor generation checks plus stale-rejection and fresh-full assertions;
no weak form, solver threshold, geometry, quadrature or physical option changed.

R5 targeted 46169775 terminated in 103 s exactly at the fixture's stale-J capture
after one-step residual-only convergence; stderr is the new production guard
message. This confirms the guard and identifies harness-only failure. Pending R5
46169825/46169826/46169828 were cancelled; R5 no-sim build 46169774 and preflight
46169827 passed and remain archived. P3ValidationR6 SHA
192352c89fcb33bfd115beec0912defada78df28f8b6f7a60994e287655cada7 fixes
fixture semantics: same-generation residual-only may retain current J, while a
SetTrialState generation change followed by residual-only must reject stale J;
targeted comparison obtains J only after full rebuild. P3VerifierR3 SHA
245e4a5a06a04d17256359185ba76c43c5f4f374d9b72be588f22cabbd40bcc5.
R6 targeted 46170114 PASSED in 1174 s: all eight branch pairs, retry/operator
checks, full-versus-residual results, and generation contracts passed bitwise.
No-simulation build 46170113 and W2 staging preflight 46170195 also passed.
Full original suite 46170193 is resource-pending (no scheduler start estimate),
then small FSI 46170194 and W2 46170196 remain dependency-gated. Adaptive
conversation sleep 1800 s, then scheduler-only query of those three handles.

Scheduler update 2026-09-16T20:50:22Z: full suite 46170193 is RUNNING on r224,
elapsed 46:16. Historical runtime is about 95.5 minutes, so the next adaptive
conversation sleep is 3300 s (55 minutes). This is a conversation wait, not an
SSH/login-node sleep process. Small FSI 46170194 and W2 46170196 remain gated.

Terminal update: full original suite 46170193 PASSED in 3865 s. Matched small
FSI 46170194 PASSED in 433 s with every native field and coupling/conservation
summary bitwise exact; full assemblies 12 -> 8 plus 8 residual-only calls, and
one screening pair elapsed 196.923 -> 151.988 s (~22.8%). W2 46170196 is now
RUNNING; at 7:30 elapsed, use a 1800 s conversation sleep based on the historical
35--42 minute W2 duration, then query only that Job ID.

P3 PROMOTED: original-Y W2 46170196 PASSED in 2022 s. Original first-two
history, every native field, coupling, formal short gates, state-aligned logical
work sequence, and current R/J generation checks are exact. Full assembly calls
40 -> 30 plus 30 residual-only calls; one screening pair elapsed 1035.563 ->
965.625 s (-6.754%). Accepted variant is P3ValidationR6_batch32. No live job;
next action is P4 selection from the measured accepted-P3 hot-kernel profile.

P4 selected only the measured hottest remaining worker kernel: full-assembly
volume work ~36.5 of ~52.3 aggregate worker-seconds per assembly, versus wall
~12.0 and trace ~3.75. P4Validation SHA
28d5f7946d9c30f8e7d6c441c9c7aede8a8c21c07c9c2f49945cb7a37d03471d
adds an immutable geometry-epoch volume basis/gradient/Hessian/physical-point
cache, default-on with a disabled control and visible hits/misses/bytes. Targeted
Slurm compile/contract Job 46184469 is resource-pending. Conversation sleep
900 s, then query only this Job ID; matched W1 is gated on targeted success.

While 46184469 remained resource-pending, the matched W1 harness was prepared
but not frozen/submitted: three independent process pairs in off-on/on-off/off-on
order, one warmup plus three timed full assemblies per process, exact residual/Jv
hashes, constructor cost and cache counters/bytes. No result is assumed. With no
other independent work, next conversation sleep is 1800 s.

Scheduler update 2026-09-16T23:01:08Z: 46184469 remains PENDING for Priority;
`squeue --start` reports N/A. After repeated pending checks and with W1 preparation
complete, the next adaptive conversation sleep is 3600 s.

46184469 failed only while copying evidence because the submitted root was
relative and the runner had changed directory; no durable numerical evidence is
claimed. P4Validation source is unchanged. R2 Job 46192353 resubmitted with
absolute root/card paths and is resource-pending. Based on the prior ~50-minute
queue wait and ~1-minute execution, next conversation sleep is 2700 s.

Scheduler update 2026-09-17T00:48:55Z: 46192353 remains PENDING for Priority
with no start estimate. No independent work remains; next adaptive conversation
sleep is 5400 s.

P4 targeted R2 46192353 PASSED in 64 s. Cache on/off R and Jv were bitwise
exact; same-epoch hits, different-epoch invalidation, coefficient/port safety,
fault retry atomicity, and visible counters passed. Depth2 cache: 25752 points,
176968392 bytes; process MaxRSS 552892 KiB. P4ValidationR2 SHA
bf43dca26b8d77f5a73870e61153915be0eb36b33d0cc3fb2a6bff72c54b25cb
adds only W1 harness files; production cache files hash-identical to targeted.
Matched W1 Job 46201812 is RUNNING: three independent off/on process pairs.
Next conversation sleep 1200 s.

P4 matched W1 46201812 PASSED in 878 s. All three independent off/on process
pairs have bitwise-identical residual/Jv hashes and candidate reductions
12.860%, 12.751%, and 13.150% (median 12.860%). Depth3 cache: 180480 points,
1240259208 bytes. P4Verifier SHA
96fb6c015c4a2b4ebd21722650434acada9ce5a9c1e84ff582b75ad1b4668a49.
W2 chain submitted: no-simulation build 46203125 -> preflight 46203289 ->
full original-Y P4 W2 46203291. Build is priority-pending; later jobs are
afterok-gated. P4 remains unpromoted until W2 exactness/cache/timing gates pass.
Next adaptive conversation sleep 1800 s; query only these three job IDs.

After that wait, 46203125 remains pending because matching nodes are down,
drained, or reserved for higher-priority partitions; 46203289 and 46203291
remain dependency-pending. No job log was read. Next conversation sleep is
3600 s before one scheduler-only check.

P4 **PROMOTED**. Build 46203125 and preflight 46203289 passed without running
a simulation; original-Y W2 46203291 passed in 2328 s. History, every native
field, coupling, state-aligned assembly sequence, cache epoch/key/counter
contract, and RSS gate are exact/passed. W2 elapsed 1157.074 -> 1157.644 s
(-0.049%, within the 3% no-regression gate); this one screening pair is not a
final speed claim. Combined with three independent W1 reductions of
12.751--13.150%, accepted variant is P4ValidationR2_batch32. P5 begins by
profiling the actual PETSc backend and factorization split on Slurm.

P5 profile 46208944 failed during staging because a patch transport artifact
inserted literal plus signs at shell continuation boundaries;
simulation_executed=false. The runner was corrected and syntax-checked.
Replacement 46208974 profiles one original-Y step with PETSc log_view and KSP
view; next adaptive conversation sleep is 1200 s.

46208974 reached only demo preflight, then failed because the runner pre-created
the protected results directory and used an invalid PETSc viewer string; no
physics step ran. Both harness issues were corrected. R3 Job 46211622 is running
on r190; next conversation sleep is 900 s.

P5 profile R3 46211622 PASSED in 483 s. Actual solver is FGMRES with right
PCLU using PETSc's factor package. Five runtimes performed 15 Newton solves:
MatLUFactorSym 5 calls / 5.535 s and MatGetOrdering 5 / 0.239 s prove
same-pattern symbolic/ordering reuse is already active; MatLUFactorNum
15 / 86.189 s dominates PCSetUp 91.976 s. The non-redundant candidate reuses
the previous same-runtime LU only as a right preconditioner for the current
Jacobian, with iteration/explicit-true-residual gates and one current-J rebuild.
P5Validation SHA b60ff956609a7bb092b72e13183894be23fac5dd1df72aeb376543813ebf3c8b;
targeted Job 46213444 pending. Next conversation sleep 1200 s.

P5Validation targeted 46213444 compiled, then exposed an incorrect test
assumption: strict-tolerance telemetry showed five stale-PC attempts correctly
triggered current-J rebuilds and one was safely accepted, while the test had
required zero rebuilds. R2 verifies the actual accounting invariant and adds a
relative iteration-increase gate against the latest fresh factorization.
P5ValidationR2 SHA
9ba7ef56189f19a29126acb5f461dcfee94cd1204e3e7490035eaeeb01b86c5e;
Job 46214496 pending. Next conversation sleep 900 s.

P5 targeted R2 46214496 PASSED in 130 s. Disabled built 7 current LUs;
candidate attempted reuse 6 times, accepted 1, rebuilt 5, and built 6 current
LUs. Forced fallback rebuilt all six stale attempts and was bitwise identical
to disabled; candidate relative state difference was 1.18e-18. Immutable
heartbeat build 46215836 passed in 145 s. Matched original-Y one-step W1
46215911 is dependency-releasing: three independent AB/BA/AB process pairs,
frozen per-quantity rtol 1e-8 and median total-time reduction gate 5%.
Next conversation sleep 3600 s.

P5 matched W1 46215911 PASSED all 94 gates in 3015 s. Three independent
original-Y one-step AB/BA/AB pairs reduced elapsed by 11.450%, 10.593%, and
10.884% (median 10.884%); setup fell from 97--101 s to 35.6--38.2 s. Each
candidate process built five current LUs and accepted ten reuse attempts with
zero rebuilds. Baseline repeats were bitwise deterministic and all candidate
native physical fields passed the frozen per-quantity rtol 1e-8. P5Verifier
SHA d4ed66aaf1fe81fb9c7365db016dc0a14c0fe3896146e5ce93908bb442309130.
Preflight 46219622 is running; full W2 46219623 is afterok-gated. Next
conversation sleep 2700 s.

P5 preflight 46219622 PASSED. W2 46219623 completed both simulations but the
initial verifier marked one check failed: near-zero dimensionless continuity
was incorrectly compared with a relative scale of 1e-8. Its baseline norm was
4.12e-14 and candidate difference 6.00e-15; all 26 native physical field gates,
coupling, formal gates, logical path, true residual and reuse accounting passed.
Elapsed 1061.255 -> 937.780 s (-11.635%); setup 206.206 -> 74.199 s; candidate
built 10 current LUs and accepted 20 reuses with zero rebuilds. P5VerifierR2
uses only 256*epsilon*sqrt(rows) at natural dimensionless scale 1 for this
near-zero metric, SHA
ad7f7ec989d82acff70dcb5a404f29745523457ab81facfa33ca8b2a921de74c.
Verification-only 46224110 and an independent full W2 46224111 are queued.
Next conversation sleep 3000 s; no promotion until both pass.

R5/R6 source hash audit: every production include/src file is byte-identical;
only test_accepted_line_search_assembly.cpp and test_immersed_transient_flow.cpp
changed. Thus the R6 rerun is solely a corrected validation-harness revision.

MAIN 46083059 TERMINAL FAILED 124:0, elapsed 02:30:49; targeted-regression
process budget expired (9000 s), before matched small FSI. Status/stage checked
after user interrupted conversation sleep to request Slurm/patch teaching.
Durable terminal evidence exists in job-46083059; numerical diagnosis pending,
do not rerun blindly or promote P1B. No primary job currently live.
Original allocation 8 CPUs/15200 MB/4 h. Frozen P1B SHA
a588354c3a9b53049097b2fbe9da544312a8f951287dd1c801cada8570779d34.
Baseline P1BBefore SHA 32b17336ad3dd025bea354c0f2df9100b68cfc6071d2bfcbcc7162c4dcb95a9f.
Independent P1B-runtime-delta.patch. Work-only option can disable local reuse.
Slurm first proves prior actual W2 line-search -> next-Newton identical inputs,
then candidate eight off/on additional depth-3 branch pairs: zero/one/multistep,
first reject/half damping/NaN rejection/backtracking failure/ready exception,
fresh R/J action + first assembly failure/retry for EVERY pair. Original depth-4
suite untouched. After targeted gate complete, matched original small FSI pair.
Earlier sleep interrupted by user; terminal evidence is now diagnosed above.
Next action is four-hour adaptive conversation sleep for the new gate chain.
P1B unpromoted until original Y W2 also passes; do not infer global completion.

MAIN 46070334 COMPLETED 0:0, 8 CPUs/15200 MB/4 h.
Runs FULL original immersed unit suite with original depth 4/branch depth 3,
quadrature, controller/gauge/wall/compact/Newton/transaction assertions unchanged.
GNU process timeout 12600 s leaves copy-back margin. Frozen numeric P1A source.

NEXT 46071831 COMPLETED 0:0, 8 CPUs/15200 MB/2.5 h:
original Y W2 2-step matched before/after, dt .05, 16x16x5, pressure 5 + pulse
amplitude 30/period 1, cut depth 1/rescue 3, inertial gamma 1, formal mode.
Real staging preflight 46071756 passed in 9 s, simulation_executed=false.
Reuses hashed binaries from passed no-simulation native/export compile 46071195;
its native/baseline prefix 46070816 was rehashed/reused, not recompiled.
Native reference is DEFAULT OFF and NOT checkpoint; optional post-commit CSV
uses sorted global IDs, 17-digit scientific values, SHA and true per-step gates.
Workload W2Tools-source.tar; runner/validator W2Runner-source.tar, all digests saved.
Historical small input historical-heartbeat-history.csv SHA verified, compare
original first 2 rows, not whole-cycle shrink. All field/coupling/assembly inputs
must match; real full assembly count decrease and no W2 timing regression.

## Exact next action / constraints

Execute independent P1-B local accepted-line-search FULL R/J reuse, preserve
all convergence checks and invalidate before every candidate update. Validate
normal/zero/multistep/damping/failure/NaN/retry/operator exactness before promotion.
After submission use adaptive conversation sleep, scheduler-only checks and
terminal-only evidence. No background monitor/external handoff.
P1A promoted after W2; P1-B submitted, P2-P7 pending, full 20-step best candidate,
repeated fair final comparison, native VTK expansion/contraction/archive pending.
Original 17/20 conservation failures unchanged; formal_physics_acceptance=false.
All computation/large builds in Slurm LOCAL; no VTU restart or old scratch input.
User dirty/untracked edits preserved; no reset/clean/commit/push, old full solver
or old VTK repair rerun. Detailed event ledger: docs/progress/
BIFURCATION_FSI_TIME_OPTIMIZATION_RESULTS.md and state.json.

P5 verification-only 46224110 failed without simulation because it checked saved
8-thread assembly evidence against the 1-thread verification process environment;
dependent 46224111 was therefore cancelled. P5VerifierR3 checks saved evidence
against the case-recorded thread count and records verification threads separately.
Frozen SHA 24db08eeac25334f72f8b33dc6cbc9aed72e1fa23f94de2408577f9899d267c5.
Replacement verification-only 46225784 is RUNNING; independent 8-thread W2
46225785 is afterok-gated. P5 remains unpromoted. Next adaptive conversation
sleep is 3300 s; query only those two jobs and read evidence only when terminal.

P5 verification-only 46225784 **PASSED** in 10 s with
`simulation_executed=false`: all corrected history/coupling/native/reuse gates
passed, and its manifest distinguishes evidence threads 8 from verification
threads 1. Independent W2 46225785 is RUNNING on r211. Next conversation sleep
is 2700 s based on the prior 34--39 minute W2 runtimes; query only 46225785.

P5 **PROMOTED**. Independent original-Y W2 46225785 passed all 214 checks in
2162 s. The fresh matched sample improved 1149.077 -> 997.231 s (13.215%);
solver setup improved 239.788 -> 89.881 s. Candidate accounting was 10 current
LU builds, 20 accepted same-runtime reuses, and zero rebuilds. Together with
three W1 reductions of 10.593--11.450% and the earlier independent W2 sample
of 11.635%, accepted variant is P5ValidationR2_batch32. Full-cycle formal
physics acceptance remains false. Next action is the plan's conditional P6
selection from measured P5 evidence; no Slurm job is currently live.

P6 conditional gate is active: accepted-P5 W2 geometry was 352.675 s / 35.37%
of elapsed. Existing details attribute ten moving BeginTrial calls to 229.761 s
geometry build, 52.272 s runtime/layout/preallocation, and 48.996 s velocity +
scalar extension; spatial-index construction was only 0.009 s. P6Profile adds
timing only, preserving operation order, and is frozen at SHA
0619efb045869debda2ffecd0e328ebe5550e2432081f4be0aa1e10e52f902f6.
Slurm Job 46230573 is RUNNING on r186. Next conversation sleep is 900 s.

P6 profile 46230573 **PASSED** in 503 s. Six runtime builds shared one exact
active-layout pattern, but active-layout construction itself was only 0.013 s.
Jacobian preallocation cost 12.172 s total; excluding the unavoidable initial
epoch gives an ideal one-step saving ceiling of about 10.1 / 452.8 s = 2.24%,
below the plan's 5% uncertainty threshold. Geometry-dependent volume-basis cache
(13.669 s), gauge weights (3.748 s), quadrature and extension cannot be reused.
P6 cross-epoch cache is therefore skipped; profiling-only production changes
were removed and the runtime is byte-identical to P5ValidationR2. P7 is also
skipped because existing Aitken coupling (~4 iterations/step, 11.27% of P5 W2)
does not dominate geometry or assembly. Next is R6 thread scaling 1/2/4/8,
fixed batch32, followed by W2 for the fastest two configurations.

R6Scaling archive SHA
5506098e7efa796b243382a6855523d0d174bc02703952b305adaebd199ad522;
all 153 production include/src/Makefile files are byte-identical to accepted P5.
W1 Job 46232121 measures 1/2/4/8 threads in forward and reverse process order,
fixed batch32, exact two-state R/Jv hashes and zero matrix mallocs. It is pending.
Next adaptive conversation sleep is 1800 s; W2 waits for the frozen top-two result.

R6 W1 46232121 **PASSED** in 1762 s. Paired medians for 1/2/4/8 threads
were 53.267/28.817/16.741/11.306 s; 1->8 speedup 4.711, efficiency 58.9%.
All two-state operator hashes were exact and matrix mallocs zero. Frozen W2
top two are 8 and 4; W1's within-5% choice is 8. R6Verifier source SHA
565e1db4a2dc7f0153846fec022f7329b530e005c004a6e66741111406206022,
again 153 production files byte-identical to P5. Original-Y W2 Job 46235778
is pending. Next adaptive conversation sleep is 3000 s.

R6 W2 46235778 completed both simulations successfully: 8 threads 14:50.58,
4 threads 17:46.25. The verifier then failed before reading evidence because
the runner saved `cases-r6-w2.json` while the verifier opened `case.json`.
No numerical gate failed. R6VerifierR2 fixes only that filename and is frozen
at SHA 24c1aea54a6bc1beda0f515b54c21bce03e3bd4df310a920d220ba3dee6e4d75.
Verification-only Job 46238720 is pending; next conversation sleep is 600 s.

R6 **PASSED**. Verification-only 46238720 exposed only a second harness
assumption: it required exactly one cumulative `hpc_profile` record, while the
driver intentionally emits step-1 and final step-2 records. R6VerifierR3 takes
the final cumulative record and is frozen at SHA
4b5ab21428f78197777abfec9bb5b9ece034cfcf7320ff8e0c7ecf69ff6da5ab.
Verification-only Job 46239782 passed every saved numerical and configuration
gate without rerunning a simulation. Cumulative solver elapsed was 890.242 s at
8 threads and 1065.908 s at 4 threads, so production selects 8 threads with
batch32. History, coupling, all native fields, assembly and linear paths are
exact; cache/reuse/true-residual gates pass.

R7Final source is frozen at SHA
0aeb8a2d12cc09035514b27c838defd84c55916a52d2b8e6e693f2d9b96d4993;
all 224 production files compared are byte-identical to accepted
P5ValidationR2. Final 20-step Job 46243889 uses 8 CPUs, 15200 MB, batch32 and a
5-hour wall limit, and performs the full CSV, expansion/contraction, native VTK
9.4.1 and ParaView package/hash gates. It is resource-pending with no start
estimate. Use one 10800 s conversation wait, then query only this Job ID and
read evidence only if terminal.

R7 and the optimization-plan goal are **COMPLETE**. Job 46243889 completed 0:0
on r030 in 1:58:34; solver elapsed was 6968.953 s. All 19 final checks passed:
20 finite accepted outputs, frozen-tolerance history/coupling comparison,
expansion and contraction, native VTK reads and exact frame topology/times,
assembly/cache/preconditioner/true-residual contracts, and archive integrity.
The historical 13603 s comparison is 48.769% lower but remains explicitly
cross-job/cross-node, not a paired speedup claim. The original 17/20
conservation failures are preserved, so formal physics acceptance remains false.
Final archive: `job-46243889/heartbeat-paraview.tar.gz` (93174492 bytes), SHA256
5936dd2905e2dbf1dcf5a0bc3314bf1ef3180b44ed3a7e64985e14932613dce5.
No Slurm job is live and no further calculation is required by this plan.
