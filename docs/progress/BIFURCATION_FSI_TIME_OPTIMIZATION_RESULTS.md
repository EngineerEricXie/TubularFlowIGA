# PSC FSI time optimization execution

Execution started 2026-09-15 UTC after reading the complete handoff and linked plan.
Existing uncommitted changes are preserved; no commit/push or long full-cycle run.

Evidence: `benchmarks/fsi-time-opt/20260915T035654Z/STATE.md` and `state.json`.
B0 (dirty/untracked source) and P0D (diagnostics only) are immutable source archives.
Baseline audit establishes exact match to the tested parallel flow header and all
other original include/CPU headers; valid-case driver changes are non-arithmetic.

P0 W0/W1 validation submitted as **46014291**, 4 CPUs, 7600 MB, RM-shared,
account mch260002p, 1-hour limit (at most 4 allocated core-hours). All compilation
and tests use Slurm compute `$LOCAL`; binaries, input/workload/source hashes,
logs, environment, correctness/timing manifests are preserved on exit.

Implemented default-off assembly detail profiling, private worker kernel counters,
caller/batch timing, geometry/extension stage timing, call reasons and input/state
identities. W1 uses the existing archived driver with its required 0.05 s evaluated
geometry; no runtime preflight was relaxed. Same-binary on/off controls measure
profiling overhead, separately from baseline/candidate comparison.

**P0 gate passed in 46014291.** All 225 checks passed; B0/P0D R/J action hashes
match the archived tested operator for zero and nonzero Y-vessel states. Expanded/
compact, caller/FUNNELED/SINGLE negative checks and worker/wall fault retry passed.
Saved binaries were independently rehashed and match their node-local manifests.
Current numerical source file hashes still match the frozen accepted P0D source.
Detail/coarse wall accounting passed. Median same-binary profiling overhead was
-0.158% (noise-level, not a speedup); all three overhead pairs were below 3%.
Worker hot paths across 21 detailed assemblies: volume 760.530 work seconds,
wall 248.229, trace 78.070. These overlapping worker times are not additive wall
time. Caller ghost 39.516 and scatter 27.004 wall seconds are also measured;
scatter/ports are nested in consume. Baseline assembly process medians:
19.695, 19.765, 19.772 s; P0D disabled: 19.778, 19.774, 19.776 s (no meaningful
performance change claimed for diagnostic-only work).

Job elapsed 1602 s, 4 allocated cores: **1.78 allocated core-hours**.
Slurm TotalCPU 3910 s; MaxRSS 666544 KiB (~651 MiB). This accounting includes
builds and tests; CPU utilization is not parallel efficiency or solver speedup.
P0D is accepted as the diagnostic baseline. P1-A is now at before-change trace
proof completed; P1-B/P2-P7 remain pending. No full-FSI speedup claim yet.
The original formal conservation failures are unchanged; formal acceptance is false.

Latest user instruction: monitor with this conversation's adaptive sleep (20 minutes
or longer when remaining-work estimates justify it), not a background script or
external-event handoff. Read only scheduler state after sleep, then new evidence
on terminal. The former script monitor was stopped without touching the job.

P1-A before-change proof **46021328 passed**: four preassembly/initial pairs had
identical complete input identities, with 130.313880 s unused preassembly wall
time. Job 719 s / 0.799 allocated core-hours; TotalCPU 1995.368 s,
MaxRSS 747176 KiB. Only the redundant adapter Assemble call is now removed;
the public API and original SolveTrial initial assembly remain. Candidate is not
promoted yet. Test-only initial-assembly fault seam and zero-residual/retry checks
are frozen with identical test source in before/after variants.

V0/V1/V2 verification **46022197 submitted**, 8 CPUs / 15200 MB / 2-hour limit.
Includes unchanged original operator/gauge/port assertions, adapter zero-residual,
first-assembly failure/retry, late-publication and prepare failure lifecycle checks,
six expanded/compact volume fault/retry cases, and same-allocation complete small
FSI comparison with exact native reference hashes and assembly trajectories.
Original Y two-step W2 and final cycle remain pending. Initial state PENDING;
next conversation sleep 45 minutes, adapt upward if warranted after scheduler check.

46022197 reached the hard Slurm time limit (7215 s, 16.033 allocated core-hours,
TotalCPU 41925 s, MaxRSS 858768 KiB); EXIT copy did not execute. No gate inferred.
Diagnostic recovery 46031486 on r282 found no retained scratch. This is an
inconclusive execution-budget/evidence failure, not a numerical pass or failure.

Durable bounded gate runner now saves completed-stage evidence and stage cursor
before the next expensive stage, with per-process timeouts and cleanup margin.
No background monitoring was introduced. Numerical before/after archives unchanged.
**46031776 RUNNING r026**, 8 CPUs / 15200 MB / 2.5-hour limit: full adapter tests
and matched complete small FSI. Full immersed suite and six volume tests are
isolated for separate budgets. Next conversation sleep 1 hour; P1A not promoted.

46031776 failed in the new zero-residual test (170 s / 0.378 allocated core-hours,
TotalCPU 234.702 s, MaxRSS 858728 KiB). New assertion used the wrong transaction
boundary: Diagnostics().ports retains committed data until finalize. Original
abort/solve/prepare unit assertions establish this contract. No production port
code was changed. Corrected test checks unchanged ports through solve/prepare,
independent exact zero trial flux, measured positive-area zero-flow ports after
finalize and no trial leakage. This strengthens the transaction-aware regression.

**46057889 submitted**, paired complete adapter tests against both unchanged
numerical archives plus complete small-FSI same-allocation comparison. Corrected
common tests/validator frozen as P1AValidationR2, SHA256
0435f2794d6f97ce2b2fa7bc4b8c8a1bb8250a90834814c294b4f5feb36f54b9.
8 CPUs / 15200 MB / 3.5-hour hard limit; per-process bounds and durable stage saves.
P1A is still unpromoted; broad immersed/volume, original Y W2 and later cards pending.

46057889 failed in the baseline new zero test: reject/abort helper also rejects
conservation, while committed transition conservation is legitimately readable
after finalize. No before/after numerical divergence was observed. Original
helper semantics preserved; post-finalize test now rejects five trial-only
accessors, requires idle/no active trial, and compares exact retained conservation
and committed traction. All original lifecycle assertions remain.

**46060472 submitted**, same complete paired gates with P1AValidationR3 common
tests, SHA256 8e991d78a2ceea28fd41a5c2970d0731bd8e25fa52604cec211e5f18e06ceddd.
No production algorithm changes beyond the proven redundant adapter call removal.
Continue with one-hour adaptive conversation sleep and terminal-only evidence.

**46060472 both complete adapter processes passed (exit 0).** Zero-residual/
first-assembly fault retry markers and original lifecycle assertions passed in
both variants; physical summary identical. The overall runner then failed before
FSI execution because binary and sample directory shared a name. This is an
execution-path defect, not numerical failure. Elapsed 2063 s / 4.584 allocated
core-hours, TotalCPU 7723 s, MaxRSS 798604 KiB; durable adapter evidence retained.

Runner now separates binaries/ and samples, correctly handles empty regression
lists, and classifies process timeout separately. **46066555 submitted**, ONLY
missing W1 exact operator 3 process pairs AB/BA/AB and complete small FSI pair.
Prior successful adapter stages are checked by exact test/numerical source and
saved binary hashes plus process exit 0, not rerun. P1AValidationR4 SHA256
b58a245dace296e09ebcb067165ab93efb465e18d226fed8dcdc102b797bd354.
8 CPUs / 15200 MB / 2.5-hour hard limit; conversation sleep 1 hour.
P1A promotion still requires remaining gates, including original Y W2.

P1-A V0/V1/small V2 **passed**, numerical job 46066555, verification-only
46069954 (17 s): prior complete paired adapter source/binary hashes and process
exits reverified; three independent W1 process pairs AB/BA/AB matched archived
R/J hashes with zero Mat mallocs and correct teams/generations. All native small
FSI fields, coupling/conservation and remaining assembly inputs were identical.
Full assembly calls 20 -> 16. Single screening complete-FSI elapsed
533.813350 -> 438.588160 s (~17.8% reduction), not a repeated final speedup claim.
46066555 cost 2049 s / 4.553 allocated core-hours, TotalCPU 7118 s,
MaxRSS 829308 KiB. Original verifier parent-path error corrected; saved results
reanalyzed in Slurm, no solver rerun. P1A remains unpromoted pending W2.

Full original immersed suite **46070334 RUNNING**, 8 CPUs / 15200 MB / 4 h,
per-process 12600 s. All original geometry depths, quadrature and assertions
remain; separate budget avoids combined-suite hard timeout/evidence loss.

Default-off post-commit native Y reference export added for W2/final comparison:
ID-aligned CSV, 17-digit scientific values and SHA; includes fluid coefficients,
membrane DOFs, traction/force, material geometry/velocity, port measurements and
Newton residual/damping/KSP diagnostics. Per-step gates are truthful; no checkpoint.
Native finite/sorted/duplicate/negative ID and partial-failure manifest tests passed.
Both heartbeat binaries compiled in Slurm; utility prefix 46070816 only failed
because rg absent on compute node. Passed prefix hashes revalidated/reused,
candidate compile completed in **46071195** (55 s, 1 CPU, no simulation).
Workload W2Tools SHA256 698fe1ea19f3e7d53a0bcf4639f4ddcd162545eb7a7b74628a898142166b09df.

Real W2 runner staging-only preflight **46071756 passed**, 9 s, explicitly no
simulation. W2Runner SHA256 32b17336ad3dd025bea354c0f2df9100b68cfc6071d2bfcbcc7162c4dcb95a9f.
Original Y two-step matched regression **46071831 queued afterok:46070334:46071756**,
8 CPUs / 15200 MB / 2.5 h, per-case 4000 s; only one main calculation at a time.
Reuses verified binaries rather than recompiling under 8-core allocation. dt .05,
16x16x5, base pressure 5 + pulse amplitude 30/period 1, cut 1/rescue 3, inertial
gamma 1, formal mode (not visualization-only), native reference enabled both.
Checks original first-two history rows, full native fields, coupling, assembly
input identities, exact redundant-call reduction and no timing regression.
Whole-cycle shrink verifier intentionally not applied to still-expanding first 2.
Adaptive conversation sleep next 3 h; no live-log reads or background monitor.

P1-A **promoted after ALL gates passed**. Full original immersed suite 46070334
COMPLETED 0:0: 7234 s / 16.076 allocated core-hours, TotalCPU 42595 s,
MaxRSS 673380 KiB. Original Y W2 46071831 COMPLETED 0:0: 3234 s / 7.187
allocated core-hours, TotalCPU 9022 s, MaxRSS 500408 KiB. Exact original first-two
history, all native fields/manifest hashes, coupling and remaining assembly inputs;
formal first-two gates unchanged. Full calls 70 -> 60, elapsed 1682.416 ->
1541.453 s (8.379% reduction, ONE screening pair, not final repeated speed claim).
Whole-cycle formal_physics_acceptance remains false. Next independent card P1-B.

P1-B local accepted full-R/J reuse implemented as independent runtime delta;
first assembly and every original convergence check retained. Ready flag consumed
once, invalidated before candidate state changes, stack-local across return/throw/
rollback/retry. Work-only disable control does not enter physical input hash.
Testing-only candidate-norm and ready-boundary hooks compile out of production.
Eight off/on additional original depth-3 branch fixtures cover zero/one/multistep,
half damping, NaN first/all, backtracking failure and exception while ready;
every fixture checks fresh R/J action and first-assembly failure + exact retry.
Original depth-4 suite unchanged. **46083059 submitted**, 8 CPUs/15200 MB/4 h,
includes prior actual W2 input-pair proof and matched original small FSI pair.
P1B SHA a588354c3a9b53049097b2fbe9da544312a8f951287dd1c801cada8570779d34;
P1BBefore SHA 32b17336ad3dd025bea354c0f2df9100b68cfc6071d2bfcbcc7162c4dcb95a9f.
Not promoted; W2 still required. Adaptive conversation sleep next 3 h.

46083059 TERMINAL FAILED 124:0, targeted process budget 9000 s expired;
9049 s / 20.109 allocated core-hours, TotalCPU 53237 s, MaxRSS 661556 KiB.
NOT numerical failure or full-suite pass. Actual prior W2 proof passed: 20
identical accepted-line-search/next-Newton pairs, redundant full wall 315.454 s.
Seven complete assertion-bearing off/on pairs passed before timeout: zero,
one-step, normal, half, nan-first, failure, nan-all. Normal 7 Newton steps exact,
full calls 14 -> 8; half/NaN-first 15 -> 9. All exact fresh R/J actions and retries.
Preserve prefix and validate exact source/binary SHA + CLI-selector-only test
delta; do not rerun seven pairs or call timed-out whole process passed.

46096311 submitted ONLY missing throw-ready pair and matched small FSI, 8 CPU/2 h.
Production P1B unchanged. ValidationR1 SHA e43f2fda138e99c0d0107588eeacd2c47725d43afd1919083b83b427d8f49442.
46096441 1 CPU/20 min native-export/heartbeat build (no simulation).
46096555 FULL original unchanged depth4/depth3 immersed suite afterok:46096311,
8 CPU/4 h. 46096556 actual W2 staging-only preflight afterok:46096441, 1 CPU/10 min.
46096574 original Y W2 afterok:46096555:46096556, 8 CPU/2.5 h; all formal inputs,
historical first-two history, native fields, Newton/KSP and conservation exact;
before next-Newton proof + candidate local full readiness + real count decrease
and no timing regression. W2 runner SHA d32662dd92575b0c8cdcf1f35db5d907d05aa056a94f3f5b506fc7546dc332fe.
One primary calculation at a time; tiny no-simulation utilities independent.
All five handles confirmed PENDING. Adaptive conversation sleep next 4 h.
P1A accepted; P1B not promoted; full plan and final whole-cycle delivery pending.

46096311 **PASSED** in 2247 s: all eight targeted pairs and matched small FSI,
exact native fields/coupling/full-input reuse; full calls removed 4; one-screen
elapsed 370.210 -> 282.967 s. 46096441 native build/export passed (137 s, no sim),
46096556 real W2 staging preflight passed (7 s, no sim).

46096555 failed at exactly one obsolete original assertion: the old gate required
at least 2 full assemblies per Newton update. P1-B correctly has one initial full
assembly plus one full candidate per update, while preserving every convergence,
Newton/KSP and physical check. Updated only this pure-work lower bound; production
P1B is unchanged. ValidationR2 SHA b5f8449b483ab99bc5234a18f5b2dd4f07c399dd8c1f7f6bd939ad1d2cd0f73e.
46112307 submitted for the entire original suite; 46113038 original Y W2 afterok.
Old dependency-never 46096574 cancelled. Adaptive conversation sleep next 3 h.

P1-B **PROMOTED**. Full original immersed suite 46112307 passed (5730 s,
TotalCPU 33394 s, MaxRSS 667932 KiB). Original Y W2 46113038 passed (3294 s,
TotalCPU 8455 s, MaxRSS 494620 KiB): exact history/native fields/coupling/
Newton/KSP and full/reuse inputs; assemblies 60 -> 40. One-screen elapsed
1816.384 -> 1463.325 s (-19.437%); not a repeated final speed claim.

P2 default-off profiling now records worker completion spread and exact dynamic
result payload capacity; disabled unit path asserts no profile work. 46126551
submitted, 8 CPU/3 h, immutable binary and fixed 8 threads: batches 8/16/32 in
forward/reverse process orders, exact W1 + zero Mat mallocs, RSS<70%. Select the
smallest within 2% of fastest only with improvement in both processes, then run
all W0 serial/OpenMP/FUNNELED/SINGLE-negative/expanded/compact/ordered/fault-retry
gates. P2Validation SHA 9352ab9e2aa132e17a33dc3a8b79601685e7c8ed3add5e7e7cf38a67027d1ba0.
P1B remains accepted; P2 unpromoted and still requires W2. Sleep next 3 h.

P2 **PROMOTED** after build 46149676, preflight 46149881 and original Y W2
46149882 passed. Exact native/history/coupling and full/reuse sequence/call count;
batch8 1326.014 -> batch32 1163.091 s (-12.287%, one screening pair).

P3 true ResidualOnly implemented across volume/trace/wall/ghost/ports/gauge;
dense J work skipped and separate R/J generations guard every Newton solve.
Preliminary 46159278 exposed a test assumption that J remains current after
residual-only accepted convergence; corrected test fresh-builds J at the same
state and requires exact R/J. Production request stays internal. P3ValidationR3
SHA 5596d6009419bea4dc820a8894d8f621902b16541e94e7f7412720ae37e59c6c;
46161562 running eight residual/full branch comparisons. Full suite/W2 pending.

46126551 failed in 9 s before runner/build/simulation: P2 runner initially did
not parse common sbatch named flags. Fixed only harness argument parsing and
added frozen case manifest; numerical candidate unchanged. P2ValidationR2 SHA
e476464325c9ec258c004b2c80b431435c3badb3a4c68aec4c1a6ed9412e4b51.
46138318 submitted/running with the same P2 gate; sleep next 3 h.

P2 W0/W1 46138318 **PASSED** (868 s, TotalCPU 2914 s, MaxRSS 661716 KiB).
Exact W1 operators/zero Mat mallocs. Median batch8/16/32: 16.821/13.453/12.283 s;
batch32 reproducibly faster in both orders. Batch32 MaxRSS 265560 KiB, payload
50528256 bytes, worker spread 0.4395 s, prepare 0.0291%, consume 11.805%; all W0
serial/OpenMP/caller/FUNNELED/SINGLE-negative/expanded/compact/ordered/fault-retry
passed. 46149676 no-sim build -> 46149881 no-sim preflight -> 46149882 original Y
W2 queued. Baseline batch8 vs candidate batch32 saved per case; exact full/reuse
sequence and native results required. P2 unpromoted. Sleep next 3 h.

P3 targeted R4 46163904 **PASSED** in 1179 s. All zero/one-step/normal/half/
nan-first/failure/nan-all/throw-ready full-versus-residual branch pairs were
bitwise exact. The candidate records separate full/residual generations and
skips dense-J result payload in line search. P3VerifierR1 SHA
e59f327fb1983b22ee6602e2fe69e79f1714f3f1f48ea940d93fd6061d45fb30.
Complete original immersed suite 46168275 is running; no-simulation W2 build
46168602 is running. afterok matched small FSI 46168976, preflight 46168977, and
original-Y W2 46168978 are dependency-gated. P3 remains unpromoted pending all.

P3 no-simulation native/export build 46168602 PASSED in 107 s and staging
preflight 46168977 PASSED in 10 s. Both explicitly recorded
simulation_executed=false. Full original suite 46168275 continues; downstream
small FSI/W2 remain dependency-gated.

P3 static contract audit found R4's public assembled R/J accessors did not reject
storage from an older state generation. Newton was guarded, but public use still
violated the explicit current-R/J invariant. R4 full job 46168275 and pending
downstream 46168976/46168978 were cancelled before promotion. P3ValidationR5
SHA f083218030094ec8816bbb1401a606fa92451aff2c3d180ec38f074ee227620e
adds accessor guards plus stale-rejection/fresh-rebuild assertions. Corrected
targeted 46169775 -> full 46169825 -> small FSI 46169826 and no-sim build
46169774 -> preflight 46169827 -> W2 46169828 are dependency-gated.

R5 targeted 46169775 failed in 103 s only because its fixture captured J before
fresh rebuild; terminal stderr is exactly the new stale-generation guard. R5
dependent simulations were cancelled. R6 keeps the production guard and fixes
test semantics: same-state residual-only preserves a current J, while generation-
changing residual-only rejects stale J; targeted captures candidate J only after
full rebuild. P3ValidationR6 SHA
192352c89fcb33bfd115beec0912defada78df28f8b6f7a60994e287655cada7;
P3VerifierR3 SHA 245e4a5a06a04d17256359185ba76c43c5f4f374d9b72be588f22cabbd40bcc5.
R6 chain: 46170114 -> 46170193 -> 46170194 and 46170113 -> 46170195;
46170196 W2 depends on small FSI and preflight.

R6 targeted 46170114 **PASSED** in 1174 s. All zero/one-step/normal/half/
nan-first/failure/nan-all/throw-ready pairs passed bitwise, including fresh
operator/retry checks and the corrected R/J generation contract. No-simulation
build 46170113 (110 s) and staging preflight 46170195 (6 s) passed. Full original
suite 46170193 is waiting for resources with no scheduler start estimate; small
FSI 46170194 and original-Y W2 46170196 remain dependency-gated. Next adaptive
conversation sleep is 1800 s; only scheduler state is read until a job terminates.

P3 R6 full original immersed suite 46170193 **PASSED** in 3865 s, with the
unchanged mathematical regression gates. Matched complete small FSI 46170194
**PASSED** in 433 s: all native fields and coupling/conservation are bitwise
exact; full assemblies 12 -> 8 with 8 residual-only calls. One screening pair
elapsed 196.923 -> 151.988 s (~22.8%, not a final repeated performance claim).
Original-Y W2 46170196 is running; P3 remains unpromoted until its terminal
correctness and timing evidence pass.

P3 **PROMOTED**. Original-Y W2 46170196 passed in 2022 s. Original first-two
history, all native fields, coupling, formal short gates, state-aligned logical
assembly sequence, and R/J generation validity are exact. Full assemblies
40 -> 30 with 30 residual-only calls; one screening pair elapsed 1035.563 ->
965.625 s (-6.754%, not an independent repeated final claim). The accepted
baseline is now P3ValidationR6 with batch32. P4 begins by selecting only the
hottest remaining measured kernel; no cache type is assumed before profiling.

P4 measured full-assembly worker work at about 52.3 aggregate worker-seconds:
volume ~36.5, wall ~12.0, trace ~3.75; ghost was ~1.88 wall-seconds and ports
~0.04. The single selected candidate is therefore an immutable geometry-epoch
volume basis/gradient/Hessian/physical-point cache. Its key binds the complete
MovingCutGeometry identity and storage mode; state, coefficients and ports are
not cached. Workers only read entries built by the constructor/caller thread.
P4Validation SHA 28d5f7946d9c30f8e7d6c441c9c7aede8a8c21c07c9c2f49945cb7a37d03471d;
targeted Slurm Job 46184469 tests enabled/disabled bitwise R/Jv, same/different
epoch, coefficient/port changes, fault retry, and visible hits/misses/bytes.

P4 targeted 46184469 produced no durable numerical evidence: the submission
used a relative evidence root, so its final copy failed after the runner changed
directory. This is a submission-path failure, not a candidate result. The exact
same P4Validation archive was resubmitted as 46192353 with absolute root/card
paths; no source, numerical option, or gate changed.

P4 targeted R2 46192353 **PASSED** in 64 s. Enabled/disabled residual and Jv
are bitwise exact; same-epoch reuse, different-epoch invalidation, coefficient/
port safety, fault retry atomicity, and hits/misses/bytes gates passed. The depth2
cache held 25752 points / 176968392 bytes; MaxRSS was 552892 KiB. P4ValidationR2
SHA bf43dca26b8d77f5a73870e61153915be0eb36b33d0cc3fb2a6bff72c54b25cb
changes only W1 harness files from the targeted archive; both production file
hashes are exact. Matched W1 46201812 runs three independent process pairs in
off-on/on-off/off-on order; P4 remains unpromoted and W2 remains gated.

P4 matched W1 46201812 **PASSED** in 878 s. Every residual and Jv hash is
bitwise exact across cache-off/cache-on, and all three independent process
pairs favor cache-on: 12.860%, 12.751%, and 13.150% reductions (median
12.860%). The depth3 cache contains 180480 points / 1240259208 bytes. Frozen
P4Verifier SHA is
96fb6c015c4a2b4ebd21722650434acada9ce5a9c1e84ff582b75ad1b4668a49.
The P4 W2 chain is no-simulation build 46203125, afterok preflight 46203289,
then full original-Y W2 46203291. P4 remains unpromoted pending its terminal
exactness, cache-invariant, RSS, and no-regression timing gates.

P4 **PROMOTED** after build 46203125, preflight 46203289, and original-Y W2
46203291 passed. W2 preserved the original first-two history, all native fields,
coupling, state-aligned assembly sequence, cache epoch/key/counter invariants,
and RSS bound. Elapsed was 1157.074 -> 1157.644 s (-0.049%); this is within the
3% no-regression gate and remains a single screening pair, not a final speed
claim. The three independent matched W1 pairs all favored cache-on by
12.751--13.150%, so P4ValidationR2_batch32 is the accepted baseline. P5 starts
from its measured 219.665 s solver setup versus 3.189 s linear solve in W2.

P5 backend profile 46211622 **PASSED**. One original-Y step used five runtime
instances and 15 Newton solves. The actual configuration is FGMRES with a right
PCLU using PETSc's built-in factor package. MatLUFactorSym ran 5 times
(5.535 s) and MatGetOrdering 5 times (0.239 s), while MatLUFactorNum ran all
15 times (86.189 s) and dominated PCSetUp (91.976 s). Thus same-pattern
symbolic/ordering reuse is already active; enabling its flags cannot remove the
measured numeric factorizations. P5Validation instead tests same-runtime old LU
strictly as a preconditioner for the current Jacobian, with iteration and
explicit true-residual gates plus one current-Jacobian rebuild fallback.

P5 targeted R2 46214496 **PASSED**. The disabled path built seven current LUs.
The candidate made six reuse attempts, accepted one, rejected/rebuilt five, and
ended within 1.18e-18 relative state difference. A forced zero-iteration reuse
limit rebuilt all six stale attempts and recovered disabled state/ports
bitwise. The matched W1 uses three independent original-Y one-step AB/BA/AB
process pairs, native physical fields with a frozen per-quantity rtol of 1e-8,
and requires every setup/elapsed pair to improve with median elapsed reduction
at least 5%.

P5 matched W1 46215911 **PASSED** all 94 checks. Independent pair reductions
were 11.450%, 10.593%, and 10.884% (median 10.884%); solver setup fell from
97--101 s to 35.6--38.2 s. Each candidate one-step process built five current
LU factors and accepted ten same-runtime reuse attempts without fallback.
Baseline repeats were bitwise deterministic; candidate native physical fields,
history/coupling quantities, conservation and explicit true residuals passed
the frozen per-quantity rtol 1e-8. P5 remains unpromoted pending two-step W2.

P5 W2 46219623 completed both simulations: 1061.255 -> 937.780 s (-11.635%),
setup 206.206 -> 74.199 s, with 10 current-LU builds, 20 accepted reuses and no
rebuilds. All native physical fields at both steps, coupling, formal gates,
logical path, true residual and reuse accounting passed. The initial verifier
failed only its aggregate history check because it applied relative tolerance
to near-zero dimensionless continuity (difference 6.00e-15 versus an invalid
1.00e-16 limit). The corrected check is roundoff-only:
256*epsilon*sqrt(rows) at natural dimensionless scale 1; it does not use the
candidate magnitude. Initial failure evidence is retained, and P5 will not be
promoted unless both verification-only and an independent full W2 rerun pass
the frozen corrected verifier.

The first verification-only rerun, 46224110, exposed a separate harness-only
error: it compared the saved 8-thread assembly records with the re-verifier's
1-thread `OMP_NUM_THREADS`. No simulation was executed, and dependent full W2
46224111 was cancelled by Slurm. P5VerifierR3 instead uses the thread count in
the frozen case as the evidence invariant and records the verifier process's
thread count separately. Its source archive SHA is
24db08eeac25334f72f8b33dc6cbc9aed72e1fa23f94de2408577f9899d267c5.
Replacement verification-only Job 46225784 and afterok-dependent independent
8-thread W2 Job 46225785 are submitted; P5 remains unpromoted.

Replacement verification-only Job 46225784 **PASSED** in 10 s without executing
a simulation. Every saved-evidence gate passed; the manifest records 8 evidence
threads and 1 verification thread. The dependency released independent W2 Job
46225785, now running on r211. P5 remains unpromoted until that fresh run passes.

Clarification for W1: the 94 gates directly compared the candidate's native
physical fields at frozen per-quantity rtol 1e-8, formal history/conservation
conditions, explicit true residual, and reuse accounting. They did not perform
a separate candidate-versus-baseline numeric comparison of every history and
coupling column; that broader comparison belongs to W2.

P5 **PROMOTED** after independent original-Y W2 Job 46225785 passed all 214
checks. The fresh pair improved elapsed time from 1149.077 to 997.231 s
(13.215%) and solver setup from 239.788 to 89.881 s. The candidate built 10
current LU factors, accepted all 20 same-runtime reuse attempts, and required no
fallback rebuild. This independently repeats Job 46219623's 11.635% reduction,
while the three W1 process pairs improved by 10.593--11.450%. The accepted
baseline is now P5ValidationR2_batch32. Formal full-cycle physics acceptance
remains false; the next step is conditional P6 selection from measured evidence.

P6 is justified by the accepted-P5 profile: geometry consumed 352.675 s, or
35.37% of elapsed time. Across ten moving trials, geometry construction consumed
229.761 s, layout/runtime/preallocation 52.272 s, and velocity plus scalar
extension 48.996 s. The candidate spatial-index construction itself consumed
only 0.009 s, so caching it alone is rejected by measurement. Profiling-only
P6Profile splits fixed-runtime construction without changing operation order;
source SHA 0619efb045869debda2ffecd0e328ebe5550e2432081f4be0aa1e10e52f902f6,
Slurm Job 46230573. No P6 optimization is promoted yet.

P6 profile Job 46230573 **PASSED** in 503 s. Six runtime builds had one common
active-layout pattern, but layout construction itself took only 0.013 s.
Jacobian preallocation took 12.172 s total; after retaining the unavoidable
initial epoch, even eliminating every later preallocation would save only about
10.1 s of a 452.8 s one-step run (2.24%). That is below the plan's 5% timing
uncertainty threshold, while a cross-epoch pattern cache would add exact-pattern,
active-set transition, and fault-cleanup contracts. Geometry-dependent volume
basis cache (13.669 s), gauge weights (3.748 s), quadrature, and extension are
not legally reusable. P6 is therefore skipped after measurement, and its
profiling-only production instrumentation was removed byte-exactly back to
P5ValidationR2. P7 is also skipped: existing Aitken coupling averages about
four iterations per step and, at 11.27% of P5 W2, does not dominate geometry or
assembly. Work proceeds to R6 1/2/4/8-thread production configuration selection.

R6Scaling is frozen at SHA
5506098e7efa796b243382a6855523d0d174bc02703952b305adaebd199ad522;
all 153 production include/src/Makefile files are byte-identical to
P5ValidationR2. W1 Job 46232121 compares 1/2/4/8 threads in forward/reverse
process order with batch32, exact two-state operator hashes, zero matrix mallocs,
and recorded team/RSS/payload evidence. The fastest two configurations will
advance to W2; no thread selection is assumed before this job passes.

R6 W1 Job 46232121 **PASSED**. Paired process medians for 1/2/4/8 threads
were 53.267, 28.817, 16.741, and 11.306 s. The 1-to-8 speedup is 4.711 with
58.9% parallel efficiency; every two-state R/Jv hash is exact and every matrix
malloc count is zero. W2 top two are therefore 8 and 4 threads. Frozen
R6Verifier SHA is
565e1db4a2dc7f0153846fec022f7329b530e005c004a6e66741111406206022;
production remains byte-identical to accepted P5. Original-Y top-two W2 Job
46235778 is submitted; production thread selection remains pending its result.

R6 W2 Job 46235778 ran both simulations to exit zero: 8 threads took 14:50.58
and 4 threads 17:46.25. Verification then raised `FileNotFoundError` before
reading numerical evidence because the runner saved `cases-r6-w2.json` instead
of the verifier's fixed `case.json`. The simulations are retained. Filename-only
R6VerifierR2 SHA is
24c1aea54a6bc1beda0f515b54c21bce03e3bd4df310a920d220ba3dee6e4d75;
verification-only Job 46238720 is submitted, with no FSI recomputation.

R6 production configuration **PASSED**. Verification-only Job 46238720 found a
harness-only exact-one-profile assumption; the driver correctly writes a step-1
and a final cumulative step-2 profile. Frozen R6VerifierR3 SHA
4b5ab21428f78197777abfec9bb5b9ece034cfcf7320ff8e0c7ecf69ff6da5ab
selects the final cumulative record. Job 46239782 passed all saved gates without
executing a simulation: both thread configurations preserve history, coupling,
all native fields, assembly/linear paths, cache and preconditioner contracts.
Eight threads took 890.242 s versus 1065.908 s at four threads (16.48% faster),
so the production configuration is 8 threads plus batch32.

R7Final SHA
0aeb8a2d12cc09035514b27c838defd84c55916a52d2b8e6e693f2d9b96d4993
contains the final 20-step runner/verifier and is byte-identical to accepted P5
across all 224 compared production files. Slurm Job 46243889 is submitted with
8 CPUs, 15200 MB, batch32 and 5 hours. It runs the solver once, then CSV and
expansion/contraction checks, native VTK 9.4.1 reads, and a verified ParaView
archive/hash. The 13603 s historical result remains explicitly a cross-job,
cross-node reference rather than a paired full-cycle speedup claim. Initial
state is resource-pending; next action is one 10800 s conversation wait and one
named-job scheduler query.

## Final R7 result and handoff

R7 Job 46243889 **PASSED** all 19 checks and completed 0:0 on r030. The solver
ran all 20 steps in 6968.953 s. History and coupling remain within the frozen
per-quantity tolerance, every recorded value is finite, expansion and
contraction are observed, and native VTK 9.4.1 reads all 20 frames with frame
times and point topology exactly matching the historical delivery. The volume
basis cache, residual-only/full paths, stale-LU preconditioner reuse, explicit
true linear residual, and archive checksum contracts all pass.

| Evidence | Result |
|---|---:|
| Queue wait | 21:42 |
| Slurm elapsed | 1:58:34 |
| Solver elapsed | 6968.953 s |
| Total CPU | 4:55:39 |
| Allocated core-hours | 15.809 |
| MaxRSS | 1134052 KiB |
| Production configuration | 8 threads, batch32, 1 MPI rank |
| ParaView archive size | 93174492 bytes |

The historical solver elapsed was 13603 s, so the observed final elapsed is
48.769% lower. This is a cross-job/cross-node historical comparison, not a
paired full-cycle timing claim. The independently repeated P5 W2 screening pair
showed 1149.077 -> 997.231 s (13.215%), while R6 measured 8 threads at
890.242 s versus 1065.908 s for 4 threads (16.481% faster) on identical
two-step output. These measurements are separate experiments and are not
multiplied together.

The final exclusive profile remains dominated by assembly (2974.686 s) and
geometry (2463.835 s), followed by coupling-exclusive work (864.772 s) and
solver setup (577.579 s). P6 did not introduce a cross-epoch pattern cache
because its measured ideal ceiling was only 2.24%, below the plan's 5%
uncertainty threshold and not worth the lifecycle risk. P7 was skipped because
the existing Aitken iteration count and coupling share did not justify a new
coupling algorithm.

The final run deliberately retains the historical model's 17 conservation
failures out of 20 steps. Therefore this is a verified visualization delivery,
not formal physics acceptance; the optimization introduced no new failure and
did not hide the pre-existing limitation.

The download artifact is
`benchmarks/fsi-time-opt/20260915T035654Z/job-46243889/heartbeat-paraview.tar.gz`.
Its verified SHA256 is
`5936dd2905e2dbf1dcf5a0bc3314bf1ef3180b44ed3a7e64985e14932613dce5`.
The immutable final source SHA256 is
`0aeb8a2d12cc09035514b27c838defd84c55916a52d2b8e6e693f2d9b96d4993`.
No Slurm job remains live. Existing uncommitted user changes were preserved;
no reset, clean, commit, or push was performed.
