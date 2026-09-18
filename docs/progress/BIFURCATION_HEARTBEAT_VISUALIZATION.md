# Heartbeat FSI visualization

User requested a functional pulsating-vessel ParaView demonstration and explicitly
allowed conservation-gate violations for this visualization. Formal acceptance
is unchanged; the new mode always records `visualization_only`.

Job 45971412: 8 CPU cores, 15200 MB, RM-shared, 4-hour limit. Builds from an
archived source copy and runs in Slurm. Results:
`.nsvms_diagnostics/heartbeat/job-45971412/results/`.

Configuration: 20 steps × 0.05 s, 1 s cosine pressure cycle, inlet 5–35 Pa,
outlets 0 Pa, grid 16×16×5, cut depth 1, Nitsche inertial gamma 1. Membrane
material and coupling convergence thresholds retain their existing values.
Actual-scale `fsi.pvd` and 100× displacement `wall-display.pvd` are exported
incrementally. Only the latter changes displayed coordinates; physical arrays
retain their original values. No synthetic temporal frames are inserted.

Status: submitted and build passed; full cycle and independent expansion /
contraction verification pending. See the example README for reproduction and
ParaView instructions. The original constant-pressure job 45961543 is separate.

Initial job 45971160 failed before stepping: grid 6×6×4 had no required ghost
coverage for positive cut cell 12. Logs are retained. The rerun restores the
original grid16 while keeping cut depth 1 for a faster visualization.

Job 45971317 failed before stepping because depth-1 quadrature left an unresolved
empty cut-cell rule. The next revision enables the existing bounded empty-rule
rescue through depth 3 only for low-depth demonstration quadrature. It does not
disable geometry checks.

Current job 45971412 includes bounded empty-rule rescue through depth 3. The
worktree subsequently splits one same-line statement to remove a compiler
indentation warning and adds finite checks for CLI pressure/time step.

Native VTK verification and packaging are queued as job **45972918**, depending
on successful completion of 45971412. The generated archive will be
`.nsvms_diagnostics/heartbeat/job-45971412/heartbeat-paraview.tar.gz`.
The verifier checks every real fluid/membrane/display frame with VTK 9.4.1 from
the existing research Conda environment. Verification status remains specific to
visualization and does not close the formal conservation acceptance case.

The running snapshot's verifier was supplemented before execution with coupling
history and membrane backward-Euler velocity checks; `verify_heartbeat_vtk.py`
was copied into its scripts directory for the dependent postprocessing job.
The current reproduction wrapper additionally enables IGA_PROFILE for iteration
logging. Job 45971412 itself logs at completed-step boundaries.

At 7/20 completed steps (0.35 s), job 45971412 had used 2h49m of its 4h
allocation. Slurm refused an extension to 12h; evidence is in
`extend-time.stderr`. Recovery job 45980201 depends on afterany:45971412 and
has a 12h allocation. If the original run completed, it only verifies/renders/
packs it; otherwise it reruns with the exact archived binary, surface and
parameters. It includes native VTK verification and packaging. Pending standalone
postprocessing job 45972918 was cancelled to prevent duplicate postprocessing.

Latest: original 45971412 timed out after 4h. Old queued retry 45980201 was
cancelled and replaced with **45982619** (12h, same physical case, tested
parallel wall assembly). New outputs are in
`benchmarks/heartbeat-optimized-45982619/results/`. Numerical regression and
performance evidence are in BIFURCATION_FSI_ASSEMBLY_PERFORMANCE.md.
