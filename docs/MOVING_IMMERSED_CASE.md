# Prescribed moving immersed flow cases

The MPI case owner accepts `prescribed_motion` in `immersed_geometry.json` for
backward-Euler flow. The ordinary fixed-geometry configuration remains valid.
Full moving FSI and moving checkpoint/repartition restart are outside this case
owner's supported scope.

Keep the Cartesian grid fixed. Every VTP frame must preserve the reference
file's vertex order, directed triangle connectivity, and boundary labels.
Coordinates are in metres; scaling and welding are rejected by the material
reader. Each surface must be closed and lie within the configured grid bounds.
Frame files must be regular files contained within the case directory.

Add this object alongside the existing `surface`, grid, quadrature, wall-label,
and runtime settings:

```json
"prescribed_motion": {
  "frames": [
    {"time_s": 0, "surface": "surface.vtp"},
    {"time_s": 0.25, "surface": "motion.vtp"}
  ],
  "extension_layers": 2,
  "conservation_limits": {
    "divergence_theorem": 0.01,
    "reynolds": 1e-8,
    "moving_mass": 0.01,
    "wall_relative_leakage": 0.01,
    "discrete_continuity": 1e-8
  },
  "volume_fitting": {"support_expansion": 3}
}
```

This example uses two solver steps of 0.125 seconds. The first frame must have
time zero and name the configured reference `surface` file. Supply 2–4096
frames with strictly increasing times covering the entire simulation. Interior
frame boundaries must coincide with configured timesteps; a solver step cannot
cross a change in the prescribed piecewise-linear velocity. Boundary comparisons
allow bounded floating-point roundoff: interpolation queries may snap to a
nearby frame, while input frame times and published clocks retain their original
values. The factory derives the allowance from the configured number of clock
additions, with a small endpoint-comparison margin. It is capped by neighbouring
frame durations so that close
physical knots cannot merge. The extension
runtime checks the actual step displacement against its configured band.

All five conservation limits are required, finite, and nonnegative. The values
above belong to the translating unit-cube regression; choose the runtime's
`flow_controller_reference_flow_m3_s` and conservation limits from the physical
scale and accuracy requirements of the actual case. Failed physical gates
reject the trial before publishing accepted fields or graph measurements.

`volume_fitting` is optional. It can also specify `candidate_orders`,
`max_columns`, `max_workspace_bytes`, and `max_point_queries`. The positive-weight
and complete degree-six moment audits remain mandatory. Increasing candidate
support does not automatically increase these resource limits.

Use the ordinary graph runner and matching domain/graph `dt` and `steps`:

```bash
make -C solvers/coupling iga_multidomain_flow PETSC_DIR=/path/to/petsc
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 IGA_PROFILE=1 \
  mpiexec -n 2 solvers/coupling/iga_multidomain_flow \
  --graph-case /path/to/graph-case --output-dir /path/to/new-results
```

The runner emits `hpc_immersed_moving_step` diagnostics when profiling is enabled,
including accepted clock/index, assembly and linear-solve times, and normalized
Reynolds, moving-mass, and wall-relative-leakage defects. Numerical fields use
owned MPI vectors and required-state communication. Geometry and prescribed
frames are currently replicated; this is not a distributed geometry-storage
implementation.

To reproduce the translating-cube graph used for local validation, generate a
new case directory from the repository's existing 1D/immersed chain:

```bash
python3 scripts/hpc_make_moving_graph_case.py --output-dir outputs/moving-demo
mkdir outputs/moving-demo/mpi-2
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 IGA_PROFILE=1 \
  mpiexec -n 2 python3 scripts/hpc_rank_run.py \
  --output-dir outputs/moving-demo/mpi-2 --expected-ranks 2 --timeout 14400 -- \
  solvers/coupling/iga_multidomain_flow \
  --graph-case outputs/moving-demo/fixture --output-dir outputs/moving-demo/mpi-2/result \
  -domain_immersed_flow_ksp_type gmres \
  -domain_immersed_flow_pc_factor_mat_solver_type mumps \
  -domain_immersed_flow_mat_mumps_icntl_14 100
```

The generator only writes input files and their SHA-256 manifest; it does not
build or run a solver. Existing output directories are rejected. The default
is four steps of 0.0625 s, ending at 0.25 s, with density 1, viscosity 0.1 and
wall inertial penalty 1. Its files reproduce the validated two-rank input
byte-for-byte. The one- and four-rank full graph comparisons remain pending.
`--steps 2` reproduces the documented coarse-step probe that fails Newton
backtracking on its second step; it is provided for reproducing that failure,
not as an accepted baseline.

Validate the saved rank reports and accepted history after the run finishes:

```bash
python3 scripts/hpc_check_moving_graph.py \
  --case-dir outputs/moving-demo/fixture --run outputs/moving-demo/mpi-2
```

The checker verifies every rank's successful exit and log hashes, ownership
coverage, accepted clocks, the three logged moving conservation diagnostics,
accepted edge residual gates and complete port coverage. It emits JSON and
returns nonzero on missing, failed or inconsistent results. Repeat `--run` for
additional runs and pass `--reference-run /path/to/mpi-1` to compare accepted
ports using the existing graph tolerance `1e-12 + 1e-6 * abs(reference)`.
Without a reference, it performs no cross-rank numerical comparison. Full-field
norms, source provenance, the physical gates absent from CLI diagnostics, and
performance acceptance remain separate checks in the progress reports.
