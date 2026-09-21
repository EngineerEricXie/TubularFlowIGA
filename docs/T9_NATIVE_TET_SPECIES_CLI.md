# Native tetra prescribed-velocity species CLI

`solvers/cpu/native_tet_species_transport` is a user-facing, project-owned
P1 tetrahedral backward-Euler species solver. It reads an external labelled
ASCII Gmsh 4.1 mesh, assembles its own ALE weak form and face fluxes, and uses
PETSc only for sparse algebra and MPI. fTetWild or Gmsh may create the input
tetrahedra; neither supplies the FEM discretization or solver.

This route prescribes a spatially constant fluid velocity and a rigid mesh
translation. It does **not** solve the fluid equations, connect 0D/1D graph
ports, provide multi-species chemistry or a 3D tissue field, or establish
physiological accuracy. The optional reaction is linear decay. Wall exchange
may use either a prescribed external concentration or one solved, fixed-volume
0D tissue reservoir, or multiple label-specific 0D reservoirs.

For surface or volume inputs, the shared workflow entry point now accepts
`backend=native_tet_p1_prescribed_species`. It runs the same native input
preflight and can generate a labelled tetra mesh with fTetWild or Gmsh:

```json
{
  "schema_version": 1,
  "backend": "native_tet_p1_prescribed_species",
  "input_route": "volume",
  "input_file": "labelled_pipe.msh",
  "case_file": "species.json",
  "output_directory": "species_workflow_001",
  "mpi_ranks": 2
}
```

Run `python3 scripts/run_native_tet_workflow.py workflow.json`. For a surface,
change `input_route` to `surface`, supply a labelled VTP/STL `input_file`, and
add the fTetWild or Gmsh `mesher` object documented in
[the native tetra workflow](T9_NATIVE_TET_WORKFLOW.md). The workflow retains
the derived `volume.msh`, effective `case.json`, mesh manifest when applicable,
preflight and simulation logs, SHA-256 provenance, stage timings, and
`workflow_status.json`. Its `fields/` directory contains this CLI's P1 VTK,
checkpoint, state and final summary. `--stop-after-step N` and `--resume` use
the species checkpoint without remeshing; altered inputs, ranks, solver or
incomplete publication are rejected. The status distinguishes functional-only
numerical acceptance, unestablished physical acceptance, and visualization
availability. This is an additional *prescribed-velocity species* backend,
not coupled fluid/species execution.

Build and use a version-1 case:

```bash
make -C solvers/cpu native_tet_species_transport PETSC_DIR=/path/to/petsc
mpiexec -np 2 solvers/cpu/native_tet_species_transport species.json --check-input
mpiexec -np 2 solvers/cpu/native_tet_species_transport species.json --output-dir species_run_001
```

To pause after an accepted step and continue in a new process, use the exact
same case, mesh, rank count, and output directory:

```bash
mpiexec -np 2 solvers/cpu/native_tet_species_transport species.json \
  --output-dir species_run_001 --stop-after-step 1
mpiexec -np 2 solvers/cpu/native_tet_species_transport species.json \
  --resume species_run_001
```

```json
{
  "schema_version": 1,
  "mesh_file": "labelled_pipe.msh",
  "species_id": "tracer",
  "initial_concentration_mol_m3": 0.0,
  "diffusivity_m2_s": 0.00001,
  "source_mol_m3_s": 0.0,
  "first_order_decay_rate_s_inv": 0.0,
  "fluid_velocity_m_s": [0.01, 0.0, 0.0],
  "mesh_velocity_m_s": [0.0, 0.0, 0.0],
  "inflow_concentration_by_label_mol_m3": {"1": 1.0},
  "wall_exchange_by_label": {
    "0": {"transfer_coefficient_m_s": 0.0001,
          "external_concentration_mol_m3": 0.0}
  },
  "monotone": true,
  "time": {"dt_s": 0.01, "steps": 2}
}
```

`mesh_file` is relative to the case JSON. All values are SI. Inflow labels
must exist in the mesh, and **every** face with inward relative velocity
`u−w` needs a supplied concentration; other labels may have natural outward
transport. A positive constant body source adds mol/(m³·s). Set `monotone`
explicitly: the optional low-order graph diffusion helps a sharp front stay
nonnegative but can obscure small physical diffusion on a coarse mesh.
Optional `first_order_decay_rate_s_inv` (default zero) removes species at
`k c` mol/(m³·s). Its consistent P1 reaction mass matrix is implicit; each
step reports `reaction_sink_mol_s = k × current_inventory_mol`, and the
balance includes that sink. With uniform concentration and zero relative flow,
the backward-Euler reference is
`c_next=(c_previous+dt×source)/(1+dt×k)`. This linear decay is not a
multi-species chemistry or physiological metabolism model.
Optional `wall_exchange_by_label` prescribes a linear exchange on labelled
zero-relative-flow faces: outward flux density is
`h × (c − c_external)` in mol/(m²·s), with `h` in m/s and both concentrations
in mol/m³. Positive flux removes species from the FEM volume; negative flux
supplies it from the prescribed external reservoir. The native P1 face mass
matrix and donor RHS are integrated on the **current** mesh, and the total
outward wall flux enters the accepted-step mass balance separately from
advective flux and body reaction. Missing labels, negative parameters and
nonzero relative flow on an exchange label fail input preflight. This option
does not update any external reservoir. Alternatively, omit
`wall_exchange_by_label` and specify one finite reservoir:

```json
"finite_wall_reservoir": {
  "boundary_label": 0,
  "transfer_coefficient_m_s": 0.0001,
  "volume_m3": 0.000001,
  "initial_concentration_mol_m3": 0.0
}
```

For multiple regions use `finite_wall_reservoirs_by_label` instead:

```json
"finite_wall_reservoirs_by_label": {
  "0": {"transfer_coefficient_m_s": 0.0001,
        "volume_m3": 0.000001, "initial_concentration_mol_m3": 0.0},
  "2": {"transfer_coefficient_m_s": 0.0002,
        "volume_m3": 0.000002, "initial_concentration_mol_m3": 1.0}
}
```

Both finite options use the same native P1 wall term, but solve each reservoir's
backward-Euler amount `M=V C` together with the vessel step. Positive outward
wall flux is added to the reservoir; the calculation checks each region and
the combined species budget, while the CLI reports final amounts and the
combined defect. All three wall options are mutually exclusive.
The finite options accept one species, and each specified label must have
its own boundary faces; the boundary
must have zero relative through-flow. This is storage/exchange, not a 3D
tissue perfusion model; see the [coupled numerical contract](T6_NATIVE_WALL_RESERVOIR_EXCHANGE.md).
On a coarse zero-diffusion mesh, ordinary Galerkin transport can conserve
total species while overshooting local concentrations; use the explicit
`monotone` mode when a bounded front is required, without assuming it proves
accuracy on arbitrary geometry.
The mesh velocity translates every vertex by `mesh_velocity_m_s × dt_s` each
step. The native assembler checks that displacement and mesh velocity agree.

`--check-input` uses the same parser, native mesh reader, Jacobian checks and
first-step species assembly without creating output. A run requires a new
output directory and will not overwrite an existing one. Each accepted step
publishes an immutable rank-owned linear-tetra `step_N/snapshot.pvtu` with
`concentration_mol_m3`, `reference_position_m` and `displacement_m`; VTK
`Points` are current coordinates. `run_summary.json` records case/mesh
SHA-256, ranks, accepted steps, final time, inventory, reaction sink, outward
wall exchange and balance defect.
With `finite_wall_reservoir`, it also records the final tissue amount and
concentration and combined budget defect. With
`finite_wall_reservoirs_by_label`, `tissue_regions_by_label` records each
label's volume, final amount and concentration, plus the combined defect.
Each accepted step also publishes a checksummed native concentration/ALE
checkpoint and, for finite reservoirs, `tissue_step_N.bin`. The original
single-reservoir sidecar retains its v1 format; multi-region sidecars use a
checksummed v2 format with all sorted labels and amounts. `run_state.json`
is updated only after the VTU and both applicable checkpoints exist. Resume
checks exact case/mesh hashes, MPI count, clock,
current ALE translation, both latest checkpoint identities, and all previous
snapshot/checkpoint paths including every rank piece. It rejects changed inputs, completed runs, missing
or corrupted files, and orphan next-step output. An OS-held output lock
rejects concurrent writers. A crash between snapshot/checkpoint/manifest
publication may leave an orphan; the CLI rejects that state rather than
silently repairing it. `run_summary.json` appears only after all requested
steps finish. The finite-reservoir mode has a paired vessel/tissue restart;
neither mode is a coupled fluid/0D graph checkpoint.

Local functional regression:

```bash
make t9-native-tet-species-cli-test PETSC_DIR=/path/to/petsc
make t9-native-tet-species-cli-test PETSC_DIR=/path/to/petsc \
  FTETWILD_BIN=/path/to/FloatTetwild_bin
make t9-native-tet-species-workflow-test PETSC_DIR=/path/to/petsc \
  FTETWILD_BIN=/path/to/FloatTetwild_bin
```

The test covers fixed and rigid ALE fields, nonzero source, a monotone inflow
front, implicit source/decay, labelled prescribed wall exchange, one- and
two-region finite exchange with paired restart and rejected corrupt/missing sidecars,
resolvable diffusion, missing donors/labels, fail-closed case keys,
no-overwrite behavior, 1/2/4 MPI ranks, source and moving-ALE cross-process
restart field parity within `1e-10 mol/m³` absolute difference, bad restart
rejection, and optionally a newly generated fTetWild pipe with two labelled
finite reservoirs exchanging in opposite directions. On the four-tetra
fixture `D=0.05 m²/s` did not change the
monotone front beyond floating-point noise because low-order graph diffusion
dominated; `D=1 m²/s` changed the field by `0.110226 mol/m³`. This is a
functional material-sensitivity observation, not a calibrated diffusion study.
