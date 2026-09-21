# Shared solver entry

Run `python3 scripts/solver.py CONFIG.json`. Use `--dry-run` to print the
selected backend command. The entry reads a `solver-v1` outer JSON, or accepts
an existing `simulation_config.json`, native tetra case JSON, or native tetra
workflow JSON directly. Direct IGA imports also need `--database DB.ntiga`;
direct native tetra cases need `--backend native_tet_hydraulic_graph`,
`native_tet_species_transport`, `native_tet_darcy`,
`native_tet_flow_cuda`, `native_tet_darcy_cuda`, or
`native_tet_species_cuda`. A direct
`native_tet_flow_darcy_chain` case is recognized from its three `meshes`;
`native_tet_fsi_steps` is recognized from its fluid and solid mesh paths.

`examples/solver/` contains editable cases. `fem_source_vessel_rcr.json`
shows the supported named-port connection: one 0D source, one P2/P1 tetra
vessel, and one RCR terminal. Add another RCR domain, a labelled vessel
outlet port, and a connection to use multiple outlets. The entry writes a
native hydraulic graph case, runs its input check and solver, and writes
`series.pvd` for ParaView. The native graph uses sequential pressure and
flow-rate exchange. `--stop-after-step N` and `--resume` use the graph's
checkpoint directory; the same options pass through to a native tetra
workflow. The solver still rejects incompatible restart identities.

`one_d_rcr_terminal.json` connects a compliant 1D network to two named 0D
RCR terminals. A 1D outlet port names its centerline `node_id`; each terminal
has an inlet port and `materials` with resistance, capacitance and pressure
values. The entry writes a runtime copy of the referenced 1D case, replaces
the selected outlet boundary conditions with `resistance`, `windkessel_rc`,
or `windkessel_rcr` according to each terminal's `solver: r|rc|rcr`, then runs
the existing 1D solver. The original case and its geometry stay intact.
Connections use `exchange: ["pressure", "flow_rate"]` and
`scheme: "sequential"`; outer `time`, `waveforms` and `output.every_steps`
apply to the generated case. This route supports 1D→0D terminal coupling;
the 1D solver advances each terminal state with the network flow.

`one_d_fem_rcr.json` connects a rigid 1D upstream vessel to the P2/P1 tetra
vessel and a 0D RCR outlet. The upstream domain references an existing 1D
case directory and gives its outlet port a centerline `node_id`. The graph
runner advances the 1D network, tetra flow, and RCR state together using
the pressure/flow port adapters. Its 1D case must use the same `dt_s` as the
outer graph; this route currently uses a 1D steady or explicit flow scheme.
It writes the tetra field series and coupled flow ledger.
`one_d_fem_rcr_species.json` uses the staged T7 flow/transport executor for
the complete 1D→FEM→0D chain. `species_bindings` maps the logical species to
the transported 1D field, the FEM P1 solve consumes the 1D outlet donor, and
the terminal RCR stores the arriving amount in its well-mixed
`species_volume_m3`. The terminal also sets
`species_initial_concentration_mol_m3`. Flow and species amounts are checked
at both graph edges and across all three domains before a step is committed.

`fem_source_vessel_rcr_species.json` adds two independent P1 species to the
same moving FEM vessel. At each accepted step, the P2 fluid velocity and ALE
mesh velocity drive both species solves. Each species has its own initial
concentration, diffusivity, source, optional decay and wall exchange, and a
concentration donor for every boundary label where flow can enter. The solver
writes `fields/species_ID.pvd` alongside the flow `series.pvd`. Coupled
flow/species runs do not yet have a paired checkpoint, so `--resume` is
unavailable for this case.

`fem_flow_darcy_wall_species.json` imports the existing idealized dual-tree
reference through the same command. It generates the labelled artery, fixed
Darcy tissue, venous and wall tetra meshes, runs the conservative sequential
artery→tissue→vein hydraulic transfer with a quasi-steady wall response, then
uses the frozen resulting fields for transient P1 species propagation. The
`fem_moving_wall.json` example runs just the quasi-steady flow, Darcy and
wall response on the smaller reference geometry. The
referenced flow and species case files must match. Results are under
`flow/` and `species/` in the configured output directory. This reference has
four matched arterial and venous terminal caps; it is a single functional
case, not a general vessel–Darcy graph.

For user supplied labelled tetra meshes, a three-domain graph can instead
declare `role: artery`, `role: tissue`, and `role: vein`. Set their `physics`
to `["flow"]`, `["darcy"]`, and `["flow"]`, and connect any number of artery
exchange ports to tissue inlet ports and tissue outlet ports to vein inlet
ports with `exchange: ["flow_rate"]` and `scheme: "sequential"`. The artery
domain supplies `materials.density_kg_m3`,
`materials.dynamic_viscosity_pa_s`, `boundaries.wall_label`,
`boundaries.inlet_velocity_m_s`, and an inlet port with a surface `label`.
The tissue domain supplies `materials.darcy_mobility_m2_pa_s`; the vein
domain supplies `boundaries.wall_label`, `boundaries.outlet_pressure_pa`, and
an outlet port. The solver matches artery and tissue facets geometrically,
maps arterial flux into conservative Darcy cell sources, and applies each
Darcy outlet flux to its named vein inlet. It writes `artery.pvd`,
`tissue.pvd`, and `vein.pvd` with corresponding `.pvtu` snapshots.
The fluid and Darcy fields are steady. Add `species` to the physics lists of
all three domains, a `species` array to the artery domain, and outer
`time.dt_s`/`time.steps` to propagate any number of species through those
solved fields. Each species defines `id`, `initial_concentration_mol_m3`,
`diffusivity_m2_s`, and `inlet_concentration_mol_m3`. The route writes one
species `.pvd` series per region. Wall motion is not coupled in this route.
`fem_flow_darcy_graph.json` and `fem_flow_darcy_graph_species.json` are
editable named-port examples; their geometry
paths refer to the meshes produced by `fem_flow_darcy_wall_species.json`.

For multi-step fluid–solid ALE, generate a small editable case with
`python3 scripts/generate_fsi_channel_fixture.py DIRECTORY`, then run
`python3 scripts/solver.py DIRECTORY/case.json`. The direct case selects
`native_tet_fsi_steps`; an outer single FEM domain may reference it using
`backend: native_tet_fsi_steps` and `case_file`. The case takes separate
labelled fluid and wall tetra meshes, `interface_label`, fluid density and
viscosity, optional `boundary_velocity_by_label_m_s` and
`natural_boundary_labels`, solid Young modulus, Poisson ratio, density and
`fixed_nodes`, plus `time.dt_s` and `time.steps`. Matching interface facets
are aligned by coordinates even when the mesh node IDs differ. This current
dense FSI runtime uses one CPU rank. Each step writes flow and wall `.pvtu`
snapshots and updates `flow.pvd` and `wall.pvd`. Its step ledger reports
interface force and power, wall displacement, moving volume and flux, and ALE
GCL residual.

The generator also writes `fsi_species_solver.json` and
`fsi_species_case.json`. This variant advances two independent P1 species
after each accepted fluid–solid ALE step, using the current P2 fluid field,
mesh velocity, and deformed tetra geometry. The case's `species` array gives
each `id`, initial concentration, diffusivity, source, optional decay,
labelled inflow donors, and optional wall exchange. The runner writes
`species_ID.pvd` alongside `flow.pvd` and `wall.pvd` and reports the
per-step species inventory and balance defect.

For the single-GPU P2/P1 tetra flow route, build with
`conda run -n tubularflow-cuda make cuda CUDA_ARCHS=89` on this workstation.
Generate the labelled channel with
`python3 scripts/generate_fsi_channel_fixture.py DIRECTORY`, then run
`python3 scripts/solver.py DIRECTORY/gpu_flow_solver.json`. The generated
`gpu_flow_case.json` sets density, viscosity, initial velocity, Newton
controls, labelled velocity and pressure boundaries, `dt_s`, and step count.
`examples/solver/zero_d_fem_gpu_rcr.json` connects the same GPU flow backend
to a CPU 0D source and RCR through named pressure/flow ports. The first route
uses explicit sequential coupling: the source turns the preceding inlet
pressure into the next prescribed inlet flow, while each measured GPU outlet
flow updates its RCR state and next outlet pressure. It writes the usual GPU
VTK series and reports port pressure, flow, peer flow, and 0D storage residual
for every step.
The GPU computes the 64-point tetra element residual and Jacobian, assembles
the global system, solves Newton updates, and advances backward-Euler steps.
It writes `flow_step_N.vtu` and `flow.pvd`; the final ledger reports assembly
time, solve time, host peak RSS, and requested CUDA peak allocation. The
current GPU linear solve is dense and targets small single-GPU cases.
Local assembly, solve, memory, and CPU comparison measurements are in
[`benchmarks/fem_gpu_native_evidence.json`](../benchmarks/fem_gpu_native_evidence.json).

For steady P1 Darcy on the same labelled tetra format, run
`python3 scripts/solver.py examples/solver/fem_darcy_gpu.json`. The
`native_tet_darcy_cuda` backend assembles the P1 matrix and source on the
GPU, solves pressure with cuSOLVER, and evaluates cell flux on the GPU. It
then recovers conservative face flows and RT0 fluxes and writes `tissue.vtu`
and `tissue.pvd`. The case accepts the same mobility, source, pressure, flux,
and per-cell override fields as the CPU Darcy CLI. The GPU solve is dense and
targets small single-GPU meshes.

For explicit sequential GPU fluid–structure coupling, generate the channel
fixture and run `python3 scripts/solver.py DIRECTORY/gpu_fsi_solver.json`.
The `native_tet_fsi_cuda` backend solves each P2/P1 flow step on the GPU,
evaluates matching-interface traction, solves the tetra solid wall, extends
its displacement through the fluid mesh with harmonic ALE, and uses the new
wall and mesh velocity on the next flow step. It writes `flow.pvd` and
`wall.pvd` and reports GCL residual, interface force and power, wall
displacement, solid residual, timings, host RSS, and CUDA allocation. The
solid and ALE algebra run on the host.

For fixed-mesh P1 species transport, run
`python3 scripts/solver.py examples/solver/fem_species_gpu.json`. The
`native_tet_species_cuda` backend assembles the backward-Euler P1 mass,
diffusion, conservative advection, source, decay, and labelled boundary
terms on the GPU, then solves with cuSOLVER. It accepts constant fluid
velocity, prescribed inflow concentrations, and wall exchange from the
native species case. It writes `species_step_N.vtu` and `species.pvd`.
The current GPU case excludes moving meshes, monotone correction, and finite
wall reservoirs.

To drive species from the GPU P2 flow field, generate the channel fixture
and run `python3 scripts/solver.py DIRECTORY/gpu_flow_species_solver.json`.
Its two GPU FEM domains have a sequential `velocity` connection. The entry
solves flow, passes each step's complete P2 velocity field to species, and
writes separate `flow/flow.pvd` and `species/species.pvd` series. Both cases
must use the same labelled tetra mesh and time grid. The exchange carries
velocity on a fixed mesh; moving-wall species remains on the CPU FSI route.
`gpu_flow_multispecies_solver.json` uses the same flow field for two separate
species domains with independent initial values, diffusivities, decay rates,
and inlet concentrations.

`examples/solver/fem_flow_darcy_species_gpu.json` connects a GPU P2/P1 flow
boundary to a GPU Darcy source port and then to GPU P1 species. Each target
Darcy port gives explicit `cell_id_weights`; the host coupling layer integrates
the named boundary of the GPU velocity field and conservatively distributes
that signed flow as source or sink on the selected tissue cells. The tissue
mesh may differ from the vessel mesh. Darcy assembly and solve remain on the
GPU, and its recovered RT0 face flow directly drives species. The example
maps `0.1` m3/s with zero map defect; its GPU Darcy pressure matches the CPU
reference exactly.

`examples/solver/fem_darcy_species_gpu.json` pairs steady GPU Darcy with
transient GPU species on the same labelled tetra mesh. Its sequential
`flow_rate` connection passes recovered conservative RT0 face flows into
the species advection assembly. The entry writes `darcy/tissue.pvd` and
`species/species.pvd`. This does not yet transfer material across separate
vessel and tissue meshes.

For a single domain, `case_dir`, `case_file`, or `workflow_file` references an
existing backend case. Paths in `solver-v1` are relative to that outer JSON.
The existing file retains its full equation, material, boundary, port, and
initial-value syntax. On a referenced 0D, 1D, or IGA `case_dir`, outer
`time.dt_s`, `time.steps`, `waveforms`, and `output.every_steps` create a
runtime `simulation_config.json` and leave the source case intact. The
`waveforms` array uses the existing `temporal_functions` objects. A native
FEM `case_file` accepts outer `dt_s` and `steps`; a FEM workflow's outer time
must match its referenced case. `options` on a domain passes existing CLI
options through to the selected backend, including its restart and solver
controls. `resources.mpi_ranks` selects local MPI ranks. GPU is selectable
for IGA standalone flow or transport, FEM flow, Darcy, fixed-mesh species,
moving-wall flow, sequential flow/Darcy→species, and the explicit 0D
source→FEM→RCR port graph.
The referenced equation system selects steady or transient integration;
`time.mode` describes that choice and does not change it. The current
entry limits `resources.mpi_ranks` to 12.

For an IGA domain, provide a packed `database` and `case_dir`, or set
`prepare: true` on a source case directory. Preparation runs
`scripts/generate_case.sh`, then the chosen IGA solver. `template_directory`
selects the existing template directory with the required filenames.
Schema-v4 `mesh.cross_section` in the referenced case accepts
`{"kind":"default"}` or `{"kind":"circle","target_size":0.25}`.
Preparation preserves `cross_section_template.vtk` and
`merge_template.vtk` previews in `generated/preprocessing/`.
Multi-rank IGA flow writes a ParaView `.pvd` series and rank-owned
`snapshot.pvtu` files.

## Existing controls and their shared-entry route

This inventory was checked against remote `main` at `dab5f25` on
2026-09-21. The shared entry also lists the native FEM work in this tree.

| Function | Existing case or CLI control | Shared entry |
| --- | --- | --- |
| 0D R, RC, RCR network | schema-v3 `equation_systems`, `lumped_parameters`, `boundaries` | `dimension: 0d`, `case_dir`; edit the referenced case |
| 1D rigid/compliant A/Q and transport | schema-v3 `equation_systems`, `fields`, `boundaries` | `dimension: 1d`, `case_dir` |
| 1D→0D R/RC/RCR terminal graph | schema-v3 outlet `node_ids` and `resistance`/`windkessel_rc`/`windkessel_rcr` | named 1D outlet `node_id` ports connected to 0D terminal domains; generated runtime case |
| 1D→FEM→0D flow/species graph | 1D case directory, named outlet node and species binding, labelled tetra inlet/outlet, and RCR species reservoir | `one_d_fem_rcr.json` for flow or `one_d_fem_rcr_species.json` for staged T7 pressure/flow and species ports |
| 0D/1D geometry, multiple inlet/outlet BCs | `geometry.file`, node IDs, boundary names and conditions | referenced `case_dir` |
| 0D/1D time functions and species | `temporal_functions`, `fields`, transport systems | referenced `case_dir`; outer `waveforms` can replace the time functions |
| 0D/1D output and restart | `time.output_every`; `--output-dir`, `--system`, `--checkpoint`, `--restart`, `--visualization-format` | outer `output`, `system`, `options` |
| IGA control mesh and template | schema-v4 `geometry`, `mesh.cross_section`; `tubular_mesh pipeline CASE_DIR TEMPLATE_DIR` | `case_dir`, `prepare`, `template_directory` |
| IGA template previews and auxiliary conversion | `cross_section_template.vtk`, `merge_template.vtk`, `examples/vtk2obj.py` | `prepare` keeps both VTK previews; the existing conversion script remains directly available |
| IGA packed database | `spline`, METIS, `iga_pack`, `.ntiga` | `prepare` or `database` |
| IGA CPU/GPU flow | `iga_navier_stokes` or `iga_cuda navier-stokes`; schema-v4 `equation_systems`, `boundaries` | `method: IGA`, `physics: ["flow"]`, `device` |
| IGA CPU/GPU species | `iga_solve` or `iga_cuda transport`; `velocity_sources`, species systems, wall exchange | `method: IGA`, `physics: ["species"]`, referenced case and `options` |
| IGA output, diagnostics, restart | `--output`, `--output-every`, `--checkpoint`, `--restart`, `--visualization-format`, nonlinear and mass tolerances | outer `output`; backend `options` |
| FEM labelled tetra or surface meshing | `input_route`, `input_file`, `mesher` in native workflow | `method: FEM`, `workflow_file` |
| FEM P2/P1 flow with 0D source/RCR | native hydraulic case `fluid`, `source`, `terminal_rcrs`, `boundaries`, `motion`, `coupling`, `time` | referenced case/workflow or named-port source–vessel–RCR graph |
| FEM GPU P2/P1 flow and 0D ports | `native_tet_flow_cuda` labelled tetra, velocity/pressure labels, fluid and time controls; optional generated `t7_ports` | `method: FEM`, `device: GPU`, `backend: native_tet_flow_cuda`, `case_file`; standalone or `zero_d_fem_gpu_rcr.json` on one GPU |
| FEM GPU moving-wall flow | `native_tet_fsi_cuda` labelled fluid and wall tetra meshes, matching interface, material, BC and time controls | `method: FEM`, `device: GPU`, `backend: native_tet_fsi_cuda`, `physics: ["flow", "solid"]`; explicit sequential coupling |
| FEM P1 Darcy and RT0 flux | native Darcy `mobility`, source and boundary-label maps | `backend: native_tet_darcy`, `case_file` or workflow |
| FEM GPU P1 Darcy and RT0 flux | same native Darcy case; GPU P1 assembly, solve, and cell flux, followed by conservative RT0 recovery; optional named GPU-flow boundary to explicit Darcy cell source map | `method: FEM`, `device: GPU`, `backend: native_tet_darcy_cuda`, `case_file`; standalone or flow→Darcy on one GPU |
| FEM P1 species | native species concentration, velocity, diffusivity, source, decay, inflow and wall-exchange maps | `backend: native_tet_species_transport`, `case_file` or workflow; or `physics: ["flow", "species"]` and `species` on the FEM vessel in a source–vessel–RCR graph |
| FEM GPU P1 species | native species case with fixed mesh, constant or connected P2 velocity or Darcy RT0 flux, inflow and wall exchange | `method: FEM`, `device: GPU`, `backend: native_tet_species_cuda`, `case_file`; one GPU or sequential flow→species or Darcy→species graph |
| Idealized FEM artery–Darcy–vein–wall–species reference | `cases/idealized_cube_dual_tree_refined.json`, matching oxygen case, four terminal pairs | `backend: idealized_cube_dual_tree`, `case_file`, optional `species_case_file`; one `solver-v1` command |
| FEM artery–Darcy–vein transfer with optional species | labelled tetra meshes; arbitrary arrays of matching arterial and venous port labels | three `role` domains with named `flow_rate` connections, or direct `native_tet_flow_darcy_chain` case; optional three-region species time loop |
| FEM multi-step fluid–solid ALE with optional species | `native_tet_fsi_steps` case with labelled fluid and wall meshes, material, BC, time, and optional `species` array | direct case import or single FEM domain with `backend: native_tet_fsi_steps`; one CPU rank |
| FEM fields and restart | native `--output-dir`, `--checkpoint-dir`, `--restart-dir`, `--resume`, `--stop-after-step` | `output`, `options`, or workflow `--resume`/`--stop-after-step`; ParaView `series.pvd`; coupled flow/species has no restart |

Native tetra flow, Darcy, and species cases use the labelled tetra mesh
described in [the FEM workflow guide](T9_NATIVE_TET_WORKFLOW.md). A surface
route uses Gmsh or fTetWild through the same workflow. The editable
`fem_darcy.json`, `fem_species.json`, and `fem_flow.json` examples use the
small four-tetra functional fixture.

For a surface that needs boundary labels, run
`python3 scripts/label_surface_bc_gui.py INPUT.vtp LABELLED.vtp` on a desktop
with VTK and Tk. STL, PLY and OBJ inputs are also accepted. Click a triangle
to select it; Ctrl-click toggles one triangle and Shift-drag selects visible
triangles in a box. Selected triangles are yellow. Press `i` to assign the
inflow label 1, `o` to assign the current outflow label, `n` to choose the
next outflow label, and `w` to return selected triangles to wall label 0.
Press `u` to undo, `s` to save, `x` to confirm cutting selected triangles,
and `c` to choose and confirm a cap for an open loop. The saved VTP contains
the `boundary_id` cell array used by `surface_to_fem_volume.py` and
`ftetwild_to_fem_volume.py`; a companion `.labels.json` lists the labels.
The cap tool creates a triangular fan, so inspect its preview before meshing
on a nonplanar opening.

Run `make solver-test` to check every editable shared-entry config, the direct
import routes, and the display-independent GUI state and surface-output path.

The connected graphs handle 1D→0D R/RC/RCR terminals, staged T7
1D→FEM→RCR flow and species, source→FEM→RCR flow with multiple vessel
species, and steady FEM artery→Darcy→vein transfer with multiple species.
A separate one-rank FEM fluid–solid ALE runner supports multi-step wall
motion and multiple species. GPU flow supports explicit 0D source/RCR ports;
GPU flow or Darcy can drive one or more GPU species. GPU flow can also feed a
separate Darcy mesh through explicit named boundary and tissue-cell weights.
The shared entry keeps the existing IGA core unchanged.
