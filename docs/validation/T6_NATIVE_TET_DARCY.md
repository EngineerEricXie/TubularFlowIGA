# Native tetra Darcy tissue-flow foundation

`NativeTetDarcyPetsc.hpp` implements a project-owned P1 tetrahedral FEM for
stationary, single-compartment Darcy pressure on a labelled 3D volume mesh:

```text
q = −K ∇p                         [m/s]
−∇·(K ∇p) = s                    [1/s]
K                                 [m²/(Pa·s)]
p                                 [Pa]
```

The caller supplies a positive scalar mobility and volumetric source per
tetrahedron, constant pressure on at least one boundary label, and optional
outward-positive normal flux density on other labels. Omitted labels have
natural zero flux. The native code computes element gradients, volume and
face terms, Dirichlet elimination, cell fluxes and boundary-flow diagnostics;
PETSc is used only for sparse algebra and MPI. fTetWild or Gmsh can provide
the labelled tetrahedra, not the FEM solver.

This is a fixed-domain, single-compartment, single-physics calculation API
with a standalone version-1 JSON CLI. It now also recovers one integrated
normal flow per tetra face: internal-face values are equal and opposite,
each cell's outward sum matches its prescribed volume source, and natural or
specified Neumann boundary flows are held fixed. A graph-Laplacian
least-squares correction changes only free internal and pressure-boundary
face flows; it does not change the FEM pressure. From the four corrected
flows, the native RT0 basis `q(x)=Σ_i F_i (x−x_i)/(3V)` gives a vector field
with continuous normal trace across conforming tetra faces and cellwise
`div(q)=s`. Its tangential component may jump. This is an **H(div)-conforming
postprocessed field**, not a mixed Darcy pressure/flux discretization. The
original discontinuous P1-gradient cell flux and boundary-flow diagnostics
remain available and should not be confused with the RT0 values. Variable
porosity/storage, multi-compartment exchange and vessel-to-tissue mapping
remain absent. No liver physiology or patient parameters are inferred.

Build and run a labelled Gmsh 4.1 tetra case:

```bash
make -C solvers/cpu native_tet_darcy PETSC_DIR=/path/to/petsc
mpiexec -np 2 solvers/cpu/native_tet_darcy tissue.json --check-input
mpiexec -np 2 solvers/cpu/native_tet_darcy tissue.json --output-dir tissue_run_001
```

```json
{
  "schema_version": 1,
  "mesh_file": "labelled_tissue.msh",
  "mobility_m2_pa_s": 0.0001,
  "source_s_inv": 0.0,
  "pressure_by_boundary_label_pa": {"1": 1.0, "2": 0.0},
  "outward_flux_by_boundary_label_m_s": {},
  "mobility_by_cell_id_m2_pa_s": {},
  "source_by_cell_id_s_inv": {}
}
```

`mesh_file` is relative to the JSON file. Constant `mobility_m2_pa_s` and
`source_s_inv` provide defaults; the optional maps override values on
individual Gmsh tetra cell IDs. A pressure and prescribed outward flux may
not occupy the same label. Labels without either condition are zero-flux;
at least one pressure label is required. `--check-input` uses the same
parser, mesh reader and solver but writes nothing. A run requires a new
output directory, publishes `fields/snapshot.pvtu` with rank-owned tetra
`pressure_pa` point values and raw `darcy_flux_m_s`,
`darcy_rt0_centroid_flux_m_s`, `mobility_m2_pa_s`,
`source_s_inv` cell values, plus
`conservative_outward_face_flow_m3_s` cell values with four components in
opposite-local-node order. `run_summary.json` includes case/mesh SHA-256,
both boundary-flow definitions, RT0 field kind, recovery iterations and
maximum cell-balance defect. The C++ API also evaluates RT0 flux at any
position in a tetra; the VTU centroid vector alone does not encode its
within-cell linear variation. No restart or coupled graph state is implied
by this steady one-shot CLI.

The native P1 species assembly/PETSc step can optionally consume the RT0
four-face flow array directly at volume and boundary quadrature points,
without projecting it to shared P2 velocity nodes. A unit-cube Darcy
`q=(2,0,0) m/s` case gives the same one-step tracer field and labelled
boundary flows as an independently specified constant P2 velocity on
1/2/4 MPI ranks. This is a **same-mesh transport API functional test**;
the Darcy and species CLIs are still separate, and no blood-vessel-to-liver
tissue transfer or organ case is implied.

A separate [generic signed port-to-cell source map](T6_NATIVE_VESSEL_TISSUE_SOURCE_MAP.md)
can turn explicitly supplied vessel `m³/s` flows and normalized tissue-cell
weights into Darcy `1/s` sources. It is tested with forward/reverse unit-cube
flows but does not yet connect a measured vessel FEM port or establish a
geometric liver vessel/tissue interface.

The shared surface/volume workflow also accepts the explicit
`native_tet_p1_darcy_steady` backend. It copies a labelled Gmsh volume or
meshes a labelled VTP/STL surface with fTetWild/Gmsh, rewrites `mesh_file` in
the effective case, runs Darcy preflight, and records input/mesh/case/solver
hashes, stage timings and functional-only status. For example:

```json
{
  "schema_version": 1,
  "backend": "native_tet_p1_darcy_steady",
  "input_route": "volume",
  "input_file": "labelled_tissue.msh",
  "case_file": "tissue.json",
  "output_directory": "tissue_workflow_001",
  "mpi_ranks": 2
}
```

Run `python3 scripts/run_native_tet_workflow.py workflow.json`.
For `surface`, supply the same `mesher` object used by the native tetra
hydraulic/species workflow. The effective Darcy CLI output is under
`output_directory/fields/`, including `fields/snapshot.pvtu` and
`run_summary.json`. This steady route rejects `--stop-after-step` and
`--resume` before creating output; its workflow `accepted_steps=1` means one
completed steady solve, not a time step.

Local functional gate:

```bash
make t6-native-tet-darcy-test PETSC_DIR=/path/to/petsc
make t6-native-tet-darcy-test PETSC_DIR=/path/to/petsc \
  FTETWILD_BIN=/path/to/FloatTetwild_bin
```

On a unit cube, `p=1−x Pa`, `K=2 m²/(Pa·s)`, `s=0` yields exactly
`q=(2,0,0) m/s`, inlet `−2 m³/s` and outlet `+2 m³/s`; a prescribed
outlet-flux version checks the Neumann sign. Doubling `K` doubles flow while
leaving pressure unchanged. The source manufactured solution
`p=x(1−x)/2 Pa`, `K=1`, `s=1 s⁻¹` has volume-weighted tetra-centroid RMS
pressure errors `0.0263008` and `0.0065752 Pa` on 2³ and 4³ cube grids,
respectively; nodal samples happen to be exact on these symmetric grids and
would conceal the interior interpolation error. The cube tests pass at
1/2/4 MPI ranks, check every cell balance, internal-face continuity, RT0
normal-face integrals and divergence with nonzero source, and reject bad
labels, contradictory boundary conditions
and nonpositive mobility. The CLI regression checks 1/2/4-rank VTU
ownership, provenance, no-overwrite and invalid cases. The workflow
regression covers volume, fTetWild and Gmsh surface routes, failed preflight,
and unsupported steady restart.

The optional fTetWild regression generates a newly labelled pipe. It checks
flow direction, near-zero side-wall flow and geometry-sensitive `p(x)`/cell
flux differences, with functional gates `max |Δp| ≤ 0.005 Pa` and
`max |Δq| ≤ 2e−4 m/s`. The generated mesh is not an exact prismatic mesh,
so this is not a spatial-convergence or physiological validation result.
