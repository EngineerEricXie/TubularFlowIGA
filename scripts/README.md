# Scripts

Run every script from the repository root. Python scripts print their options
with `--help`.

```text
scripts/
├── solver.py, run_cases.sh, ...   user entry points and case preparation (below)
├── integration/                   end-to-end tests that drive built native Tet binaries
├── hpc/                           cluster deployment, scaling, and provenance tools
└── tests/                         unit tests: python3 -m unittest discover -s scripts/tests
```

## Entry points

| Script | Purpose |
|---|---|
| `solver.py` | Run 0D, 1D, IGA, and native Tet FEM cases from one JSON configuration ([guide](../docs/SOLVER_ENTRY.md)) |
| `run_cases.sh` | Run the cases under `Input/` with the profile in `execution.conf` |
| `prepare_example.sh` | Build and validate the complete body-fitted IGA preprocessing pipeline for an example |
| `generate_case.sh` | Generate a control mesh, spline, partition, and `.ntiga` database for one case |
| `check_dependencies.sh` | Check preprocessing, CPU, 1D, or CUDA build dependencies |
| `run_native_tet_workflow.py` | Labelled surface or volume to native Tet hydraulic, species, or Darcy FEM ([guide](../docs/NATIVE_TET_WORKFLOW.md)) |

## Geometry and meshing

| Script | Purpose |
|---|---|
| `surface_to_fem_volume.py` | Closed labelled surface to first-order tetrahedral Gmsh mesh |
| `ftetwild_to_fem_volume.py` | Run fTetWild and normalize the result to the native labelled Gmsh 4.1 contract |
| `label_surface_bc_gui.py` | Interactive inlet/outlet/wall labelling, cutting, and capping of a vessel surface |
| `generate_circular_pipe_surface.py`, `generate_bent_pipe_surface.py`, `generate_y_pipe_surface.py` | Labelled analytic test surfaces |
| `generate_fsi_channel_fixture.py` | Small labelled fluid/wall Tet fixture for the native FSI runner |
| `vtk_centerline_to_obj.py` | VTK PolyData centerline to radius-annotated OBJ |
| `dicom_seg_to_surface.py` | One DICOM SEG label to a provenance-tracked voxel isosurface |
| `dicom_seg_to_multiregion_tet.py` | Same-grid DICOM SEG labels to a conforming multiregion Tet ROI |
| `split_multiregion_tet_for_native.py` | Split an audited multiregion mesh into native FEM submeshes |
| `audit_dicom_seg_region_overlap.py`, `audit_seg_exposed_vessel_faces.py`, `audit_multiregion_tet_mesh.py` | Segmentation and multiregion mesh audits |

## Showcase and validation cases

| Script | Purpose |
|---|---|
| `generate_idealized_cube_dual_tree.py` | Dual-tree cube geometry ([case](../docs/cases/CUBE_DUAL_TREE_FSI.md)) |
| `run_idealized_cube_dual_tree_fixed_flow.py`, `run_idealized_cube_dual_tree_oxygen.py` | Rerun the dual-tree FSI and passive-tracer cases |
| `prepare_idealized_cube_oxygen_fields.py` | Extract frozen FSI fields for the tracer case |
| `compare_idealized_cube_dual_tree.py`, `compare_idealized_cube_oxygen.py` | Compare independent dual-tree runs by global ID |
| `run_liver_roi_from_raw.py`, `run_liver_roi_functional_case.py`, `run_native_matching_roi_functional.py` | Liver ROI vessel/Darcy functional case |
| `prepare_t2_poiseuille_meshes.py`, `prepare_t2_bifurcation_meshes.py`, `run_native_t2_level.py`, `collect_native_t2_results.py`, `t2_contract.py` | T2 fixed-wall flow validation |
| `run_t7_native_0d_cross_process.py` | Native Tet/0D restart across two MPI process lifetimes |
| `validate_vca_bifurcation.sh` | Build and validate the 3D VCA bifurcation example |

## Evidence validators

`validate_*.py` scripts check the checked-in evidence under
[`benchmarks/`](../benchmarks) against its contract and source files. The
`t<N>` milestone codes match the root `make t<N>-…` targets and the documents
under [`docs/validation/`](../docs/validation).

| Script | Make target |
|---|---|
| `validate_supported_input_contract.py` | `make contract-audit` |
| `validate_t1_*` | `make t1-audit` |
| `validate_t2_*` | `make t2-contract-audit`; `validate_t2_fixed_flow_result.py` checks a three-level solver result |
| `validate_t3_*` | `make t3-solid-contract-audit` |
| `validate_t4_*` | `make t4-contract-audit` |
| `validate_t5_*` | `make t5-matching-interface-audit` |
| `validate_t6_*`, `validate_t7_*`, `validate_t9_*` | the matching `make t6-…`, `t7-…`, and `t9-…` targets |
| `validate_liver_roi_functional_evidence.py` | run directly after the liver ROI case |
| `validate_native_vtkhdf_reader.py` | run with `pvpython` against native VTKHDF output |

## Visualization

| Script | Purpose |
|---|---|
| `render_example.py` | Batch ParaView rendering of a CPU result |
| `render_multiphysics_gif.py`, `render_transport_gif.py` | Animated GIFs from PVD collections |
| `export_idealized_cube_paraview.py`, `visualize_idealized_cube_dual_tree.py` | Dual-tree cube ParaView export and figures |
| `repair_network_vtkhdf_geometry.py` | Repair legacy midpoint-based 0D/1D VTKHDF only ([guide](../docs/NETWORK_VTKHDF_REPAIR.md)) |

## Build provenance

| Script | Purpose |
|---|---|
| `native_checkpoint_source_identity.py` | Hash native graph sources at build time for checkpoint identity |
| `time_native_mpi_rank.py` | Measure one native MPI rank with GNU time |

## Subdirectories

- [`integration/`](integration): `test_native_tet_*.py` and
  `test_native_vtkhdf_cli.py` run built native Tet executables end to end.
  They are invoked by the `make t7-…` and `make t9-…` targets, not by unit-test
  discovery. `test_*_paraview.py` read solver output with ParaView's readers
  (`pvpython script.py OUTPUT_DIRECTORY`).
- [`hpc/`](hpc): test tiers, CPU and single-GPU matrices, cross-node scaling,
  and scheduler records. `hpc_*_regression.py`, `hpc_*_prefixes.py`, and
  `hpc_native_graph_checkpoint.py` drive the MPI failure-injection, output,
  checkpoint, and solver-option regressions that no make target runs. See
  [HPC deployment](../docs/hpc/HPC_DEPLOYMENT.md) and
  [HPC benchmarks](../docs/hpc/HPC_BENCHMARKS.md).
- [`tests/`](tests): fast unit tests for the scripts above.
