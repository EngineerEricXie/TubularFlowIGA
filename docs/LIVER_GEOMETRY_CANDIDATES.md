# Liver geometry source and feasibility audit

This document records the provenance and geometry checks used to evaluate
public liver CT/SEG data for CoupledFlow. It separates three distinct results:

1. verified source data and licensing;
2. single-region surface/volume meshing experiments; and
3. a bounded, voxel-faithful multiregion ROI used for numerical integration
   tests.

No patient images or derived patient meshes are committed to the repository.
The recorded ROI uses synthetic boundary conditions and material values; it is
not a whole-liver perfusion model.

## Candidate datasets

| Source | License/access reviewed | Relevant content | Decision |
|---|---|---|---|
| [TCIA Colorectal-Liver-Metastases](https://www.cancerimagingarchive.net/collection/colorectal-liver-metastases/) | DICOM/SEG available through the [NBIA REST API](https://wiki.cancerimagingarchive.net/display/Public/NBIA+Search+REST+API+Guide); selected series metadata and archive license state CC BY 4.0 | Same-patient liver, vessel, tumor, and future-liver-remnant masks described by the [data paper](https://www.nature.com/articles/s41597-024-02981-2) | Selected for provenance and local ROI work |
| [3D-IRCADb-01](https://www.ircad.fr/research-and-development/data-sets/liver-segmentation-3d-ircadb-01/) | Official page states CC BY-NC-ND 4.0 | CT, labeled DICOM, masks, and segmented surfaces | Not selected for redistributable derived inputs without a separate rights review |
| [Medical Segmentation Decathlon](https://registry.opendata.aws/msd/) | AWS registry states CC BY-SA 4.0 | Liver and hepatic-vasculature tasks use different cohorts | Suitable for tool tests, not a same-patient combined perfusion geometry |
| [LIRCAD/Inria](https://zenodo.org/records/13897086) | Public record reviewed; no explicit license text was visible in the record during the audit | Portal/hepatic vessel labels and branch names | Not used because licensing and tissue/port contracts were unresolved |
| [HVA-CT](https://zenodo.org/records/19850108) | Public files listed; explicit dataset license was not visible during the audit | Relabeled portal/hepatic veins and generated liver masks | Retained as a segmentation candidate only |
| [VSNet relabeling](https://github.com/XXYZB/VSNet) | Repository is MIT; underlying images and derived-label rights require a separate review | Portal/hepatic annotations | Not used as a replacement for the verified TCIA route |

The existing `cases/liver_vessels_simple_implicit_pde/` directory is a
synthetic 1D network case. Its waveform, outlet pressures, fluid properties,
and wall material are not sourced from either TCIA patient and are not reused
as patient boundary conditions.

## Verified TCIA sources

### CRLM-CT-1046

The CT and SEG series were downloaded through the public NBIA API on
2026-09-20 UTC. Both series belong to StudyInstanceUID
`1.3.6.1.4.1.14519.5.2.1.9203.8273.333449295436763574579232065427`.

| Field | CT | SEG |
|---|---|---|
| SeriesInstanceUID | `1.3.6.1.4.1.14519.5.2.1.9203.8273.876266207167921740530708709916` | `1.3.6.1.4.1.14519.5.2.1.9203.8273.198419865300741306933067795121` |
| SOPInstanceUID | Individual CT instances | `1.3.6.1.4.1.14519.5.2.1.9203.8273.207996444071145243031318450804` |
| Download SHA-256 | `7c6ed88e5a3014b3916c637f511010bd1491d8ff5561b9e431dafa86d83f31a8` | `e3db49bae8841517c3cbede6afdff21594ef1369e1837aaaee67db64dce2dc76` |
| Content | 47 CT DICOM files, 512×512 | One 107-frame DICOM SEG |

The SEG references 31 CT SOP instances, all present in the CT archive. The
in-plane spacing is `0.933594×0.933594 mm`; segmented planes use 5 mm spacing.
Coordinates are DICOM LPS in millimeters and must not be treated as RAS or
meters without an explicit transform.

| Segment | Nonzero voxels | 6-neighbor components | Largest component |
|---|---:|---:|---:|
| `Liver` | 373,238 | 2 | 373,219 |
| `Liver Remnant` | 223,010 | 1 | 223,010 |
| `Hepatic` | 5,346 | 12 | 5,235 |
| `Portal` | 6,965 | 13 | 6,790 |
| `Tumor_1` | 410 | 2 | 375 |
| `Tumor_2` | 3,215 | 1 | 3,215 |

The segment names identify image labels, not flow direction, inlet/outlet
surfaces, vessel walls, or exchange interfaces.

```bash
case_dir=$(mktemp -d -p /tmp tcia-crlm-1046-XXXXXX)
curl -fsSL --output "$case_dir/seg.dcm" \
  'https://services.cancerimagingarchive.net/nbia-api/services/v1/getSingleImage?SeriesInstanceUID=1.3.6.1.4.1.14519.5.2.1.9203.8273.198419865300741306933067795121&SOPInstanceUID=1.3.6.1.4.1.14519.5.2.1.9203.8273.207996444071145243031318450804'
curl -fsSL --output "$case_dir/ct.zip" \
  'https://services.cancerimagingarchive.net/nbia-api/services/v1/getImage?SeriesInstanceUID=1.3.6.1.4.1.14519.5.2.1.9203.8273.876266207167921740530708709916'
sha256sum "$case_dir/seg.dcm" "$case_dir/ct.zip"
unzip -Z1 "$case_dir/ct.zip"
```

Verify hashes, metadata, archive members, and `LICENSE` before extraction.

### CRLM-CT-1072

This finer-slice case was downloaded from the same collection on 2026-09-20
UTC. Both series metadata and the CT archive license state CC BY 4.0.

| Field | Value |
|---|---|
| CT SeriesInstanceUID | `1.3.6.1.4.1.14519.5.2.1.9203.8273.152783424215299005378585248107` |
| SEG SeriesInstanceUID | `1.3.6.1.4.1.14519.5.2.1.9203.8273.160655374625974558520205684891` |
| SEG SOPInstanceUID | `1.3.6.1.4.1.14519.5.2.1.9203.8273.178516068093830030868820769150` |
| StudyInstanceUID | `1.3.6.1.4.1.14519.5.2.1.9203.8273.186150960588321456990090403242` |
| CT archive SHA-256 | `4ad0142f2eae3320d31a14dfe9270bd2c7e8f9d201f944c41f4ef332018cfaba` |
| SEG SHA-256 | `afe7a2bdf63cb65a5f4cea53cde8f2271bac51f5f53369b2c47a028b083c1654` |

The archive contains 240 CT slices. The 657-frame SEG references 217 CT SOP
instances, all present in the same study. In-plane spacing is `0.724609 mm` and
adjacent image positions differ by approximately `0.8 mm`. The DICOM
`SliceThickness` value is `1.25 mm`; geometry construction uses actual frame
positions rather than substituting slice thickness for z spacing.

| Segment | Nonzero voxels | 6-neighbor components | Largest component |
|---|---:|---:|---:|
| `Liver` | 2,776,613 | 2 | 2,776,612 |
| `Hepatic` | 7,972 | 1 | 7,972 |
| `Portal` | 47,823 | 116 | 47,283 |

```bash
case_dir=$(mktemp -d -p /tmp tcia-crlm-1072-XXXXXX)
curl -fsSL --output "$case_dir/seg.dcm" \
  'https://services.cancerimagingarchive.net/nbia-api/services/v1/getSingleImage?SeriesInstanceUID=1.3.6.1.4.1.14519.5.2.1.9203.8273.160655374625974558520205684891&SOPInstanceUID=1.3.6.1.4.1.14519.5.2.1.9203.8273.178516068093830030868820769150'
curl -fsSL --output "$case_dir/ct.zip" \
  'https://services.cancerimagingarchive.net/nbia-api/services/v1/getImage?SeriesInstanceUID=1.3.6.1.4.1.14519.5.2.1.9203.8273.152783424215299005378585248107'
sha256sum "$case_dir/seg.dcm" "$case_dir/ct.zip"
unzip -Z1 "$case_dir/ct.zip"
```

## Single-region geometry checks

`scripts/dicom_seg_to_surface.py` converts a regular, non-flipped axial DICOM
SEG label to a 0.5-voxel isosurface in LPS coordinates. Component selection and
smoothing are explicit options recorded in the output manifest. The surface
preflight checks connectivity, manifoldness, self-intersection, orientation,
and labeling before volume meshing.

The audit produced the following bounded results:

| Case/region | Result |
|---|---|
| 1046 `Liver` | The raw two-component mask and unsmoothed largest component were rejected. A bounded 10-iteration windowed-sinc repair moved points by at most `0.973 mm` and changed volume by `0.0549%`; the surface then passed preflight. A 0.1 mm fTetWild envelope produced 481,398 tetrahedra with minimum scaled Jacobian `0.07445`. |
| 1046 `Hepatic`/`Portal` | Largest-component surfaces passed closed/manifold checks, but isosurface-to-voxel volume differences were approximately `4.09%` and `5.01%`, reflecting the coarse 5 mm slice spacing. |
| 1072 `Liver` | The bounded repaired surface produced a repeatable single-region Gmsh mesh with 6,694,418 tetrahedra, 1,170,873 nodes, minimum scaled Jacobian `0.002007649`, and exact preservation of 409,188 source boundary facets in the audited environment. |
| 1072 `Hepatic` | The raw surface and a 24,769-tetrahedron single-region Gmsh mesh passed the geometry checks. |
| 1072 `Portal` | The largest-component repaired surface and a 181,427-tetrahedron single-region Gmsh mesh passed the geometry checks. The original label contains 116 disconnected components. |

These are independent single-label meshes. The vessel labels extend both
inside and outside the liver mask, so independently meshed surfaces cannot be
combined into a conservative conforming model without explicit intersections,
openings, region subtraction, and interface policy.

Reproduce overlap and component accounting directly from the original SEG:

```bash
python3 scripts/audit_dicom_seg_region_overlap.py \
  /path/to/seg.dcm 1 3 4 \
  --manifest /tmp/overlap-audit.json
```

The audit reports component bounding boxes, voxel counts, tissue-contact
faces, background/other-label faces, and image-grid-edge contact. Grid-edge or
exposed faces are not automatically classified as flow ports.

## Voxel-faithful multiregion ROI

`scripts/dicom_seg_to_multiregion_tet.py` constructs mutually exclusive
regions from `Liver \ (Hepatic ∪ Portal)` and the vessel labels. Each voxel is
split into six tetrahedra with shared grid nodes. The contract distinguishes
exact tissue--vessel faces from `roi_cut_not_a_flow_port` and
`unclassified_exterior_not_a_flow_port`.

Before publishing a mesh, the workflow checks:

- non-overlapping masks and a single 6-neighbor component per included region;
- face-connected tetrahedral regions;
- positive Jacobians and minimum scaled Jacobian `0.001`;
- complete classification of interior and exterior faces;
- manifold edges and vertices;
- agreement with voxel volumes and contact-face counts; and
- output nonexistence, source identity, and the configured voxel limit.

For CRLM-CT-1072, the half-open ROI
`z=10:20, y=210:226, x=101:117` contains 2,220 labeled voxels and produced:

| Quantity | Value |
|---|---:|
| Nodes | 2,811 |
| Tetrahedra | 13,320 |
| Tissue tetrahedra | 13,128 |
| Portal tetrahedra | 192 |
| Portal--tissue interface triangles | 158 |
| Minimum determinant | `4.200465623×10⁻¹⁰ m³` |
| Minimum scaled Jacobian | `0.37417548` |
| Mesh SHA-256 | `75201bd15c126d5fd8ab71162e4f7f8b242972cc65f502dfb431951ee02acd94` |

```bash
python3 scripts/dicom_seg_to_multiregion_tet.py \
  /path/to/seg.dcm 1 3 4 /tmp/roi-portal-connected.msh \
  --roi 10:20,210:226,101:117 --max-voxels 50000
python3 scripts/audit_multiregion_tet_mesh.py \
  /tmp/roi-portal-connected.msh \
  /tmp/roi-portal-connected.contract.json
```

`scripts/split_multiregion_tet_for_native.py` produces matching single-region
meshes for the native readers while preserving the same physical interface
label on both sides. The audited split contained 13,128 tissue tetrahedra and
192 Portal tetrahedra, with all 158 interface faces matched by coordinates.

## Local coupled functional case

[`cases/liver_roi_functional.json`](../cases/liver_roi_functional.json) binds
the source/mesh hashes, ROI, labels, and explicit SI test values. It solves the
Portal submesh with the native P2/P1 steady flow formulation, transfers each
matching-face flux conservatively, and solves the tissue submesh with P1 Darcy
pressure and RT0 flux recovery.

The functional inputs are:

| Parameter | Value |
|---|---:|
| Artificial inlet velocity | `(0,0,-0.01) m/s` |
| Density | `1000 kg/m³` |
| Dynamic viscosity | `0.004 Pa·s` |
| Interface pressure | `0 Pa` |
| Tissue exterior pressure | `0 Pa` |
| Darcy mobility | `1 m²/(Pa·s)` |

Independent one-, two-, and four-rank runs recorded an inlet magnitude and
interface source of approximately `1.5751746086430×10⁻⁸ m³/s`, Darcy outflow
`1.5751746086512×10⁻⁸ m³/s`, and maximum cell-balance defect
`5.020×10⁻²² m³/s`. Rank-owned VTU/PVTU fields cover every global cell exactly
once and can be compared by global ID. `interface_facets.json` records the
owner cells, coordinates, and signed face flux for all 158 matching faces.

Run the versioned case against locally acquired source data:

```bash
python3 scripts/run_liver_roi_functional_case.py \
  --data-dir /path/to/crlm-ct-1072 --check-only
python3 scripts/run_liver_roi_functional_case.py \
  --data-dir /path/to/crlm-ct-1072 --ranks 2 \
  --field-output-dir /tmp/liver-roi-fields-r2 \
  --output /tmp/liver-roi-result-r2.json
python3 scripts/validate_liver_roi_functional_evidence.py \
  --case-result /tmp/liver-roi-result-r2.json
```

Or rebuild the complete ROI workflow from the original SEG in a new directory:

```bash
python3 scripts/run_liver_roi_from_raw.py \
  --source-seg /path/to/seg.dcm \
  --output-dir /tmp/crlm-1072-roi-run
```

The runner verifies source, geometry, mesh, split, case, and solver identities
before MPI launch and refuses to overwrite existing outputs. A portable summary
of the recorded evidence is stored in
[`benchmarks/liver_roi_functional_evidence.json`](../benchmarks/liver_roi_functional_evidence.json).

## Portal-vein truth-model runs

This is a 2026-10-17 numerical and functional study of the CRLM-CT-1072
`Portal` tree for use as a per-decision 3D model; it makes no physiological
claim. All inputs and outputs stay outside version control and are
identified by SHA-256:

| Artifact | SHA-256 |
|---|---|
| `portal_raw.stl` (`dicom_seg_to_surface.py`, largest component, 47,283 voxels) | `d2a89e8c0377c60dc0865c6aaa246ac2b7f72cda4300923078f31d9d9e8b8e51` |
| `portal.swc` centerline (51 terminals) | `93e915ca3e10cc288cd5e594a1b58eb4917cdc39fc11d77f7d20b46f95a39111` |
| downstream `portal.ports.json` (terminal points and tangents) | `03ec9084065ebec08f06551051cfb90620b407286ed75fede5294226d0628087` |
| derived ports file (50 ports, labels 1–16; each terminal takes the label of the nearest cap of the downstream 15-group surface) | `157498f0f199d869b41644daa60a5f762e630a94c6da2a2d8b5cc1ac1553d614` |
| downstream 1 mm mesh `portal_g15_h1.0mm.msh` | `f25ee927e1e307b862348f59bc50041f5be66b354b6307611d644d9e027e6fe1` |
| downstream case `portal_hydraulic_g15.json` | `71f2f87f7fe7f52a5f812bd17bdf6fc7b17e5277fc13b033376a9d810eee0bde` |

**Downstream mesh.** That mesh preserves a plane-cut, fan-capped surface
exactly (minimum edge `0.18 µm`; 211 of 271 inlet-cap triangles below quality
0.1). It passed only with `--minimum-scaled-jacobian 1e-9` (29 slivers,
minimum `4.05e-7`), and about 70% of its vertices lie on the boundary. The
reported failure ran fixed-point coupling, two ranks and default MUMPS for
about 60 minutes and 10 GB per rank. Two causes compound there. The outer
fixed-point loop diverges, because the terminal `R_eff ≈ 1e9 Pa·s/m³` is about
1000× the 3D resistance. The Newton line search compares mixed-unit residual
norms, so with the inlet constraint still 87% unsatisfied it accepted only
small steps and then failed.

With the implicit predictor, the scaled line search and an exact MUMPS LU at
every Newton iteration, that mesh still stalls. After 12 iterations
(8 ranks, 23 min, 2.5 GB per rank) the residual stays near `1.5e-6`, against
the `1e-8` tolerance. The cause is therefore the discrete problem on this mesh
rather than the linear solver.

**Replacement meshes.** The replacement route uses only the raw isosurface:

```bash
python3 scripts/cap_vessel_ports.py portal_raw.stl ports.json capped.vtp \
  --length-scale-to-m 0.001 --smoothing-iterations 30 --target-size H \
  --centerline-swc portal.swc
python3 scripts/surface_to_fem_volume.py capped.vtp volume.msh \
  --length-scale-to-m 0.001 --target-size-m H/1000 --split-boundary-locked-tetrahedra
```

Smoothing moves vertices by at most `0.75 mm` (median `0.23 mm`, about one
third of the `0.72 mm` in-plane voxel). The enclosed volume ends 6.7% below the
raw isosurface, from smoothing plus the removed blind ends. All 50 ports cut on
the first passing offset except three thin terminals (`outlet_83`,
`outlet_1068`, `outlet_1086`), which use the recorded largest-valid-cap
fallback. The remesh sizes are bounded by the target `H`, the distance to the
centerline (`--radius-factor 1`) and half the wall-to-wall distance. Thin
vessels and pinches therefore dominate the element count, and the three
targets give similar meshes:

| Target `H` | Surface triangles (min quality) | Tetrahedra | Min scaled Jacobian | Gate |
|---|---|---:|---:|---|
| 2 mm (`--radius-factor 2 --thickness-factor 1`) | 44,544 (0.154) | 103,900 | `7.3e-3` | default `1e-3` |
| 2 mm | 64,092 (0.136) | 142,301 | `3.4e-3` | default `1e-3` |
| 1.5 mm | 65,418 (0.037) | 143,147 | `2.0e-3` | default `1e-3` |
| 1 mm | 67,210 (0.038) | 152,181 | `1.6e-4` | explicit `--minimum-scaled-jacobian 1e-4` |

Each of the three default-sized meshes first contained 1–2 tetrahedra whose
ten P2 nodes all lay on the wall. That gives a singular P2/P1 system: MUMPS
reported INFOG(1)=−10. `native_tet_hydraulic_graph` now rejects such meshes
before solving, and the centroid split removes them.

**Solver timings.** The case is the downstream one: 15 grouped RCRs, source
`3e-6 m³/s`, `dt=0.1 s`, one step, `coupling.predictor=implicit_0d`. The
options are those in
[`native_tet_hydraulic_mumps_lagged.opts`](../solvers/cpu/options/native_tet_hydraulic_mumps_lagged.opts).
Runs are local (16-core WSL, Open MPI 4.1.2, PETSc 3.22, MUMPS); every step
is accepted after one coupling iteration and zero executor Newton iterations.
"Process" includes mesh reading, identity hashing and setup:

| Mesh (tets, unknowns) | Ranks | Newton / GMRES | Assembly s | Linear s | Step s | Process s | Peak GiB/rank | Inlet mean p (Pa) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 mm coarse (103,900, 584,355) | 1 | 4 / 59 | 70 | 202 | 282 | 314 | 5.84 | 398.5 |
|  | 2 | 4 / 60 | 49 | 111 | 168 | 186 | 4.34 | 398.5 |
|  | 4 | 4 / 59 | 31 | 74 | 112 | 123 | 3.16 | 398.5 |
|  | 8 | 4 / 59 | 19 | 50 | 75 | 87 | 1.95 | 398.5 |
| 2 mm (142,301, 811,613) | 1 | 4 / 61 | 96 | 256 | 365 | 406 | 7.62 | 390.7 |
|  | 2 | 4 / 61 | 69 | 143 | 223 | 249 | 6.45 | 390.7 |
|  | 4 | 4 / 61 | 43 | 83 | 136 | 152 | 4.11 | 390.7 |
|  | 8 | 4 / 61 | 28 | 59 | 97 | 109 | 2.61 | 390.7 |
| 1.5 mm (143,147, 819,899) | 1 | 4 / 62 | 97 | 251 | 361 | 402 | 7.63 | 390.7 |
|  | 2 | 4 / 62 | 69 | 139 | 219 | 245 | 5.91 | 390.7 |
|  | 4 | 4 / 62 | 43 | 88 | 141 | 157 | 4.25 | 390.7 |
|  | 8 | 4 / 62 | 30 | 65 | 105 | 121 | 2.52 | 390.7 |
| 1 mm (152,181, 863,276) | 1 | 4 / 61 | 103 | 283 | 399 | 444 | 8.24 | 390.8 |
|  | 2 | 4 / 61 | 73 | 159 | 244 | 274 | 6.62 | 390.8 |
|  | 4 | 4 / 61 | 47 | 94 | 151 | 170 | 4.13 | 390.8 |
|  | 8 | 4 / 61 | 30 | 65 | 106 | 118 | 3.11 | 390.8 |

Every 1, 2, 4 and 8-rank run was accepted, with the same Newton and GMRES
counts at every rank count. At 8 ranks every replacement mesh meets the
5-minute, 4 GiB-per-rank target. The process time is at most 121 s, which is
within the 2-minute stretch except for the 1.5 mm mesh at 121 s. Below eight
ranks the replicated mesh and state dominate memory. The 2, 1.5 and 1 mm inlet
pressures agree to `0.03%`; the coarse mesh differs by `2%`.

Per-decision use was measured on the 2 mm mesh at 8 ranks with steady
momentum and `boundary_updates`:

- **From rest:** 6 Newton iterations; 163 s step, 182 s process, 2.71 GiB per
  rank.
- **Controlled restart:** source flow raised 20% and terminal label 4's
  proximal and distal resistances doubled, restarted with
  `--accept-boundary-update`. It took 2 Newton iterations, 9 GMRES
  iterations and one factorization; 59 s step, 68 s process, 2.81 GiB per
  rank.

The restart summary lists the changed domains (`source`, `terminal_4`) and the
history `[1..1, 2..2]`. Inflow rises from `2.996e-6` to `3.599e-6 m³/s`, and the
label-4 outflow falls from `1.134e-7` to `7.78e-8 m³/s`.

## Current acceptance boundary

The source license, same-study registration, coordinates, labels, and local
ROI geometry are verified. The repository also contains a reproducible local
flow-to-Darcy transfer test on those labels.

A whole-liver case still requires a documented policy for disconnected vessel
and tissue components, validated openings and boundary labels, a conforming or
explicitly nonconforming interface strategy, and sourced physical boundary and
material data. The current full-volume command rejects both audited cases
before meshing because the tissue candidate contains two disconnected
6-neighbor components. That rejection is part of the geometry contract.

Patient-derived data and meshes remain outside version control. Large runs
must use an allocated compute resource and retain the source license,
collection citation, de-identified series identifiers, hashes, coordinate
convention, and every explicit geometry repair in their provenance record.
