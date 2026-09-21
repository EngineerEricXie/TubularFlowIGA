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
