# C++ Mesh Generator Validation

## Current adaptive-mesh regression

The schema-v4 generator and all committed 3D geometries were rerun locally on
2026-08-27. Every case had zero invalid hexahedra, zero exterior-surface
intersections, and a minimum scaled Jacobian above its configured `0.1` floor.

| Case | Points | Hexes | Minimum scaled J |
|---|---:|---:|---:|
| `neuron_transport/straight_neurite` | 1,005 | 720 | `0.763084` |
| `neuron_transport/branched_neurite` | 7,731 | 6,660 | `0.533324` |
| `neuron_transport/nmo_06840_bifurcation` | 41,097 | 36,540 | `0.288586` |
| `vascular_flow/straight_tube` | 1,005 | 720 | `0.763084` |
| `vascular_flow/bent_tube` | 1,809 | 1,440 | `0.745694` |
| `vascular_flow/y_bifurcation` | 7,731 | 6,660 | `0.533324` |
| `vascular_flow/vca_bifurcation` | 7,731 | 6,660 | `0.533324` |
| `vascular_flow/multispecies_pulse` | 1,809 | 1,440 | `0.745694` |
| `vascular_flow/iga_wordmark` | 41,205 | 36,720 | `0.743503` |
| `vascular_flow/liver_vein_obj_segment` | 2,412 | 1,980 | `0.763032` |

The regression suite also rejects short/thick and near-collinear Y junctions,
extreme junction radius ratios, overlapping swept centerlines, and intersecting
disconnected volumes. It verifies adaptive turn/diameter limits at very small
coordinate scale and runs cleanly under AddressSanitizer and
UndefinedBehaviorSanitizer (leak detection disabled where the sandbox cannot
support it).

The current Y and NMO meshes also completed spline extraction, METIS packing,
version-5 database inspection, boundary-label resolution, and two-rank PETSc
geometry checks. The Y reported `minimum_detJ=1.5105e-5`; NMO reported
`minimum_detJ=7.24399e-11`. Both had zero bad elements and zero bad quadrature
samples. The crossing-sensitive wordmark also passed all 36,720 packed elements
with `minimum_detJ=2.24357e-8`.

## Historical MATLAB compatibility

The following pre-adaptive generator was validated on PSC Bridges-2 compute
nodes on 2026-08-04. These results document lineage and the retained file
interfaces, not byte-for-byte equivalence of the current adaptive meshes.
Timings use optimized builds and `/usr/bin/time -v`; generated case data is not
committed.

### MATLAB compatibility

Reference meshes were regenerated with the tracked `TreeSmooth.m` and
`Hexmesh_main.m`, using TREES 1.15. Comparisons parse VTK data rather than
depending on file ordering outside the documented interfaces.

| Case | Points | Hexes | Connectivity | Labels | Maximum point difference |
| --- | ---: | ---: | --- | --- | ---: |
| Cylinder | 4,221 | 3,600 | exact | exact | `1e-6` |
| Bifurcation | 14,565 | 12,780 | exact | exact | `1e-6` |
| Three bifurcations | 34,851 | 30,780 | exact | exact | `1e-6` |

The point tolerance equals the six-decimal VTK text precision. Smoothed SWC
parent topology matched exactly and coordinates/radii agreed within `5e-9`.
The three-bifurcation velocity field matched exactly. A regression test also
checks rotations, cubic B-spline endpoints and derivatives, unit-cube
Jacobians, SWC section traversal, and a 3,600-element cylinder mesh.

### Resource results

| Case | Wall time | Peak RSS | Minimum determinant | Minimum scaled J |
| --- | ---: | ---: | ---: | ---: |
| Cylinder | 0.06 s | 4.5 MiB | `7.6168e-5` | `0.763084` |
| Bifurcation | 0.21 s | 6.6 MiB | `2.39749e-5` | `0.382533` |
| NMO_54499_new | 0.76 s | 10.8 MiB | `1.22726e-6` | `0.0315119` |

The corrected NMO case produced 35,949 control points and 31,680 elements with
zero invalid elements. The older VTK stored in the historical example folder
has 35,748 points and is stale, so it was not used as a golden mesh.

### Downstream large-case gate

The NMO C++ output was passed through spline extraction, partitioned with METIS
for eight ranks, packed into `.ntiga`, and checked with:

```bash
OMP_NUM_THREADS=8 ./preprocessing/spline/spline \
  "$CASE_DIR/" --no-legacy-text
mpmetis "$CASE_DIR/bzmeshinfo.txt" 8
./solvers/cpu/iga_pack "$CASE_DIR" 8 case-8.ntiga
mpiexec -np 8 ./solvers/cpu/iga_mesh_check case-8.ntiga
```

All 31,680 extracted elements and 2,027,520 geometry samples passed:
`min(detJ)=1.27675e-7`, with no bad elements or samples. The original spline
extractor took 62.50 s and about 2.19 GiB peak RSS. Its optimized chunked version
took 9.34 s and 158.1 MiB with eight OpenMP threads; all four text outputs and
the resulting `.ntiga` database matched byte for byte. Packing remained a
separate 26.86-second text-parsing stage with about 4.1 MiB RSS. The subsequent
sparse-cache path reduced extraction to 3.74 s and packing to 4.47 s, while
producing a byte-identical `.ntiga` database. Stale and truncated cache
regressions were rejected before packing.

Scheduler job IDs were 43011901 (bifurcation MATLAB reference), 43011942
(three-bifurcation reference), 43011981 (NMO C++ mesh), 43012005 (NMO
downstream validation), and 43012857 (sparse-cache preprocessing benchmark).
