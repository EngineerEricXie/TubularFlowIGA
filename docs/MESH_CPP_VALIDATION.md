# C++ Mesh Generator Validation

The C++ control-mesh generator was validated on PSC Bridges-2 compute nodes on
2026-08-04. Timings use optimized builds and `/usr/bin/time -v`; generated case
data is not committed.

## MATLAB Compatibility

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

## Resource Results

| Case | Wall time | Peak RSS | Minimum determinant | Minimum scaled J |
| --- | ---: | ---: | ---: | ---: |
| Cylinder | 0.06 s | 4.5 MiB | `7.6168e-5` | `0.763084` |
| Bifurcation | 0.21 s | 6.6 MiB | `2.39749e-5` | `0.382533` |
| NMO_54499_new | 0.76 s | 10.8 MiB | `1.22726e-6` | `0.0315119` |

The corrected NMO case produced 35,949 control points and 31,680 elements with
zero invalid elements. The older VTK stored in the historical example folder
has 35,748 points and is stale, so it was not used as a golden mesh.

## Downstream Large-Case Gate

The NMO C++ output was passed through the unmodified legacy spline extractor,
partitioned with METIS for eight ranks, packed into `.ntiga`, and checked with:

```bash
./preprocessing/spline/spline "$CASE_DIR/"
mpmetis "$CASE_DIR/bzmeshinfo.txt" 8
./solvers/cpu/iga_pack "$CASE_DIR" 8 case-8.ntiga
mpiexec -np 8 ./solvers/cpu/iga_mesh_check case-8.ntiga
```

All 31,680 extracted elements and 2,027,520 geometry samples passed:
`min(detJ)=1.27675e-7`, with no bad elements or samples. The original spline
extractor took 62.50 s and about 2.19 GiB peak RSS. Its optimized chunked version
took 9.34 s and 158.1 MiB with eight OpenMP threads; all four text outputs and
the resulting `.ntiga` database matched byte for byte. Packing remained a
separate 26.86-second text-parsing stage with about 4.1 MiB RSS.

Scheduler job IDs were 43011901 (bifurcation MATLAB reference), 43011942
(three-bifurcation reference), 43011981 (NMO C++ mesh), and 43012005 (NMO
downstream validation).
