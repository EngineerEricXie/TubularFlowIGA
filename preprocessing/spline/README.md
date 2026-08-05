# Bezier Extraction Preprocessor

This C++11/OpenMP program converts `controlmesh.vtk` into the Bezier extraction
files consumed by METIS and `iga_pack`. It retains the historical text
interfaces while bounding peak memory independently of the total number of
elements.

## Build and Run

From the repository root:

```bash
make spline EIGEN_DIR=/path/to/eigen3
OMP_NUM_THREADS=8 ./preprocessing/spline/spline /path/to/case/
```

The trailing slash is required because the argument is a directory prefix.
Run large cases on an allocated compute resource and choose a thread count that
matches the allocated CPU cores.

The case directory must contain `controlmesh.vtk`. The extractor writes:

- `bzmeshinfo.txt`: element basis connectivity used by METIS;
- `cmat.txt`: dense legacy extraction coefficients;
- `bzpt.txt`: 64 Bezier control points per element;
- `bzmesh.vtk`: eight-corner visualization mesh.

## Resource Design

Elements are processed in 256-element chunks. Extraction matrices are moved,
not copied, and released after their formatted records are written. Record
formatting runs in parallel, while final file writes remain ordered. A
thread-local control-point map replaces the former full-node allocation and
initialization performed for every element. Unused eager 64-point buffers were
also removed from `Element3D`.

The output format and floating-point text are unchanged. Cylinder and
`NMO_54499_new` outputs match the legacy extractor byte for byte, and the
resulting NMO `.ntiga` database is also identical.

On the 31,680-element NMO case with eight OpenMP threads, the optimized version
used 158.1 MiB peak RSS and 9.34 seconds, versus about 2.19 GiB and 62.50 seconds
for the legacy implementation. See [the benchmark report](../../docs/BENCHMARKS.md)
for validation context.
