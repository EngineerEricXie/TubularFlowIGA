# Bezier Extraction Preprocessor

This C++11/OpenMP program converts `controlmesh.vtk` into the Bezier extraction
data consumed by METIS and `iga_pack`. Its primary interface is a versioned
sparse binary cache; historical text output remains available.

## Build and Run

From the repository root:

```bash
make spline EIGEN_DIR=/path/to/eigen3
OMP_NUM_THREADS=8 ./preprocessing/spline/spline /path/to/case/
OMP_NUM_THREADS=8 ./preprocessing/spline/spline \
  /path/to/case/ --no-legacy-text
OMP_NUM_THREADS=8 ./preprocessing/spline/spline \
  /path/to/case/ --no-legacy-text --legacy-vtk
```

The trailing slash is required because the argument is a directory prefix.
Run large cases on an allocated compute resource and choose a thread count that
matches the allocated CPU cores.

The case directory must contain `controlmesh.vtk`. All invocations write:

- `bzmeshinfo.txt`: element basis connectivity used by METIS;
- `spline_cache.igacache`: sparse coefficients and 64 Bezier points per element;
- `geometry_transform.json`: source origin and normalization scale.

The first command also writes dense `cmat.txt` and `bzpt.txt` for legacy tools.
The second is recommended for the native pipeline and avoids those large files.
Add `--legacy-vtk` only when the historical eight-corner linear
`bzmesh.vtk` preview is needed. It is not consumed by METIS, `iga_pack`, the
solvers, or the cubic VTKHDF exporter.

## Resource Design

Elements are processed in 256-element chunks. Extraction matrices are moved,
not copied, and released after their formatted records are written. Record
formatting runs in parallel, while final file writes remain ordered. A
thread-local control-point map replaces the former full-node allocation and
initialization performed for every element. Unused eager 64-point buffers were
also removed from `Element3D`.

The cache includes a format version, dimensions, and control-mesh content hash.
The packer rejects stale, corrupt, or truncated input rather than silently
falling back. Use `iga_pack ... --legacy-text` to explicitly test text input.

The legacy output format and floating-point text are unchanged. Cylinder and
`NMO_54499_new` outputs match the legacy extractor byte for byte, and the
resulting NMO `.ntiga` database is also identical.

On the 31,680-element NMO case with eight OpenMP threads, the optimized version
used 158.1 MiB peak RSS and 9.34 seconds, versus about 2.19 GiB and 62.50 seconds
for the legacy implementation. Cache-only extraction took 3.74 seconds and
packing took 4.47 seconds. See [the benchmark report](../../docs/BENCHMARKS.md)
for validation context.
