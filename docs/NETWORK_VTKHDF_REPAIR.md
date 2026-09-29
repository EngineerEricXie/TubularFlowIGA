# Repair midpoint-based network visualization

`scripts/repair_network_vtkhdf_geometry.py` repairs existing 0D/1D network
VTKHDF visualization without rerunning a simulation. It is a postprocessor,
not a change to the native solver output writer.

The old visualization places vertices at segment midpoints and connects those
vertices across junctions. This cuts bends, shifts branch junctions, and omits
half of terminal segments. The repair uses the original OBJ vertices and edges
as VTK line cells. Segment-valued pressure, flow, oxygen, and other arrays move
from `PointData` to `CellData`, without interpolation or numerical changes.
In ParaView, select **Cell Data** to color the repaired geometry.

## Usage

Requires Python 3, NumPy, SciPy, and h5py. Supply the configuration used for the
original run; its geometry path is resolved relative to the configuration file.

```sh
python scripts/repair_network_vtkhdf_geometry.py \
  /path/to/results/profile_1d.vtkhdf \
  --config /path/to/case/simulation_config.json \
  --output /path/to/results/profile_1d_geometry_corrected.vtkhdf
```

Without `--output`, the corrected file is written beside the source, with
`_geometry_corrected` appended to its stem. An accompanying `.validation.json`
records source/output paths, array names, geometry counts, and the time range.
Existing output or report files are never overwritten. The source is read-only.
If conversion fails after output creation, the incomplete output must be
removed before retrying; only a completed conversion produces the report.

## Supported inputs and checks

- Original midpoint-based temporal PolyData or UnstructuredGrid output, with
  exactly one cell per OBJ segment, fixed geometry, and strictly increasing times.
- OBJ `v` and `l` records with positive, one-based indices; coordinates are scaled
  using `geometry.length_scale_to_m`. Use the exact OBJ from the original run.
- Native parent/child node IDs determine field-to-edge ordering when available.
  Otherwise the tool checks OBJ order, the solver's sorted undirected-edge order,
  then a unique midpoint matching fallback (absolute tolerance `1e-12` metres).
- All stored timesteps are retained in one temporal UnstructuredGrid VTKHDF.
  After writing, coordinates, connectivity, times, and every dynamic field slice
  are reopened and checked for exact equality with their intended values.

This does not repair solver physics, recover missing timesteps, support moving
geometry or subdivided segments, or turn cell values into nodal values. Do not
use it on an already corrected file or arbitrary multipart/volume VTKHDF files.

## Regression tests

```sh
python -m unittest discover -s scripts/tests -p test_repair_network_vtkhdf_geometry.py -v
```

Tests cover a bent/branched network with a cycle, merged and native storage,
field ordering, scalar/vector fields, static metadata, exact time/value
preservation, source immutability, overwrite protection, and invalid inputs.
Set `TMPDIR` to keep temporary test artifacts inside a desired output directory.
