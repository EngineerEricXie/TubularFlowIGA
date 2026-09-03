#!/usr/bin/env python3
"""Convert VTK PolyData centerlines to TubularFlowIGA radius-annotated OBJ."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any, Sequence


def write_radius_annotated_obj(
    mesh: Any,
    output: Path,
    *,
    scale: float,
    radius_array: str,
) -> None:
    """Write validated point radii and polyline segments using OBJ indices."""
    if not math.isfinite(scale) or scale <= 0.0:
        raise ValueError("scale must be finite and positive")
    if radius_array not in mesh.point_data:
        available = ", ".join(sorted(mesh.point_data.keys())) or "none"
        raise ValueError(
            f"point-data radius array '{radius_array}' was not found; "
            f"available arrays: {available}"
        )

    points = scale * mesh.points
    radii = scale * mesh.point_data[radius_array]
    if len(radii) != mesh.n_points:
        raise ValueError("radius array length does not match the VTK point count")
    for index, radius in enumerate(radii):
        value = float(radius)
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"radius at point {index} must be finite and positive")

    lines = mesh.lines
    if len(lines) == 0:
        raise ValueError("VTK PolyData contains no centerline polylines")
    segments = []
    offset = 0
    while offset < len(lines):
        point_count = int(lines[offset])
        if point_count < 2 or offset + point_count >= len(lines):
            raise ValueError(f"invalid VTK polyline record at offset {offset}")
        for local in range(point_count - 1):
            first = int(lines[offset + 1 + local])
            second = int(lines[offset + 2 + local])
            if not (0 <= first < mesh.n_points and 0 <= second < mesh.n_points):
                raise ValueError("VTK polyline contains an out-of-range point index")
            segments.append((first, second))
        offset += point_count + 1

    with output.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(f"# Centerline converted from VTK; points={mesh.n_points}\n")
        stream.write("# v x y z radius auxiliary0 auxiliary1\n")
        for point, radius in zip(points, radii):
            x, y, z = point
            stream.write(
                f"v {float(x):.17g} {float(y):.17g} {float(z):.17g} "
                f"{float(radius):.17g} 0 0\n"
            )
        for first, second in segments:
            stream.write(f"l {first + 1} {second + 1}\n")


def parse_args(arguments: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert a VTK centerline to TubularFlowIGA radius-annotated OBJ."
    )
    parser.add_argument("input", type=Path, help="input VTK PolyData file")
    parser.add_argument("output", type=Path, help="output radius-annotated OBJ file")
    parser.add_argument(
        "--radius-array",
        default="radius",
        help="point-data array containing centerline radii (default: radius)",
    )
    parser.add_argument(
        "--scale",
        type=float,
        default=1.0,
        help="positive scale applied to both coordinates and radii (default: 1)",
    )
    return parser.parse_args(arguments)


def main(arguments: Sequence[str] | None = None) -> int:
    options = parse_args(arguments)
    try:
        import pyvista as pv
    except ImportError as error:
        raise SystemExit(
            "PyVista is required for VTK conversion; install it with "
            "'python3 -m pip install pyvista'"
        ) from error

    mesh = pv.read(options.input)
    write_radius_annotated_obj(
        mesh,
        options.output,
        scale=options.scale,
        radius_array=options.radius_array,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
