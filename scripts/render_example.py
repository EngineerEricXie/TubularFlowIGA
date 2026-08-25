#!/usr/bin/env pvpython

"""Render a TubularFlowIGA CPU result with ParaView in batch mode."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from paraview.simple import (
    ColorBy,
    CreateView,
    GetColorTransferFunction,
    GetOpacityTransferFunction,
    GetScalarBar,
    Render,
    SaveScreenshot,
    Show,
    Slice,
    Text,
    TrivialProducer,
    _DisableFirstRenderCameraReset,
)
from vtkmodules.vtkCommonCore import vtkDoubleArray
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid
from vtkmodules.vtkIOLegacy import vtkUnstructuredGridReader


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Attach CPU coefficient output to controlmesh.vtk and render a center slice."
    )
    parser.add_argument("--kind", choices=("flow", "transport"), required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--array", default="Nplus", help="transport field to color")
    parser.add_argument("--title", required=True)
    parser.add_argument("--subtitle", default="CPU / MPI / PETSc")
    parser.add_argument("--legend", default="")
    parser.add_argument("--width", type=int, default=1400)
    parser.add_argument("--height", type=int, default=800)
    return parser.parse_args()


def read_grid(path: Path) -> vtkUnstructuredGrid:
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()
    source = reader.GetOutput()
    if source is None or source.GetNumberOfPoints() == 0:
        raise RuntimeError(f"mesh has no points: {path}")
    grid = vtkUnstructuredGrid()
    # The reader is local to this function. Own the points and connectivity so
    # later field attachment cannot outlive reader-managed buffers.
    grid.DeepCopy(source)
    return grid


def add_scalar(grid: vtkUnstructuredGrid, name: str, values: list[float]) -> tuple[float, float]:
    if len(values) != grid.GetNumberOfPoints():
        raise RuntimeError(
            f"field {name} has {len(values)} values for {grid.GetNumberOfPoints()} mesh points"
        )
    array = vtkDoubleArray()
    array.SetName(name)
    array.SetNumberOfComponents(1)
    array.SetNumberOfTuples(len(values))
    for index, value in enumerate(values):
        array.SetValue(index, value)
    grid.GetPointData().AddArray(array)
    return min(values), max(values)


def add_flow(grid: vtkUnstructuredGrid, result: Path) -> tuple[str, tuple[float, float], bool]:
    vectors: list[tuple[float, float, float]] = []
    with result.open("r", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            columns = line.split()
            if not columns:
                continue
            if len(columns) != 3:
                raise RuntimeError(f"{result}:{line_number}: expected three velocity columns")
            vectors.append(tuple(float(value) for value in columns))
    if len(vectors) != grid.GetNumberOfPoints():
        raise RuntimeError(
            f"velocity has {len(vectors)} rows for {grid.GetNumberOfPoints()} mesh points"
        )

    velocity = vtkDoubleArray()
    velocity.SetName("velocity")
    velocity.SetNumberOfComponents(3)
    velocity.SetNumberOfTuples(len(vectors))
    magnitudes: list[float] = []
    for index, vector in enumerate(vectors):
        velocity.SetTuple3(index, vector[0], vector[1], vector[2])
        magnitudes.append(math.sqrt(sum(component * component for component in vector)))
    grid.GetPointData().AddArray(velocity)
    grid.GetPointData().SetActiveVectors("velocity")

    pressure_path = Path(str(result) + ".pressure")
    if pressure_path.exists():
        pressure = [float(line) for line in pressure_path.read_text(encoding="utf-8").splitlines() if line]
        add_scalar(grid, "pressure", pressure)
    return "velocity", (min(magnitudes), max(magnitudes)), True


def add_transport(
    grid: vtkUnstructuredGrid, result: Path, selected: str
) -> tuple[str, tuple[float, float], bool]:
    fields_path = Path(str(result) + ".fields")
    fields = [line.strip() for line in fields_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not fields:
        raise RuntimeError(f"transport field list is empty: {fields_path}")
    if selected not in fields:
        raise RuntimeError(f"unknown transport array {selected}; available: {', '.join(fields)}")

    values = {name: [math.nan] * grid.GetNumberOfPoints() for name in fields}
    seen: set[int] = set()
    with result.open("r", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            columns = line.split()
            if not columns:
                continue
            if len(columns) != len(fields) + 1:
                raise RuntimeError(
                    f"{result}:{line_number}: expected node id plus {len(fields)} fields"
                )
            node = int(columns[0])
            if node < 0 or node >= grid.GetNumberOfPoints() or node in seen:
                raise RuntimeError(f"{result}:{line_number}: invalid or duplicate node id {node}")
            seen.add(node)
            for field, value in zip(fields, columns[1:]):
                values[field][node] = float(value)
    if len(seen) != grid.GetNumberOfPoints():
        raise RuntimeError(
            f"transport has {len(seen)} node rows for {grid.GetNumberOfPoints()} mesh points"
        )
    ranges = {name: add_scalar(grid, name, field_values) for name, field_values in values.items()}
    grid.GetPointData().SetActiveScalars(selected)
    return selected, ranges[selected], False


def render(
    grid: vtkUnstructuredGrid,
    array_name: str,
    value_range: tuple[float, float],
    vector_magnitude: bool,
    arguments: argparse.Namespace,
) -> None:
    _DisableFirstRenderCameraReset()
    producer = TrivialProducer(registrationName="TubularFlowIGA result")
    producer.GetClientSideObject().SetOutput(grid)
    producer.UpdatePipeline()

    bounds = grid.GetBounds()
    center = [0.5 * (bounds[2 * axis] + bounds[2 * axis + 1]) for axis in range(3)]
    spans = [bounds[2 * axis + 1] - bounds[2 * axis] for axis in range(3)]

    view = CreateView("RenderView")
    view.ViewSize = [arguments.width, arguments.height]
    view.Background = [1.0, 1.0, 1.0]
    view.UseColorPaletteForBackground = 0
    view.OrientationAxesVisibility = 0
    view.CameraParallelProjection = 1

    surface_display = Show(producer, view)
    surface_display.Representation = "Surface With Edges"
    surface_display.ColorArrayName = [None, ""]
    surface_display.DiffuseColor = [0.42, 0.47, 0.54]
    surface_display.EdgeColor = [0.22, 0.25, 0.30]
    surface_display.Opacity = 0.13
    surface_display.LineWidth = 0.6

    slice_filter = Slice(registrationName="Center plane", Input=producer)
    slice_filter.SliceType = "Plane"
    slice_filter.SliceType.Origin = center
    slice_filter.SliceType.Normal = [0.0, 0.0, 1.0]
    slice_display = Show(slice_filter, view)
    slice_display.Representation = "Surface"
    if vector_magnitude:
        ColorBy(slice_display, ("POINTS", array_name, "Magnitude"))
    else:
        ColorBy(slice_display, ("POINTS", array_name))

    lower, upper = value_range
    if not lower < upper:
        upper = lower + 1.0
    lookup = GetColorTransferFunction(array_name)
    preset = "Viridis (matplotlib)" if vector_magnitude else "Cool to Warm"
    try:
        lookup.ApplyPreset(preset, True)
    except RuntimeError:
        lookup.ApplyPreset("Cool to Warm", True)
    lookup.RescaleTransferFunction(lower, upper)
    opacity = GetOpacityTransferFunction(array_name)
    opacity.RescaleTransferFunction(lower, upper)
    slice_display.SetScalarBarVisibility(view, True)

    scalar_bar = GetScalarBar(lookup, view)
    scalar_bar.Title = arguments.legend or (
        "Velocity magnitude" if vector_magnitude else array_name
    )
    scalar_bar.ComponentTitle = ""
    scalar_bar.TitleColor = [0.08, 0.10, 0.14]
    scalar_bar.LabelColor = [0.08, 0.10, 0.14]
    scalar_bar.TitleFontSize = 18
    scalar_bar.LabelFontSize = 14
    scalar_bar.ScalarBarLength = 0.55
    scalar_bar.WindowLocation = "Lower Center"
    scalar_bar.Orientation = "Horizontal"

    title = Text(registrationName="Title")
    title.Text = arguments.title
    title_display = Show(title, view)
    title_display.WindowLocation = "Upper Center"
    title_display.FontSize = 24
    title_display.Bold = 1
    title_display.Color = [0.06, 0.08, 0.12]

    subtitle = Text(registrationName="Subtitle")
    subtitle.Text = arguments.subtitle
    subtitle_display = Show(subtitle, view)
    subtitle_display.WindowLocation = "Upper Left Corner"
    subtitle_display.FontSize = 13
    subtitle_display.Color = [0.28, 0.31, 0.36]

    # Materialize every representation before fixing the camera. Some ParaView
    # builds perform their first bounds update lazily during Render().
    Render(view)
    largest_span = max(spans)
    view.CameraFocalPoint = center
    view.CameraPosition = [center[0], center[1], center[2] + 3.0 * largest_span]
    view.CameraViewUp = [0.0, 1.0, 0.0]
    aspect = arguments.width / arguments.height
    view.CameraParallelScale = 0.70 * max(spans[1], spans[0] / aspect)

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    Render(view)
    SaveScreenshot(
        str(arguments.output),
        view,
        ImageResolution=[arguments.width, arguments.height],
        TransparentBackground=0,
    )


def main() -> None:
    arguments = parse_arguments()
    grid = read_grid(arguments.mesh)
    if arguments.kind == "flow":
        array_name, value_range, vector_magnitude = add_flow(grid, arguments.result)
    else:
        array_name, value_range, vector_magnitude = add_transport(
            grid, arguments.result, arguments.array
        )
    render(grid, array_name, value_range, vector_magnitude, arguments)
    print(
        f"rendered={arguments.output} points={grid.GetNumberOfPoints()} "
        f"cells={grid.GetNumberOfCells()} array={array_name} "
        f"range=[{value_range[0]:.17g},{value_range[1]:.17g}]"
    )


if __name__ == "__main__":
    main()
