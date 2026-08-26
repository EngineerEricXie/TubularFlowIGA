#!/usr/bin/env pvbatch

"""Render synchronized 1D or 3D multiphysics PVD files as an animated GIF."""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

from PIL import Image
from paraview.simple import (
    ColorBy,
    CreateLayout,
    CreateView,
    GetAnimationScene,
    GetColorTransferFunction,
    OpenDataFile,
    Render,
    SaveScreenshot,
    Show,
    Slice,
    Text,
    UpdatePipeline,
    _DisableFirstRenderCameraReset,
)


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dimension", choices=("1d", "3d"), required=True)
    parser.add_argument("--transport", type=Path, required=True)
    parser.add_argument("--flow", type=Path, help="3D flow PVD; omit for combined 1D PVD")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--title", required=True)
    parser.add_argument("--frames", type=int, default=12)
    parser.add_argument("--duration-ms", type=int, default=180)
    parser.add_argument("--width", type=int, default=1200)
    parser.add_argument("--height", type=int, default=720)
    return parser.parse_args()


def timesteps(source: object, maximum: int) -> list[float]:
    values = list(getattr(source, "TimestepValues", []) or [])
    if not values:
        values = [0.0]
    if maximum == 1:
        return [float(values[0])]
    if len(values) <= maximum:
        return [float(value) for value in values]
    selected = []
    for index in range(maximum):
        position = round(index * (len(values) - 1) / (maximum - 1))
        selected.append(float(values[position]))
    return selected


def display_input(source: object, dimension: str) -> object:
    if dimension == "1d":
        return source
    sliced = Slice(Input=source)
    sliced.SliceType = "Plane"
    sliced.SliceType.Normal = [0.0, 0.0, 1.0]
    return sliced


def add_panel(
    layout: object, slot: int, source: object, array: str, title: str, vector_velocity: bool
) -> object:
    view = CreateView("RenderView")
    view.Background = [0.035, 0.045, 0.07]
    view.UseColorPaletteForBackground = 0
    view.OrientationAxesVisibility = 0
    view.CameraParallelProjection = 1
    layout.AssignView(slot, view)
    shown = Show(source, view)
    shown.Representation = "Surface"
    shown.LineWidth = 9.0
    if array == "velocity" and vector_velocity:
        ColorBy(shown, ("POINTS", array, "Magnitude"))
    else:
        ColorBy(shown, ("POINTS", array))
    lookup = GetColorTransferFunction(array)
    try:
        lookup.ApplyPreset("Viridis (matplotlib)", True)
    except RuntimeError:
        lookup.ApplyPreset("Cool to Warm", True)
    try:
        shown.RescaleTransferFunctionToDataRangeOverTime()
    except AttributeError:
        shown.RescaleTransferFunctionToDataRange(True, False)
    shown.SetScalarBarVisibility(view, True)
    label = Text(Text=title)
    label_display = Show(label, view)
    label_display.Color = [0.96, 0.97, 1.0]
    label_display.FontSize = 16
    label_display.WindowLocation = "Upper Center"
    UpdatePipeline(proxy=source)
    bounds = source.GetDataInformation().GetBounds()
    center = [0.5 * (bounds[2 * axis] + bounds[2 * axis + 1]) for axis in range(3)]
    span = max(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4])
    view.CameraFocalPoint = center
    view.CameraPosition = [center[0], center[1], center[2] + 3.0 * max(span, 1.0e-12)]
    view.CameraViewUp = [0.0, 1.0, 0.0]
    return view


def main() -> None:
    args = arguments()
    if args.frames < 1 or args.duration_ms < 1:
        raise ValueError("frames and duration-ms must be positive")
    if args.dimension == "3d" and args.flow is None:
        raise ValueError("3D rendering requires --flow")
    _DisableFirstRenderCameraReset()
    transport = OpenDataFile(str(args.transport))
    flow = OpenDataFile(str(args.flow)) if args.flow else transport
    transport_input = display_input(transport, args.dimension)
    flow_input = display_input(flow, args.dimension)

    layout = CreateLayout(name="TubularFlowIGA multiphysics")
    layout.SplitHorizontal(0, 0.5)
    layout.SplitVertical(1, 0.5)
    layout.SplitVertical(2, 0.5)
    panels = [
        add_panel(layout, 3, flow_input, "velocity", "Velocity", args.dimension == "3d"),
        add_panel(layout, 4, transport_input, "oxygen", "Oxygen", False),
        add_panel(layout, 5, transport_input, "glucose", "Glucose", False),
        add_panel(layout, 6, transport_input, "lactate", "Lactate", False),
    ]
    layout.SetSize(args.width, args.height)
    scene = GetAnimationScene()
    times = timesteps(transport, args.frames)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    images: list[Image.Image] = []
    with tempfile.TemporaryDirectory(prefix="tubularflowiga-gif-") as temporary:
        temporary_path = Path(temporary)
        for index, time in enumerate(times):
            scene.TimeKeeper.Time = time
            UpdatePipeline(time=time, proxy=transport)
            if flow is not transport:
                UpdatePipeline(time=time, proxy=flow)
            for view in panels:
                view.ViewTime = time
                view.ResetCamera()
                Render(view)
            frame = temporary_path / f"frame-{index:04d}.png"
            SaveScreenshot(str(frame), layout, ImageResolution=[args.width, args.height])
            with Image.open(frame) as image:
                images.append(image.convert("P", palette=Image.ADAPTIVE))
    images[0].save(
        args.output,
        save_all=True,
        append_images=images[1:],
        duration=args.duration_ms,
        loop=0,
        optimize=True,
    )
    print(f"rendered {len(images)} frames to {args.output}")


if __name__ == "__main__":
    main()
