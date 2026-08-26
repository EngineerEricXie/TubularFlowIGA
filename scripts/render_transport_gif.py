#!/usr/bin/env pvbatch

"""Render one scalar array from a ParaView PVD collection as an animated GIF."""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

from PIL import Image
from paraview.simple import (
    ColorBy,
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
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--array", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--title", required=True)
    parser.add_argument("--frames", type=int, default=13)
    parser.add_argument("--duration-ms", type=int, default=180)
    parser.add_argument("--width", type=int, default=1000)
    parser.add_argument("--height", type=int, default=620)
    return parser.parse_args()


def selected_times(source: object, maximum: int) -> list[float]:
    values = [float(value) for value in (getattr(source, "TimestepValues", []) or [0.0])]
    if maximum == 1:
        return [values[0]]
    if len(values) <= maximum:
        return values
    return [
        values[round(index * (len(values) - 1) / (maximum - 1))]
        for index in range(maximum)
    ]


def main() -> None:
    args = arguments()
    if args.frames < 1 or args.duration_ms < 1:
        raise ValueError("frames and duration-ms must be positive")

    _DisableFirstRenderCameraReset()
    source = OpenDataFile(str(args.input))
    UpdatePipeline(proxy=source)
    bounds = source.GetDataInformation().GetBounds()
    center = [0.5 * (bounds[2 * axis] + bounds[2 * axis + 1]) for axis in range(3)]

    sliced = Slice(Input=source)
    sliced.SliceType = "Plane"
    sliced.SliceType.Origin = center
    sliced.SliceType.Normal = [0.0, 0.0, 1.0]

    view = CreateView("RenderView")
    view.ViewSize = [args.width, args.height]
    view.Background = [0.035, 0.045, 0.07]
    view.UseColorPaletteForBackground = 0
    view.OrientationAxesVisibility = 0
    view.CameraParallelProjection = 1

    shown = Show(sliced, view)
    shown.Representation = "Surface"
    ColorBy(shown, ("POINTS", args.array))
    lookup = GetColorTransferFunction(args.array)
    try:
        lookup.ApplyPreset("Viridis (matplotlib)", True)
    except RuntimeError:
        lookup.ApplyPreset("Cool to Warm", True)
    try:
        shown.RescaleTransferFunctionToDataRangeOverTime()
    except AttributeError:
        shown.RescaleTransferFunctionToDataRange(True, False)
    shown.SetScalarBarVisibility(view, True)

    title = Text(Text=args.title)
    title_display = Show(title, view)
    title_display.Color = [0.96, 0.97, 1.0]
    title_display.FontSize = 20
    title_display.WindowLocation = "Upper Center"

    view.CameraFocalPoint = center
    span = max(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4])
    view.CameraPosition = [center[0], center[1], center[2] + 3.0 * max(span, 1.0e-12)]
    view.CameraViewUp = [0.0, 1.0, 0.0]
    view.ResetCamera()

    scene = GetAnimationScene()
    times = selected_times(source, args.frames)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    images: list[Image.Image] = []
    with tempfile.TemporaryDirectory(prefix="tubularflowiga-transport-gif-") as temporary:
        temporary_path = Path(temporary)
        for index, time in enumerate(times):
            scene.TimeKeeper.Time = time
            UpdatePipeline(time=time, proxy=source)
            view.ViewTime = time
            Render(view)
            frame = temporary_path / f"frame-{index:04d}.png"
            SaveScreenshot(str(frame), view, ImageResolution=[args.width, args.height])
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
