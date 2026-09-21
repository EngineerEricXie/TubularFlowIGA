#!/usr/bin/env pvbatch

"""Render transport fields from ParaView PVD collections as an animated GIF."""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

from PIL import Image
from paraview.simple import (
    Clip,
    ColorBy,
    CreateView,
    GetAnimationScene,
    GetColorTransferFunction,
    GetScalarBar,
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
    parser.add_argument(
        "--overlay-input",
        type=Path,
        action="append",
        default=[],
        help="additional synchronized PVD collection to show in the same view",
    )
    parser.add_argument(
        "--overlay-color",
        type=float,
        action="append",
        nargs=3,
        default=[],
        metavar=("R", "G", "B"),
        help="solid RGB color for the corresponding overlay domain",
    )
    parser.add_argument("--array", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--title", required=True)
    parser.add_argument("--frames", type=int, default=13)
    parser.add_argument("--duration-ms", type=int, default=180)
    parser.add_argument("--width", type=int, default=1000)
    parser.add_argument("--height", type=int, default=620)
    parser.add_argument(
        "--cutaway",
        action="store_true",
        help="clip the primary domain at its x-axis midpoint",
    )
    parser.add_argument(
        "--range",
        type=float,
        nargs=2,
        metavar=("MIN", "MAX"),
        help="fixed scalar range shared by the primary and overlay domains",
    )
    parser.add_argument(
        "--log-scale",
        action="store_true",
        help="use logarithmic scalar mapping; requires a positive fixed range",
    )
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
    if args.range and not args.range[0] < args.range[1]:
        raise ValueError("range must be strictly increasing")
    if args.log_scale and (not args.range or args.range[0] <= 0.0):
        raise ValueError("log-scale requires a positive fixed range")
    if args.overlay_color and len(args.overlay_color) != len(args.overlay_input):
        raise ValueError("provide one overlay-color for each overlay-input")
    if any(component < 0.0 or component > 1.0
           for color in args.overlay_color for component in color):
        raise ValueError("overlay-color components must be between zero and one")

    _DisableFirstRenderCameraReset()
    source = OpenDataFile(str(args.input))
    UpdatePipeline(proxy=source)
    bounds = source.GetDataInformation().GetBounds()
    center = [0.5 * (bounds[2 * axis] + bounds[2 * axis + 1]) for axis in range(3)]

    if args.cutaway:
        visible_domain = Slice(Input=source)
        visible_domain.SliceType = "Plane"
        visible_domain.SliceType.Origin = center
        visible_domain.SliceType.Normal = [1.0, 0.0, 0.0]
        context_domain = Clip(Input=source)
        context_domain.ClipType = "Plane"
        context_domain.ClipType.Origin = center
        context_domain.ClipType.Normal = [1.0, 0.0, 0.0]
        context_domain.Invert = 1
    else:
        visible_domain = Slice(Input=source)
        visible_domain.SliceType = "Plane"
        visible_domain.SliceType.Origin = center
        visible_domain.SliceType.Normal = [0.0, 0.0, 1.0]
        context_domain = None

    overlays = [OpenDataFile(str(path)) for path in args.overlay_input]
    for overlay in overlays:
        UpdatePipeline(proxy=overlay)

    view = CreateView("RenderView")
    view.ViewSize = [args.width, args.height]
    view.Background = [0.035, 0.045, 0.07]
    view.UseColorPaletteForBackground = 0
    view.OrientationAxesVisibility = 0
    view.CameraParallelProjection = 1

    shown = Show(visible_domain, view)
    shown.Representation = "Surface"
    ColorBy(shown, ("POINTS", args.array))
    if context_domain:
        context_display = Show(context_domain, view)
        context_display.Representation = "Surface"
        ColorBy(context_display, ("POINTS", args.array))
        context_display.Opacity = 0.12
    lookup = GetColorTransferFunction(args.array)
    try:
        lookup.ApplyPreset("Viridis (matplotlib)", True)
    except RuntimeError:
        lookup.ApplyPreset("Cool to Warm", True)
    for index, overlay in enumerate(overlays):
        overlay_display = Show(overlay, view)
        overlay_display.Representation = "Surface"
        if args.overlay_color:
            overlay_display.AmbientColor = args.overlay_color[index]
            overlay_display.DiffuseColor = args.overlay_color[index]
            overlay_display.Specular = 0.25
        else:
            ColorBy(overlay_display, ("POINTS", args.array))
    if args.range:
        lookup.RescaleTransferFunction(*args.range)
    else:
        try:
            shown.RescaleTransferFunctionToDataRangeOverTime()
        except AttributeError:
            shown.RescaleTransferFunctionToDataRange(True, False)
    if args.log_scale:
        lookup.MapControlPointsToLogSpace()
        lookup.UseLogScale = 1
    if args.cutaway:
        shown.Opacity = 0.78
    shown.SetScalarBarVisibility(view, True)
    legend = GetScalarBar(lookup, view)
    legend.TitleColor = [0.96, 0.97, 1.0]
    legend.LabelColor = [0.96, 0.97, 1.0]

    title = Text(Text=args.title)
    title_display = Show(title, view)
    title_display.Color = [0.96, 0.97, 1.0]
    title_display.FontSize = 20
    title_display.WindowLocation = "Upper Center"

    time_label = Text(Text="")
    time_display = Show(time_label, view)
    time_display.Color = [0.96, 0.97, 1.0]
    time_display.FontSize = 16
    time_display.WindowLocation = "Upper Left Corner"

    view.CameraFocalPoint = center
    span = max(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4])
    if args.cutaway:
        view.CameraPosition = [
            center[0] + 2.8 * span,
            center[1] + 2.0 * span,
            center[2] + 1.8 * span,
        ]
        view.CameraViewUp = [0.0, 0.0, 1.0]
    else:
        view.CameraPosition = [
            center[0],
            center[1],
            center[2] + 3.0 * max(span, 1.0e-12),
        ]
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
            for overlay in overlays:
                UpdatePipeline(time=time, proxy=overlay)
            time_label.Text = f"t = {time:g} s"
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
