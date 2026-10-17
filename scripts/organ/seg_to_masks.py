#!/usr/bin/env python3
"""Extract named segments of a DICOM SEG into boolean masks on one LPS grid.

Writes an .npz with one boolean array per requested segment (axis order z, y, x),
`spacing_zyx_mm`, `origin_lps_mm` (centre of the first voxel) and the segment
names. Only BINARY axial SEGs whose image orientation (shared or per frame) is
the identity are accepted; anything else fails closed rather than being
resampled or thresholded silently.

    python3 scripts/organ/seg_to_masks.py seg.dcm masks.npz --segment Liver \
        --segment Portal --segment Hepatic --segment Tumor_1
"""
import argparse
import json
import sys

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("seg")
    parser.add_argument("output")
    parser.add_argument("--segment", action="append", required=True,
                        help="segment label to export (repeatable)")
    args = parser.parse_args()
    import pydicom

    data = pydicom.dcmread(args.seg)
    if getattr(data, "SegmentationType", None) != "BINARY":
        sys.exit(f"unsupported SEG type {getattr(data, 'SegmentationType', None)!r}; only BINARY is accepted")
    shared = data.SharedFunctionalGroupsSequence[0]
    # The orientation may be shared or stored per frame; every frame must be identity axial.
    orientations = []
    if "PlaneOrientationSequence" in shared:
        orientations.append(shared.PlaneOrientationSequence[0].ImageOrientationPatient)
    else:
        for frame in data.PerFrameFunctionalGroupsSequence:
            if "PlaneOrientationSequence" not in frame:
                sys.exit("SEG frame has no PlaneOrientationSequence")
            orientations.append(frame.PlaneOrientationSequence[0].ImageOrientationPatient)
    if not orientations:
        sys.exit("SEG has no plane orientation")
    for values in orientations:
        orientation = [float(v) for v in values]
        if not np.allclose(orientation, [1, 0, 0, 0, 1, 0], atol=1e-6):
            sys.exit(f"unsupported SEG orientation {orientation}; only identity axial frames are accepted")
    row_spacing, column_spacing = (float(v) for v in shared.PixelMeasuresSequence[0].PixelSpacing)
    labels = {int(item.SegmentNumber): str(item.SegmentLabel) for item in data.SegmentSequence}
    numbers = {}
    for name in args.segment:
        found = [number for number, label in labels.items() if label == name]
        if len(found) != 1:
            sys.exit(f"segment {name!r} not found exactly once; available: {sorted(labels.values())}")
        numbers[name] = found[0]
    pixels = np.asarray(data.pixel_array)
    frames = data.PerFrameFunctionalGroupsSequence
    positions = {}
    for index, frame in enumerate(frames):
        position = tuple(float(v) for v in frame.PlanePositionSequence[0].ImagePositionPatient)
        number = int(frame.SegmentIdentificationSequence[0].ReferencedSegmentNumber)
        positions.setdefault(number, {})[position] = index
    all_positions = sorted({p for per in positions.values() for p in per}, key=lambda p: p[2])
    xy = {p[:2] for p in all_positions}
    if len(xy) != 1:
        sys.exit("SEG frames do not share one in-plane origin")
    z_values = np.array(sorted({p[2] for p in all_positions}))
    steps = np.diff(z_values)
    if len(steps) and (steps.max() - steps.min()) > 1e-3 * steps.mean():
        sys.exit(f"SEG frames are not uniformly spaced in z ({steps.min()}..{steps.max()} mm)")
    z_spacing = float(steps.mean()) if len(steps) else 1.0
    shape = (len(z_values), pixels.shape[1], pixels.shape[2])
    masks = {}
    for name, number in numbers.items():
        mask = np.zeros(shape, dtype=bool)
        for position, index in positions.get(number, {}).items():
            k = int(np.argmin(np.abs(z_values - position[2])))
            mask[k] |= pixels[index].astype(bool)
        masks[name] = mask
    origin = np.array([next(iter(xy))[0], next(iter(xy))[1], z_values[0]])
    np.savez_compressed(args.output, spacing_zyx_mm=np.array([z_spacing, row_spacing, column_spacing]),
                        origin_lps_mm=origin, names=np.array(list(masks)), **masks)
    voxel = z_spacing * row_spacing * column_spacing
    print(json.dumps({"shape_zyx": shape, "spacing_zyx_mm": [z_spacing, row_spacing, column_spacing],
                      "volumes_ml": {name: float(mask.sum() * voxel / 1000.0) for name, mask in masks.items()}},
                     indent=1))


if __name__ == "__main__":
    main()
