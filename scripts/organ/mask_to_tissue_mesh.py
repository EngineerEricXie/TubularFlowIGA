#!/usr/bin/env python3
"""Turn an organ mask into a coarse, conforming tetrahedral tissue mesh.

The mask is block-averaged onto cubes of edge `--size-mm`; a cube is kept when
at least `--occupancy` of it lies inside the organ, and only the largest
face-connected set of cubes is kept. Every cube is split into the six Kuhn
tetrahedra sharing its main diagonal, which conforms across neighbours and
keeps every Jacobian positive. All boundary faces get `boundary_label_1`; the
volume is named `fluid` as the native tetra reader requires. Coordinates are
written in millimetres (use `mesh_length_scale_to_m: 0.001`).

`--fill NAME` adds other masks (e.g. a small tumour to be modelled as
parenchyma) to the organ before meshing; record that substitution.

    python3 scripts/organ/mask_to_tissue_mesh.py masks.npz liver_h5.msh \
        --mask Liver --fill Tumor_1 --size-mm 5 --manifest liver_h5.json
"""
import argparse
import hashlib
import json

import numpy as np
from scipy import ndimage

KUHN = ((0, 1, 3, 7), (0, 3, 2, 7), (0, 2, 6, 7), (0, 6, 4, 7), (0, 4, 5, 7), (0, 5, 1, 7))
CORNERS = [(dx, dy, dz) for dz in (0, 1) for dy in (0, 1) for dx in (0, 1)]  # bit order x, y, z


def block_occupancy(mask, factors):
    padded = [int(np.ceil(n / f)) * f for n, f in zip(mask.shape, factors)]
    grid = np.zeros(padded, dtype=np.float32)
    grid[:mask.shape[0], :mask.shape[1], :mask.shape[2]] = mask
    shape = (padded[0] // factors[0], factors[0], padded[1] // factors[1], factors[1],
             padded[2] // factors[2], factors[2])
    return grid.reshape(shape).mean(axis=(1, 3, 5))


def write_gmsh(path, points, triangles, tetrahedra):
    lower, upper = points.min(axis=0), points.max(axis=0)
    box = " ".join(f"{v:.9g}" for v in (*lower, *upper))
    with open(path, "w") as out:
        out.write("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n")
        out.write('$PhysicalNames\n2\n2 1 "boundary_label_1"\n3 2 "fluid"\n$EndPhysicalNames\n')
        out.write(f"$Entities\n0 0 1 1\n1 {box} 1 1 0\n1 {box} 1 2 1 1\n$EndEntities\n")
        n = len(points)
        out.write(f"$Nodes\n1 {n} 1 {n}\n3 1 0 {n}\n")
        out.write("".join(f"{i}\n" for i in range(1, n + 1)))
        out.write("".join(f"{x:.9g} {y:.9g} {z:.9g}\n" for x, y, z in points))
        out.write("$EndNodes\n")
        total = len(triangles) + len(tetrahedra)
        out.write(f"$Elements\n2 {total} 1 {total}\n2 1 2 {len(triangles)}\n")
        out.write("".join(f"{i} {a + 1} {b + 1} {c + 1}\n" for i, (a, b, c) in enumerate(triangles, start=1)))
        out.write(f"3 1 4 {len(tetrahedra)}\n")
        start = len(triangles) + 1
        out.write("".join(f"{i} {a + 1} {b + 1} {c + 1} {d + 1}\n"
                          for i, (a, b, c, d) in enumerate(tetrahedra, start=start)))
        out.write("$EndElements\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("masks")
    parser.add_argument("output")
    parser.add_argument("--mask", required=True)
    parser.add_argument("--fill", action="append", default=[])
    parser.add_argument("--size-mm", type=float, required=True)
    parser.add_argument("--occupancy", type=float, default=0.5)
    parser.add_argument("--manifest")
    args = parser.parse_args()
    data = np.load(args.masks)
    mask = data[args.mask].astype(bool)
    for name in args.fill:
        mask |= data[name].astype(bool)
    spacing = data["spacing_zyx_mm"].astype(float)
    origin = data["origin_lps_mm"].astype(float)  # x, y, z of voxel (0, 0, 0)
    factors = [max(int(round(args.size_mm / s)), 1) for s in spacing]
    block = np.array(factors) * spacing  # z, y, x edge lengths in mm
    occupancy = block_occupancy(mask, factors)
    keep = occupancy >= args.occupancy
    labels, count = ndimage.label(keep)  # 6-connectivity
    if count == 0:
        raise SystemExit("no block passes the occupancy threshold")
    sizes = np.bincount(labels.ravel())[1:]
    keep = labels == int(np.argmax(sizes)) + 1
    cubes = np.argwhere(keep)  # (k, j, i)
    # Block corners in voxel-corner coordinates: voxel (0,0,0) centre is the origin.
    corner_origin = origin - 0.5 * spacing[::-1]
    node_id = {}
    points = []
    tetrahedra = []
    for k, j, i in cubes:
        ids = []
        for dx, dy, dz in CORNERS:
            key = (i + dx, j + dy, k + dz)
            if key not in node_id:
                node_id[key] = len(points)
                points.append(corner_origin + np.array(key) * block[::-1])
            ids.append(node_id[key])
        tetrahedra.extend(tuple(ids[c] for c in pattern) for pattern in KUHN)
    points = np.array(points)
    tetrahedra = np.array(tetrahedra)
    faces = {}
    for tet in tetrahedra:
        for opposite in range(4):
            face = tuple(sorted(int(tet[n]) for n in range(4) if n != opposite))
            faces[face] = faces.get(face, 0) + 1
    triangles = [face for face, uses in faces.items() if uses == 1]
    write_gmsh(args.output, points, triangles, tetrahedra)
    volume_ml = len(cubes) * float(np.prod(block)) / 1000.0
    source_ml = float(mask.sum() * np.prod(spacing)) / 1000.0
    summary = {
        "mask": args.mask, "fill": args.fill, "size_mm": args.size_mm,
        "block_edges_xyz_mm": block[::-1].tolist(), "occupancy_threshold": args.occupancy,
        "cubes": int(len(cubes)), "nodes": int(len(points)), "tetrahedra": int(len(tetrahedra)),
        "boundary_triangles": int(len(triangles)), "mesh_volume_ml": volume_ml, "mask_volume_ml": source_ml,
        "discarded_components": int(count - 1),
        "mesh_sha256": hashlib.sha256(open(args.output, "rb").read()).hexdigest(),
        "masks_sha256": hashlib.sha256(open(args.masks, "rb").read()).hexdigest(),
        "length_unit": "mm",
    }
    if args.manifest:
        json.dump(summary, open(args.manifest, "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
