#!/usr/bin/env python3
"""Assign tissue cells to the nearest terminal of each vascular tree.

For every tree, each tetrahedron goes to the terminal (SWC leaf other than the
root) closest to its centroid; weights within a territory are cell volume over
territory volume. Terminals that win no cell get no territory and carry no
flow. This Voronoi rule is organ-agnostic; it ignores tissue geodesics, so for
strongly non-convex organs check the reported distances.

    python3 scripts/organ/build_territories.py liver_h5.msh territories.json \
        --tree portal=portal.swc --tree hepatic=hepatic.swc
"""
import argparse
import json

import numpy as np
from scipy.spatial import cKDTree


def read_gmsh_tets(path):
    """Minimal ASCII Gmsh 4.1 reader: node coordinates and tetrahedra (type 4)."""
    with open(path) as handle:
        lines = iter(handle.read().splitlines())
    nodes, tets, tet_ids = {}, [], []
    for line in lines:
        if line == "$Nodes":
            blocks, _, _, _ = map(int, next(lines).split())
            for _ in range(blocks):
                _, _, parametric, count = map(int, next(lines).split())
                if parametric:
                    raise SystemExit("parametric Gmsh nodes are not supported")
                tags = [int(next(lines)) for _ in range(count)]
                for tag in tags:
                    nodes[tag] = [float(v) for v in next(lines).split()[:3]]
        elif line == "$Elements":
            blocks, _, _, _ = map(int, next(lines).split())
            for _ in range(blocks):
                _, _, kind, count = map(int, next(lines).split())
                for _ in range(count):
                    values = list(map(int, next(lines).split()))
                    if kind == 4:
                        tet_ids.append(values[0])
                        tets.append(values[1:5])
    tags = np.array(sorted(nodes))
    index = {tag: i for i, tag in enumerate(tags)}
    points = np.array([nodes[t] for t in tags])
    tets = np.array([[index[n] for n in tet] for tet in tets])
    return points, tets, np.array(tet_ids)


def read_swc(path):
    rows = [line.split() for line in open(path) if line.strip() and not line.startswith("#")]
    ids = [int(r[0]) for r in rows]
    xyz = {int(r[0]): np.array([float(v) for v in r[2:5]]) for r in rows}
    parents = {int(r[0]): int(r[6]) for r in rows}
    has_child = {p for p in parents.values() if p != -1}
    roots = [i for i in ids if parents[i] == -1]
    if len(roots) != 1:
        raise SystemExit(f"{path}: expected one SWC root")
    terminals = [i for i in ids if i not in has_child and i != roots[0]]
    return terminals, np.array([xyz[i] for i in terminals])


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("mesh")
    parser.add_argument("output")
    parser.add_argument("--tree", action="append", required=True, help="ID=SWC (same length unit as the mesh)")
    parser.add_argument("--summary")
    args = parser.parse_args()
    points, tets, tet_ids = read_gmsh_tets(args.mesh)
    corners = points[tets]
    centroids = corners.mean(axis=1)
    edges = corners[:, 1:] - corners[:, :1]
    volumes = np.abs(np.einsum("ij,ij->i", edges[:, 0], np.cross(edges[:, 1], edges[:, 2]))) / 6.0
    result = {"schema_version": 1, "rule": "nearest_terminal_volume_weighted", "trees": {}}
    summary = {"cells": int(len(tets)), "trees": {}}
    for item in args.tree:
        name, path = item.split("=", 1)
        terminals, positions = read_swc(path)
        distance, nearest = cKDTree(positions).query(centroids)
        territories = {}
        sizes = []
        for k, terminal in enumerate(terminals):
            members = np.flatnonzero(nearest == k)
            if len(members) == 0:
                continue
            weights = volumes[members] / volumes[members].sum()
            territories[str(terminal)] = {"cell_ids": tet_ids[members].tolist(),
                                          "weights": [float(w) for w in weights]}
            sizes.append(float(volumes[members].sum()))
        result["trees"][name] = territories
        summary["trees"][name] = {
            "terminals": len(terminals), "terminals_with_cells": len(territories),
            "territory_volume_percentiles": np.percentile(sizes, [0, 50, 100]).tolist(),
            "cell_to_terminal_distance_percentiles": np.percentile(distance, [50, 95, 100]).tolist(),
        }
    with open(args.output, "w") as out:
        json.dump(result, out)
    if args.summary:
        json.dump(summary, open(args.summary, "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
