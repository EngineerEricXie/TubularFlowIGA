#!/usr/bin/env python3
"""Vessel mask -> rooted, radius-annotated, smoothed centreline tree (SWC, LPS mm).

Steps: keep the largest 26-connected component; skeletonize; connect skeleton
voxels into a graph and keep the maximum-radius spanning tree; prune terminal
spurs shorter than `--prune-factor` times their junction radius; pick the root
as the terminal with the largest radius near it (the trunk: portal hilum or
hepatic-vein ostium); collapse degree-2 chains into branches, smooth each
branch, resample it every `--spacing-mm`, and drop terminal branches whose
median radius is below `--min-radius-mm` (unreliable at voxel resolution).
Radii are the Euclidean distance to the mask boundary. Requires numpy, scipy,
scikit-image and networkx.

    python3 scripts/organ/vessel_mask_to_swc.py masks.npz portal.swc --mask Portal \
        --spacing-mm 2 --min-radius-mm 1.0 --summary portal.json
"""
import argparse
import hashlib
import json

import networkx as nx
import numpy as np
from scipy import ndimage
from skimage.morphology import skeletonize


def voxel_tree(mask, spacing, origin, prune_factor):
    labels, count = ndimage.label(mask, structure=np.ones((3, 3, 3)))
    sizes = np.bincount(labels.ravel())[1:]
    mask = labels == int(np.argmax(sizes)) + 1
    distance = ndimage.distance_transform_edt(mask, sampling=spacing)
    skeleton = skeletonize(mask)
    voxels = [tuple(v) for v in np.argwhere(skeleton)]
    position = {v: np.array([origin[0] + v[2] * spacing[2], origin[1] + v[1] * spacing[1],
                             origin[2] + v[0] * spacing[0]]) for v in voxels}
    graph = nx.Graph()
    for v in voxels:
        graph.add_node(v, xyz=position[v], r=float(distance[v]))
    offsets = [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1) if (a, b, c) != (0, 0, 0)]
    for v in voxels:
        for o in offsets:
            w = (v[0] + o[0], v[1] + o[1], v[2] + o[2])
            if w in position:
                graph.add_edge(v, w, length=float(np.linalg.norm(position[v] - position[w])),
                               weight=-(float(distance[v]) + float(distance[w])))
    graph = graph.subgraph(max(nx.connected_components(graph), key=len)).copy()
    tree = nx.minimum_spanning_tree(graph, weight="weight")
    changed = True
    while changed:
        changed = False
        for leaf in [n for n in tree if tree.degree(n) == 1]:
            if leaf not in tree:
                continue
            path, node, length = [leaf], leaf, 0.0
            while True:
                nbrs = [n for n in tree.neighbors(node) if n not in path]
                if len(nbrs) != 1:
                    break
                nxt = nbrs[0]
                length += tree.edges[node, nxt]["length"]
                path.append(nxt)
                node = nxt
                if tree.degree(node) != 2:
                    break
            junction = path[-1]
            if tree.degree(junction) >= 3 and length < prune_factor * tree.nodes[junction]["r"]:
                tree.remove_nodes_from(path[:-1])
                changed = True
    return tree, int(count)


def root_terminal(tree, steps=6):
    def mean_radius(n):
        seen, current, radii = {n}, n, [tree.nodes[n]["r"]]
        for _ in range(steps):
            nxt = [m for m in tree.neighbors(current) if m not in seen]
            if not nxt:
                break
            current = nxt[0]
            seen.add(current)
            radii.append(tree.nodes[current]["r"])
        return float(np.mean(radii))
    return max((n for n in tree if tree.degree(n) == 1), key=mean_radius)


def branches_from(tree, root):
    """Root-first list of branches (node paths between root/junction/terminal nodes)."""
    key = {n for n in tree if n == root or tree.degree(n) != 2}
    parent = {root: None}
    for a, b in nx.bfs_edges(tree, root):
        parent[b] = a
    children = {n: [] for n in tree}
    for n, p in parent.items():
        if p is not None:
            children[p].append(n)
    branches = []
    queue = [root]
    while queue:
        start = queue.pop(0)
        for child in children[start]:
            path = [start, child]
            while path[-1] not in key:
                path.append(children[path[-1]][0])
            branches.append(path)
            queue.append(path[-1])
    return branches, children


def smooth(values, width):
    if len(values) < 3 or width < 2:
        return values
    kernel = np.ones(width) / width
    padded = np.pad(values, ((width // 2, width // 2), (0, 0)), mode="edge")
    result = np.vstack([np.convolve(padded[:, j], kernel, mode="valid")[:len(values)]
                        for j in range(values.shape[1])]).T
    result[0], result[-1] = values[0], values[-1]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("masks")
    parser.add_argument("output")
    parser.add_argument("--mask", required=True)
    parser.add_argument("--prune-factor", type=float, default=2.0)
    parser.add_argument("--spacing-mm", type=float, default=2.0)
    parser.add_argument("--min-radius-mm", type=float, default=0.0)
    parser.add_argument("--radius-floor-mm", type=float, default=0.5)
    parser.add_argument("--summary")
    args = parser.parse_args()
    data = np.load(args.masks)
    spacing = data["spacing_zyx_mm"].astype(float)
    origin = data["origin_lps_mm"].astype(float)
    tree, components = voxel_tree(data[args.mask].astype(bool), spacing, origin, args.prune_factor)
    root = root_terminal(tree)
    branches, _ = branches_from(tree, root)
    # Drop thin terminal branches repeatedly, then merge pass-through junctions.
    while args.min_radius_mm > 0:
        ends = {b[-1] for b in branches}
        starts = {b[0] for b in branches}
        thin = [b for b in branches if b[-1] not in starts
                and np.median([tree.nodes[n]["r"] for n in b[1:]]) < args.min_radius_mm]
        if not thin:
            break
        for b in thin:
            branches.remove(b)
        merged = True
        while merged:
            merged = False
            starts_count = {}
            for b in branches:
                starts_count[b[0]] = starts_count.get(b[0], 0) + 1
            for b in list(branches):
                end = b[-1]
                if end != root and starts_count.get(end) == 1:
                    follower = next(c for c in branches if c[0] == end)
                    branches.remove(b)
                    branches.remove(follower)
                    branches.append(b + follower[1:])
                    merged = True
                    break
        del ends
    order = []
    queue = [root]
    while queue:
        start = queue.pop(0)
        for b in [b for b in branches if b[0] == start]:
            order.append(b)
            queue.append(b[-1])
    rows = [(1, tree.nodes[root]["xyz"], tree.nodes[root]["r"], -1)]
    swc_id = {root: 1}
    for b in order:
        points = smooth(np.array([tree.nodes[n]["xyz"] for n in b]), 5)
        radii = smooth(np.array([[tree.nodes[n]["r"]] for n in b]), 7)[:, 0]
        arc = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
        samples = max(int(round(arc[-1] / args.spacing_mm)), 1)
        s = np.linspace(0.0, arc[-1], samples + 1)
        xyz = np.vstack([np.interp(s, arc, points[:, j]) for j in range(3)]).T
        rr = np.interp(s, arc, radii)
        parent = swc_id[b[0]]
        for q in range(1, len(s)):
            rows.append((len(rows) + 1, xyz[q], rr[q], parent))
            parent = len(rows)
        swc_id[b[-1]] = parent
    with open(args.output, "w") as out:
        out.write(f"# {args.mask} centreline from {args.masks}; LPS mm; radius mm (distance to mask boundary)\n")
        for index, xyz, radius, parent in rows:
            kind = 1 if parent == -1 else 3
            out.write(f"{index} {kind} {xyz[0]:.4f} {xyz[1]:.4f} {xyz[2]:.4f} "
                      f"{max(radius, args.radius_floor_mm):.4f} {parent}\n")
    has_child = {parent for _, _, _, parent in rows}
    terminals = [index for index, _, _, _ in rows if index not in has_child and index != 1]
    summary = {
        "mask": args.mask, "mask_components": components, "swc_points": len(rows),
        "branches": len(order), "terminals": len(terminals),
        "root_xyz_mm": rows[0][1].tolist(), "root_radius_mm": float(rows[0][2]),
        "total_length_mm": float(sum(np.linalg.norm(rows[p - 1][1] - xyz) for _, xyz, _, p in rows if p > 0)),
        "prune_factor": args.prune_factor, "spacing_mm": args.spacing_mm, "min_radius_mm": args.min_radius_mm,
        "swc_sha256": hashlib.sha256(open(args.output, "rb").read()).hexdigest(),
    }
    if args.summary:
        json.dump(summary, open(args.summary, "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
