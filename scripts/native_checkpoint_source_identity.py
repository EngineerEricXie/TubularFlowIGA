#!/usr/bin/env python3
"""Hash native graph implementation at build time, including dirty sources.

Paths are repository-relative; timestamps, checkout paths, and Git metadata do
not affect the result. The conservative header superset deliberately rejects
unverified source changes rather than claiming cross-version compatibility.
"""
import hashlib
from pathlib import Path


def source_identity(root):
    paths = [root / 'solvers/coupling/src/iga_1d_3d_bifurcation.cpp',
             root / 'solvers/coupling/Makefile',
             root / 'scripts/native_checkpoint_source_identity.py']
    for directory in ('include', 'solvers/cpu/include', 'solvers/one_d/include'):
        paths.extend((root / directory).glob('*.hpp'))
    digest = hashlib.sha256(b'IGA_NATIVE_GRAPH_SOURCE/1\n')
    for path in sorted(paths, key=lambda p: p.relative_to(root).as_posix()):
        for value in (path.relative_to(root).as_posix().encode(), path.read_bytes()):
            digest.update(len(value).to_bytes(8, 'little'))
            digest.update(value)
    return digest.hexdigest()


if __name__ == '__main__':
    print(source_identity(Path(__file__).resolve().parents[1]))
