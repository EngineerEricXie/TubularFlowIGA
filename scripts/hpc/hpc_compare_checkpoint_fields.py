#!/usr/bin/env python3
"""Compare a body-fitted field across rank layouts on the same discretization."""
import argparse
import json
import hashlib
from pathlib import Path
import struct

from hpc_native_graph_checkpoint import generations
from hpc_solver_prefixes import norms


def owned_field(paths):
    chunks = []
    for path in paths:
        envelope = path.read_bytes().split(b'\n', 6)
        if len(envelope) != 7 or len(envelope[6]) < 24:
            raise ValueError('truncated owned field: '+str(path))
        payload = envelope[6]
        total, begin, end = struct.unpack('<QQQ', payload[:24])
        if total == 0 or not begin <= end <= total or len(payload) != 24+8*(end-begin):
            raise ValueError('invalid owned field extent: '+str(path))
        chunks.append((begin, end, total, struct.unpack('<'+'d'*(end-begin), payload[24:])))
    if not chunks:
        raise ValueError('no owned field shards')
    values, cursor, total = [], 0, chunks[0][2]
    for begin, end, count, local in sorted(chunks):
        if count != total or begin != cursor:
            raise ValueError('owned fields overlap, have gaps, or disagree on global size')
        values.extend(local)
        cursor = end
    if cursor != total:
        raise ValueError('owned fields do not cover the global vector')
    return values


def compare_bundles(reference, candidate, domain, role='flow'):
    if not domain or any(c in domain for c in '/\\*?[]'):
        raise ValueError('domain must be a literal shard identifier')
    shard_prefix = hashlib.sha256(('native-graph-domain-v1:'+domain).encode()).hexdigest()
    bundles = []
    for root in (reference, candidate):
        accepted = generations(root)
        if not accepted or len(accepted) != len(list(root.glob('*/manifest'))):
            raise ValueError('bundle has no accepted generations or contains invalid manifests')
        bundles.append(accepted)
    if set(bundles[0]) != set(bundles[1]):
        raise ValueError('accepted step coverage differs')
    rows = []
    for step in sorted(bundles[0]):
        paths = [list(epochs[step].glob(shard_prefix+'.'+role+'.rank-*.shard')) for epochs in bundles]
        first, second = [owned_field(parts) for parts in paths]
        metric = norms(first, second)
        row = dict(step=step, rows=len(first), reference_shards=len(paths[0]), candidate_shards=len(paths[1]), **metric)
        if role == 'flow':
            if len(first) % 4:
                raise ValueError('body-fitted flow vector is not interleaved velocity/pressure')
            row['velocity'] = norms([x for i,x in enumerate(first) if i%4 != 3],
                                    [x for i,x in enumerate(second) if i%4 != 3])
            row['pressure'] = norms(first[3::4], second[3::4])
        rows.append(row)
    return dict(status='passed', domain=domain, role=role, steps=rows,
                scope='Same discretization, coefficient ordering and pressure reference required; physical gates are separate.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference', type=Path)
    parser.add_argument('candidate', type=Path)
    parser.add_argument('--domain', required=True)
    parser.add_argument('--role', choices=['flow', 'transport'], default='flow')
    args = parser.parse_args()
    print(json.dumps(compare_bundles(args.reference, args.candidate, args.domain, args.role), indent=2))


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print(json.dumps(dict(status='failed', error=repr(error))))
        raise SystemExit(1)
