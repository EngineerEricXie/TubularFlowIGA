#!/usr/bin/env python3
"""Validate and compare saved two-step strong-FSI regression output."""
import argparse
import json
import math
from pathlib import Path

from hpc_check_moving_fsi import compare as compare_fields, digest


def require(value, message):
    if not value:
        raise ValueError(message)


def numbers(tokens, count):
    result = tuple(map(float, tokens))
    require(len(result) == count and all(map(math.isfinite, result)), 'invalid numeric record')
    return result


def read_run(directory, steps=2):
    reports = sorted(directory.glob('rank-*/run.json'))
    require(reports, 'missing rank reports')
    seen, hashes, snapshots, common = set(), {}, {}, {}
    expected_steps = set(range(1, steps+1))
    for report in reports:
        record = json.loads(report.read_text())
        rank = record['rank']
        require(rank not in seen and record['ranks'] == len(reports) and record['returncode'] == 0
                and not record['timed_out'] and record['status'] == 'process_passed', 'invalid rank report')
        seen.add(rank)
        hashes[str(report)] = digest(report)
        for name in ('stdout.log', 'stderr.log'):
            path = report.parent/name
            require(digest(path) == record['logs'][name], 'log digest mismatch')
            hashes[str(path)] = digest(path)
        local = {}
        for line in (report.parent/'stdout.log').read_text().splitlines():
            tokens = line.split()
            if not tokens or tokens[0] not in ('strong_fsi_step', 'strong_layout', 'strong_field',
                                               'strong_surface', 'strong_iteration', 'strong_port', 'strong_conservation'):
                continue
            kind = tokens[0]
            if kind == 'strong_fsi_step':
                values = dict(token.split('=', 1) for token in tokens[1:] if '=' in token)
                step = int(values['step'])
                require(int(values['rank']) == rank and int(values['ranks']) == len(reports)
                        and tokens[-1] == 'passed', 'invalid completion record')
                value = (int(values['iterations']), *numbers([values['rms'], values['velocity_residual_times_dt'], values['threshold']], 3))
                key = (step, 'completion')
            else:
                step = int(tokens[1])
                snapshot = snapshots.setdefault(step, {'fields': {}, 'surface': {}, 'mapping': {}})
                if kind == 'strong_field':
                    require(len(tokens) == 6, 'invalid field record')
                    row, node, component = map(int, tokens[2:5])
                    require(row >= 0 and 0 <= node <= 2**64-1 and row not in snapshot['fields'], 'duplicate or invalid field ID')
                    snapshot['fields'][row] = numbers(tokens[5:], 1)
                    snapshot['mapping'][row] = (node, component)
                    continue
                if kind == 'strong_surface':
                    node = int(tokens[2])
                    require(0 <= node <= 2**64-1 and node not in snapshot['surface'], 'duplicate or invalid surface ID')
                    snapshot['surface'][node] = numbers(tokens[3:], 12)
                    continue
                if kind == 'strong_layout':
                    value = tuple(map(int, tokens[2:]))
                    require(len(value) == 3, 'invalid layout record')
                    key = (step, 'layout')
                elif kind == 'strong_conservation':
                    value = numbers(tokens[2:], 10)
                    require(all(x >= 0 for x in value) and all(a <= b for a, b in zip(value[:5], value[5:])), 'conservation gate failed')
                    key = (step, 'conservation')
                elif kind == 'strong_iteration':
                    iteration = int(tokens[2])
                    require(iteration >= 0, 'negative iteration')
                    value = numbers(tokens[3:], 6)
                    require(all(x >= 0 for x in value) and value[2] > 0, 'invalid convergence diagnostic')
                    key = (step, 'iteration', iteration)
                else:
                    value = numbers(tokens[3:], 4)
                    require(value[0] > 0, 'invalid port area')
                    key = (step, 'port', tokens[2])
            require(step in expected_steps and key not in local, 'unexpected step or duplicate shared record')
            local[key] = value
        require({key[0] for key in local} == expected_steps, 'missing accepted step')
        if common:
            require(local == common, 'shared histories differ across ranks')
        else:
            common = local
    require(seen == set(range(len(reports))) and snapshots.keys() == expected_steps, 'missing rank or snapshot coverage')
    for step, snapshot in snapshots.items():
        layout = common[(step, 'layout')]
        node_rows, rows, nodes = layout
        require(node_rows > 0 and node_rows % 4 == 0 and rows >= node_rows and nodes > 0
                and set(snapshot['fields']) == set(range(rows)) and len(snapshot['surface']) == nodes, 'incomplete snapshot')
        node_ids = []
        for row, (node, component) in sorted(snapshot['mapping'].items()):
            require(component == (row % 4 if row < node_rows else 4), 'field component mapping differs')
            if row < node_rows:
                if row % 4 == 0:
                    node_ids.append(node)
                require(node == node_ids[-1], 'inconsistent node component IDs')
            else:
                require(node == 2**64-1, 'invalid auxiliary row ID')
        require(len(set(node_ids)) == len(node_ids), 'duplicate physical node mapping')
        count, rms, velocity, threshold = common[(step, 'completion')]
        history = sorted(key[2] for key in common if key[:2] == (step, 'iteration'))
        require(0 < count <= 16 and history == list(range(count)), 'incomplete iteration history')
        final = common[(step, 'iteration', count-1)]
        require((rms, velocity, threshold) == (final[0], final[4], final[3])
                and rms <= threshold and velocity <= threshold and final[5] == 0, 'invalid convergence decision')
        for iteration in range(count-1):
            value = common[(step, 'iteration', iteration)]
            require((value[0] > value[3] or value[4] > value[3]) and 0 < value[5] <= 1, 'invalid unconverged iteration')
        require({key[2] for key in common if key[:2] == (step, 'port')} == {'inlet', 'outlet'}, 'missing port history')
        require((step, 'conservation') in common, 'missing conservation record')
        snapshot['layout'] = layout
    return snapshots, common, hashes


def compare(reference, candidate):
    ref_steps, ref_common, _ = reference
    steps, common, hashes = candidate
    require(ref_steps.keys() == steps.keys() and ref_common.keys() == common.keys(), 'step, port or iteration coverage differs')
    field_results = {}
    for step, ref in ref_steps.items():
        other = steps[step]
        require(ref['mapping'] == other['mapping'], 'stable node/component identities differ')
        field_results[step] = compare_fields((ref['layout'], ref['fields'], ref['surface'], {}),
                                             (other['layout'], other['fields'], other['surface'], {}))['scaled_l2_by_component']
    maximum = 0.
    for key, values in ref_common.items():
        if key[1] == 'layout':
            require(values == common[key], 'layout differs')
            continue
        for a, b in zip(values, common[key]):
            ratio = abs(a-b)/(1e-12+1e-6*abs(a))
            require(math.isfinite(ratio) and ratio <= 1, 'history or port numerical mismatch')
            maximum = max(maximum, ratio)
        if key[1] == 'conservation':
            require(values[5:] == common[key][5:], 'physical policy differs')
    return {'field_scaled_l2_by_step': field_results, 'maximum_history_error_over_tolerance': maximum, 'sha256': hashes}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-run', type=Path, required=True)
    parser.add_argument('--run', type=Path, action='append', required=True)
    args = parser.parse_args()
    reference = read_run(args.reference_run)
    print(json.dumps({'status': 'passed', 'reference_sha256': reference[2],
                      'runs': {str(path): compare(reference, read_run(path)) for path in args.run},
                      'limitations': ['Validates saved two-step regression output; not source provenance, independent discretization, checkpoint restart or scaling.']}, indent=2))


if __name__ == '__main__':
    main()
