#!/usr/bin/env python3
"""Compare owned regression output from the deforming moving-FSI paired test."""
import argparse
import hashlib
import json
import math
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_run(directory):
    reports = sorted(directory.glob('rank-*/run.json'))
    if not reports:
        raise ValueError('missing rank reports')
    ranks = set()
    fields, surface, hashes = {}, {}, {}
    layout = None
    for report in reports:
        record = json.loads(report.read_text())
        rank = record['rank']
        if (rank in ranks or record['ranks'] != len(reports)
                or record['returncode'] != 0 or record['timed_out']
                or record['status'] != 'process_passed'):
            raise ValueError('failed or inconsistent rank report')
        ranks.add(rank)
        hashes[str(report)] = digest(report)
        for name in ('stdout.log', 'stderr.log'):
            path = report.parent/name
            if digest(path) != record['logs'][name]:
                raise ValueError('log digest mismatch')
            hashes[str(path)] = digest(path)
        has_layout = False
        passed = False
        for line in (report.parent/'stdout.log').read_text().splitlines():
            tokens = line.split()
            if not tokens:
                continue
            if tokens[0] == 'paired_layout':
                value = tuple(map(int, tokens[1:]))
                if has_layout or len(value) != 3 or (layout is not None and value != layout):
                    raise ValueError('inconsistent layout metadata')
                layout = value
                has_layout = True
            elif tokens[0] in ('paired_field', 'paired_surface'):
                rows = fields if tokens[0] == 'paired_field' else surface
                key = int(tokens[1])
                values = tuple(map(float, tokens[2:]))
                expected = 1 if rows is fields else 12
                if key < 0 or key in rows or len(values) != expected or not all(map(math.isfinite, values)):
                    raise ValueError('duplicate, nonfinite, or malformed owned row')
                rows[key] = values
            elif tokens[0] == 'moving_fsi_paired':
                if passed or f'rank={rank}' not in tokens or f'ranks={len(reports)}' not in tokens:
                    raise ValueError('inconsistent completion marker')
                passed = all(marker in tokens for marker in (
                    'deforming=1', 'prepare_failure=1', 'context_failure=1',
                    'exact_rollback=1', 'paired_commit=1', 'passed'))
        if not has_layout or not passed:
            raise ValueError('missing deforming regression completion')
    if ranks != set(range(len(reports))):
        raise ValueError('missing rank coverage')
    node_rows, total_rows, nodes = layout
    if (node_rows <= 0 or node_rows % 4 or total_rows < node_rows or nodes <= 0
            or set(fields) != set(range(total_rows)) or len(surface) != nodes):
        raise ValueError('incomplete owned field or surface coverage')
    return layout, fields, surface, hashes


def compare(reference, candidate):
    layout, fields, surface, _ = reference
    other_layout, other_fields, other_surface, hashes = candidate
    if layout != other_layout or surface.keys() != other_surface.keys():
        raise ValueError('reference layout or surface IDs differ')
    errors = [[] for _ in range(17)]
    values = [[] for _ in range(17)]
    for row, value in fields.items():
        group = row % 4 if row < layout[0] else 4
        errors[group].append((other_fields[row][0]-value[0])**2)
        values[group].append(value[0]**2)
    for node, row in surface.items():
        for component, value in enumerate(row):
            errors[5+component].append((other_surface[node][component]-value)**2)
            values[5+component].append(value**2)
    norms = [math.sqrt(math.fsum(error))/max(1., math.sqrt(math.fsum(value)))
             for error, value in zip(errors, values)]
    if not all(math.isfinite(value) and value < 1e-8 for value in norms):
        raise ValueError('field or surface scaled L2 exceeds 1e-8')
    return {'scaled_l2_by_component': norms, 'sha256': hashes}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-run', type=Path, required=True)
    parser.add_argument('--run', type=Path, action='append', required=True)
    args = parser.parse_args()
    reference = read_run(args.reference_run)
    results = {str(path): compare(reference, read_run(path)) for path in args.run}
    print(json.dumps({'status': 'passed', 'reference_sha256': reference[3], 'runs': results,
                      'component_order': ['u', 'v', 'w', 'p', 'auxiliary'] +
                      [f'{quantity}_{axis}' for quantity in ('displacement', 'velocity', 'traction', 'force')
                       for axis in ('x', 'y', 'z')],
                      'limitations': ['Compares test-only owned output; does not certify source provenance, strong-coupling convergence, or scaling.']}, indent=2))


if __name__ == '__main__':
    main()
