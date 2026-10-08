#!/usr/bin/env python3
"""Validate a complete paired-FSI file restart and its immutable input bundle."""
import argparse
import json
import math
from pathlib import Path

from hpc_check_strong_fsi import read_run
from hpc_inventory import digest


def require(value, message):
    if not value:
        raise ValueError(message)


def fields(line):
    tokens = line.split()
    require(tokens and tokens[-1] == 'passed', 'invalid paired restart marker')
    return {token.split('=', 1)[0]: token.split('=', 1)[1]
            for token in tokens[1:-1] if '=' in token}


def marker(stdout, prefix):
    lines = [line for line in stdout.splitlines() if line.startswith(prefix+' ')]
    require(len(lines) == 1, f'missing or duplicate {prefix} marker')
    return fields(lines[0])


def finite_below(value, limit, name):
    number = float(value)
    require(math.isfinite(number) and 0 <= number < limit, f'{name} exceeds limit')
    return number


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--source-ranks', type=int, required=True)
    parser.add_argument('--input-hashes', type=Path, required=True)
    parser.add_argument('--field-limit', type=float, default=1e-8)
    args = parser.parse_args()
    require(args.source_ranks > 0 and math.isfinite(args.field_limit) and args.field_limit > 0,
            'invalid source ranks or field limit')

    _, common, hashes = read_run(args.run)
    reports = sorted(args.run.glob('rank-*/run.json'))
    target_ranks = len(reports)
    require(target_ranks > 0, 'missing rank reports')
    accepted, continued, iterations = [], [], set()
    for report in reports:
        record = json.loads(report.read_text())
        rank = record['rank']
        stdout = (report.parent/'stdout.log').read_text()
        restored = marker(stdout, 'moving_fsi_file_checkpoint')
        next_step = marker(stdout, 'moving_fsi_file_continuation')
        for values in (restored, next_step):
            require(int(values['rank']) == rank and int(values['ranks']) == target_ranks,
                    'paired marker rank differs from report')
        require(int(restored['source_ranks']) == args.source_ranks
                and int(restored['read_only']) == 1
                and int(restored['exact_file_publications']) == (args.source_ranks == target_ranks)
                and int(restored['surface_repartition']) == (args.source_ranks != target_ranks)
                and int(restored['checksum_cleanup_retry']) == 1,
                'paired restore mode or failure coverage differs')
        accepted.append(finite_below(restored['accepted_field_scaled_l2'], args.field_limit,
                                     'accepted field error'))
        continued.append(finite_below(next_step['field_scaled_l2'], args.field_limit,
                                      'continued field error'))
        require(int(next_step['surface_history_ports_conservation']) == 1,
                'continued coupled fields were not certified')
        iterations.add(int(next_step['iterations']))
    require(len(iterations) == 1 and next(iter(iterations)) > 0,
            'continuation iteration counts differ or are empty')
    require(common[(2, 'completion')][0] == next(iter(iterations)),
            'restart marker differs from complete strong history')

    saved = json.loads(args.input_hashes.read_text())
    require(saved, 'empty input hash catalog')
    current = {name: digest(Path(name)) for name in saved}
    require(saved == current, 'checkpoint input bundle changed during reader execution')
    hashes[str(args.input_hashes)] = digest(args.input_hashes)
    print(json.dumps({'status': 'passed', 'source_ranks': args.source_ranks,
                      'target_ranks': target_ranks, 'accepted_field_scaled_l2': max(accepted),
                      'next_field_scaled_l2': max(continued), 'iterations': next(iter(iterations)),
                      'sha256': hashes, 'input_bundle_sha256': current,
                      'coverage': 'rank reports, log hashes, complete strong histories, restart markers, immutable input bundle',
                      'limitations': ['Functional validation only; wall time is not a scaling result.',
                                      'Does not establish cross-node execution.']}, indent=2))


if __name__ == '__main__':
    main()
