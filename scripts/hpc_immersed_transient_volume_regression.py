#!/usr/bin/env python3
"""Check owned backward-Euler volume assembly, frozen forces, and failure retry."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys

from hpc_inventory import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--ranks', type=int, nargs='+', choices=[1, 2, 4], default=[1, 2, 4])
    parser.add_argument('--split', action='store_true', help='Also validate separate 1+2 groups on three ranks')
    args = parser.parse_args()
    if len(set(args.ranks)) != len(args.ranks):
        parser.error('rank counts must be distinct')
    repo = Path(__file__).resolve().parents[1]
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    binary = repo/'solvers/cpu/immersed_distributed_transient_volume_test'
    summary = dict(status='running', binary=str(binary), binary_sha256=digest(binary), cases=[])
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    env.pop('PETSC_OPTIONS', None)
    try:
        cases = [(ranks, False) for ranks in args.ranks]
        if args.split:
            cases.append((3, True))
        for ranks, split in cases:
            directory = root/(str(ranks)+('-split' if split else ''))
            directory.mkdir()
            command = ['timeout', '--kill-after=5s', '210s', 'mpiexec', '--oversubscribe', '-np', str(ranks),
                       sys.executable, str(repo/'scripts/hpc_rank_run.py'), '--expected-ranks', str(ranks),
                       '--timeout', '180', '--output-dir', str(directory), '--', str(binary)]
            if split:
                command.append('split')
            case = dict(ranks=ranks, split=split, argv=command, rank_reports=[], observations=[])
            summary['cases'].append(case)
            with (directory/'launcher.log').open('w') as log:
                result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT)
            case['returncode'] = result.returncode
            if result.returncode:
                raise RuntimeError(f'transient-volume/{ranks}: nonzero launcher exit')
            for rank in range(ranks):
                rd = directory/f'rank-{rank}'
                report = json.loads((rd/'run.json').read_text())
                if report['returncode'] or report['timed_out'] or not report['resource']:
                    raise RuntimeError('failed or unmeasured transient volume rank')
                case['rank_reports'].append(report)
                lines = [line for line in (rd/'stdout.log').read_text().splitlines()
                         if line.startswith('immersed_transient_volume_mpi ')]
                if len(lines) != 6 or any(not line.endswith(' passed') for line in lines):
                    raise RuntimeError('missing transient volume MPI completion observation')
                expected_size = (1 if rank == 0 else 2) if split else ranks
                expected_rank = (0 if rank == 0 else rank-1) if split else rank
                for ordinal, line in enumerate(lines):
                    values = {key: float(value) for key, value in
                              (token.split('=') for token in line.split()[1:-1])}
                    if (any(not math.isfinite(value) for value in values.values())
                            or values['rank'] != expected_rank or values['ranks'] != expected_size
                            or values['mode'] != ordinal//2 or values['compact'] != ordinal%2
                            or not 0 <= values['maximum_error'] < 1e-12
                            or values['inertia_effect'] <= 1e-6 or values['global_rows'] != 451
                            or not 0 <= values['owned_cells'] <= 4
                            or values['frozen_points'] < 0
                            or (values['frozen_points'] == 0) != (values['owned_cells'] == 0)
                            or min(values['integration_s'], values['halo_s'], values['stash_s']) < 0):
                        raise RuntimeError('invalid transient volume MPI numerical observations')
                    case['observations'].append(values)
            case['status'] = 'passed'
            print('transient-volume', directory.name, 'passed', flush=True)
        if digest(binary) != summary['binary_sha256']:
            raise RuntimeError('transient volume MPI binary changed during validation')
        summary['status'] = 'passed'
    except BaseException:
        summary['status'] = 'failed'
        raise
    finally:
        (root/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
