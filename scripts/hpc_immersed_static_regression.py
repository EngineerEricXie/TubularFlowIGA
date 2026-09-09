#!/usr/bin/env python3
"""Check distributed steady immersed Newton solves and transactional recovery."""
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
    parser.add_argument('--mode', choices=['closed', 'flow', 'pressure', 'empty-work'], default='closed')
    parser.add_argument('--split', action='store_true', help='Also validate separate 1+2 groups on three ranks')
    parser.add_argument('--partition', choices=['cell-count', 'weighted'], default='cell-count')
    args = parser.parse_args()
    if len(set(args.ranks)) != len(args.ranks):
        parser.error('rank counts must be distinct')
    repo = Path(__file__).resolve().parents[1]
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    binary = repo/'solvers/cpu/immersed_distributed_static_flow_test'
    summary = dict(status='running', mode=args.mode, partition=args.partition, binary=str(binary), binary_sha256=digest(binary), cases=[])
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    env.pop('PETSC_OPTIONS', None)
    try:
        cases = [(ranks, False) for ranks in args.ranks]
        if args.split:
            cases.append((3, True))
        for ranks, split in cases:
            directory = root/(str(ranks)+('-split' if split else ''))
            directory.mkdir()
            command = ['timeout', '--kill-after=5s', '660s', 'mpiexec', '--oversubscribe', '-np', str(ranks),
                       sys.executable, str(repo/'scripts/hpc_rank_run.py'), '--expected-ranks', str(ranks),
                       '--timeout', '630', '--output-dir', str(directory), '--', str(binary), args.mode]
            if args.partition == 'weighted':
                command.append('weighted')
            if split:
                command.append('split')
            case = dict(ranks=ranks, split=split, argv=command, rank_reports=[], observations=[])
            summary['cases'].append(case)
            with (directory/'launcher.log').open('w') as log:
                result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT)
            case['returncode'] = result.returncode
            if result.returncode:
                raise RuntimeError(f'{args.mode}/{ranks}: nonzero launcher exit')
            for rank in range(ranks):
                rd = directory/f'rank-{rank}'
                report = json.loads((rd/'run.json').read_text())
                if report['returncode'] or report['timed_out'] or not report['resource']:
                    raise RuntimeError('failed or unmeasured static rank')
                case['rank_reports'].append(report)
                lines = [line for line in (rd/'stdout.log').read_text().splitlines()
                         if line.startswith('immersed_static_mpi ')]
                if len(lines) != 1 or not lines[0].endswith(' passed'):
                    raise RuntimeError('missing static MPI completion observation')
                values = {key: float(value) for key, value in (token.split('=') for token in lines[0].split()[1:-1])}
                expected_size = (1 if rank == 0 else 2) if split else ranks
                expected_rank = (0 if rank == 0 else rank-1) if split else rank
                if (any(not math.isfinite(value) for value in values.values())
                        or values['rank'] != expected_rank or values['ranks'] != expected_size
                        or not 0 <= values['field_relative_l2'] <= 1e-6
                        or values['nonlinear_iterations'] < 1 or values['ksp_iterations'] < 1
                        or values['assembly_s'] <= 0 or values['linear_s'] <= 0):
                    raise RuntimeError('invalid static MPI numerical/solve observations')
                case['observations'].append(values)
            case['status'] = 'passed'
            print(args.mode, directory.name, 'passed', flush=True)
        if digest(binary) != summary['binary_sha256']:
            raise RuntimeError('static MPI binary changed during validation')
        summary['status'] = 'passed'
    except BaseException:
        summary['status'] = 'failed'
        raise
    finally:
        (root/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
