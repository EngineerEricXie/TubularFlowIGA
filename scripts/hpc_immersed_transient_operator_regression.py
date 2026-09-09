#!/usr/bin/env python3
"""Check fixed-geometry transient MPI operators against the serial runtime."""
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
    parser.add_argument('--modes', nargs='+', choices=['flow','pressure','traction','closed','inertial'], default=['flow','pressure','traction','closed','inertial'])
    args = parser.parse_args()
    if len(set(args.ranks)) != len(args.ranks):
        parser.error('rank counts must be distinct')
    repo = Path(__file__).resolve().parents[1]
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    binary = repo/'solvers/cpu/immersed_transient_distributed_operator_test'
    summary = dict(status='running', binary=str(binary), binary_sha256=digest(binary), cases=[])
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    env.pop('PETSC_OPTIONS', None)
    try:
        cases = [(mode, ranks, False) for mode in args.modes for ranks in args.ranks]
        if args.split:
            cases.append(('flow', 3, True))
        for mode, ranks, split in cases:
            directory = root/(mode+'-'+str(ranks)+('-split' if split else ''))
            directory.mkdir()
            command = ['timeout', '--kill-after=5s', '240s', 'mpiexec', '--oversubscribe', '-np', str(ranks),
                       sys.executable, str(repo/'scripts/hpc_rank_run.py'), '--expected-ranks', str(ranks),
                       '--timeout', '210', '--output-dir', str(directory), '--', str(binary), mode]
            if split:
                command.append('split')
            case = dict(mode=mode, ranks=ranks, split=split, argv=command, rank_reports=[], observations=[])
            summary['cases'].append(case)
            with (directory/'launcher.log').open('w') as log:
                result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT)
            case['returncode'] = result.returncode
            if result.returncode:
                raise RuntimeError(f'transient/{mode}/{ranks}: nonzero launcher exit')
            for rank in range(ranks):
                rd = directory/f'rank-{rank}'
                report = json.loads((rd/'run.json').read_text())
                if report['returncode'] or report['timed_out'] or not report['resource']:
                    raise RuntimeError('failed or unmeasured transient operator rank')
                case['rank_reports'].append(report)
                lines = [line for line in (rd/'stdout.log').read_text().splitlines()
                         if line.startswith('immersed_transient_operator_mpi ')]
                if len(lines) != 1 or not lines[0].endswith(' passed'):
                    raise RuntimeError('missing transient operator MPI completion observation')
                values = {key: float(value) for key, value in (token.split('=') for token in lines[0].split()[1:-1])}
                expected_size = (1 if rank == 0 else 2) if split else ranks
                expected_rank = (0 if rank == 0 else rank-1) if split else rank
                if (any(not math.isfinite(value) for value in values.values())
                        or values['rank'] != expected_rank or values['ranks'] != expected_size
                        or not 0 <= values['residual_relative_l2'] < 1e-9
                        or not 0 <= values['action_relative_l2'] < 1e-9
                        or not 0 <= values['physical_error'] < 1e-12
                        or values['global_rows'] != (451 if mode in ['flow','inertial'] else 449)
                        or values['owned_rows'] < 0 or not 0 < values['owned_cells'] <= 4):
                    raise RuntimeError('invalid transient operator MPI numerical/solve observations')
                case['observations'].append(values)
            case['status'] = 'passed'
            print('transient-operator', directory.name, 'passed', flush=True)
        if digest(binary) != summary['binary_sha256']:
            raise RuntimeError('transient operator MPI binary changed during validation')
        summary['status'] = 'passed'
    except BaseException:
        summary['status'] = 'failed'
        raise
    finally:
        (root/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
