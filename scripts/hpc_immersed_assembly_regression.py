#!/usr/bin/env python3
"""Validate owned immersed assembly, halo exchange, and serial physical parity."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

from hpc_inventory import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--ranks', nargs='+', type=int, choices=[1, 2, 4], default=[1, 2, 4])
    parser.add_argument('--kind', choices=['unit', 'physics', 'all'], default='all')
    parser.add_argument('--split', action='store_true', help='Also run the unit checks in 1+2 groups on three ranks')
    args = parser.parse_args()
    if any(r < 1 for r in args.ranks) or len(set(args.ranks)) != len(args.ranks):
        parser.error('require distinct positive rank counts')
    repo = Path(__file__).resolve().parents[1]
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    binaries = {'unit': repo/'solvers/cpu/immersed_distributed_assembly_test',
                'physics': repo/'solvers/cpu/immersed_distributed_physics_test'}
    selected = ['unit', 'physics'] if args.kind == 'all' else [args.kind]
    cases = [(kind, ranks, False) for kind in selected for ranks in args.ranks]
    if args.split:
        cases.append(('unit', 3, True))
    summary = dict(status='running', cases=[], binaries={str(binaries[k]): digest(binaries[k]) for k in {c[0] for c in cases}})
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    env.pop('PETSC_OPTIONS', None)
    try:
        for kind, ranks, split in cases:
            directory = root/(kind+'-'+str(ranks)+('-split' if split else ''))
            directory.mkdir()
            command = ['timeout', '--kill-after=5s', '240s', 'mpiexec', '--oversubscribe', '-np', str(ranks),
                       sys.executable, str(repo/'scripts/hpc_rank_run.py'), '--expected-ranks', str(ranks),
                       '--timeout', '210', '--output-dir', str(directory), '--', str(binaries[kind])]
            if split:
                command.append('split')
            record = dict(kind=kind, ranks=ranks, split=split, argv=command, rank_reports=[], observations=[])
            summary['cases'].append(record)
            with (directory/'launcher.log').open('w') as log:
                result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT)
            record['returncode'] = result.returncode
            if result.returncode:
                raise RuntimeError(directory.name+': nonzero launcher exit')
            for rank in range(ranks):
                rd = directory/f'rank-{rank}'
                report = json.loads((rd/'run.json').read_text())
                record['rank_reports'].append(report)
                if report['returncode'] or report['timed_out'] or not report['resource']:
                    raise RuntimeError(directory.name+': failed or unmeasured rank')
                prefix = 'immersed_physics ' if kind == 'physics' else 'immersed_distributed_test '
                lines = [line for line in (rd/'stdout.log').read_text().splitlines() if line.startswith(prefix)]
                if len(lines) != (1 if kind == 'physics' else 2):
                    raise RuntimeError(directory.name+': missing completion observations')
                for line in lines:
                    if not line.endswith(' passed'):
                        raise RuntimeError(directory.name+': incomplete observation')
                    values = {key:float(value) for key, value in (token.split('=') for token in line.split()[1:-1])}
                    values['world_rank'] = rank
                    record['observations'].append(values)
            if kind == 'physics':
                values = record['observations']
                if (sum(v['owned_rows'] for v in values) != values[0]['global_rows']
                        or any(v['cells'] != 27 or v['ghost_faces'] != 54 for v in values)
                        or any(not (v['residual_relative_l2'] <= 1e-6 and v['action_relative_l2'] <= 1e-6) for v in values)):
                    raise RuntimeError(directory.name+': physical parity/ownership failed')
                if ranks > 1 and any(v['halo_rows'] >= v['global_rows'] for v in values):
                    raise RuntimeError(directory.name+': physical halo replicates full state')
            else:
                grouped = {}
                for value in record['observations']:
                    grouped.setdefault((value['ranks'], value['empty']), []).append(value)
                for (size, empty), values in grouped.items():
                    if (len(values) != int(size) or sum(v['owned_stencils'] for v in values) != 63
                            or sum(v['local_nz'] for v in values) != 1008):
                        raise RuntimeError(directory.name+': unit ownership/pattern failed')
                    if not empty and size > 1 and any(v['halo_rows'] >= v['global_rows'] for v in values):
                        raise RuntimeError(directory.name+': unit halo replicates state')
            record['status'] = 'passed'
            print(directory.name, 'passed', flush=True)
        for name, sha in summary['binaries'].items():
            if digest(Path(name)) != sha:
                raise RuntimeError('binary changed during validation')
        summary['status'] = 'passed'
    except BaseException:
        summary['status'] = 'failed'
        raise
    finally:
        (root/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
