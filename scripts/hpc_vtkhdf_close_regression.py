#!/usr/bin/env python3
"""Inject final H5Fclose failures in native CPU CLIs, using staged regression inputs."""

import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

from hpc_compare_fields import compare
from hpc_inventory import digest

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--flow-fixtures', type=Path, required=True,
                        help='Output directory of hpc_flow_output_regression.py')
    parser.add_argument('--transport-fixtures', type=Path, required=True,
                        help='Output directory of hpc_transport_cli_regression.py')
    parser.add_argument('--reference-flow', type=Path, required=True)
    parser.add_argument('--reference-transport', type=Path, required=True)
    parser.add_argument('--preload', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--launcher', default='mpiexec --map-by core --bind-to core')
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    preload = args.preload.resolve()
    binaries = {'flow': ROOT/'solvers/cpu/iga_navier_stokes', 'transport': ROOT/'solvers/cpu/iga_solve'}
    references = {'flow': args.reference_flow.resolve(), 'transport': args.reference_transport.resolve()}
    summary = dict(status='running', cases=[], binaries={str(p): digest(p)
                   for p in [*binaries.values(), *references.values(), preload]})
    env = dict(os.environ, OMP_NUM_THREADS='1', IGA_ASSEMBLY_THREADS='1',
               OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', BLIS_NUM_THREADS='1',
               PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -ksp_rtol 1e-12')
    try:
        for kind in ('flow', 'transport'):
            env['PETSC_OPTIONS'] = ('-ksp_type '+('preonly' if kind == 'flow' else 'gmres')
                                    +' -pc_type lu -pc_factor_mat_solver_type mumps -ksp_rtol 1e-12')
            for mode, ranks, failed, old in [('healthy-1', 1, False, False),
                                           ('healthy-2', 2, False, False),
                                           ('before-failure', 1, True, True),
                                           ('failure-1', 1, True, False),
                                           ('failure-2', 2, True, False)]:
                name = kind+'-'+mode
                case = output/name
                case.mkdir()
                source = (args.flow_fixtures/f'vtkhdf-{ranks}'/'input' if kind == 'flow' else
                          args.transport_fixtures/('vtkhdf-one' if ranks == 1 else 'vtkhdf-two')/'input-0')
                local = case/'input'
                shutil.copytree(source, local)
                inputs = {str(p): digest(p) for p in local.iterdir() if p.is_file()}
                results, records = case/'result', case/'ranks'
                results.mkdir()
                records.mkdir()
                target = results/'field.txt'
                child = [str(references[kind] if old else binaries[kind]), str(local/'database.ntiga'),
                         str(local), '--output', str(target), '--output-every', '1',
                         '--visualization-format', 'vtkhdf']
                child += ['--max-newton', '12'] if kind == 'flow' else ['--system', 'transport']
                if failed:
                    child = ['env', 'LD_PRELOAD='+str(preload),
                             'TUBULARFLOWIGA_TEST_HDF_CLOSE='+str(target.with_suffix('.vtkhdf')), *child]
                command = ['timeout', '--kill-after=5s', '90s', *shlex.split(args.launcher), '-np', str(ranks),
                           sys.executable, str(ROOT/'scripts/hpc_rank_run.py'), '--output-dir', str(records),
                           '--expected-ranks', str(ranks), '--timeout', '60', '--', *child]
                expected = 1 if failed and not old else 0
                row = dict(case=name, command_argv=command, expected_returncode=expected,
                           input_sha256=inputs, ranks=[], fields=[], status='running')
                summary['cases'].append(row)
                with (case/'launcher.log').open('x') as log:
                    result = subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
                row['launcher_returncode'] = result.returncode
                if result.returncode != expected:
                    raise RuntimeError(f'{name}: unexpected launcher exit {result.returncode}')
                success_marker = 'navier_stokes_v2 seconds=' if kind == 'flow' else 'iga_solve system='
                for rank in range(ranks):
                    directory = records/f'rank-{rank}'
                    report = json.loads((directory/'run.json').read_text())
                    row['ranks'].append(report)
                    if (report['returncode'] != expected or report['timed_out']
                            or report['rank'] != rank or report['ranks'] != ranks):
                        raise RuntimeError(f'{name}: bad rank result')
                    for filename, checksum in report['logs'].items():
                        if digest(directory/filename) != checksum:
                            raise RuntimeError('rank log changed')
                    text = (directory/'stdout.log').read_text()+(directory/'stderr.log').read_text()
                    injected = text.count('[vtkhdf-close-injection]')
                    if injected != int(failed and rank == 0):
                        raise RuntimeError(f'{name}: wrong injection count on rank {rank}')
                    if failed and not old:
                        if (f'{kind} visualization close: rank 0: cannot close VTKHDF output' not in text
                                or success_marker in text):
                            raise RuntimeError(f'{name}: missing common error or false success')
                    elif rank == 0 and success_marker not in text:
                        raise RuntimeError(f'{name}: missing expected success summary')
                if old:
                    row['observed_false_success'] = True
                if not target.with_suffix('.vtkhdf').is_file():
                    raise RuntimeError('missing HDF5 output before close failure')
                baseline = output/(kind+'-healthy-1')/'result'
                for field in sorted(results.glob('*.txt')) + sorted(results.glob('*.pressure')):
                    observation = compare(baseline/field.name, field, 1e-6, 1e-12, node_ids=kind == 'transport')
                    row['fields'].append(dict(file=field.name, **observation))
                    if not observation['passed']:
                        raise RuntimeError(f'{name}: field changed')
                if any(digest(Path(p)) != checksum for p, checksum in inputs.items()):
                    raise RuntimeError('inputs changed')
                row['status'] = 'passed'
                (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
                print(name, 'passed', flush=True)
        if any(digest(Path(p)) != checksum for p, checksum in summary['binaries'].items()):
            raise RuntimeError('binaries changed')
        summary['status'] = 'passed'
    finally:
        (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
