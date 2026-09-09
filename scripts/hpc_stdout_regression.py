#!/usr/bin/env python3
"""Exercise native MPI tool output failures and verify every rank's exit status."""

import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import time

from hpc_inventory import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture', type=Path, required=True)
    parser.add_argument('--baseline-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--launcher', default='mpiexec --oversubscribe')
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    summary = dict(status='running', cases=[], comparisons=[], binaries={})

    def run(name, binary, ranks, argv, expected=0, target=None, mode=None):
        directory = output / name
        directory.mkdir()
        binary = binary.resolve()
        summary['binaries'][str(binary)] = digest(binary)
        command = ['timeout', '--kill-after=5s', '100s', *shlex.split(args.launcher),
                   '-np', str(ranks), 'python3', str(repo/'scripts/hpc_rank_run.py'),
                   '--output-dir', str(directory), '--expected-ranks', str(ranks),
                   '--timeout', '90', '--', str(binary), *map(str, argv)]
        settings = dict(env)
        if target:
            settings.update(IGA_TEST_STDOUT_TARGET=target, IGA_TEST_STDOUT_MODE=mode)
        record = dict(name=name, command_argv=command, expected=expected,
                      fault_target=target, fault_mode=mode)
        summary['cases'].append(record)
        started = time.monotonic()
        with (directory/'launcher.log').open('x') as log:
            result = subprocess.run(command, env=settings, cwd=repo, stdout=log, stderr=subprocess.STDOUT)
        record.update(returncode=result.returncode, elapsed_s=time.monotonic()-started)
        if result.returncode != expected:
            raise RuntimeError(f'{name}: unexpected launcher exit {result.returncode}')
        reports = []
        texts = []
        for rank in range(ranks):
            path = directory/f'rank-{rank}'
            report = json.loads((path/'run.json').read_text())
            if report['returncode'] != expected or report['timed_out'] or report['resource'] is None:
                raise RuntimeError(f'{name}: rank {rank} failed exit/resource checks')
            stdout = (path/'stdout.log').read_text()
            stderr = (path/'stderr.log').read_text()
            if expected:
                if f'stdout_test rank={rank} injected={int(rank == 0)} native_status=1' not in stderr:
                    raise RuntimeError(f'{name}: missing injection or common native failure')
                stage = ('execution resource report' if target.startswith('execution_resources')
                         else 'mesh check result logging' if 'mesh' in binary.name
                         else 'assembly smoke result logging')
                if stage+': rank 0:' not in stderr:
                    raise RuntimeError(f'{name}: wrong collective failure boundary')
            reports.append(report)
            texts.append(stdout)
        record.update(status='passed', rank_reports=reports)
        print(name, 'passed', flush=True)
        return texts[0]

    try:
        for tool, wrapper, prefix in (
                ('iga_mesh_check', 'mesh_stdout_failure_test', 'elements='),
                ('iga_assembly_smoke', 'assembly_stdout_failure_test', 'global_rows=')):
            for ranks, filename in ((1, 'serial.ntiga'), (2, 'group.ntiga')):
                database = args.fixture.resolve()/filename
                argv = [database] + (['2'] if tool == 'iga_assembly_smoke' else [])
                old = run(f'{tool}-{ranks}-before', args.baseline_dir/tool, ranks, argv)
                reference = [line for line in old.splitlines() if line.startswith(prefix)]
                if len(reference) != 1:
                    raise RuntimeError('baseline numerical summary missing')
                if tool == 'iga_mesh_check':
                    values = dict(token.split('=') for token in reference[0].split())
                    if float(values['minimum_detJ']) <= 0 or values['bad_elements'] != '0':
                        raise RuntimeError('baseline geometry is invalid')
                for target in ('execution_resources ranks=', prefix):
                    for mode in ('short', 'flush', 'exception', 'allocation', 'nonstandard'):
                        label = f'{tool}-{ranks}-{target.split()[0].rstrip("=")}-{mode}'
                        run(label, repo/'solvers/cpu'/wrapper, ranks, argv, 1, target, mode)
                        text = run(label+'-retry', repo/'solvers/cpu'/tool, ranks, argv)
                        actual = [line for line in text.splitlines() if line.startswith(prefix)]
                        if actual != reference:
                            raise RuntimeError(f'{label}: healthy geometry/matrix summary changed')
                        summary['comparisons'].append(dict(name=label, exact=True, summary=actual[0]))
        summary['status'] = 'passed'
    except BaseException:
        summary['status'] = 'failed'
        raise
    finally:
        (output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
