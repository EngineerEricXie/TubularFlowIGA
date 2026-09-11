#!/usr/bin/env python3
"""Validate independent field, scalar-diagnostic, and checkpoint frequencies."""

import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import time

from hpc_inventory import digest


def profile(text):
    rows = [json.loads(line.split(' ', 1)[1]) for line in text.splitlines()
            if line.startswith('hpc_profile ')]
    if len(rows) != 1 or rows[0]['status'] != 0:
        raise RuntimeError('missing successful rank-0 profile')
    return rows[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--flow-fixture-root', type=Path, required=True,
                        help='root containing input-1 and input-2 flow fixtures')
    parser.add_argument('--transport-fixture', type=Path, required=True,
                        help='transport fixture containing serial.ntiga and group.ntiga')
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--launcher', default='mpiexec --oversubscribe')
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    flow_root = args.flow_fixture_root.resolve()
    transport = args.transport_fixture.resolve()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    binaries = {'flow': repo/'solvers/cpu/iga_navier_stokes',
                'transport': repo/'solvers/cpu/iga_solve'}
    inputs = []
    for ranks in (1, 2):
        folder = flow_root/f'input-{ranks}'
        inputs += [folder/name for name in ('database.ntiga', 'simulation_config.json',
                                            'controlmesh.vtk', 'initial_velocityfield.txt')]
    inputs += [transport/name for name in ('serial.ntiga', 'group.ntiga',
                                            'simulation_config.json', 'controlmesh.vtk',
                                            'initial_velocityfield.txt')]
    for path in inputs:
        if not path.is_file():
            parser.error(f'missing fixture {path}')
    original_inputs = {str(path): digest(path) for path in inputs}
    summary = {'status': 'running', 'cases': [], 'inputs': original_inputs,
               'binaries': {name: digest(path) for name, path in binaries.items()}}
    environment = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                       IGA_PROFILE='1',
                       PETSC_OPTIONS='-ksp_type gmres -pc_type lu '
                                     '-pc_factor_mat_solver_type mumps -ksp_rtol 1e-12')

    def run(kind, ranks, frequency, name=None):
        name = name or f'{kind}-ranks-{ranks}-frequency-{frequency}'
        case = out/name
        result = case/'result'
        result.mkdir(parents=True)
        if kind == 'flow':
            fixture = flow_root/f'input-{ranks}'
            arguments = [fixture/'database.ntiga', fixture,
                         '--output', result/'field.txt', '--visualization-format', 'pvtu']
        else:
            database = transport/('serial.ntiga' if ranks == 1 else 'group.ntiga')
            arguments = [database, transport, '--system', 'transport',
                         '--output', result/'field.txt', '--visualization-format', 'pvtu']
        arguments += ['--output-every', str(frequency),
                      '--checkpoint', result/'checkpoint',
                      '--checkpoint-every', str(frequency),
                      '--diagnostic-every', str(frequency)]
        command = ['timeout', '--kill-after=5s', '180s', *shlex.split(args.launcher),
                   '-np', str(ranks), 'python3', repo/'scripts/hpc_rank_run.py',
                   '--output-dir', case, '--expected-ranks', str(ranks), '--timeout', '150',
                   '--', binaries[kind], *arguments]
        record = {'name': name, 'kind': kind, 'ranks': ranks, 'frequency': frequency,
                  'command_argv': list(map(str, command)), 'status': 'running'}
        summary['cases'].append(record)
        started = time.monotonic()
        launch = case/'launcher.log'
        with launch.open('x') as stream:
            completed = subprocess.run(list(map(str, command)), cwd=repo, env=environment,
                                       stdout=stream, stderr=subprocess.STDOUT)
        if completed.returncode != 0:
            raise RuntimeError(f'{name} returned {completed.returncode}')
        reports = [json.loads((case/f'rank-{rank}/run.json').read_text())
                   for rank in range(ranks)]
        if any(row['returncode'] or row['timed_out'] or not row['resource'] for row in reports):
            raise RuntimeError(f'{name} has a failed rank report')
        stdout = (case/'rank-0/stdout.log').read_text()
        diagnostics = sum(line.startswith('solver_configuration ') for line in stdout.splitlines())
        checkpoints = sum(line.startswith('checkpoint=') for line in stdout.splitlines())
        expected = 2 if frequency == 1 else 1
        if diagnostics != expected or checkpoints != expected:
            raise RuntimeError(f'{name} frequency mismatch: diagnostics={diagnostics}, checkpoints={checkpoints}')
        output_files = [path for path in result.rglob('*') if path.is_file()
                        and (path.suffix in ('.vtu', '.pvtu', '.pvd', '.vtkhdf')
                             or '.step' in path.name)]
        row_profile = profile(stdout)
        output_seconds = row_profile['phases']['output']['exclusive_s']
        output_bytes = sum(path.stat().st_size for path in output_files)
        final_field = result/'field.txt'
        record.update(status='passed', returncode=0, elapsed_s=time.monotonic()-started,
                      rank_reports=len(reports), diagnostic_records=diagnostics,
                      checkpoint_records=checkpoints, visualization_files=len(output_files),
                      visualization_bytes=output_bytes, output_exclusive_s=output_seconds,
                      observed_bytes_per_output_second=(output_bytes/output_seconds
                                                        if output_seconds > 0 else None),
                      final_state_sha256=digest(result/'checkpoint.state'),
                      final_field_sha256=digest(final_field) if final_field.is_file() else None,
                      launcher_log_sha256=digest(launch))
        print(name, 'passed', flush=True)
        return record

    def run_blocked_transport_index():
        name = 'transport-ranks-2-blocked-index'
        case = out/name
        result = case/'result'
        result.mkdir(parents=True)
        (result/'field.pvd.pending').mkdir()
        arguments = [transport/'group.ntiga', transport, '--system', 'transport',
                     '--output', result/'field.txt', '--visualization-format', 'pvtu',
                     '--output-every', '1', '--diagnostic-every', '1']
        command = ['timeout', '--kill-after=5s', '180s', *shlex.split(args.launcher),
                   '-np', '2', 'python3', repo/'scripts/hpc_rank_run.py',
                   '--output-dir', case, '--expected-ranks', '2', '--timeout', '150',
                   '--', binaries['transport'], *arguments]
        record = {'name': name, 'kind': 'transport', 'ranks': 2,
                  'frequency': 1, 'expected': 'collective rejection',
                  'command_argv': list(map(str, command)), 'status': 'running'}
        summary['cases'].append(record)
        started = time.monotonic()
        launch = case/'launcher.log'
        with launch.open('x') as stream:
            completed = subprocess.run(list(map(str, command)), cwd=repo, env=environment,
                                       stdout=stream, stderr=subprocess.STDOUT)
        reports = [json.loads((case/f'rank-{rank}/run.json').read_text())
                   for rank in range(2)]
        if completed.returncode == 0:
            raise RuntimeError(f'{name} unexpectedly succeeded')
        if any(row['returncode'] == 0 or row['timed_out'] for row in reports):
            raise RuntimeError(f'{name} did not reject collectively in finite time')
        if (result/'field.pvd').exists():
            raise RuntimeError(f'{name} published a final PVD index')
        if any('iga_solve system=' in (case/f'rank-{rank}/stdout.log').read_text()
               for rank in range(2)):
            raise RuntimeError(f'{name} printed a successful summary')
        record.update(status='passed', returncode=completed.returncode,
                      elapsed_s=time.monotonic()-started, rank_reports=len(reports),
                      final_index_published=False,
                      launcher_log_sha256=digest(launch))
        print(name, 'expected rejection', flush=True)
        return record

    try:
        by_layout = {}
        for kind in ('flow', 'transport'):
            for ranks in (1, 2):
                dense = run(kind, ranks, 1)
                sparse = run(kind, ranks, 2)
                if dense['final_state_sha256'] != sparse['final_state_sha256']:
                    raise RuntimeError(f'{kind} ranks={ranks} checkpoint differs by frequency')
                if dense['final_field_sha256'] != sparse['final_field_sha256']:
                    raise RuntimeError(f'{kind} ranks={ranks} final field differs by frequency')
                if dense['visualization_files'] <= sparse['visualization_files']:
                    raise RuntimeError(f'{kind} ranks={ranks} field frequency did not reduce file count')
                by_layout[f'{kind}-{ranks}'] = {
                    'dense_files': dense['visualization_files'],
                    'sparse_files': sparse['visualization_files'],
                    'dense_bytes': dense['visualization_bytes'],
                    'sparse_bytes': sparse['visualization_bytes']}
        run_blocked_transport_index()
        retry = run('transport', 2, 1, name='transport-retry')
        reference = next(row for row in summary['cases']
                         if row['name'] == 'transport-ranks-2-frequency-1')
        if retry['final_state_sha256'] != reference['final_state_sha256']:
            raise RuntimeError('healthy transport retry changed the final checkpoint')
        if any(digest(Path(path)) != expected for path, expected in original_inputs.items()):
            raise RuntimeError('fixture changed during regression')
        summary.update(status='passed', comparisons=by_layout)
    except BaseException:
        summary['status'] = 'failed'
        raise
    finally:
        (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
