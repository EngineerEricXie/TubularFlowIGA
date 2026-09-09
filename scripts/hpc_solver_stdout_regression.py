#!/usr/bin/env python3
"""Inject native CPU solver stdout failures and check collective exits and retries."""
import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import time
from hpc_compare_fields import compare
from hpc_inventory import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ('flow-case', 'transport-case', 'legacy-case', 'baseline-dir', 'output-dir'):
        parser.add_argument('--'+flag, type=Path, required=True)
    parser.add_argument('--launcher', default='mpiexec --oversubscribe')
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', IGA_PROFILE='1',
               PETSC_OPTIONS='-ksp_type gmres -pc_type lu -pc_factor_mat_solver_type mumps -ksp_rtol 1e-12')
    summary = dict(status='running', cases=[], comparisons=[], binaries={}, fixtures={})
    fixtures = {}
    for kind, source in [('flow', args.flow_case), ('transport', args.transport_case), ('legacy', args.legacy_case)]:
        target = output/(kind+'-input')
        target.mkdir()
        names = ['controlmesh.vtk', 'initial_velocityfield.txt', 'serial.ntiga', 'group.ntiga',
                 'simulation_parameter.txt' if kind == 'legacy' else 'simulation_config.json']
        for name in names:
            origin = source/name
            if kind == 'flow' and not origin.exists() and name in ('serial.ntiga', 'group.ntiga'):
                origin = source/('fixture.ntiga' if name == 'serial.ntiga' else 'fixture-2.ntiga')
            shutil.copy2(origin, target/name)
        if kind != 'legacy':
            config = json.loads((target/'simulation_config.json').read_text())
            config.pop('coupling', None)
            config.pop('external_circuit', None)
            config['simulation_scope'] = {'mode': 'flow_only'}
            # The adapter fixture starts at rest. Drive the CLI fixture so it
            # actually exercises Newton iterations and evolving scalar fields.
            velocity = target/'initial_velocityfield.txt'
            rows = len(velocity.read_text().splitlines())
            velocity.write_text(('1 0 0\n' if kind == 'flow' else '0.01 0 0\n')*rows)
            if kind == 'flow':
                inlet = next(b for b in config['boundaries'] if b['label'] == 1)
                inlet['conditions'][0]['scale'] = 0.001
            else:
                for field in config['fields']:
                    if field['name'] == 'three_red': field['initial_value'] = 0.5
                    if field['name'] == 'three_blue': field['initial_value'] = 2.0
            (target/'simulation_config.json').write_text(json.dumps(config)+'\n')
        fixtures[kind] = target
    outlet = output/'outlet-input'
    shutil.copytree(fixtures['flow'], outlet)
    config = json.loads((outlet/'simulation_config.json').read_text())
    for boundary in config['boundaries']:
        if boundary['label'] in (2, 3):
            boundary['conditions'] = [dict(field='pressure', type='resistance', resistance=0.01)]
    (outlet/'simulation_config.json').write_text(json.dumps(config)+'\n')
    fixtures['outlet'] = outlet
    for fixture in fixtures.values():
        for path in fixture.iterdir(): summary['fixtures'][str(path)] = digest(path)
    tools = {'flow': ('iga_navier_stokes', 'flow_stdout_failure_test'),
             'transport': ('iga_solve', 'transport_stdout_failure_test'),
             'legacy': ('iga_transport', 'legacy_stdout_failure_test')}

    def run(name, kind, ranks, variant='normal', old=False, target=None, mode=None, stage=None, restart=None, stop=False):
        directory = output/name
        result = directory/'result'
        result.mkdir(parents=True)
        fixture = fixtures['outlet' if variant == 'outlet' else kind]
        database = fixture/('serial.ntiga' if ranks == 1 else 'group.ntiga')
        binary = (args.baseline_dir.resolve() if old else repo/'solvers/cpu')/tools[kind][bool(target)]
        command_args = [str(database), str(fixture)]
        if kind == 'legacy': command_args += ['2', str(result/'field.txt')]
        else:
            command_args += ['--output', str(result/'field.txt'), '--output-every', '1',
                             '--checkpoint', str(result/'checkpoint'), '--checkpoint-every', '1',
                             '--visualization-format', 'vtkhdf' if variant == 'hdf' else 'vtu']
            if kind == 'flow': command_args += ['--max-newton', '12']
            else: command_args += ['--system', 'transport']
            if restart: command_args += ['--restart', str(restart)]
            if stop: command_args += ['--stop-after-step', '1']
        command = ['timeout', '--kill-after=5s', '100s', *shlex.split(args.launcher), '-np', str(ranks),
                   'python3', str(repo/'scripts/hpc_rank_run.py'), '--output-dir', str(directory),
                   '--expected-ranks', str(ranks), '--timeout', '90', '--', str(binary), *command_args]
        settings = dict(env)
        if target: settings.update(IGA_TEST_STDOUT_TARGET=target, IGA_TEST_STDOUT_MODE=mode)
        expected = 1 if target else 0
        record = dict(name=name, argv=command, expected=expected, target=target, mode=mode, stage=stage, variant=variant)
        summary['cases'].append(record)
        summary['binaries'][str(binary)] = digest(binary)
        started = time.monotonic()
        with (directory/'launcher.log').open('w') as log:
            process = subprocess.run(command, env=settings, stdout=log, stderr=subprocess.STDOUT)
        record.update(returncode=process.returncode, elapsed_s=time.monotonic()-started)
        if process.returncode != expected: raise RuntimeError(f'{name}: exit {process.returncode}, expected {expected}')
        reports = []
        for rank in range(ranks):
            report = json.loads((directory/f'rank-{rank}/run.json').read_text())
            if report['returncode'] != expected or report['timed_out'] or report['resource'] is None:
                raise RuntimeError(f'{name}: invalid rank {rank} report')
            if target:
                text = (directory/f'rank-{rank}/stderr.log').read_text()
                if f'stdout_test rank={rank} injected={int(rank == 0)} native_status=1' not in text:
                    raise RuntimeError(f'{name}: missed fault or different native status')
                if stage+': rank 0:' not in text: raise RuntimeError(f'{name}: wrong failure boundary')
            reports.append(report)
        record.update(status='passed', rank_reports=reports)
        if old and kind == 'flow':
            text = (directory/'rank-0/stdout.log').read_text()
            if ' newton=' not in text or 'total_linear_iterations=0 ' in text:
                raise RuntimeError(name+': fixture must exercise Newton linear solves')
        print(name, 'passed', flush=True)
        return result

    targets = {
        'flow': [('input', 'boundary_config=', 'flow boundary input', 'normal'),
                 ('assembly', 'body_fitted_element_assembly rank=', 'flow element assembly', 'normal'),
                 ('iteration', ' newton=', 'flow iteration logging', 'normal'),
                 ('convergence', ' converged newton=', 'flow convergence logging', 'normal'),
                 ('checkpoint', 'checkpoint=', 'checkpoint write logging', 'normal'),
                 ('restart', 'restart=', 'checkpoint restart logging', 'restart'),
                 ('visualization', 'bezier_geometry_points=', 'flow visualization initialization', 'hdf'),
                 ('outlet', ' outlet_iteration=', 'flow outlet evaluation', 'outlet'),
                 ('outlet-state', 'outlet label=', 'flow outlet evaluation', 'outlet'),
                 ('completion', 'navier_stokes_v2 seconds=', 'flow completion logging', 'normal'),
                 ('profile', 'hpc_profile ', 'flow profile logging', 'normal')],
        'transport': [('input', 'configuration=simulation_config.json system=', 'transport input logging', 'normal'),
                      ('fields', 'field[', 'transport input logging', 'normal'),
                      ('checkpoint', 'checkpoint=', 'transport checkpoint logging', 'normal'),
                      ('restart', 'restart=', 'transport restart logging', 'restart'),
                      ('visualization', 'bezier_geometry_points=', 'transport visualization initialization', 'hdf'),
                      ('completion', 'iga_solve system=', 'transport final summary', 'normal'),
                      ('profile', 'hpc_profile ', 'transport profile logging', 'normal')],
        'legacy': [('input', 'boundary_config=', 'legacy transport input logging', 'normal'),
                   ('completion', 'transport_v2 nodes=', 'legacy transport final summary', 'normal')]}
    try:
        for kind in tools:
            for ranks in (1, 2):
                references = {}
                variants = ('normal', 'outlet') if kind == 'flow' else ('normal',)
                for variant in variants:
                    references[variant] = run(f'{kind}-{ranks}-{variant}-before', kind, ranks, variant, old=True)
                checkpoint = None
                if kind != 'legacy':
                    seed = run(f'{kind}-{ranks}-restart-seed', kind, ranks, stop=True)
                    checkpoint = seed/'checkpoint'
                for label, target, stage, variant in targets[kind]:
                    for mode in ('short', 'flush', 'exception', 'allocation', 'nonstandard'):
                        name = f'{kind}-{ranks}-{label}-{mode}'
                        restart = checkpoint if variant == 'restart' else None
                        run(name, kind, ranks, variant, target=target, mode=mode, stage=stage, restart=restart)
                        retry = run(name+'-retry', kind, ranks, variant, restart=restart)
                        reference = references['outlet' if variant == 'outlet' else 'normal']
                        comparison = compare(reference/'field.txt', retry/'field.txt', 1e-6, 1e-12, kind != 'flow')
                        if not comparison['passed']: raise RuntimeError(name+': field parity failed')
                        summary['comparisons'].append(dict(name=name, **comparison))
                        if kind == 'flow':
                            pressure = compare(reference/'field.txt.pressure', retry/'field.txt.pressure', 1e-6, 1e-12)
                            if not pressure['passed']: raise RuntimeError(name+': pressure parity failed')
                            summary['comparisons'].append(dict(name=name+'-pressure', **pressure))
        for name, sha in summary['fixtures'].items():
            if digest(Path(name)) != sha: raise RuntimeError('fixture changed')
        summary['status'] = 'passed'
    except BaseException:
        summary['status'] = 'failed'
        raise
    finally:
        (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
