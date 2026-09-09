#!/usr/bin/env python3
"""Validate terminal runtime cleanup using native MPI CLI entry points."""
import argparse
import json
import os
import shlex
import subprocess
from pathlib import Path
from hpc_compare_fields import compare
from hpc_inventory import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-dir', type=Path, required=True)
    parser.add_argument('--fixture-root', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--launcher', default='mpiexec --oversubscribe')
    parser.add_argument('--scope', choices=('all', 'coupling'), default='all')
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    fixtures = args.fixture_root.resolve()
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               PETSC_OPTIONS='-ksp_type gmres -pc_type lu -pc_factor_mat_solver_type mumps -ksp_rtol 1e-12')
    summary = dict(status='running', cases=[], comparisons=[], binaries={}, inputs={})
    specs = []
    for ranks, group in ((1, 'single'), (2, 'pair')):
        flow = fixtures/'hpc01/solver-stdout/native-driven/flow-input'
        vca = fixtures/'hpc01/checkpoint-write/vca/vca_3d_smoke_test/tubularflowiga-vca-3d-smoke'
        for kind, case, database in (('flow', flow, flow/('serial.ntiga' if ranks == 1 else 'group.ntiga')),
                                    ('vca', vca, vca/('fixture.ntiga' if ranks == 1 else 'fixture-2.ntiga'))):
            specs.append(dict(name=f'{kind}-{ranks}', ranks=ranks, folder='cpu', native='iga_navier_stokes',
                wrapper='runtime_cleanup_cli_test', argv=[database, case, '--max-newton', '30',
                    '--nonlinear-rtol', '1e-8', '--nonlinear-atol', '1e-12', '--mass-rtol', '1e-6',
                    '--visualization-format', 'vtu'], output='field', transport=kind == 'vca',
                success='navier_stokes_v2 seconds='))
        for kind in ('flow', 'species', 'bifurcation'):
            bif = kind == 'bifurcation'
            case = fixtures/'hpc01/registry/native-final'/group/('flow' if bif else kind)
            specs.append(dict(name=f'graph-{kind}-{ranks}', ranks=ranks, folder='coupling',
                native='iga_1d_3d_bifurcation' if bif else 'iga_multidomain_flow',
                wrapper='bifurcation_cleanup_cli_test' if bif else 'graph_cleanup_cli_test',
                argv=['--graph-case', case], output='directory', transport=kind == 'species', success='completed schema-v'))
        for scheme in ('explicit', 'fixed', 'aitken'):
            case = fixtures/'hpc01/sequential-output/accepted'/f'input-{ranks}-{scheme}'
            specs.append(dict(name=f'sequential-{scheme}-{ranks}', ranks=ranks, folder='coupling',
                native='iga_1d_3d_explicit', wrapper='sequential_cleanup_cli_test', argv=['--graph-case', case],
                output='directory', transport=False, success='completed '))
    if args.scope == 'coupling': specs = [spec for spec in specs if spec['folder'] == 'coupling']
    summary['scope'] = args.scope
    for spec in specs:
        for argument in spec['argv']:
            if isinstance(argument, Path):
                paths = argument.rglob('*') if argument.is_dir() else [argument]
                for path in paths:
                    if path.is_file(): summary['inputs'][str(path)] = digest(path)

    def execute(spec, label, old=False, target=None):
        directory = root/(spec['name']+'-'+label)
        directory.mkdir()
        (directory/'ranks').mkdir()
        output = directory/'result'
        if spec['output'] == 'field': output.mkdir()
        binary = (args.baseline_dir.resolve() if old else repo/'solvers'/spec['folder'])/(spec['wrapper'] if target else spec['native'])
        argv = [str(binary), *map(str, spec['argv'])]
        argv += ['--output', str(output/'field.txt')] if spec['output'] == 'field' else ['--output-dir', str(output)]
        command = ['timeout', '--kill-after=5s', '120s', *shlex.split(args.launcher), '-np', str(spec['ranks']),
            'python3', str(repo/'scripts/hpc_sequential_output_regression.py'), '--rank-worker',
            '--output-dir', str(directory/'ranks'), '--', *argv]
        settings = dict(env)
        if target: settings.update(IGA_TEST_CLEANUP_TARGET=target, IGA_TEST_CLEANUP_RANK=str(spec['ranks']-1))
        record = dict(name=directory.name, command=command, target=target, ranks=spec['ranks'])
        summary['cases'].append(record)
        summary['binaries'][str(binary)] = digest(binary)
        with (directory/'launcher.log').open('w') as log:
            process = subprocess.run(command, env=settings, stdout=log, stderr=subprocess.STDOUT)
        expected = int(target is not None)
        record['returncode'] = process.returncode
        if process.returncode != expected: raise RuntimeError(directory.name+': wrong launcher exit')
        if target and spec['output'] == 'directory' and output.exists():
            raise RuntimeError(directory.name+': cleanup failure published coupling output')
        reports = []
        for rank in range(spec['ranks']):
            path = directory/'ranks'/f'rank-{rank}'
            report = json.loads((path/'run.json').read_text())
            if report['returncode'] != expected or report['timed_out'] or report['resource'] is None:
                raise RuntimeError(directory.name+': wrong rank outcome')
            stdout = (path/'stdout.log').read_text()
            stderr = (path/'stderr.log').read_text()
            if target:
                if spec['success'] in stdout: raise RuntimeError(directory.name+': premature success summary')
                if f'cleanup_test rank={rank} injected={int(rank == spec["ranks"]-1)} native_status=1' not in stderr:
                    raise RuntimeError(directory.name+': missing native injection result')
                if rank == 0 and f'{target.split()[0]} runtime cleanup: rank {spec["ranks"]-1}:' not in stderr:
                    raise RuntimeError(directory.name+': wrong common failure stage')
            elif rank == 0 and spec['success'] not in stdout: raise RuntimeError(directory.name+': missing completion')
            reports.append(report)
        record.update(status='passed', rank_reports=reports)
        print(directory.name, 'passed', flush=True)
        return output

    def compare_outputs(spec, before, after):
        files = {str(p.relative_to(before)) for p in before.rglob('*') if p.is_file()}
        if not files or files != {str(p.relative_to(after)) for p in after.rglob('*') if p.is_file()}:
            raise RuntimeError(str(after)+': different output set')
        if spec['output'] == 'field':
            for name in ('field.txt', 'field.txt.pressure'):
                result = compare(before/name, after/name, 1e-6, 1e-12)
                if not result['passed']: raise RuntimeError(str(after)+': field mismatch')
                summary['comparisons'].append(dict(case=str(after), file=name, result=result))
            files -= {'field.txt', 'field.txt.pressure'}
        for name in sorted(files):
            if (before/name).read_bytes() != (after/name).read_bytes(): raise RuntimeError(str(after)+': changed '+name)
            summary['comparisons'].append(dict(case=str(after), file=name, exact=True))

    try:
        for spec in specs:
            before = execute(spec, 'before', old=True)
            compare_outputs(spec, before, execute(spec, 'healthy'))
            targets = ['flow solver', 'flow jacobian']
            if spec['transport']: targets += ['transport solver', 'transport left']
            for target in targets:
                execute(spec, target.replace(' ', '-'), target=target)
                compare_outputs(spec, before, execute(spec, target.replace(' ', '-')+'-retry'))
        for name, sha in summary['inputs'].items():
            if digest(Path(name)) != sha: raise RuntimeError('input changed: '+name)
        summary['status'] = 'passed'
    except BaseException:
        summary['status'] = 'failed'
        raise
    finally:
        (root/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
