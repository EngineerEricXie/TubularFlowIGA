#!/usr/bin/env python3
"""Check native 1D and coupling stdout failures with per-rank exit evidence."""
import argparse
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import time
from hpc_inventory import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('baseline-dir', 'graph-fixtures', 'sequential-fixtures', 'output-dir'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--launcher', default='mpiexec --oversubscribe')
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    one_d = root/'one-d-input'
    shutil.copytree(repo/'examples/one_d/rigid_straight', one_d)
    config = json.loads((one_d/'simulation_config.json').read_text())
    config['time'].update(steps=2, output_every=1)
    (one_d/'simulation_config.json').write_text(json.dumps(config)+'\n')
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps')
    summary = dict(status='running', cases=[], comparisons=[], binaries={}, inputs={})
    specs = []
    for ranks in (1, 3):
        for check in (True, False):
            specs.append(dict(name=f'one-d-{ranks}-'+('check' if check else 'solve'), ranks=ranks,
                folder='one_d', native='iga_1d', wrapper='one_d_stdout_failure_test', fixture=one_d,
                check=check, target='schema_version=3 dimension=' if check else 'completed 1d system=',
                stage='1d check output' if check else '1d final output'))
    for ranks, group in ((1, 'single'), (2, 'pair')):
        for kind in ('flow', 'species', 'zero-d', 'bifurcation'):
            bif = kind == 'bifurcation'
            specs.append(dict(name=f'graph-{ranks}-{kind}', ranks=ranks, folder='coupling',
                native='iga_1d_3d_bifurcation' if bif else 'iga_multidomain_flow',
                wrapper='bifurcation_stdout_failure_test' if bif else 'graph_stdout_failure_test',
                fixture=args.graph_fixtures.resolve()/group/('flow' if bif else kind), check=False,
                target='completed schema-v5 1D--3D bifurcation steps=' if bif else 'completed schema-v',
                stage='graph completion logging'))
        for scheme in ('explicit', 'fixed', 'aitken'):
            specs.append(dict(name=f'sequential-{ranks}-{scheme}', ranks=ranks, folder='coupling',
                native='iga_1d_3d_explicit', wrapper='sequential_stdout_failure_test',
                fixture=args.sequential_fixtures.resolve()/f'input-{ranks}-{scheme}', check=False,
                target='completed ', stage='sequential completion logging'))
    for spec in specs:
        for path in spec['fixture'].rglob('*'):
            if path.is_file(): summary['inputs'][str(path)] = digest(path)

    def run(spec, label, old=False, mode=None):
        directory = root/(spec['name']+'-'+label)
        directory.mkdir()
        (directory/'ranks').mkdir()
        output = directory/'result'
        binary = (args.baseline_dir.resolve() if old else repo/'solvers'/spec['folder'])/(spec['wrapper'] if mode else spec['native'])
        argv = ([str(spec['fixture'])] if spec['folder'] == 'one_d' else ['--graph-case', str(spec['fixture'])])
        argv += ['--output-dir', str(output)]
        if spec['check']: argv += ['--check']
        command = ['timeout', '--kill-after=5s', '120s', *shlex.split(args.launcher), '-np', str(spec['ranks']),
                   'python3', str(repo/'scripts/hpc_sequential_output_regression.py'), '--rank-worker',
                   '--output-dir', str(directory/'ranks'), '--', str(binary), *argv]
        settings = dict(env)
        if mode: settings.update(IGA_TEST_STDOUT_TARGET=spec['target'], IGA_TEST_STDOUT_MODE=mode)
        expected = 1 if mode else 0
        record = dict(name=directory.name, argv=command, expected=expected, mode=mode,
                      target=spec['target'], stage=spec['stage'], ranks=spec['ranks'])
        summary['cases'].append(record)
        summary['binaries'][str(binary)] = digest(binary)
        start = time.monotonic()
        with (directory/'launcher.log').open('w') as log:
            result = subprocess.run(command, env=settings, stdout=log, stderr=subprocess.STDOUT)
        record.update(returncode=result.returncode, elapsed_s=time.monotonic()-start)
        if result.returncode != expected: raise RuntimeError(f'{directory.name}: exit {result.returncode}, expected {expected}')
        reports = []
        for rank in range(spec['ranks']):
            base = directory/'ranks'/f'rank-{rank}'
            report = json.loads((base/'run.json').read_text())
            if report['returncode'] != expected or report['timed_out'] or report['resource'] is None:
                raise RuntimeError(f'{directory.name}: invalid rank {rank} outcome')
            if mode:
                text = (base/'stderr.log').read_text()
                if f'stdout_test rank={rank} injected={int(rank == 0)} native_status=1' not in text:
                    raise RuntimeError(directory.name+': injection/common exit failed')
                if rank == 0 and spec['stage']+': rank 0:' not in text:
                    raise RuntimeError(directory.name+': wrong failure stage')
            reports.append(report)
        record.update(status='passed', rank_reports=reports)
        print(directory.name, 'passed', flush=True)
        return directory

    def compare(spec, before, after):
        if spec['check']:
            old = (before/'ranks/rank-0/stdout.log').read_text().splitlines()
            new = (after/'ranks/rank-0/stdout.log').read_text().splitlines()
            a = [x for x in old if x.startswith(spec['target'])]
            b = [x for x in new if x.startswith(spec['target'])]
            if len(a) != 1 or a != b: raise RuntimeError('1D check summary changed')
            summary['comparisons'].append(dict(case=after.name, kind='check-summary', exact=True))
            return
        left, right = before/'result', after/'result'
        files = {str(p.relative_to(left)) for p in left.rglob('*') if p.is_file()}
        if not files or files != {str(p.relative_to(right)) for p in right.rglob('*') if p.is_file()}:
            raise RuntimeError(after.name+': output set changed')
        for name in sorted(files):
            a, b = (left/name).read_bytes(), (right/name).read_bytes()
            # Only the 1D summary contains timing/RSS observations. Preserve its
            # other numerical/convergence fields in the comparison.
            if spec['folder'] == 'one_d' and name == 'summary.json':
                a, b = json.loads(a), json.loads(b)
                for key in ('setup_seconds', 'solve_seconds', 'output_seconds', 'peak_rss_kib'):
                    for value in (a.pop(key), b.pop(key)):
                        if not math.isfinite(value) or value < 0 or (key == 'peak_rss_kib' and value == 0):
                            raise RuntimeError('invalid timing/RSS observation')
            if a != b: raise RuntimeError(after.name+': changed '+name)
            summary['comparisons'].append(dict(case=after.name, file=name, exact=True,
                scope='excluding validated 1D timing/RSS observations' if name == 'summary.json' else 'all bytes'))

    try:
        for spec in specs:
            before = run(spec, 'before', old=True)
            compare(spec, before, run(spec, 'healthy'))
            for mode in ('short', 'flush', 'exception', 'allocation', 'nonstandard'):
                run(spec, mode, mode=mode)
                compare(spec, before, run(spec, mode+'-retry'))
        for name, sha in summary['inputs'].items():
            if digest(Path(name)) != sha: raise RuntimeError('fixture changed: '+name)
        summary['status'] = 'passed'
    except BaseException:
        summary['status'] = 'failed'
        raise
    finally:
        (root/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
