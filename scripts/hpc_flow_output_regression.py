#!/usr/bin/env python3
"""Verify native CPU flow output failures, final success ordering, and field parity."""

import argparse
import json
import os
from pathlib import Path
import resource
import shlex
import signal
import subprocess
import sys

from hpc_compare_fields import compare
from hpc_inventory import digest

ROOT = Path(__file__).resolve().parents[1]


def limited_child(command):
    def limit():
        signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
        resource.setrlimit(resource.RLIMIT_FSIZE, (64, resource.getrlimit(resource.RLIMIT_FSIZE)[1]))
    # Keep the report files outside the limited process; otherwise its diagnostic
    # would itself be truncated at 64 bytes. Pipes have no regular-file limit.
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=limit)
    sys.stdout.buffer.write(result.stdout)
    sys.stderr.buffer.write(result.stderr)
    return result.returncode


def main():
    if sys.argv[1:2] == ['--limited-child']:
        return limited_child(sys.argv[2:])
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case-dir', type=Path, required=True)
    parser.add_argument('--reference-binary', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--launcher', default='mpiexec --map-by core --bind-to core')
    parser.add_argument('--only', action='append', help='Run named cases only; a complete acceptance run omits this option')
    args = parser.parse_args()
    source, output, reference = args.case_dir.resolve(), args.output_dir.resolve(), args.reference_binary.resolve()
    output.mkdir(parents=True, exist_ok=False)
    binary = ROOT/'solvers/cpu/iga_navier_stokes'
    openmp = ROOT/'solvers/cpu/iga_navier_stokes_openmp'
    base = json.loads((source/'simulation_config.json').read_text())
    base['simulation_scope'] = {'mode': 'flow_only'}
    base.pop('coupling', None)
    base.pop('external_circuit', None)
    base['temporal_functions'] = [dict(name='inlet_scale', kind='periodic_table',
        units='dimensionless', period=1.0, file='scale.csv', interpolation='linear')]
    base['boundaries'][0]['conditions'][0].update(waveform='inlet_scale', scale=0.001)
    payload = {name: (source/name).read_bytes() for name in ('controlmesh.vtk', 'initial_velocityfield.txt')}
    payload['simulation_config.json'] = (json.dumps(base)+'\n').encode()
    payload['scale.csv'] = b'time,value\n0,1\n0.5,1\n'
    # name, ranks, format, output frequency, damage target/type, executable
    cases = [('vtu-before', 1, 'vtu', 1, '', '', reference)]
    cases += [(f'{kind}-{ranks}', ranks, fmt, every, '', '', binary)
        for kind, fmt, every in [('vtu', 'vtu', 1), ('final', 'vtu', 0), ('none', 'vtu', -1)]
        for ranks in (1, 2)]
    cases += [('vtkhdf-before', 1, 'vtkhdf', 1, '', '', reference)]
    cases += [(f'vtkhdf-{ranks}', ranks, 'vtkhdf', 1, '', '', binary) for ranks in (1, 2)]
    cases += [(f'openmp-{ranks}', ranks, 'vtu', 1, '', '', openmp) for ranks in (1, 2)]
    targets = {}
    for phase, stem in [('initial', 'flow.step000000'), ('step1', 'flow.step000001'),
                         ('step2', 'flow.step000002'), ('final', 'flow')]:
        for part, suffix in [('velocity', '.txt'), ('pressure', '.txt.pressure'), ('vtu', '.vtu')]:
            targets[phase+'-'+part] = stem+suffix
    targets.update(index='flow.pvd', series='flow.series.csv',
                   geometry='flow.bezier_geometry.json', hdf='flow.vtkhdf')
    for target in targets:
        for damage in ('directory', 'fifo'):
            cases.append((target+'-'+damage, 2, 'vtkhdf' if target in ('geometry', 'hdf') else 'vtu',
                          1, target, damage, binary))
    # Multi-rank file limits are injected after MPI initialization in the native writer test.
    cases += [('file-limit-1', 1, 'vtu', 1, '', 'limit', binary)]
    cases += [('old-final-failure', 1, 'vtu', 0, 'final-velocity', 'directory', reference)]
    if args.only:
        unknown = set(args.only)-{case[0] for case in cases}
        if unknown: raise ValueError(f'unknown cases: {sorted(unknown)}')
        selected = set(args.only)
        if any(case[0] in selected and not case[5] and case[3] >= 0 for case in cases):
            selected.add('vtu-before')
        cases = [case for case in cases if case[0] in selected]
    summary = dict(status='running', requested_cases=args.only, binaries={str(p): digest(p) for p in (binary, reference, openmp)}, cases=[])
    try:
        for name, ranks, fmt, every, target, damage, executable in cases:
            case = output/name
            local, results, records = case/'input', case/'result', case/'ranks'
            local.mkdir(parents=True)
            results.mkdir()
            records.mkdir()
            for filename, data in payload.items(): (local/filename).write_bytes(data)
            (local/'database.ntiga').write_bytes((source/('fixture.ntiga' if ranks == 1 else 'fixture-2.ntiga')).read_bytes())
            inputs = {str(p): digest(p) for p in local.iterdir()}
            if target:
                broken = results/targets[target]
                if damage == 'directory': broken.mkdir()
                else: os.mkfifo(broken, 0o600)
            child = [str(executable), str(local/'database.ntiga'), str(local), '--max-newton', '12',
                     '--visualization-format', fmt]
            if every >= 0: child += ['--output', str(results/'flow.txt')]
            if every > 0: child += ['--output-every', str(every)]
            if damage == 'limit': child = [sys.executable, str(Path(__file__).resolve()), '--limited-child', *child]
            env = dict(os.environ, OMP_NUM_THREADS='2' if executable == openmp else '1',
                IGA_ASSEMBLY_THREADS='2' if executable == openmp else '1', OPENBLAS_NUM_THREADS='1',
                MKL_NUM_THREADS='1', BLIS_NUM_THREADS='1',
                PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -ksp_rtol 1e-12')
            worker = ([sys.executable, str(ROOT/'scripts/hpc_flow_input_regression.py'), '--rank-worker',
                       '--output-dir', str(records)] if ranks == 2 else
                      [sys.executable, str(ROOT/'scripts/hpc_rank_run.py'), '--output-dir', str(records), '--timeout', '60'])
            command = ['timeout', '--kill-after=5s', '90s', *shlex.split(args.launcher), '-np', str(ranks), *worker, '--', *child]
            stage = ('flow visualization initialization' if target in ('geometry', 'hdf') else
                     'flow output index' if target in ('index', 'series') else 'flow field output') if damage else ''
            row = dict(case=name, command_argv=command, input_sha256=inputs, expected_stage=stage,
                       expected_returncode=1 if damage else 0, status='running', ranks=[], fields=[])
            summary['cases'].append(row)
            with (case/'launcher.log').open('x') as log:
                result = subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
            row['launcher_returncode'] = result.returncode
            if result.returncode != row['expected_returncode']: raise RuntimeError(f'{name}: unexpected exit {result.returncode}')
            for rank in range(ranks):
                rd = records/f'rank-{rank}'
                report = json.loads((rd/'run.json').read_text())
                row['ranks'].append(report)
                if report['returncode'] != row['expected_returncode'] or report['timed_out']: raise RuntimeError('bad rank exit')
                if report['rank'] != rank or report['ranks'] != ranks: raise RuntimeError('wrong rank identity')
                for filename, checksum in report['logs'].items():
                    if digest(rd/filename) != checksum: raise RuntimeError('rank log changed')
                text = (rd/'stdout.log').read_text()+(rd/'stderr.log').read_text()
                success = 'navier_stokes_v2 seconds=' in text
                if name == 'old-final-failure':
                    if not success: raise RuntimeError('pre-change success-order defect was not reproduced')
                    row['observed_false_success_summary'] = True
                elif damage:
                    if stage+': rank 0:' not in text or success: raise RuntimeError('missing common failure or false success summary')
                elif rank == 0 and not success: raise RuntimeError('missing successful solve summary')
            if not damage and every >= 0:
                files = ['flow.txt', 'flow.txt.pressure']
                if every > 0:
                    files += [f'flow.step{step:06d}.txt{suffix}' for step in range(3) for suffix in ('', '.pressure')]
                for field in files:
                    comparison = compare(output/'vtu-before/result'/field, results/field, 1e-6, 1e-12)
                    row['fields'].append({'file': field, **comparison})
                    if not comparison['passed']: raise RuntimeError('flow output field changed')
            if target.startswith('step1-') or target.startswith('step2-'):
                previous = 0 if target.startswith('step1-') else 1
                if not (results/f'flow.step{previous:06d}.vtu').is_file(): raise RuntimeError('previous accepted step missing')
                if (results/'flow.txt').exists(): raise RuntimeError('final output followed a failed step')
            if damage == 'limit':
                if (results/'flow.step000000.txt').stat().st_size != 64: raise RuntimeError('file limit fault missing')
            if any(digest(Path(p)) != checksum for p, checksum in inputs.items()): raise RuntimeError('inputs changed')
            row['status'] = 'passed'
            (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
            print(name, 'passed', flush=True)
        if any(digest(Path(p)) != checksum for p, checksum in summary['binaries'].items()): raise RuntimeError('binaries changed')
        summary['status'] = 'passed'
    finally:
        (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main_result = main()
    sys.exit(main_result)
