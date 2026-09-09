#!/usr/bin/env python3
"""Exercise sequential rank-local input failures and identical local replicas."""

import argparse
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from hpc_inventory import digest
from hpc_rank_run import rank_identity, run as run_rank

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case-dir', type=Path)
    parser.add_argument('--reference-binary', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--commands', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.commands:
        rank, ranks = rank_identity(os.environ)
        entry = json.loads(args.commands.read_text())[rank]
        os.environ['PETSC_OPTIONS'] = entry['petsc_options']
        status = run_rank(entry['argv'], args.output_dir, ranks, 60)
        (args.output_dir/f'rank-{rank}/report-ready').touch(exist_ok=False)
        deadline = time.monotonic()+15
        while not all((args.output_dir/f'rank-{i}/report-ready').exists() for i in range(ranks)):
            if time.monotonic() >= deadline:
                return 2
            time.sleep(0.01)
        return status
    if args.case_dir is None or args.reference_binary is None:
        parser.error('require --case-dir and --reference-binary')
    source, out = args.case_dir.resolve(), args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    binary = ROOT/'solvers/coupling/iga_1d_3d_explicit'
    reference = args.reference_binary.resolve()
    graph = json.loads((source/'simulation_config.json').read_text())
    graph['execution'].update(kind='explicit', maximum_iterations=1)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    solver = '-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -ksp_rtol 1e-12'
    summary = dict(status='running', binaries={str(p): digest(p) for p in (binary, reference)}, cases=[])

    def save():
        (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')

    def fixture(path, ranks):
        path.mkdir()
        for name in ('three_d', 'upstream', 'downstream'):
            shutil.copytree(source/name, path/name)
        shutil.copy2(source/('one.ntiga' if ranks == 1 else 'two.ntiga'), path/'one.ntiga')
        (path/'simulation_config.json').write_text(json.dumps(graph)+'\n')
        for domain in ('upstream', 'three_d'):
            config_path = path/domain/'simulation_config.json'
            config = json.loads(config_path.read_text())
            if domain == 'upstream':
                function = config['temporal_functions'][0]
                value = function.pop('value')
                function.update(kind='periodic_table', period=1.0, file='inlet.csv', interpolation='linear')
            else:
                value = 1.0
                function = dict(name='inlet_scale', kind='periodic_table', units='dimensionless',
                                period=1.0, file='inlet.csv', interpolation='linear')
                config['temporal_functions'] = [function]
                next(b for b in config['boundaries'] if b['label'] == 1)['conditions'][0]['waveform'] = 'inlet_scale'
            config_path.write_text(json.dumps(config)+'\n')
            (path/domain/'inlet.csv').write_text(f'time,value\n0,{value:.17g}\n0.5,{value:.17g}\n')

    baselines = {}

    def run(name, ranks=2, graph_mode=False, old=False, replicas=True, target='', damage='', control=''):
        case = out/name
        case.mkdir()
        inputs = []
        for i in range(ranks if replicas else 1):
            local = case/f'input-{i}'
            fixture(local, ranks)
            inputs.append(local)
        if not replicas:
            inputs *= ranks
        if target:
            broken = inputs[-1]/target
            if damage == 'different':
                with broken.open('ab') as stream:
                    stream.write(b'\n')
            else:
                broken.unlink()
                if damage == 'directory':
                    broken.mkdir()
                elif damage == 'fifo':
                    os.mkfifo(broken, 0o600)
        if control == 'unused-table':
            for i, local in enumerate(inputs):
                path = local/'upstream/simulation_config.json'
                config = json.loads(path.read_text())
                config['temporal_functions'].append(dict(name='unused', kind='periodic_table', units='m3/s',
                    period=1.0, file='unused.csv', interpolation='linear'))
                path.write_text(json.dumps(config)+'\n')
                (local/'upstream/unused.csv').write_text(f'unparsed and rank-specific: {i}\n')
        # PETSc 3.15 normalizes these environment entries during initialization.
        # Effective rank-local options are injected after initialization in the
        # native sequential_initialization_failure_test, not inferred here.
        expected = 1 if target or (control and control not in ('unused-table', 'petsc', 'unused-petsc')) else 0
        result_dir = case/'result'
        commands = []
        for i, local in enumerate(inputs):
            use_graph = not graph_mode if control == 'input-mode' and i == ranks-1 else graph_mode
            command = [str(reference if old else binary)]
            if use_graph:
                command += ['--graph-case', str(local)]
            else:
                command += [str(local/'one.ntiga'), str(local/'three_d'), str(local/'upstream'),
                            str(local/'downstream'), '--upstream-terminal-node', '2']
            command += ['--output-dir', str(result_dir)]
            options = solver
            if i == ranks-1:
                if control == 'argument':
                    command += ['--invalid-sequential-option']
                elif control == 'stop':
                    command += ['--stop-after-step', '1']
                elif control == 'newton':
                    command += ['--three-d-max-newton', '31']
                elif control == 'petsc':
                    options = solver.replace('-pc_type lu', '-pc_type bjacobi')
                elif control == 'unused-petsc':
                    options += ' -sequential_unused_probe rank1'
            commands.append(dict(argv=command, petsc_options=options))
        (case/'commands.json').write_text(json.dumps(commands, indent=2)+'\n')
        records = case/'ranks'
        records.mkdir()
        launcher = ['timeout', '--kill-after=5s', '90s', 'mpiexec', '--map-by', 'core', '--bind-to', 'core',
                    '-np', str(ranks), sys.executable, str(Path(__file__).resolve()),
                    '--output-dir', str(records), '--commands', str(case/'commands.json')]
        row = dict(name=name, ranks=ranks, target=target, damage=damage, control=control,
                   expected_returncode=expected, command=launcher)
        summary['cases'].append(row)
        with (case/'launcher.log').open('x') as log:
            result = subprocess.run(launcher, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
        row['launcher_returncode'] = result.returncode
        reports = [json.loads((records/f'rank-{i}/run.json').read_text()) for i in range(ranks)]
        row['rank_returncodes'] = [r['returncode'] for r in reports]
        save()
        assert result.returncode == expected and row['rank_returncodes'] == [expected]*ranks, name
        assert not any(r['timed_out'] for r in reports), name
        if expected:
            assert not result_dir.exists(), name
            stderr = (records/'rank-0/stderr.log').read_text()
            assert 'rank 1:' in stderr or 'value differs within communicator' in stderr, name
        else:
            history = result_dir/'explicit_coupling_history.csv'
            with history.open() as stream:
                rows = list(csv.DictReader(stream))
            assert len(rows) == 3 and all(math.isfinite(float(v)) for row in rows for v in row.values()), name
            key = (ranks, graph_mode)
            if old:
                baselines[key] = history
            else:
                assert history.read_bytes() == baselines[key].read_bytes(), name
                row['csv_identical_to'] = str(baselines[key])
                row['csv_sha256'] = digest(history)
        print(name, expected, flush=True)

    try:
        for ranks in (1, 2):
            for graph_mode in (False, True):
                prefix = f'{ranks}-'+('graph' if graph_mode else 'positional')
                run(prefix+'-before', ranks, graph_mode, old=True, replicas=False)
                run(prefix+'-after', ranks, graph_mode, replicas=False)
                if ranks == 2:
                    run(prefix+'-replicas', ranks, graph_mode)
                    run(prefix+'-unused-table', ranks, graph_mode, control='unused-table')
        for target in ['simulation_config.json', 'upstream/simulation_config.json',
                       'downstream/simulation_config.json', 'three_d/simulation_config.json', 'one.ntiga',
                       'three_d/controlmesh.vtk', 'three_d/initial_velocityfield.txt', 'upstream/tree.swc',
                       'downstream/tree.swc', 'upstream/inlet.csv', 'three_d/inlet.csv']:
            run('different-'+target.replace('/', '-'), graph_mode=target == 'simulation_config.json',
                target=target, damage='different')
        for damage in ('missing', 'directory', 'fifo'):
            run('table-'+damage, target='upstream/inlet.csv', damage=damage)
        for control in ('argument', 'input-mode', 'stop', 'newton', 'petsc', 'unused-petsc'):
            run('control-'+control, control=control)
        summary['status'] = 'passed'
        save()
        print('sequential input regression passed', len(summary['cases']), 'cases', flush=True)
        return 0
    except Exception:
        summary['status'] = 'failed'
        save()
        raise


if __name__ == '__main__':
    sys.exit(main())
