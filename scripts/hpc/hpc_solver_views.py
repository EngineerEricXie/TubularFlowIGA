#!/usr/bin/env python3
"""Validate native KSP/SNES views without changing accepted simulation outputs."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

from hpc_inventory import digest
from hpc_native_graph_checkpoint import fixture


def main():
    repo = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, default=repo / 'solvers/coupling/iga_multidomain_flow')
    parser.add_argument('--fixture-builder', type=Path, default=repo / 'solvers/coupling/native_graph_checkpoint_fixture_test')
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    binary, builder, root = args.binary.resolve(), args.fixture_builder.resolve(), args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    (root / 'fixtures').mkdir()
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps')
    report = dict(status='running', binary_sha256=digest(binary), builder_sha256=digest(builder), jobs=[], inputs={})

    def save():
        (root / 'acceptance.json').write_text(json.dumps(report, indent=2) + '\n')

    def run(name, case, flags, expected=0):
        directory = root / name
        directory.mkdir()
        result = directory / 'result'
        command = ['timeout', '--kill-after=5s', '270s', 'mpiexec', '--oversubscribe', '-np', '3']
        if expected == 0:
            command += [sys.executable, str(repo / 'scripts/hpc/hpc_rank_run.py'), '--output-dir', str(directory),
                        '--expected-ranks', '3', '--timeout', '240', '--']
        command += [str(binary), '--graph-case', str(case), '--output-dir', str(result)] + flags
        record = dict(name=name, command=command, expected=expected)
        report['jobs'].append(record)
        with (directory / 'launcher.log').open('w') as stream:
            completed = subprocess.run(command, cwd=repo, env=env, stdout=stream, stderr=subprocess.STDOUT)
        record['returncode'] = completed.returncode
        save()
        if completed.returncode != expected:
            raise RuntimeError(name + ' returned unexpected status')
        if expected:
            if result.exists():
                raise RuntimeError('failed viewer published accepted output')
        else:
            for rank in range(3):
                measured = json.loads((directory / f'rank-{rank}/run.json').read_text())
                if measured['returncode'] or measured['timed_out'] or not measured['resource']:
                    raise RuntimeError('unsuccessful or unmeasured view rank')
        return directory, result, record

    try:
        for mode in ['nonlinear_aq', 'implicit_1d_pde', 'species']:
            case = root / 'fixtures' / mode
            fixture(builder, case, 'species' if mode == 'species' else 'flow', 3)
            if mode != 'species':
                for path in case.glob('*/simulation_config.json'):
                    config = json.loads(path.read_text())
                    if config.get('dimension') != '1d':
                        continue
                    for system in config['equation_systems']:
                        if system['kind'] == 'network_flow_1d':
                            system.update(model='compliant', scheme='implicit_petsc', formulation=mode,
                                          wall=dict(model='linear', young_modulus=1e9, thickness_ratio=0.1, reference_pressure=0))
                    path.write_text(json.dumps(config, indent=2) + '\n')
                flags = ['-domain_source_flow_snes_view', '-domain_source_flow_ksp_view', '-domain_junction_flow_ksp_view']
                markers = ['SNES Object: (domain_source_flow_)', 'KSP Object: (domain_source_flow_)', 'KSP Object: (domain_junction_flow_)']
            else:
                flags = ['-domain_junction_flow_ksp_view', '-domain_junction_transport_ksp_view']
                markers = ['KSP Object: (domain_junction_flow_)', 'KSP Object: (domain_junction_transport_)']
            report['inputs'].update({str(path): digest(path) for path in case.rglob('*') if path.is_file()})
            _, baseline, _ = run(mode + '-baseline', case, [])
            directory, actual, record = run(mode + '-view', case, flags)
            log = (directory / 'rank-0/stdout.log').read_text()
            if any(marker not in log for marker in markers):
                raise RuntimeError(mode + ' is missing scoped solver view')
            record['view_markers'] = markers
            record['equal_files'] = {str(path.relative_to(baseline)): path.read_bytes() == (actual / path.relative_to(baseline)).read_bytes()
                                     for path in baseline.rglob('*') if path.is_file()}
            if {str(path.relative_to(actual)) for path in actual.rglob('*') if path.is_file()} != set(record['equal_files']):
                raise RuntimeError('view changed output file coverage')
            if not all(record['equal_files'].values()):
                raise RuntimeError(mode + ' view changed accepted output')
            save()
        failure = root / 'bad-view'
        directory, _, record = run('bad-view', root / 'fixtures/nonlinear_aq',
                                  ['-domain_junction_flow_ksp_view', 'ascii:' + str(failure / 'missing/view.txt')], 1)
        failure_log = (directory / 'launcher.log').read_text()
        if 'PETSc returned error' not in failure_log or 'MPI_ABORT was invoked' in failure_log:
            raise RuntimeError('viewer failure did not return through the PETSc scope')
        if digest(binary) != report['binary_sha256'] or digest(builder) != report['builder_sha256']:
            raise RuntimeError('view binary changed during validation')
        if any(digest(Path(path)) != value for path, value in report['inputs'].items()):
            raise RuntimeError('view fixture changed during validation')
        report['status'] = 'passed'
    except Exception as error:
        report['status'] = 'failed'
        report['error'] = str(error)
        raise
    finally:
        save()


if __name__ == '__main__':
    main()
