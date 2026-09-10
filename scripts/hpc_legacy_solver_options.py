#!/usr/bin/env python3
"""Compare legacy transport solver snapshots with an archived binary."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

from hpc_compare_fields import compare
from hpc_inventory import digest


def main():
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture', type=Path, required=True)
    parser.add_argument('--reference-binary', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    fixture, old, root = args.fixture.resolve(), args.reference_binary.resolve(), args.output_dir.resolve()
    binary = repo / 'solvers/cpu/iga_transport'
    root.mkdir(parents=True, exist_ok=False)
    report = dict(status='running', binaries={str(p): digest(p) for p in (old, binary)},
                  inputs={str(p): digest(p) for p in fixture.iterdir() if p.is_file()}, jobs=[])
    prefix = 'domain_neuron_transport_transport_'
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')

    def run(name, ranks, executable, options, negative=False):
        directory = root / name
        directory.mkdir()
        field = directory / 'field.txt'
        command = ['timeout', '--kill-after=5s', '90s', 'mpiexec', '-np', str(ranks)]
        if not negative:
            command += [sys.executable, str(repo / 'scripts/hpc_rank_run.py'), '--output-dir', str(directory),
                        '--expected-ranks', str(ranks), '--timeout', '60', '--']
        command += [str(executable), str(fixture / ('serial.ntiga' if ranks == 1 else 'group.ntiga')),
                    str(fixture), '2', str(field)]
        with (directory / 'launcher.log').open('w') as log:
            result = subprocess.run(command, cwd=repo, env=dict(env, PETSC_OPTIONS=options),
                                    stdout=log, stderr=subprocess.STDOUT)
        row = dict(name=name, ranks=ranks, command=command, options=options, returncode=result.returncode)
        report['jobs'].append(row)
        if result.returncode != int(negative):
            raise RuntimeError(name + ': unexpected exit status')
        if negative:
            log = (directory / 'launcher.log').read_text()
            if field.exists() or 'MPI_ABORT was invoked' in log or 'transport_v2 nodes=' in log:
                raise RuntimeError(name + ': uncontrolled failure or output publication')
            if 'PETSc returned error' not in log and 'backend' not in log:
                raise RuntimeError(name + ': missing failure diagnostic')
            return field
        for rank in range(ranks):
            measured = json.loads((directory / f'rank-{rank}/run.json').read_text())
            if measured['returncode'] or measured['timed_out'] or not measured['resource']:
                raise RuntimeError(name + ': invalid rank report')
        log = (directory / 'rank-0/stdout.log').read_text()
        if 'transport_v2 nodes=' not in log:
            raise RuntimeError(name + ': missing success summary')
        if executable == binary:
            rows = [dict(part.split('=', 1) for part in line.split()[1:]) for line in log.splitlines()
                    if line.startswith('solver_configuration ')]
            expected_ksp = 'fgmres' if '-'+prefix+'ksp_type fgmres' in options else 'gmres'
            if len(rows) != 2 or any(r['prefix'] != prefix or r['ksp'] != expected_ksp or int(r['reason']) <= 0 for r in rows):
                raise RuntimeError(name + ': incorrect effective solver configuration')
            row['solver_configuration'] = rows
        return field

    try:
        for ranks in (1, 2):
            for strategy, common in [('default', ''), ('lu', '-ksp_type gmres -pc_type lu -pc_factor_mat_solver_type mumps')]:
                reference = run(f'{strategy}-before-{ranks}', ranks, old, common)
                actual = run(f'{strategy}-after-{ranks}', ranks, binary, common)
                if reference.read_bytes() != actual.read_bytes():
                    raise RuntimeError('unmodified options changed default fields')
                report['jobs'][-1]['byte_equal'] = True
                actual = run(f'{strategy}-scoped-{ranks}', ranks, binary,
                             common + ' -'+prefix+'ksp_type fgmres -'+prefix+'ksp_rtol 1e-12 -'+prefix+'ksp_view')
                observation = compare(reference, actual, 1e-6, 1e-12, node_ids=True)
                report['jobs'][-1]['field_comparison'] = observation
                if not observation['passed']:
                    raise RuntimeError('scoped fields differ')
            for mode, options in [('bad-ksp', '-'+prefix+'ksp_type unavailable_ksp'),
                                  ('bad-backend', '-'+prefix+'pc_type lu -'+prefix+'pc_factor_mat_solver_type unavailable_backend'),
                                  ('bad-view', '-'+prefix+'ksp_view ascii:'+str(root/'missing'/'view.txt'))]:
                run(f'{mode}-{ranks}', ranks, binary, options, True)
            run(f'healthy-retry-{ranks}', ranks, binary, '')
        for path, expected in {**report['binaries'], **report['inputs']}.items():
            if digest(Path(path)) != expected:
                raise RuntimeError('binary or fixture changed')
        report['status'] = 'passed'
    except BaseException:
        report['status'] = 'failed'
        raise
    finally:
        (root / 'acceptance.json').write_text(json.dumps(report, indent=2) + '\n')
    print('legacy solver options:', report['status'], len(report['jobs']), 'jobs', flush=True)


if __name__ == '__main__':
    main()
