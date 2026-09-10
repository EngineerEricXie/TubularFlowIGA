#!/usr/bin/env python3
"""Check immersed defaults against a saved executable and reject invalid solvers."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

from hpc_inventory import digest
from hpc_immersed_graph_regression import fixture, table, compare_ports


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    binary, reference, root = args.binary.resolve(), args.reference.resolve(), args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    frozen = {str(path): digest(path) for path in [binary, reference]}
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps')
    report = dict(status='running', binaries=frozen, jobs=[])

    def save():
        (root / 'acceptance.json').write_text(json.dumps(report, indent=2) + '\n')

    def run(name, ranks, executable, case, flags=None):
        directory = root / name
        directory.mkdir()
        result = directory / 'result'
        app = [str(executable), '--graph-case', str(case), '--output-dir', str(result)] + (flags or [])
        command = ['timeout', '--kill-after=5s', '210s', 'mpiexec', '--oversubscribe', '-np', str(ranks)]
        if flags is None:
            command += [sys.executable, str(repo / 'scripts/hpc_rank_run.py'),
                        '--output-dir', str(directory), '--expected-ranks', str(ranks), '--timeout', '180', '--']
        command += app
        record = dict(name=name, ranks=ranks, argv=command, expected_exit=1 if flags else 0)
        report['jobs'].append(record)
        with (directory / 'launcher.log').open('w') as stream:
            completed = subprocess.run(command, cwd=repo, env=env, stdout=stream, stderr=subprocess.STDOUT)
        record['returncode'] = completed.returncode
        save()
        if completed.returncode != record['expected_exit']:
            raise RuntimeError(name + ' returned unexpected exit status')
        if flags:
            log = (directory / 'launcher.log').read_text()
            if result.exists() or not ('PETSc returned error' in log or 'factor' in log):
                raise RuntimeError(name + ' did not reject cleanly before output publication')
        else:
            for rank in range(ranks):
                measured = json.loads((directory / f'rank-{rank}' / 'run.json').read_text())
                if measured['returncode'] or measured['timed_out'] or not measured['resource']:
                    raise RuntimeError(name + ' has an unsuccessful or unmeasured rank')
        print(name + ' passed', flush=True)
        return result, record

    try:
        for transient in [False, True]:
            mode = 'transient' if transient else 'static'
            folder = root / (mode + '-fixture')
            folder.mkdir()
            case = fixture(repo, folder, 'fixed', transient)
            inputs = {str(path): digest(path) for path in case.rglob('*') if path.is_file()}
            for ranks in [1, 2]:
                old, _ = run(f'{mode}-{ranks}-old', ranks, reference, case)
                new, record = run(f'{mode}-{ranks}-new', ranks, binary, case)
                record['port_error_fraction'] = compare_ports(table(old / 'pressure_flow_ports.csv'), table(new / 'pressure_flow_ports.csv'))
                record['byte_equal_outputs'] = {str(path.relative_to(old)): path.read_bytes() == (new / path.relative_to(old)).read_bytes()
                                                for path in old.rglob('*') if path.is_file()}
                if not all(record['byte_equal_outputs'].values()):
                    raise RuntimeError('default accepted outputs changed')
                for kind, key in [('ksp', 'ksp_type'), ('backend', 'pc_factor_mat_solver_type')]:
                    run(f'{mode}-{ranks}-{kind}', ranks, binary, case,
                        ['-domain_immersed_flow_' + key, 'unavailable_' + kind])
                # A new healthy process follows both failed launches.
                run(f'{mode}-{ranks}-retry', ranks, binary, case)
            if any(digest(Path(path)) != value for path, value in inputs.items()):
                raise RuntimeError('fixture changed during validation')
            report.setdefault('inputs', {}).update(inputs)
        if any(digest(Path(path)) != value for path, value in frozen.items()):
            raise RuntimeError('binary changed during validation')
        report['status'] = 'passed'
    except Exception as error:
        report['status'] = 'failed'
        report['error'] = str(error)
        raise
    finally:
        save()


if __name__ == '__main__':
    main()
