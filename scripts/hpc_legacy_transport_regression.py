#!/usr/bin/env python3
"""Validate legacy transport replicas, controlled input errors and output compatibility."""

import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

from hpc_compare_fields import compare
from hpc_inventory import digest

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture', type=Path, required=True,
                        help='Legacy-compatible fixture with serial.ntiga/group.ntiga and labels 0/1/2')
    parser.add_argument('--reference-binary', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--launcher', default='mpiexec --map-by core --bind-to core')
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    fixture, old = args.fixture.resolve(), args.reference_binary.resolve()
    binary = ROOT/'solvers/cpu/iga_transport'
    payload = {name: (fixture/name).read_bytes() for name in
               ('simulation_parameter.txt', 'controlmesh.vtk', 'initial_velocityfield.txt')}
    config = json.dumps({'schema_version': 1, 'boundaries': {'inherit_legacy': True, 'conditions': []}})+'\n'
    targets = dict(database='database.ntiga', parameters='simulation_parameter.txt', mesh='controlmesh.vtk',
                   velocity='initial_velocityfield.txt', configuration='case_config.json')
    # name, ranks, archived executable, expected diagnostic.
    cases = [('before-1', 1, True, ''), ('before-2', 2, True, ''),
             ('healthy-1', 1, False, ''), ('healthy-2', 2, False, ''),
             ('configuration-before', 1, True, ''), ('configuration', 2, False, ''),
             ('override', 2, False, ''), ('unused-velocity', 2, False, ''),
             ('no-output', 2, False, ''), ('zero-before', 1, True, ''), ('zero', 2, False, '')]
    for asset in targets:
        for damage in ('different', 'missing', 'fifo', 'directory'):
            stage = ('asset catalog agreement' if asset == 'configuration' and damage == 'missing' else
                     'legacy transport database input' if asset == 'database' and damage != 'different' else
                     'legacy transport asset '+asset if damage == 'different' else 'asset content read')
            cases.append((asset+'-'+damage, 2, False, stage+': rank 1:'))
    cases += [('parameters-different-before', 2, True, ''),
              ('petsc-env-different', 2, False, ''),
              ('petsc-invalid', 2, False, 'legacy transport KSPSetFromOptions: rank 0:'),
              ('output-fifo', 2, False, 'legacy transport output: rank 0:'),
              ('output-directory', 2, False, 'legacy transport output: rank 0:')]
    for asset in ('parameters', 'mesh', 'velocity', 'configuration'):
        cases.append((asset+'-malformed', 2, False, 'legacy transport case input: rank 0:'))
    for value in ('-1', '2junk', '2147483648'):
        cases.append(('steps-'+value, 2, False, 'legacy transport case input: rank 1:'))
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               BLIS_NUM_THREADS='1', PETSC_OPTIONS='-ksp_type gmres -pc_type lu -pc_factor_mat_solver_type mumps')
    summary = dict(status='running', binaries={str(p): digest(p) for p in (binary, old)}, cases=[])
    try:
        for name, ranks, previous, stage in cases:
            case = output/name
            results, records = case/'result', case/'ranks'
            results.mkdir(parents=True)
            records.mkdir()
            inputs = {}
            command = ['timeout', '--kill-after=5s', '90s', *shlex.split(args.launcher)]
            for rank in range(ranks):
                local = case/f'input-{rank}'
                local.mkdir()
                files = dict(payload)
                files['database.ntiga'] = (fixture/('serial.ntiga' if ranks == 1 else 'group.ntiga')).read_bytes()
                if name.startswith('configuration'):
                    files['case_config.json'] = config.encode()
                asset, _, damage = name.partition('-')
                if damage == 'malformed': files[targets[asset]] = b'invalid\n'
                if rank == 1 and damage.startswith('different') and asset in targets:
                    if asset == 'database':
                        # Preserve the binary parser's inputs but change the asset identity.
                        files[targets[asset]] += b'\n'
                    elif asset == 'parameters':
                        files[targets[asset]] = files[targets[asset]].replace(b'D 0.1\n', b'D 0.2\n')
                    else: files[targets[asset]] += b'\n'
                for filename, contents in files.items(): (local/filename).write_bytes(contents)
                if name in ('override', 'unused-velocity'):
                    shutil.copyfile(local/'initial_velocityfield.txt', local/'override.txt')
                    if name == 'unused-velocity':
                        (local/'initial_velocityfield.txt').unlink()
                        os.mkfifo(local/'initial_velocityfield.txt', 0o600)
                if rank == 1 and asset in targets and damage in ('missing', 'fifo', 'directory'):
                    target = local/targets[asset]
                    target.unlink()
                    if damage == 'fifo': os.mkfifo(target, 0o600)
                    if damage == 'directory': target.mkdir()
                for path in local.iterdir():
                    if path.is_file(): inputs[str(path)] = digest(path)
                child = [str(old if previous else binary), str(local/'database.ntiga'), str(local),
                         '0' if name.startswith('zero') else name[6:] if rank == 1 and name.startswith('steps-') else '2']
                if name != 'no-output': child += [str(results/'field.txt')]
                if name in ('override', 'unused-velocity'): child += [str(local/'override.txt')]
                if rank == 0 and name.startswith('output-'):
                    if name == 'output-fifo': os.mkfifo(results/'field.txt', 0o600)
                    else: (results/'field.txt').mkdir()
                settings = []
                if name == 'petsc-env-different' and rank == 1:
                    settings = ['PETSC_OPTIONS='+env['PETSC_OPTIONS']+' -ksp_rtol 1e-7']
                if name == 'petsc-invalid': settings = ['PETSC_OPTIONS=-ksp_type nonexistent']
                if rank: command += [':']
                command += ['-np', '1', 'env', *settings, sys.executable, str(ROOT/'scripts/hpc_rank_run.py'),
                            '--output-dir', str(records), '--expected-ranks', str(ranks), '--timeout', '60', '--', *child]
            expected = int(bool(stage))
            row = dict(case=name, command_argv=command, input_sha256=inputs, expected_returncode=expected,
                       expected_stage=stage, ranks=[], fields=[], status='running')
            summary['cases'].append(row)
            with (case/'launcher.log').open('x') as log:
                result = subprocess.run(command, env=env, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
            row['launcher_returncode'] = result.returncode
            if result.returncode != expected: raise RuntimeError(f'{name}: unexpected exit {result.returncode}')
            for rank in range(ranks):
                directory = records/f'rank-{rank}'
                report = json.loads((directory/'run.json').read_text())
                row['ranks'].append(report)
                if (report['returncode'] != expected or report['timed_out']
                        or report['rank'] != rank or report['ranks'] != ranks):
                    raise RuntimeError(f'{name}: bad rank report')
                text = (directory/'stdout.log').read_text()+(directory/'stderr.log').read_text()
                success = 'transport_v2 nodes=' in text
                if stage and (stage not in text or success): raise RuntimeError(f'{name}: failure not coordinated')
                if not stage and rank == 0 and not success: raise RuntimeError(f'{name}: missing summary')
                if any(digest(directory/f) != h for f, h in report['logs'].items()): raise RuntimeError('log changed')
            if not stage and name != 'no-output':
                if name.startswith('zero'): reference = output/'zero-before/result/field.txt'
                elif name.startswith('configuration'): reference = output/'configuration-before/result/field.txt'
                else: reference = output/f'before-{ranks}/result/field.txt'
                observation = compare(reference, results/'field.txt', 1e-6, 1e-12, node_ids=True)
                row['fields'].append(observation)
                if name == 'parameters-different-before':
                    if observation['passed']: raise RuntimeError('old unequal-parameter effect was not reproduced')
                    row['observed_inconsistent_input_success'] = True
                elif not observation['passed']: raise RuntimeError(f'{name}: numerical regression')
            if any(digest(Path(f)) != h for f, h in inputs.items()): raise RuntimeError('input changed')
            row['status'] = 'passed'
            (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
            print(name, 'passed', flush=True)
        if any(digest(Path(f)) != h for f, h in summary['binaries'].items()): raise RuntimeError('binary changed')
        summary['status'] = 'passed'
    except BaseException:
        summary['status'] = 'failed'
        raise
    finally:
        (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
