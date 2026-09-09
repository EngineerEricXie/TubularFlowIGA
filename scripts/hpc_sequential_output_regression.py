#!/usr/bin/env python3
"""Check sequential coupling close failures and unchanged native outputs."""

import argparse
import copy
import json
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
    parser.add_argument('--rank-worker', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.rank_worker:
        rank, ranks = rank_identity(os.environ)
        status = run_rank(args.command[1:], args.output_dir, ranks, 90)
        (args.output_dir/f'rank-{rank}/report-ready').touch(exist_ok=False)
        deadline = time.monotonic()+15
        while not all((args.output_dir/f'rank-{i}/report-ready').exists() for i in range(ranks)):
            if time.monotonic() >= deadline:
                return 2
            time.sleep(0.01)
        return status
    if args.case_dir is None or args.reference_binary is None or args.command:
        parser.error('require --case-dir and --reference-binary without positional arguments')
    source, out = args.case_dir.resolve(), args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    binary = ROOT/'solvers/coupling/iga_1d_3d_explicit'
    reference = args.reference_binary.resolve()
    preload = ROOT/'solvers/coupling/text_close_preload.so'
    graph = json.loads((source/'simulation_config.json').read_text())
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    env.pop('PETSC_OPTIONS', None)
    env.pop('TUBULARFLOWIGA_TEST_TEXT_CLOSE', None)
    env.pop('LD_PRELOAD', None)
    report = dict(status='running', binaries={str(p): digest(p) for p in (binary, reference, preload)},
                  cases=[], comparison=[])

    def save():
        (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')

    def run(name, ranks, fixture, executable, target='', old=False):
        case = out/name
        records = case/'ranks'
        records.mkdir(parents=True)
        result_dir = case/'result'
        local_env = dict(env)
        if target:
            local_env.update(LD_PRELOAD=str(preload), TUBULARFLOWIGA_TEST_TEXT_CLOSE=str(result_dir/target))
        command = ['timeout', '--kill-after=5s', '120s', 'mpiexec', '--map-by', 'core', '--bind-to', 'core',
                   '-np', str(ranks), sys.executable, str(Path(__file__).resolve()), '--rank-worker',
                   '--output-dir', str(records), '--', str(executable), '--graph-case', str(fixture),
                   '--output-dir', str(result_dir), '-ksp_type', 'preonly', '-pc_type', 'lu']
        expected = 1 if target and (not old or target == 'graph_binding_manifest.json.tmp') else 0
        row = dict(name=name, command=command, ranks=ranks, target=target, expected_returncode=expected)
        report['cases'].append(row)
        with (case/'launcher.log').open('x') as log:
            result = subprocess.run(command, cwd=ROOT, env=local_env, stdout=log, stderr=subprocess.STDOUT)
        row['launcher_returncode'] = result.returncode
        rank_reports = [json.loads((records/f'rank-{i}/run.json').read_text()) for i in range(ranks)]
        row['rank_returncodes'] = [r['returncode'] for r in rank_reports]
        save()
        assert result.returncode == expected and row['rank_returncodes'] == [expected]*ranks, name
        assert not any(r['timed_out'] for r in rank_reports), name
        stderr = (records/'rank-0/stderr.log').read_text()
        stdout = (records/'rank-0/stdout.log').read_text()
        if target:
            assert '[text-close-injection] fclose returned failure' in stderr, name
        if expected:
            assert 'completed ' not in stdout, name
            assert not (result_dir/'graph_binding_manifest.json').exists(), name
            if not old:
                assert 'sequential output: rank 0:' in stderr, name
        else:
            assert 'completed ' in stdout and (result_dir/'graph_binding_manifest.json').is_file(), name
        print(name, expected, flush=True)
        return result_dir

    try:
        for ranks in (1, 2):
            for scheme in ('explicit', 'fixed', 'aitken'):
                fixture = out/f'input-{ranks}-{scheme}'
                fixture.mkdir()
                config = copy.deepcopy(graph)
                config['execution']['kind'] = scheme
                config['execution']['maximum_iterations'] = 1 if scheme == 'explicit' else 50
                for domain in config['domains']:
                    name = domain['case']
                    shutil.copytree(source/name, fixture/name)
                    if 'database' in domain:
                        domain['database'] = 'database.ntiga'
                shutil.copy2(source/('one.ntiga' if ranks == 1 else 'two.ntiga'), fixture/'database.ntiga')
                (fixture/'simulation_config.json').write_text(json.dumps(config)+'\n')
                prefix = f'{ranks}-{scheme}'
                before = run(prefix+'-before', ranks, fixture, reference)
                after = run(prefix+'-after', ranks, fixture, binary)
                files = sorted(p.name for p in before.iterdir() if p.suffix in ('.csv', '.json'))
                assert files and files == sorted(p.name for p in after.iterdir() if p.suffix in ('.csv', '.json'))
                for name in files:
                    identical = (before/name).read_bytes() == (after/name).read_bytes()
                    report['comparison'].append(dict(case=prefix, file=name, identical=identical,
                                                     sha256=digest(after/name)))
                    assert identical, (prefix, name)
                for name in files:
                    target = name+'.tmp' if name == 'graph_binding_manifest.json' else name
                    if ranks == 2:
                        run(prefix+'-old-fault-'+name, ranks, fixture, reference, target, old=True)
                    run(prefix+'-fault-'+name, ranks, fixture, binary, target)
                retry = run(prefix+'-retry', ranks, fixture, binary)
                assert all((retry/name).read_bytes() == (after/name).read_bytes() for name in files)
        report['status'] = 'passed'
        save()
        print('sequential close regression passed', len(report['cases']), 'cases;',
              len(report['comparison']), 'identical files', flush=True)
        return 0
    except Exception:
        report['status'] = 'failed'
        save()
        raise


if __name__ == '__main__':
    sys.exit(main())
