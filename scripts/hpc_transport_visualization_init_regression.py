#!/usr/bin/env python3
"""Verify collective transport visualization initialization with actual CLI failures."""

import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import time

from hpc_inventory import digest
from hpc_rank_run import rank_identity, run as run_rank

ROOT = Path(__file__).resolve().parents[1]
STAGE = 'transport visualization initialization: rank 0:'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixtures', type=Path, help='Completed hpc_transport_cli_regression.py directory')
    parser.add_argument('--reference-binary', type=Path)
    parser.add_argument('--preload', type=Path)
    parser.add_argument('--hdf-reader', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--launcher', default='mpiexec --map-by core --bind-to core')
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--ranks', type=int, default=1, help=argparse.SUPPRESS)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.worker:
        command = args.command[1:] if args.command[:1] == ['--'] else args.command
        status = run_rank(command, args.output_dir, args.ranks, 20)
        rank, _ = rank_identity(os.environ)
        (args.output_dir/f'rank-{rank}/ready').touch(exist_ok=False)
        deadline = time.monotonic()+30
        while not all((args.output_dir/f'rank-{i}/ready').exists() for i in range(args.ranks)):
            if time.monotonic() > deadline:
                return 2
            time.sleep(0.01)
        return status
    if not all((args.fixtures, args.reference_binary, args.preload, args.hdf_reader)) or args.command:
        parser.error('require fixtures, reference binary, preload and HDF reader')
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    reference = args.reference_binary.resolve()
    binary = ROOT/'solvers/cpu/iga_solve'
    preload, reader = args.preload.resolve(), args.hdf_reader.resolve()
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               PETSC_OPTIONS='-ksp_type gmres -pc_type lu -pc_factor_mat_solver_type mumps -ksp_rtol 1e-12')
    summary = dict(status='running', binaries={str(p): digest(p) for p in (reference, binary, preload, reader)}, cases=[])
    faults = ['error', 'bad-alloc', 'nonstandard', 'report-directory', 'report-fifo', 'hdf-directory', 'hdf-fifo']
    try:
        for ranks in (1, 2):
            source = args.fixtures.resolve()/('vtkhdf-one' if ranks == 1 else 'vtkhdf-two')
            modes = [('before', '', True), ('healthy', '', False), ('old-nonstandard', 'nonstandard', True)]
            for fault in faults:
                modes += [(fault, fault, False), ('retry-'+fault, '', False)]
            for mode, fault, old in modes:
                name = f'{ranks}-{mode}'
                case = output/name
                records, result = case/'ranks', case/'result'
                records.mkdir(parents=True)
                result.mkdir()
                if fault.endswith(('directory', 'fifo')):
                    target = result/('field.bezier_geometry.json' if fault.startswith('report') else 'field.vtkhdf')
                    if fault.endswith('fifo'):
                        os.mkfifo(target, 0o600)
                    else:
                        target.mkdir()
                command = ['timeout', '--kill-after=5s', '60s', *shlex.split(args.launcher)]
                inputs = {}
                for rank in range(ranks):
                    local = case/f'input-{rank}'
                    shutil.copytree(source/f'input-{rank}', local)
                    inputs.update({str(p): digest(p) for p in local.iterdir() if p.is_file()})
                    child = [str(reference if old else binary), str(local/'database.ntiga'), str(local),
                             '--system', 'transport', '--output', str(result/'field.txt'),
                             '--output-every', '1', '--visualization-format', 'vtkhdf']
                    if fault in ('error', 'bad-alloc', 'nonstandard'):
                        child = ['env', 'LD_PRELOAD='+str(preload),
                                 'TUBULARFLOWIGA_TEST_HDF_CREATE_PATH='+str(result/'field.vtkhdf'),
                                 'TUBULARFLOWIGA_TEST_HDF_CREATE_MODE='+fault, *child]
                    if rank:
                        command.append(':')
                    command += ['-np', '1', sys.executable, str(Path(__file__).resolve()), '--worker',
                                '--ranks', str(ranks), '--output-dir', str(records), '--', *child]
                entry = dict(name=name, rank_count=ranks, fault=fault, old=old, command=command,
                             input_sha256=inputs, ranks=[], status='running')
                summary['cases'].append(entry)
                with (case/'launcher.log').open('x') as log:
                    run = subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
                entry['launcher_returncode'] = run.returncode
                for rank in range(ranks):
                    folder = records/f'rank-{rank}'
                    report = json.loads((folder/'run.json').read_text())
                    assert report['rank'] == rank and report['ranks'] == ranks, name
                    assert all(digest(folder/p) == h for p, h in report['logs'].items()), name
                    entry['ranks'].append(report)
                    stderr = (folder/'stderr.log').read_text()
                    if mode != 'old-nonstandard':
                        expected = 1 if fault else 0
                        assert report['returncode'] == expected and not report['timed_out'], (name, rank)
                        assert run.returncode == expected, name
                        if fault:
                            assert STAGE in stderr, (name, rank)
                    if fault:
                        assert 'iga_solve system=' not in (folder/'stdout.log').read_text(), name
                        if rank == 0 and fault in ('error', 'bad-alloc', 'nonstandard'):
                            assert '[hdf-create-injection] '+fault in stderr, name
                if mode == 'old-nonstandard':
                    assert run.returncode != 0 and entry['ranks'][0]['returncode'] != 1, name
                if fault:
                    assert not list(result.glob('field*.txt')), name
                    assert not (result/'physiology_fields.json').exists(), name
                elif mode != 'before':
                    baseline = output/f'{ranks}-before'/'result'
                    names = sorted(p.name for p in baseline.glob('field*.txt*'))
                    assert names and names == sorted(p.name for p in result.glob('field*.txt*')), name
                    for filename in names:
                        assert (baseline/filename).read_bytes() == (result/filename).read_bytes(), (name, filename)
                    entry['identical_text_files'] = names
                    with (case/'hdf-comparison.log').open('x') as log:
                        checked = subprocess.run([str(reader), str(baseline/'field.vtkhdf'), str(result/'field.vtkhdf')],
                                                 stdout=log, stderr=subprocess.STDOUT)
                    assert checked.returncode == 0, name
                    entry['hdf_comparison'] = 'passed'
                assert all(digest(Path(p)) == h for p, h in inputs.items()), name
                entry['status'] = 'passed'
                print(name, 'passed', flush=True)
        assert all(digest(Path(p)) == h for p, h in summary['binaries'].items())
        summary['status'] = 'passed'
    except BaseException:
        summary['status'] = 'failed'
        raise
    finally:
        (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
