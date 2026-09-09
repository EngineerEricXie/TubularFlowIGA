#!/usr/bin/env python3
"""Check MPI tool input identity, geometry coverage and returned PETSc failures."""

import argparse
import json
import os
from pathlib import Path
import shlex
import struct
import subprocess
import sys

from hpc_inventory import digest

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture', type=Path, required=True)
    parser.add_argument('--baseline-dir', type=Path, required=True)
    parser.add_argument('--preload', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--launcher', default='mpiexec --map-by core --bind-to core')
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    source, baseline, preload = args.fixture.resolve(), args.baseline_dir.resolve(), args.preload.resolve()
    tools = ('iga_mesh_check', 'iga_assembly_smoke')
    paths = [ROOT/'solvers/cpu'/name for name in tools]+[baseline/name for name in tools]+[preload]
    payload = {r: (source/('serial.ntiga' if r == 1 else 'group.ntiga')).read_bytes() for r in (1, 2)}
    for data in payload.values():
        if struct.unpack_from('<IQQ', data, 12)[1:] != (1, 64) or struct.unpack_from('<I', data, 8)[0] != 5:
            raise ValueError('require the one-element, 64-node v5 fixture')
    # tool, name, ranks, use archived executable, expected exit, expected stage
    cases = []
    for tool in tools:
        for rank_count in (1, 2):
            cases += [(tool, f'before-{rank_count}', rank_count, True, 0, ''),
                      (tool, f'healthy-{rank_count}', rank_count, False, 0, '')]
        prefix = 'mesh check' if tool == 'iga_mesh_check' else 'assembly smoke'
        input_stage = 'mesh check database input' if tool == 'iga_mesh_check' else 'assembly smoke input'
        for damage in ('different', 'geometry', 'missing', 'fifo', 'directory', 'malformed'):
            stage = prefix+' asset database' if damage in ('different', 'geometry') else input_stage
            cases.append((tool, damage, 2, False, 1, stage+': rank 1:'))
        cases += [(tool, 'different-before', 2, True, 0, '')]
    for name, stage in [('owner-gap', 'mesh check asset database: rank 1:'),
                        ('owner-range', 'mesh check database input: rank 0:'),
                        ('owner-record', 'mesh check element input: rank 0:')]:
        cases += [('iga_mesh_check', name+'-before', 2, True, 0, ''),
                  ('iga_mesh_check', name, 2, False, 1, stage)]
    stages = dict(assembly='matrix assembly end', info='assembly smoke matrix info',
                  destroy='assembly smoke destroy', **{'mixed-insertion': 'assembly smoke element insertion'})
    for ranks in (1, 2):
        for fault, stage in stages.items():
            cases.append(('iga_assembly_smoke', f'fault-{fault}-{ranks}', ranks, False, 1, stage+': rank 0:'))
        cases.append(('iga_assembly_smoke', f'fault-destroy-before-{ranks}', ranks, True, 1,
                      'assembly smoke destroy: rank 0:'))
    summary = dict(status='running', binaries={str(p): digest(p) for p in paths}, cases=[])
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               BLIS_NUM_THREADS='1', PETSC_OPTIONS='')
    try:
        for tool, name, ranks, old, expected, stage in cases:
            tag = tool+'-'+name
            case = output/tag
            records = case/'ranks'
            records.mkdir(parents=True)
            inputs = {}
            command = ['timeout', '--kill-after=5s', '90s', *shlex.split(args.launcher)]
            fault = next((key for key in stages if name.startswith('fault-'+key+'-')), '')
            for rank in range(ranks):
                local = case/f'input-{rank}'
                local.mkdir()
                path = local/'database.ntiga'
                data = bytearray(payload[ranks])
                record = struct.unpack_from('<Q', data, 88)[0]
                if name.startswith('different') and rank == 1: data += b'\n'
                if name == 'geometry' and rank == 1:
                    end = struct.unpack_from('<Q', data, 40)[0]
                    for point in range(64):
                        position = end-64*24+point*24
                        struct.pack_into('<d', data, position, 2*struct.unpack_from('<d', data, position)[0])
                if name.startswith('owner-gap') and rank == 0:
                    struct.pack_into('<i', data, 104, 1)
                    struct.pack_into('<i', data, record+12, 1)
                if name.startswith('owner-range'): struct.pack_into('<i', data, 104, ranks)
                if name.startswith('owner-record'): struct.pack_into('<i', data, record+12, 1)
                if name == 'malformed' and rank == 1: data[:8] = b'INVALID!'
                path.write_bytes(data)
                if rank == 1 and name in ('missing', 'fifo', 'directory'):
                    path.unlink()
                    if name == 'fifo': os.mkfifo(path, 0o600)
                    if name == 'directory': path.mkdir()
                if path.is_file(): inputs[str(path)] = digest(path)
                child = [str(baseline/tool if old else ROOT/'solvers/cpu'/tool), str(path)]
                if tool == 'iga_assembly_smoke': child += ['2']
                if fault:
                    child = ['env', 'LD_PRELOAD='+str(preload), 'TUBULARFLOWIGA_TEST_SMOKE_FAULT='+fault, *child]
                if rank: command += [':']
                command += ['-np', '1', sys.executable, str(ROOT/'scripts/hpc_rank_run.py'),
                            '--output-dir', str(records), '--expected-ranks', str(ranks), '--timeout', '60', '--', *child]
            row = dict(case=tag, command_argv=command, expected_returncode=expected, expected_stage=stage,
                       input_sha256=inputs, ranks=[], status='running')
            summary['cases'].append(row)
            with (case/'launcher.log').open('x') as log:
                result = subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
            row['launcher_returncode'] = result.returncode
            if result.returncode != expected: raise RuntimeError(f'{tag}: unexpected exit {result.returncode}')
            marker = 'elements=' if tool == 'iga_mesh_check' else 'global_rows='
            for rank in range(ranks):
                directory = records/f'rank-{rank}'
                report = json.loads((directory/'run.json').read_text())
                row['ranks'].append(report)
                if report['returncode'] != expected or report['timed_out'] or report['rank'] != rank or report['ranks'] != ranks:
                    raise RuntimeError(f'{tag}: bad rank result')
                text = (directory/'stdout.log').read_text()+(directory/'stderr.log').read_text()
                lines = [line for line in text.splitlines() if line.startswith(marker)]
                false_summary = old and name.startswith('fault-destroy-')
                if stage and (stage not in text or (lines and not false_summary)):
                    raise RuntimeError(f'{tag}: missing common error or unexpected summary')
                if rank == 0 and (expected == 0 or false_summary):
                    if len(lines) != 1: raise RuntimeError(f'{tag}: missing result')
                    row['result_line'] = lines[0]
                    if old and name.startswith('owner-'):
                        if 'minimum_detJ=inf bad_elements=0' not in lines[0]:
                            raise RuntimeError(f'{tag}: empty-coverage defect was not reproduced')
                        row['observed_empty_geometry_success'] = True
                    elif expected == 0:
                        reference = next(r for r in summary['cases'] if r['case'] == tool+f'-before-{ranks}')
                        if lines[0] != reference['result_line']:
                            raise RuntimeError(f'{tag}: healthy result changed')
                if text.count('[smoke-failure-injection]') != int(bool(fault) and rank == 0):
                    raise RuntimeError(f'{tag}: incorrect fault injection count')
                if any(digest(directory/f) != h for f, h in report['logs'].items()): raise RuntimeError('log changed')
            if old and name.startswith('fault-destroy-'): row['observed_false_success_summary'] = True
            if any(digest(Path(f)) != h for f, h in inputs.items()): raise RuntimeError('input changed')
            row['status'] = 'passed'
            (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
            print(tag, 'passed', flush=True)
        if any(digest(Path(f)) != h for f, h in summary['binaries'].items()): raise RuntimeError('binary changed')
        summary['status'] = 'passed'
    except BaseException:
        summary['status'] = 'failed'
        raise
    finally:
        (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
