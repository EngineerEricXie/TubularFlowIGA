#!/usr/bin/env python3
"""Exercise configured CPU transport input, stepping, output and restart contracts."""

import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import struct
import subprocess
import sys

from hpc_compare_fields import compare
from hpc_inventory import digest

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case-dir', type=Path, required=True)
    parser.add_argument('--reference-binary', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--launcher', default='mpiexec --map-by core --bind-to core')
    args = parser.parse_args()
    source = args.case_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    binary = ROOT/'solvers/cpu/iga_solve'
    reference = args.reference_binary.resolve()
    binaries = {str(p): digest(p) for p in (binary, reference)}
    base = json.loads((source/'simulation_config.json').read_text())
    base['temporal_functions'] = [dict(name='red_scale', kind='periodic_table', units='dimensionless',
        period=1, file='red.csv', interpolation='linear')]
    inlet = next(b for b in base['boundaries'] if b['label'] == 1)
    next(c for c in inlet['conditions'] if c['field'] == 'three_red')['waveform'] = 'red_scale'
    files = {p: (source/p).read_bytes() for p in ('controlmesh.vtk', 'initial_velocityfield.txt')}
    files['red.csv'] = b'time,value\n0,1\n0.04,1.2\n'
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', BLIS_NUM_THREADS='1',
        PETSC_OPTIONS='-ksp_type gmres -pc_type lu -pc_factor_mat_solver_type mumps -ksp_rtol 1e-12')
    # name, fixture kind, MPI size, diagnostic, use archived binary
    cases = [(kind+'-'+mode, kind, ranks, '', old)
        for kind in ('prescribed', 'series', 'vtkhdf')
        for mode, ranks, old in (('before', 1, True), ('one', 1, False), ('two', 2, False))]
    cases += [(name, 'prescribed', 2, '', False) for name in ('override', 'memory', 'unused-table')]
    cases.append(('unused-snapshot', 'series', 2, '', False))
    for name, kind, stage in (
        ('arguments', 'prescribed', 'transport arguments: rank 1:'),
        ('system-control', 'prescribed', 'transport execution controls: rank 1:'),
        ('stop-control', 'prescribed', 'transport execution controls: rank 1:'),
        ('memory-control', 'prescribed', 'transport execution controls: rank 1:'),
        ('output-control', 'prescribed', 'transport execution controls: rank 1:'),
        ('bad-system', 'prescribed', 'transport case input: rank 0:'),
        ('series-override', 'series', 'transport case input: rank 0:'),
        ('configuration-malformed', 'prescribed', 'transport case input: rank 0:'),
        ('mesh-malformed', 'prescribed', 'transport case input: rank 0:'),
        ('velocity-malformed', 'prescribed', 'transport forcing input: rank 0:'),
        ('table-malformed', 'prescribed', 'transport forcing input: rank 0:'),
        ('manifest-malformed', 'series', 'transport forcing input: rank 0:'),
        ('snapshot-malformed', 'series', 'transport velocity input: rank 0:'),
        ('snapshot-late', 'series', 'transport asset lower snapshot: rank 1:'),
        ('series-range', 'series', 'transport velocity selection: rank 0:'),
        ('step-overflow', 'prescribed', 'transport step input: rank 0:'),
        ('solver-preonly', 'prescribed', 'transport linear solve: rank 0:'),
        ('solver-unavailable', 'prescribed', 'factor backend availability: rank 0:'),
        ('memory-directory', 'prescribed', 'transport memory report setup: rank 0:'),
        ('memory-fifo', 'prescribed', 'transport memory report setup: rank 0:'),
        ('output-initial', 'prescribed', 'transport field output: rank 0:'),
        ('output-step', 'prescribed', 'transport field output: rank 0:'),
        ('output-final', 'prescribed', 'transport field output: rank 0:'),
        ('output-fifo', 'prescribed', 'transport field output: rank 0:'),
        ('output-index', 'prescribed', 'transport output index: rank 0:'),
        ('output-physiology', 'prescribed', 'transport output index: rank 0:'),
        ('checkpoint-metadata-output', 'prescribed', 'transport checkpoint metadata write: rank 0:')):
        cases.append((name, kind, 2, stage, False))
    for asset in ('database', 'configuration', 'mesh', 'velocity', 'table', 'manifest', 'snapshot'):
        for damage in ('different', 'missing', 'fifo', 'directory'):
            kind = 'series' if asset in ('manifest', 'snapshot') else 'prescribed'
            if asset == 'database': stage = 'transport asset database' if damage == 'different' else 'transport database preflight'
            elif damage != 'different': stage = 'asset content read'
            else: stage = 'transport asset ' + dict(configuration='configuration', mesh='mesh', velocity='velocity',
                table='temporal red_scale', manifest='velocity manifest', snapshot='upper snapshot')[asset]
            cases.append((asset+'-'+damage, kind, 2, stage+': rank 1:', False))
    for ranks in (1, 2):
        cases += [(f'split-{ranks}', 'prescribed', ranks, '', False), (f'resume-{ranks}', 'prescribed', ranks, '', False)]
    cases.append(('resume-cross-rank', 'prescribed', 2, '', False))
    for name, stage in (
        ('metadata-missing', 'transport checkpoint metadata: rank 1:'),
        ('metadata-fifo', 'transport checkpoint metadata: rank 1:'),
        ('metadata-malformed', 'transport checkpoint metadata: rank 1:'),
        ('metadata-different', 'transport checkpoint metadata agreement: rank 1:'),
        ('state-missing', 'asset content read: rank 1:'),
        ('state-fifo', 'asset content read: rank 1:'),
        ('state-different', 'checkpoint asset state: rank 1:'),
        ('state-truncated', 'checkpoint local read: rank 0:')):
        cases.append(('restart-'+name, 'prescribed', 2, stage, False))
    summary = dict(status='running', binaries=binaries, cases=[],
        notes=['Small correctness fixtures; no scaling claim. Baseline observations are not independent comparisons.',
               'Snapshot inputs are validated when selected; unused files are deliberately not opened.',
               'Checkpoint checks retain the existing format, not a complete crash-safe publication protocol.'])
    targets = dict(database='database.ntiga', configuration='simulation_config.json', mesh='controlmesh.vtk',
        velocity='initial_velocityfield.txt', table='red.csv', manifest='velocity.csv', snapshot='v2.txt')
    try:
        for name, kind, ranks, stage, old in cases:
            case = output/name
            records, results = case/'ranks', case/'result'
            records.mkdir(parents=True); results.mkdir()
            if name.startswith('output-') and name != 'output-control':
                target = results / dict(initial='field.step000000.txt', step='field.step000001.txt',
                    final='field.txt', fifo='field.step000000.txt', index='field.pvd', physiology='physiology_fields.json')[name[7:]]
                if name == 'output-fifo': os.mkfifo(target, 0o600)
                else: target.mkdir()
            if name.startswith('memory-') and name != 'memory-control':
                if name == 'memory-fifo': os.mkfifo(results/'memory.jsonl', 0o600)
                else: (results/'memory.jsonl').mkdir()
            if name == 'checkpoint-metadata-output': (results/'checkpoint.json').mkdir()
            command = ['timeout', '--kill-after=5s', '90s', *shlex.split(args.launcher)]
            inputs = {}
            for rank in range(ranks):
                local = case/f'input-{rank}'; local.mkdir()
                config = json.loads(json.dumps(base)); payload = dict(files)
                payload['database.ntiga'] = (source/('serial.ntiga' if ranks == 1 else 'group.ntiga')).read_bytes()
                if kind == 'series':
                    config['velocity_sources'] = [dict(name='series', kind='snapshot_series', manifest='velocity.csv',
                        interpolation='linear', out_of_range='error')]
                    for system in config['equation_systems']:
                        for term in system.get('terms', []):
                            if term.get('velocity') == 'prescribed': term['velocity'] = 'series'
                    payload['velocity.csv'] = b'time,file\n0,v0.txt\n0.02,v2.txt\n'
                    payload['v0.txt'] = files['initial_velocityfield.txt']
                    payload['v1.txt'] = files['initial_velocityfield.txt']
                    # Exercise interpolation between unequal velocity fields.
                    payload['v2.txt'] = ''.join(' '.join(f'{1.5*float(v):.17g}' for v in line.split())+'\n'
                        for line in files['initial_velocityfield.txt'].decode().splitlines()).encode()
                    assert payload['v2.txt'] != payload['v0.txt']
                    if name == 'unused-snapshot': payload['velocity.csv'] += b'0.03,absent.txt\n'
                    if name == 'snapshot-late': payload['velocity.csv'] = b'time,file\n0,v0.txt\n0.01,v1.txt\n0.02,v2.txt\n'
                    if name == 'series-range': payload['velocity.csv'] = b'time,file\n0,v0.txt\n0.01,v1.txt\n'
                if name == 'unused-table':
                    config['temporal_functions'].append(dict(name='unused', kind='periodic_table', units='dimensionless',
                        period=1, file='absent.csv', interpolation='linear'))
                payload['simulation_config.json'] = (json.dumps(config)+'\n').encode()
                asset, _, damage = name.partition('-')
                if damage == 'malformed' and asset in targets: payload[targets[asset]] = b'invalid\n'
                if name == 'step-overflow': payload['red.csv'] = b'time,value\n0,1\n0.01,1.7e308\n'
                if rank == 1 and (damage == 'different' or name == 'snapshot-late'):
                    target = targets[asset]; previous = payload[target]
                    if asset == 'database':
                        data = bytearray(previous); struct.pack_into('<d', data, 48, 0.125); payload[target] = bytes(data)
                    else: payload[target] += b'\n'  # Exact input identity also rejects valid formatting differences.
                    assert payload[target] != previous
                for filename, data in payload.items(): (local/filename).write_bytes(data)
                if rank == 1 and asset in targets and damage in ('missing', 'fifo', 'directory'):
                    target = local/targets[asset]; target.unlink()
                    if damage == 'fifo': os.mkfifo(target, 0o600)
                    if damage == 'directory': target.mkdir()
                if name == 'override': shutil.copyfile(local/'initial_velocityfield.txt', local/'override.txt')
                restart = None
                if name.startswith(('resume-', 'restart-')):
                    from_ranks = 1 if name == 'resume-cross-rank' else ranks
                    restart = local/'restart'
                    for suffix in ('.json', '.state'):
                        shutil.copyfile(output/f'split-{from_ranks}/result'/('checkpoint'+suffix), Path(str(restart)+suffix))
                    # Preserve the relative metadata filename in each local replica.
                    shutil.copyfile(Path(str(restart)+'.state'), local/'checkpoint.state')
                    if name.startswith('restart-') and (rank == 1 or name == 'restart-state-truncated'):
                        mode = name[len('restart-'):]
                        target = Path(str(restart)+'.json') if mode.startswith('metadata') else local/'checkpoint.state'
                        operation = mode.split('-')[-1]
                        if operation in ('missing', 'fifo'):
                            target.unlink()
                            if operation == 'fifo': os.mkfifo(target, 0o600)
                        elif operation == 'malformed': target.write_text('{invalid\n')
                        elif operation == 'truncated': target.write_bytes(target.read_bytes()[:-1])
                        elif operation == 'different' and mode.startswith('metadata'):
                            meta = json.loads(target.read_text()); meta['completed_step'] = 0; meta['physical_time'] = 0
                            target.write_text(json.dumps(meta)+'\n')
                        else:
                            data = bytearray(target.read_bytes()); data[-1] ^= 1; target.write_bytes(data)
                for p in local.iterdir():
                    if p.is_file(): inputs[str(p)] = digest(p)
                child = [str(reference if old else binary), str(local/'database.ntiga'), str(local), '--system', 'transport',
                    '--output', str(results/'field.txt'), '--output-every', '1', '--checkpoint', str(results/'checkpoint'),
                    '--visualization-format', 'vtkhdf' if kind == 'vtkhdf' else 'vtu']
                if name.startswith('split-'): child += ['--stop-after-step', '1']
                if restart: child += ['--restart', str(restart)]
                if name in ('memory', 'memory-directory', 'memory-fifo') or (name == 'memory-control' and rank == 1):
                    child += ['--memory-report', str(results/'memory.jsonl')]
                if name in ('override', 'series-override'):
                    child += ['--velocity', str(local/('override.txt' if name == 'override' else 'initial_velocityfield.txt'))]
                if name == 'bad-system' or (name == 'system-control' and rank == 1): child[child.index('--system')+1] = 'missing'
                if name == 'arguments' and rank == 1: child += ['--unknown', '1']
                if name == 'stop-control' and rank == 1: child += ['--stop-after-step', '1']
                if name == 'output-control' and rank == 1: child[child.index('--output-every')+1] = '2'
                if name == 'solver-preonly': child = ['env', 'PETSC_OPTIONS=-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps', *child]
                if name == 'solver-unavailable': child = ['env', 'PETSC_OPTIONS=-ksp_type gmres -pc_type lu -pc_factor_mat_solver_type petsc', *child]
                worker = ([sys.executable, str(ROOT/'scripts/hpc_flow_input_regression.py'), '--rank-worker', '--output-dir', str(records)]
                    if ranks == 2 else [sys.executable, str(ROOT/'scripts/hpc_rank_run.py'), '--output-dir', str(records), '--timeout', '60'])
                if rank: command.append(':')
                command += ['-np', '1', *worker, '--', *child]
            entry = dict(case=name, command_argv=command, input_sha256=inputs, expected_stage=stage,
                expected_returncode=1 if stage else 0, status='running', ranks=[])
            summary['cases'].append(entry)
            with (case/'launcher.log').open('x') as log:
                run = subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
            entry['launcher_returncode'] = run.returncode
            if run.returncode != entry['expected_returncode']: raise RuntimeError(f'{name}: unexpected launcher exit {run.returncode}')
            for rank in range(ranks):
                d=records/f'rank-{rank}'; r=json.loads((d/'run.json').read_text())
                if r['rank'] != rank or r['ranks'] != ranks or r['returncode'] != entry['expected_returncode'] or r['timed_out']:
                    raise RuntimeError(f'{name}: rank status differs')
                for p,h in r['logs'].items():
                    if digest(d/p) != h: raise RuntimeError('rank log hash differs')
                if stage and stage not in (d/'stderr.log').read_text(): raise RuntimeError(f'{name}: wrong rank failure stage')
                if stage and 'iga_solve system=' in (d/'stdout.log').read_text(): raise RuntimeError('failed job printed successful summary')
                entry['ranks'].append(r)
            if not stage:
                baseline = output/(kind+'-before')/'result'
                step = 1 if name.startswith('split-') else 2
                expected_file = baseline/(f'field.step{step:06d}.txt')
                entry['fields'] = [compare(expected_file, results/'field.txt', 1e-6, 1e-12, node_ids=True)]
                for p in sorted(results.glob('field.step*.txt')):
                    entry['fields'].append(compare(baseline/p.name, p, 1e-6, 1e-12, node_ids=True))
                if not all(f['passed'] for f in entry['fields']): raise RuntimeError(f'{name}: field comparison failed')
                if (results/'field.txt.fields').read_text() != 'three_red\nthree_blue\n': raise RuntimeError('field order differs')
                meta=json.loads((results/'checkpoint.json').read_text())
                if meta['completed_step'] != step or not (results/'checkpoint.state').is_file(): raise RuntimeError('incomplete checkpoint')
                if name == 'memory':
                    data=[json.loads(line) for line in (results/'memory.jsonl').read_text().splitlines()]
                    if [r['stage'] for r in data] != ['database_open','required_elements_loaded','matrix_preallocation','operator_assembly','ksp_setup','linear_solve']:
                        raise RuntimeError('memory stages missing or reordered')
                    if any(r['mpi_ranks'] != ranks or len(r['rss_peak_bytes_per_rank']) != ranks for r in data): raise RuntimeError('memory rank data incomplete')
            if name in ('snapshot-late', 'series-range'):
                if not (results/'field.step000001.txt').is_file() or (results/'field.step000002.txt').exists(): raise RuntimeError('late input failure lost step boundary')
            if name == 'step-overflow':
                if not (results/'field.step000000.txt').is_file() or (results/'field.step000001.txt').exists(): raise RuntimeError('overflow crossed accepted step boundary')
            if any(digest(Path(p)) != h for p,h in inputs.items()): raise RuntimeError('case input changed during run')
            entry['status']='passed'; print(name, 'passed', flush=True)
        if any(digest(Path(p)) != h for p,h in binaries.items()): raise RuntimeError('binary changed during regression')
        summary['status']='passed'
    except BaseException:
        summary['status']='failed'; raise
    finally:
        (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__ == '__main__':
    main()
