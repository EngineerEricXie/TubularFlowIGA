#!/usr/bin/env python3
"""Check native CPU flow input failures and before/after local-replica fields."""

import argparse
import json
import os
from pathlib import Path
import shlex
import struct
import subprocess
import sys
import time

from hpc_compare_fields import compare
from hpc_inventory import digest
from hpc_rank_run import rank_identity, run as run_rank

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--case-dir', type=Path)
    parser.add_argument('--reference-binary', type=Path)
    parser.add_argument('--launcher', default='mpiexec --map-by core --bind-to core')
    parser.add_argument('--rank-worker', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.rank_worker:
        status = run_rank(args.command[1:], args.output_dir, 2, 60)
        rank, _ = rank_identity(os.environ)
        (args.output_dir / f'rank-{rank}/report-ready').touch(exist_ok=False)
        deadline = time.monotonic() + 15
        while not all((args.output_dir / f'rank-{i}/report-ready').exists() for i in range(2)):
            if time.monotonic() >= deadline:
                return 2
            time.sleep(0.01)
        return status
    if args.case_dir is None or args.reference_binary is None or args.command:
        parser.error('require --case-dir and --reference-binary without positional arguments')
    fixture = args.case_dir.resolve()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    binary = ROOT / 'solvers/cpu/iga_navier_stokes'
    reference = args.reference_binary.resolve()
    binaries = {str(p): digest(p) for p in (binary, reference)}
    base = json.loads((fixture / 'simulation_config.json').read_text())
    base['simulation_scope'] = {'mode': 'flow_only'}
    base.pop('coupling', None)
    base.pop('external_circuit', None)
    base['temporal_functions'] = [dict(name='inlet_scale', kind='periodic_table',
        units='dimensionless', period=1.0, file='scale.csv', interpolation='linear')]
    base['boundaries'][0]['conditions'][0]['waveform'] = 'inlet_scale'
    base['boundaries'][0]['conditions'][0]['scale'] = 0.001
    files = {name: (fixture / name).read_bytes()
             for name in ('controlmesh.vtk', 'initial_velocityfield.txt')}
    files['scale.csv'] = b'time,value\n0,1\n0.5,1\n'
    cases = [(kind + '-' + mode, kind, ranks, '', before)
             for kind in ('flow', 'vca', 'legacy', 'legacy-config')
             for mode, ranks, before in ([('before', 1, True), ('one', 1, False)]
                 + ([('before-two', 2, True)] if kind.startswith('legacy') else [])
                 + [('two', 2, False)])]
    for name, kind, stage in (
        ('arguments', 'flow', 'flow arguments: rank 1:'),
        ('stop-control', 'flow', 'flow execution controls: rank 1:'),
        ('tolerance-control', 'flow', 'flow execution controls: rank 1:'),
        ('output-control', 'flow', 'flow execution controls: rank 1:'),
        ('database-different', 'flow', 'flow asset database: rank 1:'),
        ('configuration-different', 'flow', 'flow asset configuration: rank 1:'),
        ('configuration-absent', 'flow', 'asset catalog agreement: rank 1:'),
        ('mesh-different', 'flow', 'flow asset mesh: rank 1:'),
        ('velocity-different', 'flow', 'flow asset velocity: rank 1:'),
        ('table-different', 'flow', 'flow asset temporal inlet_scale: rank 1:'),
        ('legacy-absent', 'legacy-config', 'asset catalog agreement: rank 1:'),
        ('legacy-different', 'legacy-config', 'flow asset legacy configuration: rank 1:'),
        ('parameters-different', 'legacy', 'flow asset parameters: rank 1:'),
        ('configuration-malformed', 'flow', 'flow case input: rank 0:'),
        ('table-malformed', 'flow', 'flow boundary input: rank 0:'),
        ('parameters-malformed', 'legacy', 'flow boundary input: rank 0:')):
        cases.append((name, kind, 2, stage, False))
    for asset in ('database', 'configuration', 'mesh', 'velocity', 'table', 'parameters', 'legacy'):
        for damage in ('fifo', 'directory', 'missing'):
            if asset in ('configuration', 'legacy') and damage == 'missing':
                continue  # Optional-file selection has its own catalogue test.
            stage = 'flow database preflight: rank 1:' if asset == 'database' else 'asset content read: rank 1:'
            kind = 'legacy-config' if asset == 'legacy' else 'legacy' if asset == 'parameters' else 'flow'
            cases.append((asset + '-' + damage, kind, 2, stage, False))
    # This PETSc installation normalizes these environment options during
    # initialization. An environment difference is not an effective-option fault.
    cases.extend([('unused-table', 'flow', 2, '', False), ('petsc-environment', 'flow', 2, '', False)])
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
        BLIS_NUM_THREADS='1', IGA_ASSEMBLY_THREADS='1',
        PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -ksp_rtol 1e-12')
    summary = dict(status='running', binaries=binaries, cases=[])
    try:
        for name, kind, ranks, stage, before in cases:
            case = out / name
            records = case / 'ranks'
            records.mkdir(parents=True)
            result_dir = case / 'result'
            if not stage:
                result_dir.mkdir()
            command = ['timeout', '--kill-after=5s', '90s', *shlex.split(args.launcher)]
            inputs = {}
            for rank in range(ranks):
                local = case / f'input-{rank}'
                local.mkdir()
                config = json.loads(json.dumps(base))
                if kind == 'vca':
                    config = json.loads((fixture / 'simulation_config.json').read_text())
                if name == 'unused-table':
                    config['temporal_functions'].append(dict(name='unused', kind='periodic_table',
                        units='dimensionless', period=1.0, file='absent.csv', interpolation='linear'))
                payload = dict(files)
                payload['database.ntiga'] = (fixture / ('fixture.ntiga' if ranks == 1 else 'fixture-2.ntiga')).read_bytes()
                if kind.startswith('legacy'):
                    payload['simulation_parameter.txt'] = b"D 0.1\nvplus 1\nvminus 1\nkplus 0\nkminus 0\nk'plus 0\nk'minus 0\ndt 0.01\nnstep 2\nN0bc 0\nNplusbc 0\nNminusbc 0\n"
                    if kind == 'legacy-config':
                        payload['case_config.json'] = b'{"schema_version":1,"boundaries":{"inherit_legacy":true,"conditions":[]}}\n'
                else:
                    payload['simulation_config.json'] = (json.dumps(config) + '\n').encode()
                targets = dict(database='database.ntiga', configuration='simulation_config.json',
                    mesh='controlmesh.vtk', velocity='initial_velocityfield.txt', table='scale.csv',
                    parameters='simulation_parameter.txt', legacy='case_config.json')
                asset, _, damage = name.partition('-')
                if damage == 'malformed':
                    payload[targets[asset]] = b'invalid\n'
                if rank == 1:
                    if damage == 'different':
                        target = targets[asset]
                        if asset == 'database':
                            data = bytearray(payload[target]); struct.pack_into('<d', data, 48, 0.125)
                            payload[target] = bytes(data)
                        elif asset == 'table': payload[target] = payload[target].replace(b'0.5,1', b'0.5,2')
                        elif asset == 'velocity': payload[target] = payload[target].replace(b'1 0 0', b'2 0 0', 1)
                        elif asset == 'parameters': payload[target] = payload[target].replace(b'vplus 1', b'vplus 2')
                        else: payload[target] += b'\n'  # Valid formatting difference, exact byte contract.
                    if damage == 'absent': payload.pop(targets[asset])
                for filename, data in payload.items(): (local / filename).write_bytes(data)
                if rank == 1 and damage in ('fifo', 'directory', 'missing'):
                    target = local / targets[asset]; target.unlink()
                    if damage == 'fifo': os.mkfifo(target, 0o600)
                    if damage == 'directory': target.mkdir()
                for p in local.iterdir():
                    if p.is_file(): inputs[str(p)] = digest(p)
                child = [str(reference if before else binary), str(local/'database.ntiga'), str(local),
                    '--max-newton', '12', '--output', str(result_dir/'flow.txt'), '--visualization-format', 'vtu']
                if rank == 1:
                    if name == 'arguments': child += ['--unknown', '1']
                    if name == 'stop-control': child += ['--stop-after-step', '1']
                    if name == 'tolerance-control': child += ['--nonlinear-rtol', '0.00002']
                    if name == 'output-control':
                        at = child.index('--output'); del child[at:at+2]
                    if name == 'petsc-environment': child = ['env', 'PETSC_OPTIONS=' + env['PETSC_OPTIONS'] + ' -unused_difference 1', *child]
                if rank: command.append(':')
                worker = ([sys.executable, str(Path(__file__).resolve()), '--rank-worker', '--output-dir', str(records)]
                    if ranks == 2 else [sys.executable, str(ROOT/'scripts/hpc_rank_run.py'), '--output-dir', str(records), '--timeout', '60'])
                command += ['-np', '1', *worker, '--', *child]
            entry = dict(case=name, command_argv=command, input_sha256=inputs,
                expected_returncode=1 if stage else 0, expected_stage=stage, status='running', ranks=[])
            summary['cases'].append(entry)
            with (case/'launcher.log').open('x') as log:
                result = subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
            entry['launcher_returncode'] = result.returncode
            if result.returncode != entry['expected_returncode']: raise RuntimeError(f'{name}: unexpected exit')
            for rank in range(ranks):
                d = records/f'rank-{rank}'; report = json.loads((d/'run.json').read_text())
                if report['returncode'] != entry['expected_returncode'] or report['timed_out']: raise RuntimeError('bad rank exit')
                if report['rank'] != rank or report['ranks'] != ranks: raise RuntimeError('bad rank identity')
                for filename, checksum in report['logs'].items():
                    if digest(d/filename) != checksum: raise RuntimeError('log hash changed')
                entry['ranks'].append(report)
            text = (records/'rank-0/stdout.log').read_text() + (records/'rank-0/stderr.log').read_text()
            if stage:
                if stage not in text or result_dir.exists(): raise RuntimeError(f'{name}: wrong stage or output published')
            else:
                if 'navier_stokes_v2 seconds=' not in text: raise RuntimeError('missing successful solve summary')
                # Legacy pressure is roundoff around zero in this fixture. Check
                # its before/after behavior at the same rank count; cross-rank
                # comparisons below are retained as observations, never relabeled.
                baseline_name = kind + ('-before-two' if kind.startswith('legacy') and ranks == 2 else '-before')
                entry['fields'] = [compare(out/baseline_name/'result'/f, result_dir/f, 1e-6, 1e-12)
                                   for f in ('flow.txt', 'flow.txt.pressure')]
                if not all(c['passed'] for c in entry['fields']): raise RuntimeError('field regression')
                if kind.startswith('legacy') and ranks == 2:
                    entry['cross_rank_observations'] = [compare(out/(kind+'-before/result')/f, result_dir/f, 1e-6, 1e-12)
                        for f in ('flow.txt', 'flow.txt.pressure')]
                    entry['comparison_scope'] = 'Same-rank before/after compatibility only; see failed cross-rank pressure observation.'
            if any(digest(Path(p)) != h for p, h in inputs.items()): raise RuntimeError('inputs changed during run')
            entry['status'] = 'passed'
            print(name, 'passed', flush=True)
        if any(digest(Path(p)) != h for p, h in binaries.items()): raise RuntimeError('binary changed during run')
        summary['status'] = 'passed'
    except BaseException:
        summary['status'] = 'failed'
        raise
    finally:
        (out/'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
