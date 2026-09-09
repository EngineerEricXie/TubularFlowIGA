#!/usr/bin/env python3
"""Probe real CUDA CLI startup rejection before GPU allocation or file output."""
import argparse
import json
import os
from pathlib import Path
import struct
import subprocess
from hpc_inventory import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--fixture-root', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    binary = repo/'solvers/cuda/iga_cuda'
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    fixtures = args.fixture_root.resolve()
    case = fixtures/'hpc00/matrix/transport-inputs'
    source = case/'straight_neurite-1.ntiga'
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    for name in ('OMPI_COMM_WORLD_RANK', 'OMPI_COMM_WORLD_SIZE', 'PMI_RANK', 'PMI_SIZE',
                 'SLURM_PROCID', 'SLURM_NTASKS', 'SLURM_STEP_NUM_TASKS', 'SLURM_STEP_ID', 'SLURM_STEPID'):
        env.pop(name, None)
    summary = dict(status='running', cases=[], binaries={str(binary):digest(binary),
        str(args.baseline.resolve()):digest(args.baseline)}, input_sha256=digest(source))

    def run(name, argv, expected, settings=None, ranks=None, old=False, diagnostic=None, full=False):
        directory = root/name
        directory.mkdir()
        executable = args.baseline.resolve() if old else binary
        command = [str(executable), *map(str, argv)]
        if ranks:
            (directory/'ranks').mkdir()
            command = ['mpiexec', '--oversubscribe', '-np', str(ranks), 'python3',
                str(repo/'scripts/hpc_sequential_output_regression.py'), '--rank-worker',
                '--output-dir', str(directory/'ranks'), '--', *command]
        command = ['timeout', '--kill-after=5s', '45s', *command]
        record = dict(name=name, argv=command, expected=expected, environment_overrides=settings or {})
        summary['cases'].append(record)
        stdout = Path('/dev/full') if full else directory/'stdout.log'
        with stdout.open('wb') as out, (directory/'stderr.log').open('wb') as err:
            result = subprocess.run(command, env=dict(env, **(settings or {})), stdout=out, stderr=err)
        record['returncode'] = result.returncode
        if result.returncode != expected: raise RuntimeError(name+': wrong exit')
        text = (directory/'stderr.log').read_text()
        rank_output = ''
        if ranks:
            reports = []
            for rank in range(ranks):
                rank_dir = directory/'ranks'/f'rank-{rank}'
                report = json.loads((rank_dir/'run.json').read_text())
                if report['returncode'] != expected or report['timed_out']:
                    raise RuntimeError(name+': wrong per-process exit')
                reports.append(report)
                text += (rank_dir/'stderr.log').read_text()
                rank_output += (rank_dir/'stdout.log').read_text()
            record['rank_reports'] = reports
        if diagnostic and diagnostic not in text: raise RuntimeError(name+': missing diagnostic')
        if expected == 1:
            if 'requested_peak_bytes=0 requested_live_bytes=0' not in text:
                raise RuntimeError(name+': allocated project GPU buffers before rejection')
            if (root/'forbidden-output').exists(): raise RuntimeError(name+': unexpected output')
            if ranks and text.count('iga_cuda: '+diagnostic) != ranks:
                raise RuntimeError(name+': not every launched process rejected')
        elif not old:
            output = rank_output if ranks else stdout.read_text()
            if 'processes=1 gpu_limit=1' not in output or 'cuda_capabilities runtime_build=' not in output:
                raise RuntimeError(name+': missing capabilities')
        record['status'] = 'passed'
        print(name, 'passed', flush=True)

    try:
        run('before-duplicate-device-info', ['device-info'], 0, ranks=2, old=True)
        for command in ('device-info', 'mesh-check', 'solve', 'transport', 'navier-stokes'):
            argv = [command] if command == 'device-info' else [command, 'missing.ntiga', 'missing-case']
            run('reject-duplicate-'+command, argv, 1, ranks=2,
                diagnostic='single-process executable requires one task; OMPI_COMM_WORLD_SIZE=2')
        for variable in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS'):
            run('reject-'+variable, ['device-info'], 1, settings={variable:'bad', 'CUDA_VISIBLE_DEVICES':''}, diagnostic=variable)
        run('reject-pmi', ['device-info'], 1, settings={'PMI_RANK':'0', 'PMI_SIZE':'2', 'CUDA_VISIBLE_DEVICES':''}, diagnostic='PMI_SIZE=2')
        run('reject-slurm-step', ['device-info'], 1, settings={'SLURM_PROCID':'0', 'SLURM_STEP_ID':'0',
            'SLURM_STEP_NUM_TASKS':'2', 'SLURM_NTASKS':'8', 'CUDA_VISIBLE_DEVICES':''}, diagnostic='SLURM_STEP_NUM_TASKS=2')
        run('reject-resource-stdout', ['device-info'], 1, settings={'CUDA_VISIBLE_DEVICES':''},
            diagnostic='cannot write complete text stream', full=True)
        run('direct', ['device-info'], 0)
        run('mpi-single', ['device-info'], 0, ranks=1)
        run('allocation-shell', ['device-info'], 0, settings={'SLURM_NTASKS':'8'})
        run('batch-shell', ['device-info'], 0, settings={'SLURM_PROCID':'0', 'SLURM_NTASKS':'8', 'SLURM_STEP_ID':'batch'})
        run('single-step', ['device-info'], 0, settings={'SLURM_PROCID':'0', 'SLURM_NTASKS':'8', 'SLURM_STEP_NUM_TASKS':'1'})
        for command, nodes in (('mesh-check', 2**31), ('transport', (2**31-1)//3+1),
                               ('solve', (2**31-1)//3+1), ('navier-stokes', (2**31-1)//4+1)):
            data = bytearray(source.read_bytes())
            struct.pack_into('<Q', data, 24, nodes)
            database = root/(command+'-oversized.ntiga')
            database.write_bytes(data)
            argv = [command, database]
            if command != 'mesh-check': argv += [case]
            if command in ('solve', 'navier-stokes'): argv += ['--output', root/'forbidden-output']
            run('reject-capacity-'+command, argv, 1, diagnostic='int32 capacity')
        if digest(source) != summary['input_sha256']: raise RuntimeError('input changed')
        summary['status'] = 'passed'
    except BaseException:
        summary['status'] = 'failed'
        raise
    finally:
        (root/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
