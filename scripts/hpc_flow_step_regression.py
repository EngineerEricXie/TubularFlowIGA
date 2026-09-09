#!/usr/bin/env python3
"""Exercise native VCA step failures and accepted-step output/checkpoint boundaries."""

import argparse
import json
import os
from pathlib import Path
import shlex
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
    source, output = args.case_dir.resolve(), args.output_dir.resolve()
    reference = args.reference_binary.resolve()
    binary = ROOT/'solvers/cpu/iga_navier_stokes'
    output.mkdir(parents=True, exist_ok=False)
    base = json.loads((source/'simulation_config.json').read_text())
    cases = [('healthy-before', 1, 'healthy', reference)]
    cases += [(f'{kind}-{ranks}', ranks, kind, binary)
              for kind in ('healthy', 'no-output', 'split', 'first-failure', 'second-failure') for ranks in (1, 2)]
    cases += [(kind, 2, kind, binary) for kind in ('history-directory', 'history-fifo')]
    summary = dict(status='running', cases=[], binaries={str(p): digest(p) for p in (reference, binary)})
    env = dict(os.environ, OMP_NUM_THREADS='1', IGA_ASSEMBLY_THREADS='1', OPENBLAS_NUM_THREADS='1',
               MKL_NUM_THREADS='1', BLIS_NUM_THREADS='1',
               PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -ksp_rtol 1e-12')
    try:
        for name, ranks, kind, executable in cases:
            case = output/name
            local, result, records = case/'input', case/'result', case/'ranks'
            local.mkdir(parents=True); result.mkdir(); records.mkdir()
            config = json.loads(json.dumps(base))
            if kind in ('first-failure', 'second-failure', 'split'):
                # Negative infusion is accepted by the existing schema. The
                # unmodified reservoir equation must reject the depleted mass.
                config['external_circuit']['infusion_rates'] = {'oxygen': -30 if kind == 'first-failure' else -15}
            (local/'simulation_config.json').write_text(json.dumps(config)+'\n')
            for filename in ('controlmesh.vtk', 'initial_velocityfield.txt'):
                (local/filename).write_bytes((source/filename).read_bytes())
            (local/'database.ntiga').write_bytes((source/('fixture.ntiga' if ranks == 1 else 'fixture-2.ntiga')).read_bytes())
            inputs = {str(p): digest(p) for p in local.iterdir()}
            if kind == 'history-directory': (result/'coupling_manifest.json').mkdir()
            if kind == 'history-fifo': os.mkfifo(result/'coupling_manifest.json', 0o600)
            child = [str(executable), str(local/'database.ntiga'), str(local), '--max-newton', '12',
                     '--visualization-format', 'vtu', '--checkpoint', str(result/'checkpoint'), '--checkpoint-every', '1']
            if kind == 'split': child += ['--stop-after-step', '1']
            if kind != 'no-output': child += ['--output', str(result/'flow.txt'), '--output-every', '1']
            worker = ([sys.executable, str(ROOT/'scripts/hpc_flow_input_regression.py'), '--rank-worker', '--output-dir', str(records)]
                      if ranks == 2 else [sys.executable, str(ROOT/'scripts/hpc_rank_run.py'), '--output-dir', str(records), '--timeout', '60'])
            command = ['timeout', '--kill-after=5s', '90s', *shlex.split(args.launcher), '-np', str(ranks), *worker, '--', *child]
            stage = ('flow VCA circuit advance' if kind.endswith('failure') else
                     'flow VCA history output' if kind.startswith('history-') else '')
            row = dict(case=name, command_argv=command, input_sha256=inputs, expected_stage=stage,
                       expected_returncode=1 if stage else 0, ranks=[], status='running', fields=[])
            summary['cases'].append(row)
            with (case/'launcher.log').open('x') as log:
                process = subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
            row['launcher_returncode'] = process.returncode
            if process.returncode != row['expected_returncode']: raise RuntimeError(f'{name}: unexpected launcher exit')
            for rank in range(ranks):
                rd = records/f'rank-{rank}'
                report = json.loads((rd/'run.json').read_text()); row['ranks'].append(report)
                if report['returncode'] != row['expected_returncode'] or report['timed_out']: raise RuntimeError('unexpected rank exit')
                if report['rank'] != rank or report['ranks'] != ranks: raise RuntimeError('wrong rank identity')
                for filename, checksum in report['logs'].items():
                    if digest(rd/filename) != checksum: raise RuntimeError('rank log changed')
                text = (rd/'stdout.log').read_text()+(rd/'stderr.log').read_text()
                if stage:
                    if stage+': rank 0:' not in text or 'navier_stokes_v2 seconds=' in text:
                        raise RuntimeError('wrong failure stage or false successful summary')
                elif rank == 0 and 'navier_stokes_v2 seconds=' not in text: raise RuntimeError('missing successful summary')
            completed = 0 if kind == 'first-failure' else 1 if kind in ('second-failure', 'split') else 2
            row['last_accepted_step'] = completed
            if completed:
                for filename in ('checkpoint.json', 'checkpoint.vca.json'):
                    metadata = json.loads((result/filename).read_text())
                    if metadata['completed_step'] != completed: raise RuntimeError('wrong accepted checkpoint step')
            elif (result/'checkpoint.json').exists() or (result/'checkpoint.vca.json').exists():
                raise RuntimeError('unaccepted checkpoint was published')
            if kind.endswith('failure'):
                if not (result/f'flow.step{completed:06d}.vtu').is_file(): raise RuntimeError('last accepted output missing')
                if (result/f'flow.step{completed+1:06d}.txt').exists() or (result/'flow.txt').exists():
                    raise RuntimeError('output followed a rejected circuit step')
                if (result/'coupling_manifest.json').exists(): raise RuntimeError('failed circuit wrote final history')
                if kind == 'second-failure':
                    stopped = output/f'split-{ranks}'/'result'
                    metadata = json.loads((result/'checkpoint.json').read_text())
                    vca = json.loads((result/'checkpoint.vca.json').read_text())
                    files = ['checkpoint.json', 'checkpoint.vca.json', metadata['state_file'], vca['transport_state_file']]
                    row['accepted_checkpoint_sha256'] = {}
                    for filename in files:
                        if (stopped/filename).read_bytes() != (result/filename).read_bytes():
                            raise RuntimeError('last checkpoint differs from an intentional stop at the accepted step')
                        row['accepted_checkpoint_sha256'][filename] = digest(result/filename)
            elif not stage:
                history_path = local/'results/vca_flow/coupling_manifest.json' if kind == 'no-output' else result/'coupling_manifest.json'
                history = json.loads(history_path.read_text())
                if len(history['arterial_inlet_history']) != completed: raise RuntimeError('wrong inlet history length')
                row['history_sha256'] = digest(history_path)
                if kind not in ('no-output', 'split'):
                    for field in ('flow.txt', 'flow.txt.pressure'):
                        comparison = compare(output/'healthy-before/result'/field, result/field, 1e-6, 1e-12)
                        row['fields'].append(comparison)
                        if not comparison['passed']: raise RuntimeError('healthy VCA field changed')
            if any(digest(Path(p)) != checksum for p, checksum in inputs.items()): raise RuntimeError('inputs changed')
            row['status'] = 'passed'
            (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
            print(name, 'passed', flush=True)
        if any(digest(Path(p)) != checksum for p, checksum in summary['binaries'].items()): raise RuntimeError('binary changed')
        summary['status'] = 'passed'
    finally:
        (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
