#!/usr/bin/env python3
"""Check native VCA checkpoint output failures on retained smoke-test fixtures."""

import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import time

from hpc_inventory import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--launcher', default='mpiexec --oversubscribe')
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    fixture = args.case_dir.resolve()
    binary = repo/'solvers/cpu/iga_navier_stokes'
    required = [fixture/p for p in ('fixture.ntiga', 'fixture-2.ntiga', 'simulation_config.json', 'controlmesh.vtk')]
    for p in required:
        if not p.is_file(): parser.error(f'missing retained smoke fixture: {p}')
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -ksp_rtol 1e-12')
    summary = dict(status='running', binary_sha256=digest(binary),
                   inputs={str(p): digest(p) for p in required},
                   environment={k: env[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','PETSC_OPTIONS')}, cases=[])
    modes = [('missing-parent', 'checkpoint write file preparation'),
             ('state-directory', 'checkpoint write file preparation'),
             ('state-fifo', 'checkpoint write file preparation'),
             ('flow-metadata-directory', 'flow checkpoint metadata write'),
             ('flow-metadata-full', 'flow checkpoint metadata write'),
             ('transport-directory', 'checkpoint write file preparation'),
             ('vca-metadata-directory', 'VCA checkpoint metadata write'),
             ('vca-metadata-full', 'VCA checkpoint metadata write'),
             ('different-path', 'checkpoint write agreement')]
    try:
        for ranks in (1, 2):
            for name, stage in modes:
                if name == 'different-path' and ranks == 1: continue
                case = output/f'ranks-{ranks}-{name}'
                case.mkdir()
                prefix = case/'checkpoint'
                if name == 'missing-parent': prefix = case/'missing'/'checkpoint'
                if name == 'state-directory': Path(str(prefix)+'.state').mkdir()
                if name == 'state-fifo': os.mkfifo(str(prefix)+'.state', 0o600)
                if name == 'flow-metadata-directory': Path(str(prefix)+'.json').mkdir()
                if name == 'flow-metadata-full': Path(str(prefix)+'.json').symlink_to('/dev/full')
                if name == 'transport-directory': Path(str(prefix)+'.vca_transport.state').mkdir()
                if name == 'vca-metadata-directory': Path(str(prefix)+'.vca.json').mkdir()
                if name == 'vca-metadata-full': Path(str(prefix)+'.vca.json').symlink_to('/dev/full')
                database = fixture/('fixture.ntiga' if ranks == 1 else 'fixture-2.ntiga')
                def child(selected):
                    return [str(binary),str(database),str(fixture),'--max-newton','12',
                            '--stop-after-step','1','--checkpoint',str(selected),'--output',str(case/'flow.txt')]
                command = ['timeout','--kill-after=5s','90s',*shlex.split(args.launcher)]
                if name == 'different-path':
                    command += ['-np','1',*child(prefix),':','-np','1',*child(case/'other')]
                else:
                    command += ['-np',str(ranks),*child(prefix)]
                entry = dict(case=case.name, ranks=ranks, command_argv=command,
                             required_stage=stage, timeout_s=90, status='running')
                summary['cases'].append(entry)
                log = case/'launcher.log'
                started = time.monotonic()
                with log.open('x') as stream:
                    run = subprocess.run(command,cwd=repo,env=env,stdout=stream,stderr=subprocess.STDOUT)
                text = log.read_text()
                entry.update(returncode=run.returncode,elapsed_s=time.monotonic()-started,log_sha256=digest(log))
                if run.returncode != 1 or stage not in text:
                    raise RuntimeError(f'{case.name}: expected exit 1 and {stage!r}')
                if 'checkpoint=' in text or (case/'flow.txt').exists():
                    raise RuntimeError(f'{case.name}: checkpoint failure published success/final field')
                if list(case.rglob('*.tmp.*')):
                    raise RuntimeError(f'{case.name}: temporary checkpoint not removed')
                entry['status'] = 'passed'
                print(case.name,'passed',flush=True)
        if any(digest(Path(p)) != h for p,h in summary['inputs'].items()):
            raise RuntimeError('source fixtures changed during regression')
        summary['status'] = 'passed'
    except BaseException:
        summary['status'] = 'failed'
        raise
    finally:
        (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__ == '__main__':
    main()
