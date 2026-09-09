#!/usr/bin/env python3
"""Exercise resource and packed-database preflight at actual CPU/1D/graph CLI entrypoints."""
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
    parser.add_argument('--case-dir', type=Path, required=True, help='retained 2-rank staged fixture (group.ntiga and serial.ntiga)')
    parser.add_argument('--flow-case-dir', type=Path, required=True, help='retained nonzero VCA smoke fixture')
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--launcher', default='mpiexec --oversubscribe')
    parser.add_argument('--only-final', action='store_true', help='run just the backend and heterogeneous-thread cases')
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    fixture = args.case_dir.resolve()
    flow_fixture = args.flow_case_dir.resolve()
    inputs = [fixture/p for p in ('group.ntiga','serial.ntiga','simulation_config.json','controlmesh.vtk','initial_velocityfield.txt')]
    inputs += [flow_fixture/p for p in ('fixture-2.ntiga','simulation_config.json','controlmesh.vtk','initial_velocityfield.txt')]
    for p in inputs:
        if not p.is_file(): parser.error(f'missing fixture {p}')
    out = args.output_dir.resolve(); out.mkdir(parents=True,exist_ok=False)
    binaries = dict(flow=repo/'solvers/cpu/iga_navier_stokes',transport=repo/'solvers/cpu/iga_solve',
                    one_d=repo/'solvers/one_d/iga_1d',graph=repo/'solvers/coupling/iga_multidomain_flow')
    env = dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',
               PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps')
    summary = dict(status='running',binaries={k:dict(path=str(p),sha256=digest(p)) for k,p in binaries.items()},
                   inputs={str(p):digest(p) for p in inputs},cases=[])
    def run(name,kind,children,stage='',expected=1):
        case = out/name; case.mkdir()
        command = ['timeout','--kill-after=5s','90s',*shlex.split(args.launcher)]
        for rank,(settings,arguments) in enumerate(children):
            if rank: command += [':']
            command += ['-np','1','env',*[k+'='+v for k,v in settings.items()],str(binaries[kind]),*map(str,arguments)]
        record = dict(case=name,kind=kind,command_argv=command,expected_returncode=expected,stage=stage,status='running')
        summary['cases'].append(record)
        started = time.monotonic()
        log = case/'launcher.log'
        with log.open('x') as stream:
            result = subprocess.run(command,cwd=repo,env=env,stdout=stream,stderr=subprocess.STDOUT)
        text = log.read_text()
        record.update(returncode=result.returncode,elapsed_s=time.monotonic()-started,timeout_s=90,log_sha256=digest(log))
        if result.returncode != expected or (stage and stage not in text):
            raise RuntimeError(f'{name}: expected {expected} and {stage!r}')
        if expected == 1 and any(word in text for word in ('converged newton=','completed schema-','navier_stokes_v2 seconds=')):
            raise RuntimeError(f'{name}: rejected configuration reached successful solve')
        record['status']='passed'; print(name,'passed',flush=True)
        return text
    try:
        for kind in ([] if args.only_final else binaries):
            for variable,value in (('OMP_NUM_THREADS','0'),('OPENBLAS_NUM_THREADS','bad')):
                run(kind+'-'+variable,kind,[({},[]),({variable:value},[])],'execution resource preflight')
        for kind in (() if args.only_final else ('flow','transport')):
            for mode in ('partitions','missing-database'):
                children=[]
                for rank in range(2):
                    db = fixture/('serial.ntiga' if mode=='partitions' else 'group.ntiga')
                    if mode=='missing-database' and rank==1: db = fixture/'missing.ntiga'
                    children.append(({},[db,out/'missing-case']))
                run(kind+'-'+mode,kind,children,kind+' database preflight')
            folder=out/(kind+'-healthy-output'); folder.mkdir()
            arguments=[fixture/'group.ntiga',fixture,'--output',folder/'field.txt','--visualization-format','vtu']
            if kind=='transport': arguments += ['--system','transport']
            healthy_settings = {'PETSC_OPTIONS':'-ksp_type gmres -pc_type lu -pc_factor_mat_solver_type mumps'} if kind=='transport' else {}
            run(kind+'-healthy',kind,[(healthy_settings,arguments),(healthy_settings,arguments)],'execution_resources ranks=2',0)
            if not (folder/'field.txt').is_file(): raise RuntimeError('healthy native run has no field')
        folder=out/'backend-output'; folder.mkdir()
        arguments=[flow_fixture/'fixture-2.ntiga',flow_fixture,'--output',folder/'field.txt','--visualization-format','vtu']
        settings={'PETSC_OPTIONS':'-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type petsc'}
        run('flow-unavailable-parallel-lu','flow',[(settings,arguments),(settings,arguments)],'factor backend availability')
        arguments=[repo/'examples/one_d/rigid_straight','--check']
        text=run('one-d-heterogeneous-threads','one_d',[({'OMP_NUM_THREADS':'1'},arguments),({'OMP_NUM_THREADS':'2'},arguments)],expected=0)
        if 'omp_max_threads_min=1 omp_max_threads_max=2' not in text: raise RuntimeError('resource summary lost thread differences')
        if any(digest(Path(p))!=h for p,h in summary['inputs'].items()): raise RuntimeError('fixture changed')
        summary['status']='passed'
    except BaseException:
        summary['status']='failed'; raise
    finally:
        (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':
    main()
