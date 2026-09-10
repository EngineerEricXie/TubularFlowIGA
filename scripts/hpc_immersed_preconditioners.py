#!/usr/bin/env python3
"""Evaluate immersed static candidates with existing field, gauge and rollback gates."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

from hpc_inventory import digest


def main():
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--binary', type=Path, default=repo/'solvers/cpu/immersed_distributed_static_flow_test')
    parser.add_argument('--modes', nargs='+', choices=['closed','flow','pressure'], default=['closed','flow','pressure'])
    parser.add_argument('--ranks', type=int, default=2)
    parser.add_argument('--timeout', type=int, default=600)
    args = parser.parse_args()
    if args.ranks<1 or args.timeout<1:
        parser.error('ranks and timeout must be positive')
    root, binary = args.output_dir.resolve(), args.binary.resolve()
    root.mkdir(parents=True, exist_ok=False)
    common = '-test_solver_prefix -domain_test_flow_ksp_type fgmres -domain_test_flow_ksp_rtol 1e-12 -domain_test_flow_ksp_converged_reason '
    candidates = {
        'lu': '-domain_test_flow_pc_type lu -domain_test_flow_pc_factor_mat_solver_type mumps',
        'block-ilu': '-domain_test_flow_pc_type bjacobi -domain_test_flow_sub_ksp_type preonly -domain_test_flow_sub_pc_type ilu',
        'block-lu-shift': '-domain_test_flow_pc_type bjacobi -domain_test_flow_sub_ksp_type preonly -domain_test_flow_sub_pc_type lu -domain_test_flow_sub_pc_factor_mat_solver_type mumps -domain_test_flow_sub_pc_factor_shift_type nonzero -domain_test_flow_sub_pc_factor_shift_amount 1e-12',
    }
    report = dict(status='running', binary_sha256=digest(binary), harness_sha256=digest(Path(__file__)), jobs=[])
    try:
        for mode in dict.fromkeys(args.modes):
            for name, options in candidates.items():
                folder = root/(mode+'-'+name)
                folder.mkdir()
                env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', PETSC_OPTIONS=common+options)
                command = ['timeout','--kill-after=5s',str(args.timeout+15),'mpiexec','-np',str(args.ranks),
                           sys.executable,str(repo/'scripts/hpc_rank_run.py'),'--output-dir',str(folder),
                           '--expected-ranks',str(args.ranks),'--timeout',str(args.timeout),'--',str(binary),mode]
                with (folder/'launcher.log').open('w') as log:
                    result = subprocess.run(command,cwd=repo,env=env,stdout=log,stderr=subprocess.STDOUT)
                job = dict(mode=mode,candidate=name,command=command,options=env['PETSC_OPTIONS'],
                           returncode=result.returncode,status='failed',ranks=[],diagnostics=[])
                report['jobs'].append(job)
                for rank in range(args.ranks):
                    directory = folder/f'rank-{rank}'
                    measured = directory/'run.json'
                    if measured.exists():
                        job['ranks'].append(json.loads(measured.read_text()))
                    for filename in ['stdout.log','stderr.log']:
                        path = directory/filename
                        if path.exists():
                            job['diagnostics'].append(dict(rank=rank,file=filename,sha256=digest(path),tail=path.read_text()[-8192:]))
                if result.returncode==0:
                    if len(job['ranks'])!=args.ranks or any(r['returncode'] or r['timed_out'] or not r['resource'] for r in job['ranks']):
                        raise RuntimeError('successful process lacks valid rank measurements')
                    records=[]
                    for rank in range(args.ranks):
                        log=(folder/f'rank-{rank}/stdout.log').read_text()
                        layouts=[dict(v.split('=',1) for v in line.split()[1:]) for line in log.splitlines() if line.startswith('immersed_candidate_layout ')]
                        if len(layouts)!=1 or layouts[0]['gauge']!=str(int(mode!='pressure')):
                            raise RuntimeError('fixture gauge layout differs')
                        rows=[dict(v.split('=',1) for v in line.split()[1:-1]) for line in log.splitlines() if line.startswith('immersed_static_mpi ') and line.endswith(' passed')]
                        if len(rows)!=1:
                            raise RuntimeError('missing original immersed acceptance record')
                        records.append(dict(layout=layouts[0],result=rows[0]))
                    job.update(status='passed',records=records)
                (root/'acceptance.json').write_text(json.dumps(report,indent=2)+'\n')
                print(mode,name,job['status'],result.returncode,flush=True)
                if name=='lu' and job['status']!='passed':
                    raise RuntimeError('LU reference failed; stop this evaluation')
        if digest(binary)!=report['binary_sha256']:
            raise RuntimeError('test binary changed during evaluation')
        report['status']='passed' if all(j['status']=='passed' for j in report['jobs']) else 'completed_with_failures'
    except BaseException:
        report['status']='failed'
        raise
    finally:
        (root/'acceptance.json').write_text(json.dumps(report,indent=2)+'\n')
    return 0 if report['status']=='passed' else 1


if __name__=='__main__':
    sys.exit(main())
