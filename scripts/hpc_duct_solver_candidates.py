#!/usr/bin/env python3
"""Evaluate solver candidates on a refined C2 duct with original Phase 9 gates."""
import argparse
import json
import os
import re
import shlex
from pathlib import Path
import subprocess
import sys

from hpc_inventory import digest
from hpc_solver_prefixes import field_comparison


def main():
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--fixture-builder', type=Path, default=repo/'solvers/coupling/hpc_duct_solver_fixture')
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--transverse', type=int, default=2)
    parser.add_argument('--axial', type=int, default=4)
    parser.add_argument('--ranks', type=int, default=2)
    parser.add_argument('--timeout', type=int, default=600)
    parser.add_argument('--candidates', nargs='*', choices=['bjacobi', 'schur-lu', 'schur-gamg', 'selfp-lu', 'selfp-gamg', 'selfp-gamg-lu'],
                        default=['bjacobi', 'schur-lu', 'schur-gamg'], help='LU reference is always evaluated first; an empty list runs only LU')
    parser.add_argument('--launcher', default='mpiexec', help='MPI launcher and binding arguments; -np is appended')
    parser.add_argument('--solver-view', action='store_true', help='capture nested diagnostics; not a timing configuration')
    args = parser.parse_args()
    if min(args.transverse, args.axial, args.ranks, args.timeout) < 1:
        parser.error('dimensions, ranks and timeout must be positive')
    launcher = shlex.split(args.launcher)
    if not launcher:
        parser.error('launcher must not be empty')
    binary, builder = args.binary.resolve(), args.fixture_builder.resolve()
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    common = '-domain_root3d_flow_ksp_type fgmres -domain_root3d_flow_ksp_rtol 1e-12 '
    split = common + ('-domain_root3d_flow_pc_type fieldsplit -domain_root3d_flow_pc_fieldsplit_block_size 4 '
                     '-domain_root3d_flow_pc_fieldsplit_0_fields 0,1,2 -domain_root3d_flow_pc_fieldsplit_1_fields 3 '
                     '-domain_root3d_flow_pc_fieldsplit_type schur -domain_root3d_flow_pc_fieldsplit_schur_fact_type full '
                     '-domain_root3d_flow_pc_fieldsplit_schur_precondition a11 '
                     '-domain_root3d_flow_fieldsplit_0_ksp_type preonly -domain_root3d_flow_fieldsplit_1_ksp_type preonly ')
    configurations = {
        'lu': common+'-domain_root3d_flow_pc_type lu -domain_root3d_flow_pc_factor_mat_solver_type mumps',
        'bjacobi': common+'-domain_root3d_flow_pc_type bjacobi -domain_root3d_flow_sub_ksp_type preonly -domain_root3d_flow_sub_pc_type ilu',
        'schur-lu': split+'-domain_root3d_flow_fieldsplit_0_pc_type lu -domain_root3d_flow_fieldsplit_0_pc_factor_mat_solver_type mumps -domain_root3d_flow_fieldsplit_1_pc_type lu -domain_root3d_flow_fieldsplit_1_pc_factor_mat_solver_type mumps',
        'schur-gamg': split+'-domain_root3d_flow_fieldsplit_0_pc_type gamg -domain_root3d_flow_fieldsplit_1_pc_type jacobi',
    }
    configurations['selfp-lu'] = configurations['schur-lu'].replace('schur_precondition a11', 'schur_precondition selfp')
    configurations['selfp-gamg'] = configurations['schur-gamg'].replace('schur_precondition a11', 'schur_precondition selfp')
    configurations['selfp-gamg-lu'] = configurations['selfp-lu'].replace('fieldsplit_0_pc_type lu', 'fieldsplit_0_pc_type gamg').replace('-domain_root3d_flow_fieldsplit_0_pc_factor_mat_solver_type mumps ', '')
    selected = ['lu'] + list(dict.fromkeys(args.candidates))
    report = dict(status='running', binaries={str(p): digest(p) for p in (binary, builder)},
                  mesh=dict(transverse=args.transverse, axial=args.axial, ranks=args.ranks), candidates=[],
                  harness_sha256=digest(Path(__file__)), selected=selected, solver_view=args.solver_view)
    try:
        for name in selected:
            options = configurations[name]
            if args.solver_view:
                options += ' -domain_root3d_flow_ksp_view'
            folder = root/name
            folder.mkdir()
            case = folder/'case'
            subprocess.run([str(builder), str(case), str(args.transverse), str(args.axial), str(args.ranks)], check=True)
            inputs = {str(p): digest(p) for p in case.rglob('*') if p.is_file()}
            output, bundle = case/'result', case/'bundle'
            command = ['timeout', '--kill-after=5s', str(args.timeout+15), *launcher, '-np', str(args.ranks),
                       sys.executable, str(repo/'scripts/hpc_rank_run.py'), '--output-dir', str(folder),
                       '--expected-ranks', str(args.ranks), '--timeout', str(args.timeout), '--', str(binary),
                       '--graph-case', str(case), '--output-dir', str(output), '--checkpoint-dir', str(bundle)]
            env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', IGA_PROFILE='1', PETSC_OPTIONS=options)
            with (folder/'launcher.log').open('w') as log:
                result = subprocess.run(command, cwd=repo, env=env, stdout=log, stderr=subprocess.STDOUT)
            row = dict(name=name, command=command, options=options, inputs=inputs, returncode=result.returncode, status='failed')
            report['candidates'].append(row)
            # Keep resource and failure evidence even when the native solve fails.
            row['rank_reports'], row['stderr_tails'] = [], {}
            for rank in range(args.ranks):
                directory = folder/f'rank-{rank}'
                measured_path = directory/'run.json'
                if measured_path.is_file():
                    row['rank_reports'].append(json.loads(measured_path.read_text()))
                stderr_path = directory/'stderr.log'
                if stderr_path.is_file():
                    row['stderr_tails'][str(rank)] = stderr_path.read_text()[-8192:]
            if result.returncode:
                row['error'] = 'native process failed or timed out; inspect rank reports and stderr tails'
            if result.returncode == 0:
                try:
                    validation = subprocess.run([str(builder), '--validate', str(output)], capture_output=True, text=True)
                    row['validation'] = dict(returncode=validation.returncode, stdout=validation.stdout, stderr=validation.stderr)
                    if validation.returncode:
                        raise RuntimeError('original coupled physical gates failed')
                    row['profiles'] = []
                    for rank in range(args.ranks):
                        measured = json.loads((folder/f'rank-{rank}/run.json').read_text())
                        if measured['returncode'] or measured['timed_out'] or not measured['resource']:
                            raise RuntimeError('rank failed or lacks measurement')
                        log = (folder/f'rank-{rank}/stdout.log').read_text()
                        profiles = [json.loads(line.split(' ', 1)[1]) for line in log.splitlines() if line.startswith('hpc_profile ')]
                        if len(profiles) != 1 or profiles[0]['status'] != 0:
                            raise RuntimeError('missing final phase profile')
                        row['profiles'].extend(profiles)
                        if rank == 0:
                            configured = [json.loads(line.split(' ', 1)[1]) for line in log.splitlines()
                                          if line.startswith('solver_configuration ')]
                            expected_pc = {'lu': 'lu', 'bjacobi': 'bjacobi'}.get(name, 'fieldsplit')
                            if len(configured) != 2 or any(c['domain'] != 'root3d' or c['prefix'] != 'domain_root3d_flow_'
                                    or c['ksp'] != 'fgmres' or c['pc'] != expected_pc or c['last_reason'] <= 0 for c in configured):
                                raise RuntimeError('effective solver configuration differs')
                            row['solver_configuration'] = configured
                            row['total_linear_iterations'] = sum(map(int, re.findall(r'linear_iterations=(\d+)', log)))
                            if args.solver_view:
                                if 'KSP Object: (domain_root3d_flow_)' not in log:
                                    raise RuntimeError('missing solver view')
                                if expected_pc == 'fieldsplit':
                                    for child in ('0', '1'):
                                        if 'KSP Object: (domain_root3d_flow_fieldsplit_'+child+'_)' not in log:
                                            raise RuntimeError('missing fieldsplit child view')
                    if name != 'lu':
                        row['fields'] = field_comparison(root/'lu/case/bundle', bundle, False)
                    row['status'] = 'passed'
                except Exception as error:
                    row['error'] = repr(error)
            if any(digest(Path(p)) != value for p, value in inputs.items()):
                raise RuntimeError('fixture changed during solve')
            (root/'acceptance.json').write_text(json.dumps(report, indent=2)+'\n')
            print(name, row['status'], row.get('error', ''), flush=True)
            if name == 'lu' and row['status'] != 'passed':
                raise RuntimeError('reference candidate failed; cannot compare fields')
        if any(digest(Path(p)) != value for p, value in report['binaries'].items()):
            raise RuntimeError('binary changed during evaluation')
        report['status'] = 'passed' if all(r['status']=='passed' for r in report['candidates']) else 'completed_with_failures'
    except BaseException:
        report['status'] = 'failed'
        raise
    finally:
        (root/'acceptance.json').write_text(json.dumps(report, indent=2)+'\n')
    return 0 if report['status']=='passed' else 1


if __name__ == '__main__':
    sys.exit(main())
