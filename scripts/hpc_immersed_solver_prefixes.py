#!/usr/bin/env python3
"""Exercise immersed per-domain options and numerical regression fixtures."""
import argparse, subprocess, pathlib, os, json, sys
from hpc_inventory import digest
from hpc_immersed_graph_regression import fixture, table, compare_ports

def main():
    parser = argparse.ArgumentParser(description='Exercise immersed solver family and domain prefixes.')
    parser.add_argument('--output-dir', type=pathlib.Path, required=True)
    parser.add_argument('--suite', choices=['core', 'native'], default='core')
    parser.add_argument('--only', nargs='+', help='Run only named cases, for focused follow-up checks.')
    parser.add_argument('--timeout', type=int, default=1200)
    options = parser.parse_args()
    repo = pathlib.Path(__file__).resolve().parents[1]
    os.chdir(repo)
    root = options.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps')
    jobs = [('static-default', 2, 'immersed_distributed_static_flow_test', ['flow'], None), ('static-scoped', 2, 'immersed_distributed_static_flow_test', ['flow', '-test_solver_prefix', '-domain_test_flow_ksp_type', 'fgmres'], ('domain_test_flow_', 'fgmres')), ('transient-default', 2, 'immersed_transient_distributed_runtime_test', ['flow'], None), ('transient-scoped', 2, 'immersed_transient_distributed_runtime_test', ['flow', '-test_solver_prefix', '-domain_test_flow_ksp_type', 'fgmres'], ('domain_test_flow_', 'fgmres')), ('transient-split', 3, 'immersed_transient_distributed_runtime_test', ['flow', 'split', '-test_solver_prefix', '-domain_test_flow_ksp_type', 'fgmres'], ('domain_test_flow_', 'fgmres')), ('moving-options', 1, 'moving_immersed_solver_options_test', [], None), ('moving-regression', 1, 'moving_immersed_transient_flow_test', [], None), ('fsi-regression', 1, 'moving_immersed_transient_flow_fsi_runtime_test', [], None)]
    if options.suite == 'native':
        jobs = []
        for transient in [False, True]:
            family = 'immersed_transient_' if transient else 'immersed_static_'
            mode = 'transient' if transient else 'static'
            folder = root / (mode + '-fixture')
            folder.mkdir()
            case = fixture(repo, folder, 'fixed', transient)
            for (label, ranks, flags) in [('reference', 1, []), ('parallel', 2, []), ('override', 2, ['-domain_immersed_flow_ksp_type', 'fgmres'])]:
                name = mode + '-' + label
                jobs.append((name, ranks, 'iga_multidomain_flow', ['--graph-case', str(case), '--output-dir', str(root / name / 'result'), '-' + family + 'ksp_type', 'gmres'] + flags, ('domain_immersed_flow_', 'fgmres' if label == 'override' else 'gmres')))
    if options.only:
        selected = set(options.only)
        if selected - {j[0] for j in jobs}:
            parser.error('unknown --only case')
        jobs = [j for j in jobs if j[0] in selected]
    if options.timeout < 1:
        parser.error('timeout must be positive')
    records = []
    references = {}
    manifest = {str(p.relative_to(repo)): digest(p) for base in ['include', 'solvers/cpu/include'] for p in (repo / base).glob('*.hpp')}
    for source in [pathlib.Path(__file__), repo / 'solvers/coupling/src/iga_1d_3d_bifurcation.cpp', repo / 'solvers/cpu/Makefile', repo / 'solvers/coupling/Makefile']:
        manifest[str(source.relative_to(repo))] = digest(source)
    for source in (repo / 'solvers/cpu/tests').glob('test_*immersed*.cpp'):
        manifest[str(source.relative_to(repo))] = digest(source)
    frozen = {str(repo / ('solvers/coupling' if options.suite == 'native' else 'solvers/cpu') / j[2]): digest(repo / ('solvers/coupling' if options.suite == 'native' else 'solvers/cpu') / j[2]) for j in jobs}
    (root / 'binaries.json').write_text(json.dumps(frozen, indent=2) + '\n')
    (root / 'harness.py').write_text(pathlib.Path(__file__).read_text())
    (root / 'source.json').write_text(json.dumps(manifest, indent=2) + '\n')
    for (name, n, b, args, expected) in jobs:
        d = root / name
        d.mkdir()
        job_env = dict(env)
        if options.suite == 'core' and '-test_solver_prefix' in args:
            offset = args.index('-test_solver_prefix')
            job_env['PETSC_OPTIONS'] += ' ' + ' '.join(args[offset:])
            args = args[:offset]
        cmd = ['mpiexec', '--oversubscribe', '-np', str(n), sys.executable, 'scripts/hpc_rank_run.py', '--output-dir', str(d), '--expected-ranks', str(n), '--timeout', str(options.timeout), '--', ('solvers/coupling/' if options.suite == 'native' else 'solvers/cpu/') + b] + args
        with (d / 'launcher.log').open('w') as f:
            r = subprocess.run(cmd, env=job_env, stdout=f, stderr=subprocess.STDOUT)
        record = dict(name=name, argv=cmd, returncode=r.returncode, binary_sha256=digest(repo / ('solvers/coupling' if options.suite == 'native' else 'solvers/cpu') / b))
        records.append(record)
        (root / 'summary.json').write_text(json.dumps(records, indent=2) + '\n')
        if r.returncode:
            raise RuntimeError(name + ' failed')
        if any((digest(pathlib.Path(p)) != h for (p, h) in frozen.items())):
            raise RuntimeError('binary changed during validation')
        for i in range(n):
            report = json.loads((d / f'rank-{i}/run.json').read_text())
            assert report['returncode'] == 0 and (not report['timed_out']) and report['resource']
            if expected and options.suite == 'core':
                assert f'prefix={expected[0]} ksp={expected[1]}' in (d / f'rank-{i}/stdout.log').read_text()
        if options.suite == 'native':
            configurations = [json.loads(line.split(' ', 1)[1]) for line in (d / 'rank-0/stdout.log').read_text().splitlines() if line.startswith('solver_configuration ')]
            selected = [v for v in configurations if v['domain'] == 'immersed']
            assert len(selected) == 2 and all((v['prefix'] == expected[0] and v['ksp'] == expected[1] and (v['pc'] == 'lu') and (v['last_reason'] > 0) for v in selected))
            mode = name.split('-')[0]
            ports = table(d / 'result/pressure_flow_ports.csv')
            if mode not in references:
                references[mode] = ports
            record['port_error_fraction'] = compare_ports(references[mode], ports)
            (root / 'summary.json').write_text(json.dumps(records, indent=2) + '\n')
        print(name + ' passed', flush=True)
    assert all((digest(repo / p) == h for (p, h) in manifest.items())), 'source changed during validation'

if __name__ == "__main__":
    main()
