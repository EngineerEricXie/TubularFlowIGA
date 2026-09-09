#!/usr/bin/env python3
import argparse, json, subprocess, os, hashlib, math, struct
from pathlib import Path
root = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description='Replay healthy CPU CLI fixtures with independent solver options.')
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--reference-dir', type=Path, required=True)
    parser.add_argument('--flow-manifest', type=Path, required=True)
    parser.add_argument('--transport-manifest', type=Path, required=True)
    args = parser.parse_args()
    p = args.output_root.resolve()
    p.mkdir()
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', IGA_PROFILE='1')
    report = dict(status='running', jobs=[], comparisons=[], binaries={})
    for name in ['iga_solve', 'iga_navier_stokes']:
        for binary in [root / 'solvers/cpu' / name, args.reference_dir.resolve() / (name + '-reference')]:
            report['binaries'][str(binary)] = hashlib.sha256(binary.read_bytes()).hexdigest()
    try:
        for (stage, names, binary, role) in [('flow-step', ['healthy-1', 'healthy-2'], 'iga_navier_stokes', 'flow'), ('transport-cli', ['prescribed-one', 'prescribed-two', 'series-one', 'series-two'], 'iga_solve', 'transport')]:
            d = json.loads((args.flow_manifest if role == 'flow' else args.transport_manifest).read_text())
            for name in names:
                c = next((x for x in d['cases'] if x['case'] == name))
                for (filename, digest) in c['input_sha256'].items():
                    assert hashlib.sha256(Path(filename).read_bytes()).hexdigest() == digest, filename
                saved = c['ranks'][0]
                env['PETSC_OPTIONS'] = saved['environment']['PETSC_OPTIONS']
                ranks = saved['ranks']
                argv = saved['command_argv']
                case = Path(argv[2])
                config = json.loads((case / 'simulation_config.json').read_text())
                system = next((s['name'] for s in config['equation_systems'] if (s['kind'] == 'navier_stokes' if role == 'flow' else s['name'] == 'transport')))
                prefix = 'domain_' + system + '_' + role + '_'
                outputs = {}
                for kind in ['reference', 'default', 'override', 'bad'] + (['backend'] if ranks > 1 else []):
                    job = p / (stage + '-' + name + '-' + kind)
                    job.mkdir()
                    out = job / 'result'
                    out.mkdir()
                    cmdargs = [str(args.reference_dir.resolve() / (binary + '-reference') if kind == 'reference' else root / 'solvers/cpu' / binary)] + argv[1:]
                    for flag in ['--output', '--checkpoint', '--memory-report']:
                        if flag in cmdargs:
                            i = cmdargs.index(flag) + 1
                            cmdargs[i] = str(out / Path(cmdargs[i]).name)
                    if kind in ['override', 'bad']:
                        cmdargs += ['-' + prefix + 'ksp_type', 'fgmres' if kind == 'override' else 'not_an_iga_solver']
                    if kind == 'backend':
                        cmdargs += ['-' + prefix + 'pc_factor_mat_solver_type', 'petsc']
                    expected = 1 if kind in ['bad', 'backend'] else 0
                    cmd = ['mpiexec', '--oversubscribe', '-np', str(ranks)]
                    if not expected:
                        cmd += ['python3', str(root / 'scripts/hpc_rank_run.py'), '--output-dir', str(job), '--expected-ranks', str(ranks), '--timeout', '120', '--']
                    cmd += cmdargs
                    with (job / 'launcher.log').open('w') as f:
                        code = subprocess.call(['timeout', '135'] + cmd, env=env, stdout=f, stderr=subprocess.STDOUT)
                    report['jobs'].append(dict(name=job.name, command=cmd, returncode=code, expected=expected))
                    print(job.name, code, flush=True)
                    assert code == expected
                    if expected:
                        diagnostic = 'factor backend availability:' if kind == 'backend' else ('flow solver options:' if role == 'flow' else 'transport KSPSetFromOptions:')
                        assert diagnostic in (job / 'launcher.log').read_text()
                    if not expected:
                        for rank in range(ranks):
                            rr = json.loads((job / f'rank-{rank}/run.json').read_text())
                            assert rr['returncode'] == 0 and (not rr['timed_out'])
                            assert (job / f'rank-{rank}/stderr.log').stat().st_size == 0
                        outputs[kind] = out
                        if kind == 'override':
                            assert 'prefix=' + prefix + ' ksp=fgmres pc=lu factor_backend=mumps' in (job / 'rank-0/stdout.log').read_text()
                (a, b) = (outputs['reference'], outputs['default'])
                assert {x.name for x in a.iterdir()} == {x.name for x in b.iterdir()}
                for f in a.iterdir():
                    assert f.read_bytes() == (b / f.name).read_bytes(), f.name
                def compare_metadata(first, second):
                    assert type(first) is type(second)
                    if isinstance(first, dict):
                        assert first.keys() == second.keys()
                        for key in first:
                            compare_metadata(first[key], second[key])
                    elif isinstance(first, list):
                        assert len(first) == len(second)
                        for x, y in zip(first, second):
                            compare_metadata(x, y)
                    elif isinstance(first, float):
                        assert math.isfinite(first) and math.isfinite(second)
                        assert abs(first-second) <= 1e-12+1e-6*abs(first)
                    else:
                        assert first == second
                for f in a.glob('*.json'):
                    compare_metadata(json.loads(f.read_text()), json.loads((outputs['override']/f.name).read_text()))
                for f in a.iterdir():
                    if f.suffix not in ('.txt', '.pressure', '.state'):
                        continue
                    def values(path):
                        if path.suffix != '.state':
                            return [float(v) for v in path.read_text().split()]
                        payload = path.read_bytes()
                        tag, count = struct.unpack('>ii', payload[:8])
                        assert tag == 1211214 and count >= 0 and len(payload) == 8+count*8
                        return list(struct.unpack('>'+'d'*count, payload[8:]))
                    x, y = values(f), values(outputs['override']/f.name)
                    assert len(x) == len(y)
                    norm = math.sqrt(math.fsum((v * v for v in x)))
                    err = math.sqrt(math.fsum(((v - w) ** 2 for (v, w) in zip(x, y))))
                    assert all((math.isfinite(v) for v in x + y))
                    assert err <= 1e-12 + 1e-06 * norm
                    report['comparisons'].append(dict(case=name, file=f.name, relative_l2=err / norm if norm else 0))
        report['status'] = 'passed'
    finally:
        (p / 'acceptance.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == "__main__":
    main()
