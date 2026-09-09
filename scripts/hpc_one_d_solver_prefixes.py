#!/usr/bin/env python3
"""Validate scoped 1D KSP/SNES options against an archived CLI executable."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]


def compare(first, second, exact):
    assert {p.name for p in first.iterdir()} == {p.name for p in second.iterdir()}
    metrics = {}
    for path in first.iterdir():
        other = second/path.name
        if path.name == 'summary.json':
            a, b = json.loads(path.read_text()), json.loads(other.read_text())
            for key in ('setup_seconds', 'solve_seconds', 'output_seconds', 'peak_rss_kib'):
                a.pop(key); b.pop(key)
            assert a['converged'] and b['converged']
            if exact:
                assert a == b
            else:
                for key in a:
                    if isinstance(a[key], float):
                        assert math.isfinite(b[key]) and abs(a[key]-b[key]) <= 1e-12+1e-6*abs(a[key]), key
                    else:
                        assert a[key] == b[key]
        elif exact:
            assert path.read_bytes() == other.read_bytes(), path.name
        elif path.suffix == '.csv':
            with path.open() as f, other.open() as g:
                a, b = list(csv.reader(f)), list(csv.reader(g))
            assert len(a) == len(b) and a[0] == b[0]
            for index, column in enumerate(a[0]):
                x, y = [], []
                for p, q in zip(a[1:], b[1:]):
                    assert len(p) == len(q) == len(a[0])
                    try:
                        x.append(float(p[index])); y.append(float(q[index]))
                    except ValueError:
                        assert p[index] == q[index]
                assert all(math.isfinite(v) for v in x+y)
                ref = math.sqrt(math.fsum(v*v for v in x))
                err = math.sqrt(math.fsum((v-w)**2 for v, w in zip(x, y)))
                assert err <= (1e-6*ref if ref else 1e-12), (path.name, column, ref, err)
                metrics[path.name+'/'+column] = dict(reference_l2=ref, difference_l2=err, relative_l2=err/ref if ref else 0)
        elif path.suffix == '.vtp':
            a, b = list(ET.parse(path).iter()), list(ET.parse(other).iter())
            assert len(a) == len(b)
            for x, y in zip(a, b):
                assert x.tag == y.tag and x.attrib == y.attrib
                if x.tag == 'DataArray':
                    first_values = (x.text or '').split(); second_values = (y.text or '').split()
                    assert len(first_values) == len(second_values)
                    if x.attrib.get('type', '').startswith(('Int', 'UInt')):
                        assert first_values == second_values
                    else:
                        for v, w in zip(first_values, second_values):
                            v, w = float(v), float(w)
                            assert math.isfinite(v) and math.isfinite(w) and abs(v-w) <= 1e-12+1e-6*abs(v)
        else:
            assert path.read_bytes() == other.read_bytes(), path.name
    assert exact or metrics
    return metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--ranks', type=int, default=3)
    args = parser.parse_args()
    assert args.ranks > 0
    root = args.output_root.resolve(); root.mkdir()
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps')
    report = dict(status='running', jobs=[], comparisons={}, binaries={str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest() for p in (args.binary, args.reference)})

    def job(name, case, flags=(), expected=0, reference=False):
        log = root/name; log.mkdir(); output = root/(name+'-output')
        command = ['mpiexec', '--oversubscribe', '-np', str(args.ranks)]
        if expected == 0:
            command += ['python3', str(ROOT/'scripts/hpc_rank_run.py'), '--output-dir', str(log), '--expected-ranks', str(args.ranks), '--timeout', '180', '--']
        command += [str((args.reference if reference else args.binary).resolve()), str(case), '--output-dir', str(output)] + list(flags)
        with (log/'launcher.log').open('w') as out:
            code = subprocess.call(['timeout', '--kill-after=5s', '195']+command, env=env, stdout=out, stderr=subprocess.STDOUT)
        report['jobs'].append(dict(name=name, command=command, returncode=code, expected=expected))
        print(name, code, 'expected', expected, flush=True); assert code == expected
        if expected:
            diagnostic = ('1d solver options: rank ' if name.endswith('-bad-ksp') else
                          '1d nonlinear options: rank ' if name.endswith('-bad-snes') else
                          'factor backend availability: rank ')
            assert diagnostic in (log/'launcher.log').read_text()
        if expected == 0:
            for rank in range(args.ranks):
                p = log/f'rank-{rank}'
                d = json.loads((p/'run.json').read_text()); assert d['returncode']==0 and not d['timed_out']
                assert (p/'stderr.log').stat().st_size == 0
        return output

    try:
        for mode in ('pressure_network', 'linearized_aq', 'nonlinear_aq', 'implicit_1d_pde'):
            case = root/'fixtures'/mode
            shutil.copytree(ROOT/'examples/one_d/compliant_bifurcation', case)
            p = case/'simulation_config.json'; d = json.loads(p.read_text())
            d['time'].update(steps=4, output_every=1); d['equation_systems'][0].update(scheme='implicit_petsc', formulation=mode)
            p.write_text(json.dumps(d, indent=2)+'\n')
            baseline = job(mode+'-reference', case, reference=True)
            default = job(mode+'-default', case); compare(baseline, default, True)
            flags = ['-domain_blood_flow_1d_flow_ksp_type', 'fgmres']
            nonlinear = mode in ('nonlinear_aq', 'implicit_1d_pde')
            if nonlinear:
                flags += ['-domain_blood_flow_1d_flow_snes_type', 'newtontr']
            override = job(mode+'-override', case, flags)
            report['comparisons'][mode] = compare(baseline, override, False)
            rows = [json.loads(line.split(' ', 1)[1]) for line in (root/(mode+'-override')/'rank-0/stdout.log').read_text().splitlines() if line.startswith('one_d_solver_configuration ')]
            assert len(rows)==4
            for row in rows:
                assert row['prefix']=='domain_blood_flow_1d_flow_' and row['ksp']=='fgmres' and row['pc']=='lu' and row['factor_backend']=='mumps'
                assert row['last_reason']>0
                if nonlinear:
                    assert row['snes']=='newtontr' and row['snes_reason']>0
            job(mode+'-bad-ksp', case, ['-domain_blood_flow_1d_flow_ksp_type', 'not_an_iga_solver'], expected=1)
            if nonlinear:
                job(mode+'-bad-snes', case, ['-domain_blood_flow_1d_flow_snes_type', 'not_an_iga_solver'], expected=1)
            if args.ranks>1:
                job(mode+'-bad-backend', case, ['-domain_blood_flow_1d_flow_pc_factor_mat_solver_type', 'petsc'], expected=1)
        report['status']='passed'
    except Exception as error:
        report.update(status='failed', error=str(error)); raise
    finally:
        (root/'acceptance.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
