#!/usr/bin/env python3
"""Exercise the CPU flow PVTU CLI with actual solves and checkpoint parity."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture-root', type=Path, required=True,
                        help='existing flow-output regression containing vtkhdf-1/input and vtkhdf-2/input')
    parser.add_argument('--reference-binary', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    binary = repo / 'solvers/cpu/iga_navier_stokes'
    reference = args.reference_binary.resolve()
    inputs = {}
    for ranks in (1, 2):
        folder = root / f'input-{ranks}'
        shutil.copytree(args.fixture_root / f'vtkhdf-{ranks}' / 'input', folder)
        inputs[ranks] = folder
    evidence = dict(status='running', binaries={str(p): digest(p) for p in (binary, reference)},
                    inputs={str(p): digest(p) for d in inputs.values() for p in d.iterdir()}, jobs=[])
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -ksp_rtol 1e-12')

    def run(name, ranks, fmt, every=1, executable=binary, extra=(), damage=None):
        folder = root / name
        folder.mkdir()
        output = folder / 'result'
        output.mkdir()
        if damage:
            (output / damage).mkdir()
        command = ['timeout', '--kill-after=5s', '150s', 'mpiexec', '-np', str(ranks), sys.executable,
                   str(repo / 'scripts/hpc_rank_run.py'), '--output-dir', str(folder),
                   '--expected-ranks', str(ranks), '--timeout', '120', '--', str(executable),
                   str(inputs[ranks] / 'database.ntiga'), str(inputs[ranks]), '--visualization-format', fmt,
                   '--output', str(output / 'flow.txt'), '--checkpoint', str(output / 'checkpoint')]
        if every:
            command += ['--output-every', str(every)]
        command += list(extra)
        with (folder / 'launcher.log').open('w') as log:
            result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT)
        row = dict(name=name, ranks=ranks, command=command, returncode=result.returncode)
        evidence['jobs'].append(row)
        reports = [json.loads((folder / f'rank-{rank}/run.json').read_text()) for rank in range(ranks)]
        row['rank_reports'] = reports
        assert all(not report['timed_out'] for report in reports)
        if damage:
            assert result.returncode != 0 and all(report['returncode'] != 0 for report in reports)
            assert all('navier_stokes_v2 seconds=' not in (folder / f'rank-{rank}/stdout.log').read_text() for rank in range(ranks))
            if damage.endswith('step000001'):
                times = [float(node.attrib['timestep']) for node in ET.parse(output / 'flow.pvd').iter('DataSet')]
                assert times == [0.0]
            else:
                assert not (output / 'flow.pvd').exists()
        else:
            assert result.returncode == 0 and all(report['returncode'] == 0 for report in reports)
            assert (output / 'checkpoint.state').is_file()
            if fmt == 'pvtu':
                assert not (output / 'flow.txt').exists() and not (output / 'flow.txt.pressure').exists()
                assert not (output / 'flow.series.csv').exists() and not (output / 'flow.vtkhdf').exists()
                frames = list(ET.parse(output / 'flow.pvd').iter('DataSet'))
                row['times'] = [float(frame.attrib['timestep']) for frame in frames]
                assert row['times'] == sorted(set(row['times']))
                for frame in frames:
                    index = output / frame.attrib['file']
                    pieces = list(ET.parse(index).iter('Piece'))
                    assert len(pieces) == ranks
                    assert all((index.parent / piece.attrib['Source']).is_file() for piece in pieces)
        (root / 'acceptance.json').write_text(json.dumps(evidence, indent=2)+'\n')
        print(name, 'expected rejection' if damage else 'passed', flush=True)
        return output

    try:
        for ranks in (1, 2):
            baseline = run(f'reference-{ranks}', ranks, 'vtkhdf', executable=reference)
            legacy = run(f'legacy-{ranks}', ranks, 'vtkhdf')
            parallel = run(f'parallel-{ranks}', ranks, 'pvtu')
            assert digest(baseline / 'checkpoint.state') == digest(legacy / 'checkpoint.state') == digest(parallel / 'checkpoint.state')
        final = run('final-1', 1, 'pvtu', every=0)
        sparse = run('sparse-2', 2, 'pvtu', every=2)
        stopped = run('stopped-2', 2, 'pvtu', extra=['--stop-after-step', '1'])
        restarted = run('restarted-2', 2, 'pvtu', extra=['--restart', str(stopped / 'checkpoint')])
        assert digest(restarted / 'checkpoint.state') == digest(root / 'parallel-2/result/checkpoint.state')
        assert next(row['times'] for row in evidence['jobs'] if row['name']=='final-1') == [.02]
        assert next(row['times'] for row in evidence['jobs'] if row['name']=='sparse-2') == [0., .02]
        assert next(row['times'] for row in evidence['jobs'] if row['name']=='stopped-2') == [0., .01]
        assert next(row['times'] for row in evidence['jobs'] if row['name']=='restarted-2') == [.01, .02]
        run('blocked-step-2', 2, 'pvtu', damage='flow.step000001')
        run('blocked-index-2', 2, 'pvtu', damage='flow.pvd.pending')
        assert all(digest(Path(p))==expected for p, expected in evidence['inputs'].items())
        assert all(digest(Path(p))==expected for p, expected in evidence['binaries'].items())
        evidence['status'] = 'passed'
    except BaseException:
        evidence['status'] = 'failed'
        raise
    finally:
        (root / 'acceptance.json').write_text(json.dumps(evidence, indent=2)+'\n')


if __name__ == '__main__':
    main()
