#!/usr/bin/env python3
"""Reject unsupported or inconsistent transient graph case configurations collectively."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

from hpc_immersed_graph_regression import fixture
from hpc_inventory import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    binary = repo/'solvers/coupling/iga_multidomain_flow'
    summary = dict(status='running', binary=str(binary), binary_sha256=digest(binary), cases=[])
    variants = [('dt', 'transient immersed domain and graph time grids must match'),
                ('steps', 'transient immersed domain and graph time grids must match'),
                ('negative-inertia', 'invalid transient wall impedance or flow tolerance'),
                ('steady-inertia', 'unknown immersed_geometry.json.runtime key'),
                ('moving', 'unknown immersed_geometry.json key'),
                ('species', 'immersed flow case rejects transport/species fields')]
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    for key in ['PETSC_OPTIONS', 'TUBULARFLOWIGA_INJECT_EXPLICIT_COUPLING_FAILURE_STEP',
                'TUBULARFLOWIGA_INJECT_BIFURCATION_FAILURE_STEP']:
        env.pop(key, None)
    try:
        for name, expected in variants:
            directory = root/name
            directory.mkdir()
            case = fixture(repo, directory, 'explicit', transient=True)
            simulation_path = case/'immersed/simulation_config.json'
            geometry_path = case/'immersed/immersed_geometry.json'
            simulation = json.loads(simulation_path.read_text())
            geometry = json.loads(geometry_path.read_text())
            if name == 'dt':
                simulation['time']['dt'] *= 2
            elif name == 'steps':
                simulation['time']['steps'] += 1
            elif name == 'negative-inertia':
                geometry['runtime']['wall_inertial_gamma0'] = -1
            elif name == 'steady-inertia':
                simulation['equation_systems'][0]['time_integration'] = 'steady'
                geometry['runtime']['wall_inertial_gamma0'] = .6
            elif name == 'moving':
                geometry['motion'] = {}
            elif name == 'species':
                simulation['fields'].append(dict(name='species', kind='scalar'))
            simulation_path.write_text(json.dumps(simulation, indent=2)+'\n')
            geometry_path.write_text(json.dumps(geometry, indent=2)+'\n')
            inputs = {str(p.relative_to(case)): digest(p) for p in case.rglob('*') if p.is_file()}
            command = ['timeout', '--kill-after=5s', '150s', 'mpiexec', '--oversubscribe', '-np', '2',
                       sys.executable, str(repo/'scripts/hpc_rank_run.py'), '--expected-ranks', '2',
                       '--timeout', '120', '--output-dir', str(directory), '--', str(binary),
                       '--graph-case', str(case), '--output-dir', str(directory/'result')]
            record = dict(name=name, argv=command, expected_error=expected, input_sha256=inputs, rank_reports=[])
            summary['cases'].append(record)
            with (directory/'launcher.log').open('w') as log:
                result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT)
            record['returncode'] = result.returncode
            if not result.returncode or (directory/'result').exists():
                raise RuntimeError('invalid immersed configuration succeeded or published outputs')
            for rank in range(2):
                rd = directory/f'rank-{rank}'
                report = json.loads((rd/'run.json').read_text())
                if (report['returncode'] != 1 or report['timed_out'] or not report['resource']
                        or (rank == 0 and expected not in (rd/'stderr.log').read_text())):
                    raise RuntimeError('missing coordinated case rejection: '+name)
                record['rank_reports'].append(report)
            if any(digest(case/p) != h for p, h in inputs.items()):
                raise RuntimeError('rejection fixture changed during validation')
            record['status'] = 'passed'
            print(name, 'passed', flush=True)
        if digest(binary) != summary['binary_sha256']:
            raise RuntimeError('binary changed during rejection validation')
        summary['status'] = 'passed'
    except BaseException:
        summary['status'] = 'failed'
        raise
    finally:
        (root/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
