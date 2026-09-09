#!/usr/bin/env python3
"""Exercise the production steady or fixed-transient immersed graph entry across local MPI ranks."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

from hpc_inventory import digest


def fixture(repo, root, execution, transient=False, wall_inertial_gamma0=None):
    source = repo/'examples/vascular_flow/immersed_aneurysm_chain'
    target = root/'fixture'
    shutil.copytree(source, target)
    configurations = ['simulation_config.json', 'source/simulation_config.json', 'sink/simulation_config.json']
    if transient:
        configurations.append('immersed/simulation_config.json')
    for relative in configurations:
        path = target/relative
        config = json.loads(path.read_text())
        config['time']['steps'] = 2
        if relative == 'immersed/simulation_config.json':
            config['equation_systems'][0]['time_integration'] = 'backward_euler'
        if relative == 'simulation_config.json':
            config['execution']['kind'] = execution
            config['execution']['maximum_iterations'] = 1 if execution == 'explicit' else 64
        path.write_text(json.dumps(config, indent=2)+'\n')
    path = target/'immersed/immersed_geometry.json'
    geometry = json.loads(path.read_text())
    geometry['grid'] = dict(lower_m=[0, 0, 0], upper_m=[1, 1, 1], cells=[4, 1, 1])
    geometry['volume_quadrature']['max_depth'] = 2
    if wall_inertial_gamma0 is not None:
        geometry['runtime']['wall_inertial_gamma0'] = wall_inertial_gamma0
    path.write_text(json.dumps(geometry, indent=2)+'\n')
    (target/'immersed/surface.vtp').write_text('''<?xml version="1.0"?>
<VTKFile type="PolyData" version="1.0" byte_order="LittleEndian"><PolyData><Piece NumberOfPoints="8" NumberOfPolys="12"><Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">0 0 0 1 0 0 1 1 0 0 1 0 0 0 1 1 0 1 1 1 1 0 1 1</DataArray></Points><Polys><DataArray type="Int32" Name="connectivity" format="ascii">0 2 1 0 3 2 4 5 6 4 6 7 0 1 5 0 5 4 1 2 6 1 6 5 2 3 7 2 7 6 3 0 4 3 4 7</DataArray><DataArray type="Int32" Name="offsets" format="ascii">3 6 9 12 15 18 21 24 27 30 33 36</DataArray></Polys><CellData Scalars="boundary_id"><DataArray type="UInt32" Name="boundary_id" format="ascii">1 1 2 2 0 0 0 0 0 0 0 0</DataArray></CellData></Piece></PolyData></VTKFile>
''')
    return target


def table(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def compare_ports(reference, actual):
    keys = ['step', 'domain_id', 'port_id']
    def indexed(rows):
        result = {tuple(row[k] for k in keys): row for row in rows}
        if len(result) != len(rows):
            raise RuntimeError('duplicate accepted graph port')
        return result
    expected, observed = indexed(reference), indexed(actual)
    if expected.keys() != observed.keys():
        raise RuntimeError('accepted graph port coverage differs')
    worst = 0.0
    for key, ref in expected.items():
        for field in ['time_s', 'area_m2', 'outward_flow_m3_s', 'mean_pressure_pa']:
            a, b = float(observed[key][field]), float(ref[field])
            error, tolerance = abs(a-b), 1e-12+1e-6*abs(b)
            if not math.isfinite(a) or not math.isfinite(b) or error > tolerance:
                raise RuntimeError('graph port parity failed: '+str((key, field, a, b)))
            worst = max(worst, error/tolerance)
    return worst


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--wall-inertial-gamma0', type=float, help='Optional transient wall impedance coefficient')
    parser.add_argument('--transient', action='store_true', help='Use the fixed-geometry backward-Euler immersed backend')
    parser.add_argument('--execution', choices=['explicit', 'fixed', 'aitken'], default='fixed')
    args = parser.parse_args()
    if args.wall_inertial_gamma0 is not None and (not args.transient or not math.isfinite(args.wall_inertial_gamma0) or args.wall_inertial_gamma0 < 0):
        parser.error('wall inertia requires --transient and a finite nonnegative value')
    repo = Path(__file__).resolve().parents[1]
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    case = fixture(repo, root, args.execution, args.transient, args.wall_inertial_gamma0)
    binary = repo/'solvers/coupling/iga_multidomain_flow'
    inputs = {str(p.relative_to(case)): digest(p) for p in case.rglob('*') if p.is_file()}
    result = dict(status='running', execution=args.execution, transient=args.transient, wall_inertial_gamma0=args.wall_inertial_gamma0, binary=str(binary), binary_sha256=digest(binary), input_sha256=inputs, cases=[])
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', IGA_PROFILE='1')
    for key in ['PETSC_OPTIONS', 'TUBULARFLOWIGA_INJECT_EXPLICIT_COUPLING_FAILURE_STEP', 'TUBULARFLOWIGA_INJECT_BIFURCATION_FAILURE_STEP']:
        env.pop(key, None)
    reference_ports = None
    reference_diagnostics = {}
    try:
        # A fresh successful process follows the failed 2-rank process. The
        # adapter test separately proves rollback/retry within one object.
        for ranks, failure, name in [(1, False, '1'), (2, False, '2'), (4, False, '4'), (2, True, '2-failure'), (2, False, '2-retry')]:
            directory = root/name
            directory.mkdir()
            command = ['timeout', '--kill-after=5s', '660s', 'mpiexec', '--oversubscribe', '-np', str(ranks),
                       sys.executable, str(repo/'scripts/hpc_rank_run.py'), '--expected-ranks', str(ranks),
                       '--timeout', '630', '--output-dir', str(directory), '--', str(binary),
                       '--graph-case', str(case), '--output-dir', str(directory/'result')]
            environment = dict(env)
            if failure:
                environment['TUBULARFLOWIGA_INJECT_EXPLICIT_COUPLING_FAILURE_STEP'] = '2'
            record = dict(ranks=ranks, injected_precommit_failure=failure, argv=command, rank_reports=[], distributions=[], transient_diagnostics=[])
            result['cases'].append(record)
            with (directory/'launcher.log').open('w') as log:
                completed = subprocess.run(command, env=environment, stdout=log, stderr=subprocess.STDOUT)
            record['returncode'] = completed.returncode
            if bool(completed.returncode) != failure:
                raise RuntimeError('unexpected graph launcher result: '+name)
            for rank in range(ranks):
                rd = directory/f'rank-{rank}'
                report = json.loads((rd/'run.json').read_text())
                record['rank_reports'].append(report)
                if report['timed_out'] or not report['resource'] or bool(report['returncode']) != failure:
                    raise RuntimeError('failed, timed out, or unmeasured graph rank')
                if ranks > 1 or args.transient:
                    lines = [line for line in (rd/'stdout.log').read_text().splitlines() if line.startswith('hpc_immersed_distribution ')]
                    if len(lines) != 1:
                        raise RuntimeError('missing distributed graph ownership report')
                    observation = json.loads(lines[0].split(' ', 1)[1])
                    if (observation['rank'] != rank or observation['ranks'] != ranks
                            or observation['owned_rows'] <= 0 or observation['owned_stencils'] <= 0
                            or not 0 < observation['required_rows'] <= observation['global_rows']
                            or (ranks > 1 and observation['required_rows'] == observation['global_rows'])
                            or observation['time_integration'] != ('backward_euler' if args.transient else 'steady')):
                        raise RuntimeError('invalid distributed graph ownership')
                    record['distributions'].append(observation)
                if args.transient:
                    lines = [line for line in (rd/'stdout.log').read_text().splitlines()
                             if line.startswith('hpc_immersed_transient_step ')]
                    if len(lines) != (1 if failure else 2):
                        raise RuntimeError('transient backend committed an incorrect number of steps')
                    for index, line in enumerate(lines, 1):
                        observation = json.loads(line.split(' ', 1)[1])
                        if (observation['rank'] != rank or observation['ranks'] != ranks
                                or observation['domain'] != 'immersed'
                                or observation['index'] != index or observation['commits'] != index
                                or observation['time_s'] != index*.01
                                or any(not math.isfinite(v) for k, v in observation.items() if k != 'domain')
                                or observation['assembly_s'] <= 0 or observation['solve_s'] <= 0
                                or not 0 <= observation['residual_norm'] < 1e-8):
                            raise RuntimeError('invalid accepted transient backend diagnostics')
                        if ranks == 1:
                            reference_diagnostics[index] = observation
                        for key in ['surface_flow', 'volume_divergence', 'wall_flow']:
                            expected = reference_diagnostics[index][key]
                            if abs(observation[key]-expected) > 1e-12+1e-6*abs(expected):
                                raise RuntimeError('global transient conservation differs across ranks')
                        record['transient_diagnostics'].append(observation)
            if (ranks > 1 or args.transient) and sum(o['owned_rows'] for o in record['distributions']) != record['distributions'][0]['global_rows']:
                raise RuntimeError('graph row ownership does not cover the operator')
            if failure:
                stderr = (directory/'rank-0/stderr.log').read_text()
                if 'injected bifurcation failure before commit' not in stderr or (directory/'result').exists():
                    raise RuntimeError('graph failure did not preserve unpublished outputs')
            else:
                output = directory/'result'
                steps, edges, ports = table(output/'pressure_flow_steps.csv'), table(output/'pressure_flow_edges.csv'), table(output/'pressure_flow_ports.csv')
                if len(steps) != 2 or len(edges) != 4 or len(ports) != 12:
                    raise RuntimeError('incomplete accepted graph output')
                for index, step in enumerate(steps, 1):
                    if int(step['step']) != index or float(step['time_s']) != 0.01*index or not 1 <= int(step['iterations']) <= 64:
                        raise RuntimeError('invalid graph accepted clock/iterations')
                for edge in edges:
                    if any(not math.isfinite(float(edge[key])) for key in edge if key != 'edge_id'):
                        raise RuntimeError('nonfinite graph interface observation')
                    if abs(float(edge['normalized_flow_residual'])) > 1e-10:
                        raise RuntimeError('graph flow interface conservation gate failed')
                    if args.execution != 'explicit' and abs(float(edge['normalized_pressure_residual'])) > 1e-6:
                        raise RuntimeError('graph pressure convergence gate failed')
                if reference_ports is None:
                    reference_ports = ports
                record['port_error_fraction_of_gate'] = compare_ports(reference_ports, ports)
                record['steps'] = steps
                record['output_sha256'] = {p.name: digest(p) for p in output.iterdir() if p.is_file()}
            record['status'] = 'passed'
            print(args.execution, name, 'passed', flush=True)
        if digest(binary) != result['binary_sha256'] or any(digest(case/p) != h for p, h in inputs.items()):
            raise RuntimeError('graph binary or fixture changed during validation')
        result['status'] = 'passed'
    except BaseException:
        result['status'] = 'failed'
        raise
    finally:
        (root/'summary.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
