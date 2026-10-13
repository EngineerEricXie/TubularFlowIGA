#!/usr/bin/env python3
"""Validate saved moving-graph rank reports and compare accepted port histories."""
import argparse
import json
import math
from pathlib import Path
import sys

from hpc_immersed_graph_regression import compare_ports, table
from hpc_inventory import digest


def finite(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError('nonfinite or nonnumeric '+name)
    return value


def records(path, prefix):
    return [json.loads(line[len(prefix):]) for line in path.read_text().splitlines()
            if line.startswith(prefix)]


def validate_run(directory, case):
    config = json.loads((case/'simulation_config.json').read_text())
    steps, dt = config['time']['steps'], config['time']['dt']
    limits = json.loads((case/'immersed/immersed_geometry.json').read_text())['prescribed_motion']['conservation_limits']
    reports = sorted(directory.glob('rank-*/run.json'))
    if not reports:
        raise ValueError('missing rank reports: '+str(directory))
    ranks = json.loads(reports[0].read_text())['ranks']
    if type(ranks) is not int or ranks < 1 or len(reports) != ranks:
        raise ValueError('incomplete rank report coverage')
    seen, distributions, accepted = set(), [], []
    for path in reports:
        report = json.loads(path.read_text())
        rank = report['rank']
        if (type(rank) is not int or rank in seen or not 0 <= rank < ranks
                or report['ranks'] != ranks or path.parent.name != f'rank-{rank}'
                or report['returncode'] != 0 or report['timed_out'] or report['status'] != 'process_passed'):
            raise ValueError('invalid or unsuccessful rank report: '+str(path))
        seen.add(rank)
        for name in ('stdout.log', 'stderr.log'):
            if digest(path.parent/name) != report['logs'][name]:
                raise ValueError('rank log digest mismatch: '+str(path.parent/name))
        command = report['command_argv']
        supplied = Path(command[command.index('--graph-case')+1])
        if not supplied.is_absolute():
            supplied = Path(report['cwd'])/supplied
        if supplied.resolve() != case.resolve():
            raise ValueError('rank graph input path differs from requested case')
        ownership = records(path.parent/'stdout.log', 'hpc_immersed_distribution ')
        if len(ownership) != 1:
            raise ValueError('missing or repeated ownership report')
        ownership = ownership[0]
        if (ownership['rank'] != rank or ownership['ranks'] != ranks or ownership['domain'] != 'immersed'
                or ownership['time_integration'] != 'backward_euler'
                or not 0 <= ownership['owned_rows'] <= ownership['global_rows']
                or not 0 <= ownership['required_rows'] <= ownership['global_rows']):
            raise ValueError('invalid ownership report')
        distributions.append(ownership)
        history = records(path.parent/'stdout.log', 'hpc_immersed_moving_step ')
        if len(history) != steps:
            raise ValueError('incorrect accepted moving step count')
        clock = 0.
        for index, step in enumerate(history, 1):
            clock += dt
            if step['rank'] != rank or step['ranks'] != ranks or step['domain'] != 'immersed' or step['index'] != index or step['time_s'] != clock:
                raise ValueError('accepted moving clock or rank differs')
            for key in ('residual_norm', 'assembly_s', 'solve_s', 'newton_iterations'):
                if finite(step[key], key) < 0:
                    raise ValueError('negative moving diagnostic '+key)
            damping = finite(step['minimum_accepted_damping'], 'minimum_accepted_damping')
            if not 0 < damping <= 1:
                raise ValueError('invalid accepted Newton damping')
            for field in ('reynolds', 'moving_mass', 'wall_relative_leakage'):
                if not 0 <= finite(step[field], field) <= limits[field]:
                    raise ValueError('moving physical gate failed: '+field)
        accepted.append(history)
    if len({row['global_rows'] for row in distributions}) != 1 or sum(row['owned_rows'] for row in distributions) != distributions[0]['global_rows']:
        raise ValueError('owned rows do not cover global operator')
    graph_steps = table(directory/'result/pressure_flow_steps.csv')
    edges = table(directory/'result/pressure_flow_edges.csv')
    if len(graph_steps) != steps or len(edges) != steps*len(config['couplings']):
        raise ValueError('incomplete accepted graph history')
    expected_edges = {(index, edge['id']) for index in range(1, steps+1) for edge in config['couplings']}
    observed_edges = {(int(row['step']), row['edge_id']) for row in edges}
    if observed_edges != expected_edges or len(observed_edges) != len(edges):
        raise ValueError('duplicate or missing accepted edge')
    clock, times = 0., {}
    for index, row in enumerate(graph_steps, 1):
        clock += dt
        times[index] = clock
        if int(row['step']) != index or float(row['time_s']) != clock or not 1 <= int(row['iterations']) <= config['execution']['maximum_iterations']:
            raise ValueError('invalid graph step history')
    for row in edges:
        if float(row['time_s']) != times[int(row['step'])]:
            raise ValueError('accepted edge time differs')
        for field, gate in [('normalized_pressure_residual', 'pressure_relative_tolerance'),
                            ('normalized_flow_residual', 'flow_relative_tolerance')]:
            value = float(row[field])
            if not math.isfinite(value) or not 0 <= value <= config['execution'][gate]:
                raise ValueError('accepted graph edge gate failed')
    ports = table(directory/'result/pressure_flow_ports.csv')
    expected_ports = {(str(index), domain['id'], port['id']) for index in range(1, steps+1)
                      for domain in config['domains'] for port in domain['ports']}
    if {(row['step'], row['domain_id'], row['port_id']) for row in ports} != expected_ports or len(ports) != len(expected_ports):
        raise ValueError('incomplete or duplicate accepted ports')
    for row in ports:
        if float(row['time_s']) != times[int(row['step'])] or not float(row['area_m2']) > 0:
            raise ValueError('invalid accepted port time or area')
    compare_ports(ports, ports)  # Reject nonfinite physical values even without a reference run.
    return dict(directory=str(directory), ranks=ranks, steps=steps, accepted=accepted,
                rank_reports_sha256={str(path): digest(path) for path in reports},
                result_files_sha256={str(path): digest(path) for path in sorted((directory/'result').iterdir()) if path.is_file()}), ports


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case-dir', required=True, type=Path)
    parser.add_argument('--run', required=True, action='append', type=Path)
    parser.add_argument('--reference-run', type=Path, help='optional independently validated port reference')
    args = parser.parse_args()
    try:
        reference, reference_record = None, None
        if args.reference_run:
            reference_record, reference = validate_run(args.reference_run, args.case_dir)
        output = dict(status='passed', reference_ranks=reference_record['ranks'] if reference_record else None, runs=[],
                      limitations=['Validates saved diagnostics and accepted ports; does not certify full solution fields, source provenance, all physical gates, or performance speedup.'])
        for run in args.run:
            record, ports = validate_run(run, args.case_dir)
            if reference is not None:
                record['port_error_over_tolerance'] = compare_ports(reference, ports)
            output['runs'].append(record)
        print(json.dumps(output, indent=2))
        return 0
    except (OSError, ValueError, KeyError, IndexError, TypeError, RuntimeError) as error:
        print(json.dumps(dict(status='failed', error=str(error))))
        return 1


if __name__ == '__main__':
    sys.exit(main())
