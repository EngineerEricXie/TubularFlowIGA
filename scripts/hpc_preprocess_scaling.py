#!/usr/bin/env python3
"""Measure mesh, spline, METIS, and packer scaling on generated tube cases."""

import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys

from hpc_inventory import digest


def vtk_counts(path):
    text = path.read_text()
    points = re.search(r'^POINTS\s+(\d+)\s+', text, re.MULTILINE)
    cells = re.search(r'^CELLS\s+(\d+)\s+', text, re.MULTILINE)
    if not points or not cells:
        raise RuntimeError(f'cannot parse VTK counts from {path}')
    return int(points.group(1)), int(cells.group(1))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--threads', type=int, nargs='+', default=[1, 4])
    parser.add_argument('--memory-limit-gib', type=float, default=2.0)
    parser.add_argument('--timeout', type=float, default=900.0)
    args = parser.parse_args()
    if not args.threads or any(value < 1 for value in args.threads):
        parser.error('threads must be positive')
    if args.memory_limit_gib <= 0 or args.timeout <= 0:
        parser.error('memory limit and timeout must be positive')

    repo = Path(__file__).resolve().parents[1]
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    tools = {
        'mesh': repo/'preprocessing/mesh/tubular_mesh',
        'spline': repo/'preprocessing/spline/spline',
        'partition': Path('/usr/bin/mpmetis'),
        'packer': repo/'solvers/cpu/iga_pack',
        'inspect': repo/'solvers/cpu/iga_inspect',
    }
    for name, path in tools.items():
        if not path.is_file() or not os.access(path, os.X_OK):
            parser.error(f'missing executable {name}: {path}')
    base_configuration = json.loads(
        (repo/'examples/vascular_flow/straight_tube/simulation_config.json').read_text())
    cases = [('small', 8.0), ('medium', 32.0), ('large', 96.0)]
    report = {
        'status': 'running',
        'cases': [],
        'threads': args.threads,
        'memory_limit_bytes': int(args.memory_limit_gib*1024**3),
        'tool_sha256': {name: digest(path) for name, path in tools.items()},
    }

    def run(stage, command, environment=None):
        run_dir = output/'runs'/stage
        run_dir.mkdir(parents=True)
        completed = subprocess.run(
            [sys.executable, str(repo/'scripts/hpc_rank_run.py'),
             '--output-dir', str(run_dir), '--timeout', str(args.timeout),
             '--', *map(str, command)], cwd=repo, env=environment)
        row = json.loads((run_dir/'rank-0/run.json').read_text())
        if completed.returncode != 0 or row['returncode'] != 0 \
                or row['timed_out'] or not row['resource']:
            raise RuntimeError(f'{stage} failed')
        if row['resource']['peak_rss_bytes'] > report['memory_limit_bytes']:
            raise RuntimeError(f'{stage} exceeded the memory limit')
        return {
            'wall_s': row['wall_s'],
            'user_cpu_s': row['resource']['user_cpu_s'],
            'system_cpu_s': row['resource']['system_cpu_s'],
            'peak_rss_bytes': row['resource']['peak_rss_bytes'],
            'stdout_sha256': row['logs']['stdout.log'],
            'stderr_sha256': row['logs']['stderr.log'],
        }

    try:
        previous_elements = 0
        for case_name, length in cases:
            case = output/'cases'/case_name
            case.mkdir(parents=True)
            configuration = json.loads(json.dumps(base_configuration))
            configuration['mesh']['centerline']['target_spacing'] = 0.5
            configuration['mesh']['centerline']['max_spacing_over_diameter'] = 0.5
            (case/'simulation_config.json').write_text(json.dumps(configuration, indent=2)+'\n')
            (case/'skeleton_initial.swc').write_text(
                f'1 2 0 0 0 0.5 -1\n2 2 {length/2:g} 0 0 0.5 1\n'
                f'3 2 {length:g} 0 0 0.5 2\n')
            inputs = {name: digest(case/name)
                      for name in ('simulation_config.json', 'skeleton_initial.swc')}
            mesh = run(f'{case_name}-mesh', [tools['mesh'], 'pipeline', case,
                       repo/'meshgeneration/template'])
            nodes, elements = vtk_counts(case/'controlmesh.vtk')
            if elements <= previous_elements:
                raise RuntimeError('generated mesh sizes are not increasing')
            previous_elements = elements
            quality = json.loads((case/'mesh_quality.json').read_text())
            if quality['bad_elements'] != 0 or quality['surface_intersections'] != 0 \
                    or quality['minimum_scaled_jacobian'] < 0.1:
                raise RuntimeError(f'{case_name} failed mesh quality')

            spline_runs = []
            expected_outputs = None
            for threads in args.threads:
                environment = dict(os.environ, OMP_NUM_THREADS=str(threads),
                                   OMP_DYNAMIC='FALSE', OMP_PROC_BIND='close',
                                   OMP_PLACES='cores')
                measured = run(f'{case_name}-spline-t{threads}',
                               [tools['spline'], str(case)+'/', '--no-legacy-text'], environment)
                files = {name: digest(case/name) for name in
                         ('bzmeshinfo.txt', 'spline_cache.igacache', 'bzmesh.vtk',
                          'geometry_transform.json')}
                if expected_outputs is None:
                    expected_outputs = files
                elif files != expected_outputs:
                    raise RuntimeError(f'{case_name} spline outputs differ by thread count')
                spline_runs.append({'threads': threads, **measured})
            if int((case/'bzmeshinfo.txt').read_text().split(None, 1)[0]) != elements:
                raise RuntimeError(f'{case_name} spline element count differs from control mesh')

            partition = run(f'{case_name}-partition',
                            [tools['partition'], case/'bzmeshinfo.txt', '4'])
            database = case/'database.ntiga'
            packer = run(f'{case_name}-packer', [tools['packer'], case, '4', database])
            inspection = run(f'{case_name}-inspect', [tools['inspect'], database])
            legacy = None
            if case_name == 'small':
                threads = max(args.threads)
                environment = dict(os.environ, OMP_NUM_THREADS=str(threads),
                                   OMP_DYNAMIC='FALSE', OMP_PROC_BIND='close',
                                   OMP_PLACES='cores')
                legacy_spline = run(f'{case_name}-spline-legacy-t{threads}',
                                    [tools['spline'], str(case)+'/'], environment)
                common = {name: digest(case/name) for name in
                          ('bzmeshinfo.txt', 'spline_cache.igacache', 'bzmesh.vtk',
                           'geometry_transform.json')}
                if common != expected_outputs:
                    raise RuntimeError('legacy-enabled spline changed common outputs')
                legacy_database = case/'database-legacy.ntiga'
                legacy_packer = run(f'{case_name}-packer-legacy',
                                    [tools['packer'], case, '4', legacy_database,
                                     '--legacy-text'])
                if digest(legacy_database) != digest(database):
                    raise RuntimeError('cache and legacy text packers differ')
                legacy = {
                    'spline': {'threads': threads, **legacy_spline},
                    'packer': legacy_packer,
                    'database_sha256': digest(database),
                    'cmat_bytes': (case/'cmat.txt').stat().st_size,
                    'bzpt_bytes': (case/'bzpt.txt').stat().st_size,
                }
            if any(digest(case/name) != expected for name, expected in inputs.items()):
                raise RuntimeError(f'{case_name} inputs changed')
            row = {
                'name': case_name,
                'length': length,
                'nodes': nodes,
                'elements': elements,
                'mesh': mesh,
                'spline': spline_runs,
                'partition': partition,
                'packer': packer,
                'inspect': inspection,
                'legacy_compatibility': legacy,
                'files': {name: {'bytes': (case/name).stat().st_size,
                                 'sha256': digest(case/name)} for name in
                          ('controlmesh.vtk', 'bzmeshinfo.txt', 'spline_cache.igacache',
                           'bzmesh.vtk', 'database.ntiga')},
                'input_sha256': inputs,
            }
            report['cases'].append(row)
            (output/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
            print(f"{case_name}: {nodes} nodes, {elements} elements passed", flush=True)

        large = report['cases'][-1]
        spline_by_threads = {row['threads']: row for row in large['spline']}
        one = spline_by_threads.get(1)
        most = spline_by_threads[max(args.threads)]
        decision = {
            'mesh_parallelization': 'not justified by measured workstation cost',
            'spline_parallelization': 'retain existing chunked OpenMP implementation',
            'packer_parallelization': 'not justified by measured workstation cost',
            'maximum_measured_peak_rss_bytes': max(
                stage['peak_rss_bytes']
                for case in report['cases']
                for stage in ([case['mesh'], case['partition'], case['packer'], case['inspect'],
                               *case['spline']]
                              + ([case['legacy_compatibility']['spline'],
                                  case['legacy_compatibility']['packer']]
                                 if case['legacy_compatibility'] else []))),
            'large_spline_speedup': (one['wall_s']/most['wall_s']
                                     if one and most['wall_s'] > 0 else None),
        }
        report.update(status='passed', decision=decision,
                      harness_sha256=digest(Path(__file__).resolve()))
    except BaseException:
        report['status'] = 'failed'
        raise
    finally:
        (output/'summary.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
