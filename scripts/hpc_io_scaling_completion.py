#!/usr/bin/env python3
"""Consolidate the HPC-06A/B/C memory, PVTU, and I/O acceptance evidence."""

import argparse
import json
from pathlib import Path

from hpc_inventory import digest


def read_json(path):
    if not path.is_file():
        raise RuntimeError(f'missing evidence {path}')
    return json.loads(path.read_text())


def require_rank_success(rows, label):
    if not rows or any(row['returncode'] != 0 or row['timed_out']
                       or not row.get('resource') for row in rows):
        raise RuntimeError(f'{label} rank reports are incomplete')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence-root', type=Path, default=Path('artifacts/benchmarks/hpc06'))
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    evidence = args.evidence_root.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)

    paths = {
        'partitioned_writer': evidence/'partitioned-v1/audit.json',
        'collective_publisher': evidence/'parallel-v2/audit.json',
        'point_identity': evidence/'identity-v2/audit.json',
        'point_values': evidence/'point-values-v2/audit.json',
        'parallel_bezier': evidence/'parallel-bezier-v3/audit.json',
        'petsc_bezier': evidence/'petsc-bezier-v1/audit.json',
        'flow_cli': evidence/'flow-pvtu-v1/audit.json',
        'flow_memory': evidence/'flow-memory-v2/acceptance.json',
        'duct_128': evidence/'duct-output-v1/audit.json',
        'duct_1024': evidence/'duct-output-large-v1/audit.json',
        'io_control': evidence/'io-control-v4/summary.json',
        'transport_reader': evidence/'io-control-v4/paraview.json',
    }
    documents = {name: read_json(path) for name, path in paths.items()}

    if documents['partitioned_writer']['results']['cpp_exit'] != 0 \
            or documents['partitioned_writer']['results']['paraview_final_exit'] != 0 \
            or documents['partitioned_writer']['paraview']['status'] != 'passed':
        raise RuntimeError('partitioned VTK format evidence failed')
    for name in ('collective_publisher', 'point_identity', 'point_values',
                 'parallel_bezier', 'petsc_bezier'):
        require_rank_success(documents[name]['rank_reports'], name)
    if documents['collective_publisher']['reader']['status'] != 'passed' \
            or documents['petsc_bezier']['reader']['status'] != 'passed':
        raise RuntimeError('distributed PVTU reader evidence failed')
    if any(value != 0 for key, value in documents['parallel_bezier']['results'].items()
           if key != 'legacy_hdf_geometry_error'):
        raise RuntimeError('parallel Bezier regression failed')
    if documents['flow_cli']['status'] != 'passed':
        raise RuntimeError('flow PVTU CLI evidence failed')

    memory = documents['flow_memory']
    if memory['status'] != 'passed' or len(memory['jobs']) != 5:
        raise RuntimeError('flow memory evidence is incomplete')
    successful_memory_jobs = [row for row in memory['jobs'] if row['name'] != 'blocked-2']
    for row in successful_memory_jobs:
        require_rank_success(row['rank_reports'], row['name'])
        if not row['memory_records']:
            raise RuntimeError(f"{row['name']} has no phase memory records")
    blocked_memory = next(row for row in memory['jobs'] if row['name'] == 'blocked-2')
    if blocked_memory['returncode'] == 0 or any(row['returncode'] == 0 or row['timed_out']
                                                for row in blocked_memory['rank_reports']):
        raise RuntimeError('memory-report failure did not reject collectively')

    duct_128 = documents['duct_128']
    duct_1024 = documents['duct_1024']
    if duct_128['status'] != 'passed_original_L2_gates' \
            or duct_128['reader']['status'] != 'passed' \
            or duct_128['same_state_reader']['status'] != 'passed':
        raise RuntimeError('128-element PVTU comparison failed')
    if duct_1024['status'] != 'passed' or duct_1024['reader']['status'] != 'passed' \
            or duct_1024['case'] != {'elements': 1024, 'nodes': 2299, 'ranks': 4}:
        raise RuntimeError('1024-element PVTU comparison failed')
    for label, stats in duct_1024['stats'].items():
        if len(stats['rss_peak']) != 4 or any(value <= 0 for value in stats['rss_peak']):
            raise RuntimeError(f'{label} does not report every rank RSS')
        if stats['files'] <= 0 or stats['bytes'] <= 0:
            raise RuntimeError(f'{label} output cost is incomplete')

    controls = documents['io_control']
    reader = documents['transport_reader']
    if controls['status'] != 'passed' or len(controls['cases']) != 10:
        raise RuntimeError('I/O frequency matrix is incomplete')
    positives = [row for row in controls['cases']
                 if row.get('expected') != 'collective rejection']
    rejection = next(row for row in controls['cases']
                     if row.get('expected') == 'collective rejection')
    if any(row['status'] != 'passed' or row['diagnostic_records'] not in (1, 2)
           or row['checkpoint_records'] not in (1, 2) for row in positives):
        raise RuntimeError('positive I/O frequency case failed')
    if rejection['status'] != 'passed' or rejection['returncode'] == 0 \
            or rejection['final_index_published']:
        raise RuntimeError('transport PVTU failure case did not preserve publication rules')
    for label, comparison in controls['comparisons'].items():
        if comparison['dense_files'] <= comparison['sparse_files'] \
                or comparison['dense_bytes'] <= comparison['sparse_bytes']:
            raise RuntimeError(f'{label} sparse output did not reduce I/O')
    if reader['status'] != 'passed' or reader['frames'] != 10 \
            or max(row['relative_l2'] for row in reader['comparisons']) > reader['relative_tolerance']:
        raise RuntimeError('transport ParaView comparison failed')

    sources = [repo/'solvers/cpu/src/iga_navier_stokes.cpp',
               repo/'solvers/cpu/src/iga_solve.cpp',
               repo/'solvers/cpu/Makefile',
               repo/'scripts/hpc_io_control_regression.py',
               repo/'scripts/test_transport_pvtu_paraview.py',
               Path(__file__).resolve()]
    report = {
        'status': 'passed',
        'tasks': ['HPC-06A', 'HPC-06B', 'HPC-06C'],
        'evidence_sha256': {str(path.relative_to(repo)): digest(path)
                            for path in paths.values()},
        'source_sha256': {str(path.relative_to(repo)): digest(path) for path in sources},
        'large_case': {
            **duct_1024['case'],
            'reader_frames': duct_1024['reader']['frames'],
            'maximum_field_relative_l2': max(
                row['relative_l2'] or 0.0 for row in duct_1024['reader']['comparisons']),
            'output_statistics': duct_1024['stats'],
        },
        'io_control': {
            'jobs': len(controls['cases']),
            'positive_jobs': len(positives),
            'expected_rejections': 1,
            'comparisons': controls['comparisons'],
            'transport_reader_frames': reader['frames'],
            'transport_maximum_relative_l2': max(
                row['relative_l2'] or 0.0 for row in reader['comparisons']),
        },
        'decision': 'retain rank-local PVTU pieces; no I/O aggregator is justified by local evidence',
        'scope': 'single-workstation validation; cross-node behavior remains HPC-09',
    }
    (output/'acceptance.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
