#!/usr/bin/env python3
"""Exercise CUDA CLI stdout boundaries and compare retained CPU/GPU fields."""
import argparse
import json
import os
import re
from pathlib import Path
from hpc_compare_fields import compare
from hpc_inventory import digest
from hpc_rank_run import run as run_rank


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-dir', required=True, type=Path)
    parser.add_argument('--fixture-root', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    fixtures = args.fixture_root.resolve()
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                      OMPI_COMM_WORLD_RANK='0', OMPI_COMM_WORLD_SIZE='1', IGA_PROFILE='1')
    summary = dict(status='running', cases=[], comparisons=[], binaries={}, inputs={}, cpu_references={})
    refs = {}
    for system in ('flow', 'transport'):
        directory = fixtures/'hpc00/matrix'/system
        metadata = json.loads((directory/'matrix.json').read_text())
        results = json.loads((directory/'summary.json').read_text())
        if results['status'] != 'field_comparisons_passed': raise RuntimeError('unaccepted CPU reference')
        observation, = [r for r in results['observations'] if r['ranks'] == 1 and r['repetition'] == 0]
        for name, field in observation['fields'].items():
            path = Path(observation['directory'])/name
            if digest(path) != field['candidate_sha256']: raise RuntimeError('changed CPU field')
        case = Path(metadata['case_dir'])
        for name, sha in metadata['input_and_binary_sha256'].items():
            path = Path(name)
            if path.is_relative_to(case) and digest(path) != sha: raise RuntimeError('changed CPU input')
        refs[system] = Path(observation['directory'])
        summary['cpu_references'][system] = dict(matrix_sha256=digest(directory/'matrix.json'),
            summary_sha256=digest(directory/'summary.json'), observation=observation,
            scope='Retained accepted CPU fields and unchanged case inputs; CPU binary not rerun')
    flow = fixtures/'hpc00/matrix/flow-inputs'
    transport = fixtures/'hpc00/matrix/transport-inputs'
    legacy = fixtures/'hpc01/tool-assets/tools/fixture'
    flowargs = ['navier-stokes', flow/'straight_tube-1.ntiga', flow, '--max-newton', '30',
                '--nonlinear-rtol', '1e-8', '--nonlinear-atol', '1e-12', '--mass-rtol', '1e-6']
    solveargs = ['solve', transport/'straight_neurite-1.ntiga', transport, '--system', 'neuron_transport']
    specs = [
        ('device', ['device-info'], 'device=', None, None),
        ('mesh', ['mesh-check', flow/'straight_tube-1.ntiga'], 'geometry min_detJ=', None, None),
        ('flow', flowargs, 'navier_stokes_cuda nodes=', 'flow', '--output'),
        ('solve', solveargs, 'iga_cuda_solve system=', 'transport', '--output'),
        ('legacy', ['transport', legacy/'serial.ntiga', legacy, '2'], 'transport_cuda nodes=', None, 'positional'),
        ('allocation', solveargs, 'cuda_allocations scope=', 'transport', '--output'),
        ('profile', solveargs, 'hpc_profile {', 'transport', '--output'),
    ]
    for case in (flow, transport, legacy):
        for path in case.rglob('*'):
            if path.is_file(): summary['inputs'][str(path)] = digest(path)

    def execute(spec, label, old=False, mode=None):
        name, argv, target, system, outputflag = spec
        directory = root/(name+'-'+label)
        directory.mkdir()
        binary = (args.baseline_dir.resolve() if old else repo/'solvers/cuda')/('cuda_stdout_failure_test' if mode else 'iga_cuda')
        command = [str(binary), *map(str, argv)]
        if outputflag:
            if outputflag != 'positional': command.append(outputflag)
            command.append(str(directory/'field.txt'))
        if mode: os.environ.update(IGA_TEST_STDOUT_TARGET=target, IGA_TEST_STDOUT_MODE=mode)
        else:
            os.environ.pop('IGA_TEST_STDOUT_TARGET', None)
            os.environ.pop('IGA_TEST_STDOUT_MODE', None)
        summary['binaries'][str(binary)] = digest(binary)
        code = run_rank(command, directory, 1, 120)
        report = json.loads((directory/'rank-0/run.json').read_text())
        record = dict(name=directory.name, target=target, mode=mode, report=report)
        summary['cases'].append(record)
        expected = int(mode is not None)
        if code != expected or report['timed_out'] or report['resource'] is None:
            raise RuntimeError(directory.name+': incorrect exit or resources')
        stdout = (directory/'rank-0/stdout.log').read_text()
        stderr = (directory/'rank-0/stderr.log').read_text()
        if mode:
            if ('stdout_test rank=0 injected=1 native_status=1' not in stderr or
                    'cannot write complete text stream' not in stderr):
                raise RuntimeError(directory.name+': missing native failure')
        elif target not in stdout: raise RuntimeError(directory.name+': missing summary')
        allocation = re.search(r'cuda_allocations scope=project_device_buffers requested_peak_bytes=(\d+) requested_live_bytes=(\d+)',
                               stderr if mode else stdout)
        if outputflag:
            if not allocation or int(allocation[1]) <= 0 or int(allocation[2]) != 0:
                raise RuntimeError(directory.name+': device allocation leak or missing measurement')
            record['requested_peak_bytes'] = int(allocation[1])
        if not mode and system:
            for name in (['field.txt', 'field.txt.pressure'] if system == 'flow' else ['field.txt']):
                comparison = compare(refs[system]/name, directory/name, 1e-5, 1e-12, node_ids=system == 'transport')
                summary['comparisons'].append(dict(case=directory.name, reference='CPU', file=name, result=comparison))
                if not comparison['passed']: raise RuntimeError(directory.name+': CPU field mismatch')
        record['status'] = 'passed'
        print(directory.name, 'passed', flush=True)
        return directory

    try:
        for spec in specs:
            before = execute(spec, 'before', old=True)
            healthy = execute(spec, 'healthy')
            modes = ('short', 'flush', 'exception', 'allocation', 'nonstandard') if spec[0] == 'device' else ('short', 'flush')
            for mode in modes: execute(spec, mode, mode=mode)
            retry = execute(spec, 'retry')
            if spec[4]:
                names = ['field.txt', 'field.txt.pressure'] if spec[0] == 'flow' else ['field.txt']
                for candidate in (healthy, retry):
                    for name in names:
                        result = compare(before/name, candidate/name, 1e-6, 1e-12, node_ids=spec[3] == 'transport')
                        summary['comparisons'].append(dict(case=candidate.name, reference='GPU before', file=name, result=result))
                        if not result['passed']: raise RuntimeError(candidate.name+': GPU field mismatch')
        for name, sha in summary['inputs'].items():
            if digest(Path(name)) != sha: raise RuntimeError('input changed: '+name)
        summary['status'] = 'passed'
    except BaseException:
        summary['status'] = 'failed'
        raise
    finally:
        (root/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
