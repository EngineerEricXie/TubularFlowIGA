#!/usr/bin/env python3
"""Check serial tool streams against retained binaries and numerical fixtures."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time
from hpc_inventory import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-dir', type=Path, required=True)
    parser.add_argument('--flow-case', type=Path, required=True)
    parser.add_argument('--transport-case', type=Path, required=True)
    parser.add_argument('--transport-history', type=Path, required=True)
    parser.add_argument('--read-preload', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    baseline = args.baseline_dir.resolve()
    current = repo/'solvers/cpu'
    flow, transport, history = (p.resolve() for p in (args.flow_case, args.transport_case, args.transport_history))
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    result = dict(status='running', cases=[], comparisons=[], binaries={}, inputs={})
    for directory in (flow, transport, history):
        for path in directory.iterdir():
            if path.is_file(): result['inputs'][str(path)] = digest(path)

    def run(label, binary, argv, expected=0, settings=None, diagnostic=None, full=False):
        directory = output/label
        directory.mkdir()
        command = [str(binary), *map(str, argv)]
        record = dict(name=label, argv=command, expected=expected, settings=settings or {})
        result['cases'].append(record)
        result['binaries'][str(binary)] = digest(binary)
        started = time.monotonic()
        with open('/dev/full' if full else directory/'stdout.log', 'w') as stdout, (directory/'stderr.log').open('w') as stderr:
            process = subprocess.run(command, env=dict(env, **(settings or {})), stdout=stdout, stderr=stderr, timeout=60)
        record.update(returncode=process.returncode, elapsed_s=time.monotonic()-started)
        if process.returncode != expected: raise RuntimeError(f'{label}: exit {process.returncode}, expected {expected}')
        if diagnostic and diagnostic not in (directory/'stderr.log').read_text(): raise RuntimeError(f'{label}: diagnostic missing')
        record['status'] = 'passed'
        print(label, 'passed', flush=True)
        return directory

    def equal(label, a, b):
        if a.read_bytes() != b.read_bytes(): raise RuntimeError(f'{label}: bytes changed')
        result['comparisons'].append(dict(name=label, exact=True, sha256=digest(a)))

    try:
        wom_args = [flow/'straight_tube-1.ntiga', flow/'controlmesh.vtk', repo/'examples/validation/womersley/womersley_reference.json']
        for name, binaries in (('before', baseline), ('after', current)):
            run('womersley-'+name, binaries/'iga_womersley_reference', wom_args+[output/('wom-'+name)])
        for p in (output/'wom-before').iterdir(): equal('womersley-'+p.name, p, output/'wom-after'/p.name)
        for mode in range(3):
            target = output/f'wom-fault-{mode}'
            run(f'womersley-format-{mode}', current/'womersley_stream_failure_test', wom_args+[target], 1,
                {'IGA_TEST_FORMAT_MODE': str(mode), 'IGA_TEST_FORMAT_PREFIX': 'womersley.step'}, 'serial_formatter injected=1 native_status=1')
            if (target/'womersley.step').exists(): raise RuntimeError('truncated filename escaped')
        close_lib = str(repo/'solvers/coupling/text_close_preload.so')
        for filename in ('womersley.step000001.txt', 'velocity_series.csv'):
            target = output/('close-'+filename)
            run('womersley-close-'+filename, current/'iga_womersley_reference', wom_args+[target], 1,
                {'LD_PRELOAD': close_lib, 'TUBULARFLOWIGA_TEST_TEXT_CLOSE': str(target/filename)}, 'cannot write Womersley')
            run('womersley-close-retry-'+filename, current/'iga_womersley_reference', wom_args+[target])
            for p in (output/'wom-before').iterdir(): equal('close-retry-'+filename+'-'+p.name, p, target/p.name)
        run('womersley-stdout-full', current/'iga_womersley_reference', wom_args+[output/'wom-stdout'], 1, full=True)
        budget_args = [transport/'straight_neurite-1.ntiga', transport, 'neuron_transport',
                       history/'field.step000000.txt', history/'field.step000001.txt', '--rtol', '1e-6', '--zero-atol', '1e-12']
        old = run('budget-before', baseline/'iga_transport_validate', budget_args)
        new = run('budget-after', current/'iga_transport_validate', budget_args)
        equal('budget-json', old/'stdout.log', new/'stdout.log')
        if not json.loads((new/'stdout.log').read_text())['accepted']: raise RuntimeError('budget rejected')
        for mode in range(3):
            run(f'budget-format-{mode}', current/'transport_budget_stream_failure_test', budget_args, 1,
                {'IGA_TEST_FORMAT_MODE': str(mode), 'IGA_TEST_FORMAT_PREFIX': '{'}, 'serial_formatter injected=1 native_status=1')
        run('budget-stdout-full', current/'iga_transport_validate', budget_args, 1, full=True)
        for name, binaries in (('before', baseline), ('after', current)):
            run('pack-'+name, binaries/'iga_pack', [flow, 1, output/(name+'.ntiga')])
        equal('packed-database', output/'before.ntiga', output/'after.ntiga')
        read_settings = {'LD_PRELOAD': str(args.read_preload.resolve()), 'IGA_TEST_READ_PATH': str(flow/'simulation_config.json')}
        run('pack-read-eio', current/'iga_pack', [flow, 1, output/'bad.ntiga'], 1, read_settings, 'cannot read complete text stream')
        run('pack-read-retry', current/'iga_pack', [flow, 1, output/'retry.ntiga'])
        equal('pack-read-retry', output/'before.ntiga', output/'retry.ntiga')
        run('pack-stdout-full', current/'iga_pack', [flow, 1, output/'stdout.ntiga'], 1, full=True)
        for path, sha in result['inputs'].items():
            if digest(Path(path)) != sha: raise RuntimeError('fixture changed: '+path)
        result['status'] = 'passed'
    except BaseException:
        result['status'] = 'failed'
        raise
    finally:
        (output/'summary.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
