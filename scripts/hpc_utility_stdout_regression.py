#!/usr/bin/env python3
"""Compare native inspection tools and reject successful writes to /dev/full."""
import argparse
import json
import subprocess
from pathlib import Path
from hpc_inventory import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--fixture-root', type=Path, required=True,
                        help='Repository outputs directory containing retained HPC fixtures')
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    fixtures = args.fixture_root.resolve()
    tube = fixtures/'hpc00/matrix/flow-inputs'
    db = tube/'straight_tube-1.ntiga'
    legacy = fixtures/'hpc01/tool-assets/tools/fixture'
    flow = fixtures/'hpc01/solver-stdout/native-driven'
    result = flow/'flow-1-normal-before/result'
    field = result/'field.txt'
    manifest = result/'field.series.csv'
    flowdb = flow/'flow-input/serial.ntiga'
    wom = fixtures/'hpc01/serial-tools-audit/native/wom-after'
    graph = fixtures/'hpc01/registry/native-final/single'
    specs = [
        ('inspect', 'iga_inspect', [db]),
        ('legacy-case', 'iga_case_check', [legacy/'serial.ntiga', legacy]),
        ('configured-case', 'iga_case_check', [db, tube]),
        ('config-3d', 'iga_config_check', [tube/'simulation_config.json']),
        ('config-1d', 'iga_config_check', [repo/'examples/one_d/rigid_straight/simulation_config.json']),
        ('config-flow-graph', 'iga_config_check', [graph/'flow/simulation_config.json']),
        ('config-species-graph', 'iga_config_check', [graph/'species/simulation_config.json']),
        ('field', 'iga_flow_validate', [flowdb, field]),
        ('compare', 'iga_flow_validate', [flowdb, '--compare', field, field]),
        ('manifest', 'iga_flow_validate', [flowdb, '--manifest', result, manifest]),
        ('compare-manifests', 'iga_flow_validate', [flowdb, '--compare-manifests', result, manifest, result, manifest]),
        ('womersley', 'iga_flow_validate', [db, '--womersley',
            repo/'examples/validation/womersley/womersley_reference.json', wom, wom/'velocity_series.csv']),
    ]
    summary = dict(status='running', cases=[], binaries={}, inputs={})
    for _, _, argv in specs:
        for value in argv:
            if isinstance(value, Path):
                paths = value.rglob('*') if value.is_dir() else [value]
                for path in paths:
                    if path.is_file(): summary['inputs'][str(path)] = digest(path)
    try:
        for name, tool, argv in specs:
            reference = None
            for version, folder in (('before', args.baseline_dir.resolve()), ('after', repo/'solvers/cpu')):
                binary = folder/tool
                summary['binaries'][str(binary)] = digest(binary)
                for broken in (False, True):
                    label = name+'-'+version+('-full' if broken else '-healthy')
                    command = [str(binary), *map(str, argv)]
                    stdout = Path('/dev/full') if broken else root/(label+'.stdout')
                    with stdout.open('wb') as out, (root/(label+'.stderr')).open('wb') as err:
                        process = subprocess.run(command, stdout=out, stderr=err, timeout=60)
                    expected = int(broken and version == 'after')
                    record = dict(name=label, argv=command, returncode=process.returncode, expected=expected)
                    summary['cases'].append(record)
                    if process.returncode != expected:
                        raise RuntimeError(f'{label}: exit {process.returncode}, expected {expected}')
                    if not broken:
                        data = stdout.read_bytes()
                        if not data: raise RuntimeError(label+': empty summary')
                        if reference is None: reference = data
                        elif data != reference: raise RuntimeError(label+': summary changed')
                    elif version == 'after':
                        if b'cannot write complete text stream' not in (root/(label+'.stderr')).read_bytes():
                            raise RuntimeError(label+': wrong failure')
                    record['status'] = 'passed'
                    print(label, 'passed', flush=True)
        for path, sha in summary['inputs'].items():
            if digest(Path(path)) != sha: raise RuntimeError('input changed: '+path)
        summary['status'] = 'passed'
    except BaseException:
        summary['status'] = 'failed'
        raise
    finally:
        (root/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
