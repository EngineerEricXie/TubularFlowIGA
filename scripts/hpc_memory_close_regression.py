#!/usr/bin/env python3
"""Reject memory-report close failures before CPU transport reports success."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

from hpc_compare_fields import compare
from hpc_inventory import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', required=True, type=Path)
    parser.add_argument('--fixture-root', required=True, type=Path,
                        help='Accepted transport CLI directory containing prescribed-one/two inputs')
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    binary = repo/'solvers/cpu/iga_solve'
    preload = repo/'solvers/coupling/text_close_preload.so'
    baseline = args.baseline.resolve()
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', IGA_PROFILE='1')
    for name in ('PETSC_OPTIONS', 'LD_PRELOAD', 'TUBULARFLOWIGA_TEST_TEXT_CLOSE'):
        env.pop(name, None)
    summary = dict(status='running', cases=[], comparisons=[], inputs={}, binaries={
        str(p):digest(p) for p in (binary, baseline, preload)})

    def execute(ranks, label, source, old=False, fault=False, enabled=True):
        directory = root/(str(ranks)+'-'+label)
        directory.mkdir()
        records = directory/'ranks'
        records.mkdir()
        child = [str(baseline if old else binary), str(source/'database.ntiga'), str(source),
                 '--system', 'transport', '--output', str(directory/'field.txt'), '--visualization-format', 'vtu']
        report_path = directory/'memory.jsonl'
        if enabled:
            child += ['--memory-report', str(report_path)]
        if fault:
            child = ['env', 'LD_PRELOAD='+str(preload), 'TUBULARFLOWIGA_TEST_TEXT_CLOSE='+str(report_path), *child]
        command = ['timeout', '--kill-after=5s', '120s', 'mpiexec', '--oversubscribe', '-np', str(ranks),
                   sys.executable, str(repo/'scripts/hpc_sequential_output_regression.py'), '--rank-worker',
                   '--output-dir', str(records), '--', *child]
        expected = int(fault and not old)
        entry = dict(name=directory.name, argv=command, expected=expected, ranks=[])
        summary['cases'].append(entry)
        with (directory/'launcher.log').open('w') as out:
            result = subprocess.run(command, env=env, stdout=out, stderr=subprocess.STDOUT)
        entry['returncode'] = result.returncode
        if result.returncode != expected:
            raise RuntimeError(directory.name+': wrong launcher exit')
        for rank in range(ranks):
            rd = records/f'rank-{rank}'
            record = json.loads((rd/'run.json').read_text())
            entry['ranks'].append(record)
            if record['returncode'] != expected or record['timed_out']:
                raise RuntimeError(directory.name+': wrong rank exit or timeout')
            stdout, stderr = (rd/'stdout.log').read_text(), (rd/'stderr.log').read_text()
            if expected and ('memory report close: rank 0: cannot close memory report' not in stderr
                             or 'iga_solve system=' in stdout):
                raise RuntimeError(directory.name+': missing rejection or premature success')
            if rank == 0:
                if fault and '[text-close-injection] fclose returned failure' not in stderr:
                    raise RuntimeError(directory.name+': close injection missed')
                if not expected and 'iga_solve system=' not in stdout:
                    raise RuntimeError(directory.name+': missing healthy completion')
        if enabled:
            lines = [json.loads(line) for line in report_path.read_text().splitlines()]
            expected_stages = ['database_open', 'required_elements_loaded', 'matrix_preallocation',
                               'operator_assembly', 'ksp_setup', 'linear_solve']
            if [line['stage'] for line in lines] != expected_stages:
                raise RuntimeError(directory.name+': changed report stages')
            if any(line['mpi_ranks'] != ranks or len(line['rss_peak_bytes_per_rank']) != ranks for line in lines):
                raise RuntimeError(directory.name+': incomplete report ranks')
        elif report_path.exists():
            raise RuntimeError(directory.name+': disabled report created')
        entry['status'] = 'passed'
        print(directory.name, 'passed', flush=True)
        return directory

    try:
        for ranks, suffix in ((1, 'one'), (2, 'two')):
            source = args.fixture_root.resolve()/('prescribed-'+suffix)/'input-0'
            for path in source.iterdir():
                if path.is_file():
                    summary['inputs'][str(path)] = digest(path)
            before = execute(ranks, 'before', source, old=True)
            candidates = [execute(ranks, 'old-fault', source, old=True, fault=True),
                          execute(ranks, 'healthy', source), execute(ranks, 'fault', source, fault=True),
                          execute(ranks, 'retry', source), execute(ranks, 'disabled', source, enabled=False)]
            for candidate in candidates:
                result = compare(before/'field.txt', candidate/'field.txt', 1e-6, 1e-12, node_ids=True)
                summary['comparisons'].append(dict(case=candidate.name, result=result))
                if not result['passed']:
                    raise RuntimeError(candidate.name+': field mismatch')
                for name in ('field.txt.fields', 'field.vtu', 'field.pvd', 'physiology_fields.json'):
                    if (before/name).read_bytes() != (candidate/name).read_bytes():
                        raise RuntimeError(candidate.name+': changed auxiliary output '+name)
        for name, sha in {**summary['inputs'], **summary['binaries']}.items():
            if digest(Path(name)) != sha:
                raise RuntimeError('input or executable changed: '+name)
        summary['status'] = 'passed'
    except BaseException:
        summary['status'] = 'failed'
        raise
    finally:
        (root/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
