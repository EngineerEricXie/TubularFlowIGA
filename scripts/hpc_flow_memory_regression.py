#!/usr/bin/env python3
"""Verify optional CPU flow memory records and numerical neutrality."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-root', type=Path, required=True,
                        help='completed hpc_flow_pvtu_regression output')
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    reference = args.reference_root.resolve()
    binary = repo / 'solvers/cpu/iga_navier_stokes'
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    report = dict(status='running', binary_sha256=sha(binary), jobs=[])
    try:
        for name, ranks, fmt, damage in [(f'{fmt}-{ranks}', ranks, fmt, False)
                for fmt in ('vtkhdf', 'pvtu') for ranks in (1, 2)] + [('blocked-2', 2, 'pvtu', True)]:
            folder = root / name
            folder.mkdir()
            memory = folder / 'memory.jsonl'
            if damage:
                memory.mkdir()
            inputs = reference / f'input-{ranks}'
            command = ['timeout', '--kill-after=5s', '150s', 'mpiexec', '-np', str(ranks),
                       sys.executable, str(repo / 'scripts/hpc_rank_run.py'), '--output-dir', str(folder),
                       '--expected-ranks', str(ranks), '--timeout', '120', '--', str(binary),
                       str(inputs / 'database.ntiga'), str(inputs), '--output', str(folder / 'flow.txt'),
                       '--output-every', '1', '--visualization-format', fmt,
                       '--checkpoint', str(folder / 'checkpoint'), '--memory-report', str(memory)]
            with (folder / 'launcher.log').open('w') as log:
                result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                    env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                    PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -ksp_rtol 1e-12'))
            row = dict(name=name, command=command, returncode=result.returncode)
            report['jobs'].append(row)
            records = [json.loads((folder / f'rank-{r}/run.json').read_text()) for r in range(ranks)]
            assert all(not r['timed_out'] for r in records)
            row['rank_reports'] = records
            if damage:
                assert result.returncode != 0 and all(r['returncode'] != 0 for r in records)
                assert all('navier_stokes_v2 seconds=' not in (folder / f'rank-{r}/stdout.log').read_text() for r in range(ranks))
            else:
                assert result.returncode == 0 and all(r['returncode'] == 0 for r in records)
                assert sha(folder / 'checkpoint.state') == sha(reference / f'parallel-{ranks}/result/checkpoint.state')
                rows = [json.loads(line) for line in memory.read_text().splitlines()]
                stages = [r['stage'] for r in rows]
                assert stages[0] == 'database_open' and stages[-1] == 'flow_closed'
                for stage in ('initialized_state', 'visualization_geometry', 'visualization_ready', 'output_begin', 'flow_trial_solved'):
                    assert stage in stages
                output_stage = 'parallel_piece_extracted' if fmt == 'pvtu' else 'serial_output_released'
                assert [r['step'] for r in rows if r['stage'] == output_stage] == ([0, 1, 2] if fmt == 'pvtu' else [0, 1, 2, 2])
                for r in rows:
                    assert r['mpi_ranks'] == ranks
                    assert len(r['rss_current_bytes_per_rank']) == ranks
                    assert all(v > 0 for v in r['rss_current_bytes_per_rank'])
                row['memory_records'] = rows
            print(name, 'passed', flush=True)
        assert sha(binary) == report['binary_sha256']
        report['status'] = 'passed'
    except BaseException:
        report['status'] = 'failed'
        raise
    finally:
        (root / 'acceptance.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
