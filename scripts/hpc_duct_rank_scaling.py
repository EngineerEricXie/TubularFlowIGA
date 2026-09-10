#!/usr/bin/env python3
"""Run sequential fixed-mesh duct repetitions and compare fields across ranks."""
import argparse
import json
from pathlib import Path
import statistics
import subprocess
import sys

from hpc_inventory import digest
from hpc_compare_checkpoint_fields import compare_bundles


def main():
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--transverse', type=int, default=4)
    parser.add_argument('--axial', type=int, default=8)
    parser.add_argument('--ranks', type=int, nargs='+', default=[1, 2, 4])
    parser.add_argument('--repeats', type=int, default=2)
    parser.add_argument('--timeout', type=int, default=2400)
    parser.add_argument('--launcher', default='mpiexec --map-by core --bind-to core')
    args = parser.parse_args()
    ranks = sorted(set(args.ranks))
    if min(ranks + [args.transverse, args.axial, args.repeats, args.timeout]) < 1 or ranks[0] != 1:
        parser.error('positive parameters and a one-rank reference are required')
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    binary = args.binary.resolve()
    files = [binary, repo/'solvers/coupling/hpc_duct_solver_fixture'] + [repo/'scripts'/name for name in
        ('hpc_duct_rank_scaling.py', 'hpc_duct_solver_candidates.py', 'hpc_rank_run.py',
         'hpc_compare_checkpoint_fields.py', 'hpc_solver_prefixes.py', 'hpc_inventory.py',
         'hpc_native_graph_checkpoint.py')]
    report = dict(status='running', files={str(p): digest(p) for p in files}, runs=[],
                  parameters=vars(args) | {'binary': str(binary), 'output_dir': str(root)},
                  source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip(),
                  limits=['Local sequential experiments do not imply exclusive workstation access.',
                          'Maximum rank phase times cannot be added as one rank wall time.',
                          'Rank peak RSS sums are not simultaneous aggregate peaks.'])
    def save():
        (root/'acceptance.json').write_text(json.dumps(report, indent=2)+'\n')
    try:
        save()
        for repeat in range(args.repeats):
            for count in ranks:
                folder = root/f'repeat-{repeat}-ranks-{count}'
                command = [sys.executable, str(repo/'scripts/hpc_duct_solver_candidates.py'),
                           '--binary', str(binary), '--output-dir', str(folder),
                           '--transverse', str(args.transverse), '--axial', str(args.axial),
                           '--ranks', str(count), '--timeout', str(args.timeout),
                           '--launcher', args.launcher, '--candidates', 'bjacobi']
                with (root/f'repeat-{repeat}-ranks-{count}.log').open('w') as log:
                    result = subprocess.run(command, cwd=repo, stdout=log, stderr=subprocess.STDOUT)
                data = json.loads((folder/'acceptance.json').read_text())
                row = dict(repeat=repeat, ranks=count, command=command, returncode=result.returncode,
                           acceptance_sha256=digest(folder/'acceptance.json'), candidates=[])
                report['runs'].append(row)
                if result.returncode or data['status'] != 'passed':
                    raise RuntimeError('candidate sweep failed: '+str(folder))
                for candidate in data['candidates']:
                    name = candidate['name']
                    comparison = compare_bundles(root/'repeat-0-ranks-1'/name/'case/bundle',
                                                 folder/name/'case/bundle', 'root3d')
                    row['candidates'].append(dict(name=name, fields=comparison,
                        iterations=candidate['total_linear_iterations'],
                        wall_s=max(r['wall_s'] for r in candidate['rank_reports']),
                        rss_bytes=[r['resource']['peak_rss_bytes'] for r in candidate['rank_reports']],
                        affinity=[r['affinity_cpus'] for r in candidate['rank_reports']],
                        phases={phase: max(r['phases'][phase]['exclusive_s'] for r in candidate['profiles'])
                                for phase in ('assembly', 'solver_setup', 'linear_solve', 'communication', 'output')}))
                save()
                print(folder.name, 'passed', flush=True)
        summary = []
        for name in ('lu', 'bjacobi'):
            reference = [c['wall_s'] for row in report['runs'] if row['ranks'] == 1
                         for c in row['candidates'] if c['name'] == name]
            for count in ranks:
                selected = [c for row in report['runs'] if row['ranks'] == count
                            for c in row['candidates'] if c['name'] == name]
                wall = statistics.median(c['wall_s'] for c in selected)
                speedup = statistics.median(reference)/wall
                summary.append(dict(candidate=name, ranks=count, repeats=len(selected), median_wall_s=wall,
                    observed_speedup=speedup, observed_efficiency=speedup/count,
                    min_wall_s=min(c['wall_s'] for c in selected), max_wall_s=max(c['wall_s'] for c in selected),
                    iterations=[c['iterations'] for c in selected],
                    median_phases={p: statistics.median(c['phases'][p] for c in selected) for p in selected[0]['phases']}))
        report['summary'] = summary
        assert all(digest(Path(p)) == value for p, value in report['files'].items())
        report['status'] = 'passed'
    except BaseException:
        report['status'] = 'failed'
        raise
    finally:
        save()


if __name__ == '__main__':
    main()
