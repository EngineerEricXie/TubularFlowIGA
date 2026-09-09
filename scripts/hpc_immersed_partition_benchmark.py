#!/usr/bin/env python3
"""Compare contiguous cell-count and work-weighted immersed assembly on local MPI."""
import argparse
import json
from pathlib import Path
import statistics
import subprocess
import sys

from hpc_inventory import digest


def summarize(records, ranks):
    result = {}
    for partition in ['cell-count', 'weighted']:
        selected = [r for r in records if r['partition'] == partition]
        samples = {key: [] for key in ['assembly_s', 'integration_insert_s', 'halo_s', 'stash_s', 'integration_imbalance']}
        resources = []
        for record in selected:
            observations = record['observations']
            if len(observations) != ranks:
                raise RuntimeError('incomplete partition observations')
            affinity = [set(r['affinity_cpus']) for r in record['rank_reports']]
            if any(not a for a in affinity) or any(affinity[i] & affinity[j] for i in range(ranks) for j in range(i)):
                raise RuntimeError('timing ranks do not have disjoint CPU affinity')
            resources.extend(r['resource']['peak_rss_bytes'] for r in record['rank_reports'])
            # Repeat zero warms matrix/halo caches. The remaining repetitions
            # are synchronized samples; retain per-rank raw data in the child.
            for repeat in range(1, record['repetitions']):
                group = [s for s in record['work_samples'] if s['repeat'] == repeat]
                if len(group) != ranks:
                    raise RuntimeError('incomplete partition timing sample')
                for key in samples:
                    if key == 'integration_imbalance':
                        values = [s['integration_insert_s'] for s in group]
                        samples[key].append(max(values)/statistics.mean(values))
                    else:
                        samples[key].append(max(s[key] for s in group))
        first = selected[0]['observations']
        work = [o['estimated_work'] for o in first]
        result[partition] = dict(samples=samples, medians={k: statistics.median(v) for k, v in samples.items()},
                                 estimated_work=work, estimated_work_imbalance=max(work)/statistics.mean(work),
                                 required_rows=[o['halo_rows'] for o in first],
                                 remote_rows=[o['remote_halo_rows'] for o in first],
                                 maximum_rank_peak_rss_bytes=max(resources))
    count, weighted = result['cell-count'], result['weighted']
    if sum(count['estimated_work']) != sum(weighted['estimated_work']):
        raise RuntimeError('partition changed total estimated work')
    if max(weighted['estimated_work']) > max(count['estimated_work']):
        raise RuntimeError('weighted contiguous bound exceeds count partition')
    result['assembly_ratio_weighted_to_count'] = weighted['medians']['assembly_s']/count['medians']['assembly_s']
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--ranks', type=int, nargs='+', choices=[1, 2, 4], default=[1, 2, 4])
    parser.add_argument('--repetitions', type=int, choices=range(3, 11), default=4)
    args = parser.parse_args()
    if len(set(args.ranks)) != len(args.ranks):
        parser.error('rank counts must be distinct')
    repo = Path(__file__).resolve().parents[1]
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    binary = repo/'solvers/cpu/immersed_distributed_physics_test'
    result = dict(status='running', binary=str(binary), binary_sha256=digest(binary),
                  revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip(),
                  cases=[], comparisons={})
    try:
        for ranks in args.ranks:
            records = []
            order = ['cell-count', 'weighted'] if ranks == 1 else ['cell-count', 'weighted', 'weighted', 'cell-count']
            for index, partition in enumerate(order):
                directory = root/f'{ranks}-{index}-{partition}'
                command = [sys.executable, str(repo/'scripts/hpc_immersed_assembly_regression.py'),
                           '--kind', 'physics', '--physics-mode', 'flow', '--partition', partition,
                           '--repetitions', str(args.repetitions), '--ranks', str(ranks), '--bind-to-core',
                           '--output-dir', str(directory)]
                with (root/f'{directory.name}.log').open('w') as log:
                    completed = subprocess.run(command, cwd=repo, stdout=log, stderr=subprocess.STDOUT)
                item = dict(argv=command, returncode=completed.returncode, summary=str(directory/'summary.json'))
                result['cases'].append(item)
                if completed.returncode:
                    raise RuntimeError('partition benchmark child failed: '+directory.name)
                child = json.loads((directory/'summary.json').read_text())
                if child['status'] != 'passed' or child['binaries'][str(binary)] != result['binary_sha256']:
                    raise RuntimeError('partition benchmark evidence/binary differs')
                records.extend(child['cases'])
                item['summary_sha256'] = digest(directory/'summary.json')
                print(directory.name, 'passed', flush=True)
            result['comparisons'][str(ranks)] = summarize(records, ranks)
        if digest(binary) != result['binary_sha256']:
            raise RuntimeError('partition benchmark binary changed')
        result['status'] = 'passed'
    except BaseException:
        result['status'] = 'failed'
        raise
    finally:
        (root/'summary.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
