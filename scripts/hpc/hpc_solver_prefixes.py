#!/usr/bin/env python3
"""Compare native per-domain PETSc overrides with an archived executable."""
import argparse
import csv
import json
import math
from pathlib import Path
import struct
from hpc_native_graph_checkpoint import Acceptance, fixture, generations, same_output, sha


def norms(first, second):
    assert len(first) == len(second)
    assert all(math.isfinite(x) for x in first + second)
    assert all(abs(a-b) <= 1e-12 + 1e-6*abs(a) for a, b in zip(first, second))
    reference = math.sqrt(math.fsum(x*x for x in first))
    difference = math.sqrt(math.fsum((a-b)**2 for a, b in zip(first, second)))
    assert difference <= 1e-12 + 1e-6 * reference, (difference, reference)
    return dict(reference_l2=reference, difference_l2=difference,
                relative_l2=difference/reference if reference else (0.0 if not difference else None))


def field_comparison(first, second, exact):
    old, new = generations(first), generations(second)
    assert old and set(old) == set(new)
    result = []
    def fields(epoch):
        return {p.name for p in epoch.glob('*.shard')
                if '.flow.rank-' in p.name or '.transport.rank-' in p.name or p.name.endswith('.one-d.fields.shard')}
    for step in old:
        names = fields(old[step])
        assert names and names == fields(new[step])
        for name in sorted(names):
            path = old[step]/name
            a = path.read_bytes().split(b'\n', 6)[6]
            b = (new[step]/path.name).read_bytes().split(b'\n', 6)[6]
            if exact:
                assert a == b, path.name
                result.append(dict(step=step, shard=path.name, exact=True))
            else:
                header = 24 if '.rank-' in path.name else 0
                assert a[:header] == b[:header] and len(a) == len(b)
                assert len(a) >= header and (len(a)-header) % 8 == 0
                n = (len(a)-header)//8
                metric = norms(list(struct.unpack('<'+'d'*n, a[header:])), list(struct.unpack('<'+'d'*n, b[header:])))
                result.append(dict(step=step, shard=path.name, **metric))
    return result


def histories(first, second, amount_tolerances):
    result = {}
    assert {p.name for p in first.iterdir()} == {p.name for p in second.iterdir()}
    for path in first.glob('*.csv'):
        with path.open() as source, (second/path.name).open() as target:
            a, b = list(csv.reader(source)), list(csv.reader(target))
        assert len(a) == len(b) and a[0] == b[0]
        columns = {}
        for index, name in enumerate(a[0]):
            numeric_a, numeric_b = [], []
            for row_a, row_b in zip(a[1:], b[1:]):
                assert len(row_a) == len(row_b) == len(a[0])
                x, y = row_a[index], row_b[index]
                if name in ('step', 'iteration', 'iterations', 'node_id', 'element_id', 'label'):
                    assert x == y
                try:
                    xx, yy = float(x), float(y)
                except ValueError:
                    assert x == y
                else:
                    numeric_a.append(xx); numeric_b.append(yy)
            if name in ('normalized_residual', 'normalized_amount_residual'):
                # This is already a dimensionless conservation error. Compare
                # using the configured relative conservation tolerance, not
                # a second relative error against a roundoff-sized residual.
                species_column = a[0].index('species_id')
                for row_a, row_b, x, y in zip(a[1:], b[1:], numeric_a, numeric_b):
                    assert row_a[species_column] == row_b[species_column]
                    tolerance = amount_tolerances[row_a[species_column]]
                    bound = tolerance['relative_tolerance'] + tolerance['absolute_tolerance']/tolerance['reference_amount']
                    assert math.isfinite(x) and math.isfinite(y) and 0 <= x <= bound and 0 <= y <= bound
                    assert abs(x-y) <= tolerance['relative_tolerance']
                columns[name] = dict(max_reference=max(numeric_a, default=0), max_candidate=max(numeric_b, default=0),
                                     max_absolute_difference=max((abs(x-y) for x, y in zip(numeric_a, numeric_b)), default=0),
                                     criterion='original configured dimensionless conservation tolerance')
            else:
                columns[name] = norms(numeric_a, numeric_b)
        result[path.name] = columns
    return result


def configurations(log, expected):
    rows = [json.loads(line.split(' ', 1)[1]) for line in log.read_text().splitlines()
            if line.startswith('solver_configuration ')]
    assert rows
    for row in rows:
        key = (row['domain'], row['role'])
        assert key in expected and row['ksp'] == expected[key]
        assert row['prefix'] == f'domain_{key[0]}_{key[1]}_'
        assert row['pc'] == 'lu' and row['factor_backend'] == 'mumps'
        assert row['last_reason'] > 0 and row['last_iterations'] >= 0
    assert {tuple((r['domain'], r['role'])) for r in rows} == set(expected)
    assert all(any((r['domain'], r['role']) == key and r['last_iterations'] > 0 for r in rows) for key in expected)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--fixture-builder', type=Path, default=Path('solvers/coupling/native_graph_checkpoint_fixture_test'))
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--ranks', type=int, default=3)
    parser.add_argument('--timeout', type=int, default=240)
    options = parser.parse_args()
    options.binary = options.binary.resolve(); options.reference = options.reference.resolve()
    options.fixture_builder = options.fixture_builder.resolve(); options.positive_only = True
    if options.ranks < 2 or options.timeout < 1:
        parser.error('this native backend rejection suite requires at least two ranks and a positive timeout')
    a = Acceptance(options)
    a.report['reference_binary'] = {str(options.reference): sha(options.reference)}
    a.report['comparisons'] = {}
    try:
        for mode in ('zero', 'flow', 'species', 'multi'):
            fixture(options.fixture_builder, a.root/'fixtures'/mode, mode, options.ranks)
            reference_bundle, default_bundle = a.root/(mode+'-reference-bundle'), a.root/(mode+'-default-bundle')
            reference = a.job(mode+'-reference', mode, ['--checkpoint-dir', reference_bundle], binary=options.reference)
            default = a.job(mode+'-default', mode, ['--checkpoint-dir', default_bundle])
            same_output(reference, default)
            a.report['comparisons'][mode+'-default'] = field_comparison(reference_bundle, default_bundle, True)
            if mode == 'zero':
                continue
            expected = {('junction', 'flow'): 'fgmres'}
            if mode == 'species': expected[('junction', 'transport')] = 'gmres'
            if mode == 'multi': expected = {('island_a', 'flow'): 'fgmres', ('island_b', 'flow'): 'gmres'}
            flags = []
            for (domain, role), kind in expected.items(): flags += [f'-domain_{domain}_{role}_ksp_type', kind]
            bundle = a.root/(mode+'-override-bundle')
            output = a.job(mode+'-override', mode, ['--checkpoint-dir', bundle]+flags)
            rows = configurations(a.root/(mode+'-override')/'rank-0/stdout.log', expected)
            a.report['comparisons'][mode+'-override'] = dict(fields=field_comparison(reference_bundle, bundle, False),
                                                           histories=histories(reference, output, json.loads((a.root/'fixtures'/mode/'simulation_config.json').read_text())['execution'].get('species_amount_tolerances', {})), configurations=rows)
        a.job('unavailable-scoped-backend', 'flow', ['-domain_junction_flow_pc_factor_mat_solver_type', 'petsc'], expected=1)
        a.job('unknown-scoped-transport', 'species', ['-domain_junction_transport_ksp_type', 'unavailable_solver'], expected=1)
        output = a.job('healthy-after-rejections', 'species')
        same_output(a.root/'species-reference-output', output)
        a.report.update(status='passed', scope='native body-fitted solver prefixes; inherited defaults exact; independent overrides retain numerical gates')
    except Exception as error:
        a.report.update(status='failed', error=str(error)); raise
    finally:
        (a.root/'acceptance.json').write_text(json.dumps(a.report, indent=2)+'\n')


if __name__ == '__main__':
    main()
