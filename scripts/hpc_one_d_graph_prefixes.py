#!/usr/bin/env python3
import argparse
import json
from pathlib import Path
root = Path(__file__).resolve().parents[1]
from hpc_native_graph_checkpoint import Acceptance, fixture, same_output, sha
from hpc_solver_prefixes import field_comparison, histories

def main():
    parser = argparse.ArgumentParser(description='Verify independent implicit 1D solvers and restart in native graphs.')
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--fixture-builder', type=Path, default=root / 'solvers/coupling/native_graph_checkpoint_fixture_test')
    parser.add_argument('--ranks', type=int, default=3)
    parser.add_argument('--timeout', type=int, default=240)
    args = parser.parse_args()
    args.binary = args.binary.resolve()
    args.reference = args.reference.resolve()
    args.fixture_builder = args.fixture_builder.resolve()
    args.positive_only = True
    if args.ranks < 1 or args.timeout < 1:
        parser.error('ranks and timeout must be positive')
    a = Acceptance(args)
    a.report['comparisons'] = {}
    a.report['reference_binary'] = {str(args.reference): sha(args.reference)}
    try:
        for mode in ['pressure_network', 'linearized_aq', 'nonlinear_aq', 'implicit_1d_pde']:
            case = a.root / 'fixtures' / mode
            fixture(a.options.fixture_builder, case, 'flow', args.ranks)
            for path in case.glob('*/simulation_config.json'):
                d = json.loads(path.read_text())
                if d.get('dimension') != '1d':
                    continue
                for system in d['equation_systems']:
                    if system['kind'] == 'network_flow_1d':
                        system.update(model='compliant', scheme='implicit_petsc', formulation=mode, wall=dict(model='linear', young_modulus=1000000000.0, thickness_ratio=0.1, reference_pressure=0))
                path.write_text(json.dumps(d, indent=2) + '\n')
            baseline = a.job(mode + '-reference', mode, ['--checkpoint-dir', a.root / (mode + '-reference-bundle')], binary=args.reference)
            default = a.job(mode + '-default', mode, ['--checkpoint-dir', a.root / (mode + '-default-bundle')])
            same_output(baseline, default)
            field_comparison(a.root / (mode + '-reference-bundle'), a.root / (mode + '-default-bundle'), True)
            flags = ['-domain_source_flow_ksp_type', 'fgmres', '-domain_branch_a_flow_ksp_type', 'gmres']
            if mode in ['nonlinear_aq', 'implicit_1d_pde']:
                flags += ['-domain_source_flow_snes_type', 'newtontr']
            alt = a.job(mode + '-override', mode, ['--checkpoint-dir', a.root / (mode + '-override-bundle')] + flags)
            a.report['comparisons'][mode] = dict(fields=field_comparison(a.root / (mode + '-reference-bundle'), a.root / (mode + '-override-bundle'), False), histories=histories(baseline, alt, {}))
            rows = [json.loads(line.split(' ', 1)[1]) for line in (a.root / (mode + '-override') / 'rank-0/stdout.log').read_text().splitlines() if line.startswith('one_d_solver_configuration ')]
            assert rows
            for row in rows:
                prefix = row['prefix']
                expected = 'fgmres' if prefix == 'domain_source_flow_' else 'gmres' if prefix == 'domain_branch_a_flow_' else 'preonly'
                assert row['ksp'] == expected and row['factor_backend'] == 'mumps'
                if mode in ['nonlinear_aq', 'implicit_1d_pde'] and prefix == 'domain_source_flow_':
                    assert row['snes'] == 'newtontr'
            bundle = a.root / (mode + '-restart-bundle')
            a.job(mode + '-save', mode, ['--checkpoint-dir', bundle, '--stop-after-step', '3'])
            resumed = a.job(mode + '-resume', mode, ['--checkpoint-dir', bundle, '--restart-dir', bundle])
            same_output(default, resumed)
            field_comparison(a.root / (mode + '-default-bundle'), bundle, True)
        a.report['status'] = 'passed'
    except Exception as error:
        a.report.update(status='failed', error=str(error))
        raise
    finally:
        (a.root / 'acceptance.json').write_text(json.dumps(a.report, indent=2) + '\n')
if __name__ == '__main__':
    main()
