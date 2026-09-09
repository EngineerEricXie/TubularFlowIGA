#!/usr/bin/env python3
"""Fresh MPI-job checkpoint acceptance; isolated fixtures and immutable evidence."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess

REPO = Path(__file__).resolve().parents[1]
MODES = ('zero', 'flow', 'species', 'zero-fixed', 'flow-aitken', 'species-reverse', 'multi')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generations(root):
    result = {}
    for manifest in root.glob('*/manifest'):
        try:
            lines = manifest.read_bytes().splitlines(keepends=True)
            if lines[-1].decode().strip() != 'sha256 ' + hashlib.sha256(b''.join(lines[:-1])).hexdigest():
                continue
            header = dict(line.decode().strip().split(' ', 1) for line in lines[:10])
            complete = True
            for line in lines[11:-1]:
                tag, name, owner, kind, size, digest = line.decode().split()
                shard = manifest.parent / (name + '.shard')
                if tag != 'shard' or not shard.is_file() or sha(shard) != digest:
                    complete = False
                    break
                fields = shard.read_bytes().split(b'\n', 6)
                if len(fields) != 7 or len(fields[6]) != int(size) or fields[1].decode() != 'epoch ' + header['epoch']:
                    complete = False
                    break
            if not complete:
                continue
            count = int(header['accepted_steps'])
        except (OSError, UnicodeError, ValueError, KeyError, IndexError):
            continue
        if count in result:
            raise AssertionError('ambiguous accepted generation')
        result[count] = manifest.parent
    return result


def images(root, relative=True):
    return {str(p.relative_to(root) if relative else p): sha(p)
            for p in sorted(root.rglob('*')) if p.is_file()}


def payloads(epoch):
    return {p.name: hashlib.sha256(p.read_bytes().split(b'\n', 6)[6]).hexdigest()
            for p in epoch.glob('*.shard')}


def same_output(first, second, relocated=False):
    files = {p.name for p in first.iterdir()}
    assert files == {p.name for p in second.iterdir()}
    for name in files:
        if relocated and name == 'graph_binding_manifest.json':
            continue  # This output intentionally names the new absolute case paths.
        assert (first / name).read_bytes() == (second / name).read_bytes(), name
        if name.endswith('.csv'):
            with (second / name).open() as source:
                for row in csv.DictReader(source):
                    for value in row.values():
                        try:
                            number = float(value)
                        except ValueError:
                            continue
                        assert math.isfinite(number), (name, value)


def fixture(builder, root, mode, ranks):
    base = mode.split('-')[0]
    subprocess.run([str(builder), str(root), base, str(ranks)], check=True)
    for path in root.rglob('simulation_config.json'):
        value = json.loads(path.read_text())
        value['time']['steps'] = 8 if base == 'zero' else 6
        path.write_text(json.dumps(value, indent=2) + '\n')
    path = root / 'simulation_config.json'
    value = json.loads(path.read_text())
    if mode in ('zero-fixed', 'flow-aitken'):
        value['execution'].update(kind='fixed' if mode == 'zero-fixed' else 'aitken',
                                  maximum_iterations=60, pressure_relative_tolerance=1e-6,
                                  pressure_reference_pa=1, flow_relative_tolerance=1e-10,
                                  relaxation_factor=0.5, minimum_relaxation=0.05,
                                  maximum_relaxation=0.999)
    path.write_text(json.dumps(value, indent=2) + '\n')
    if mode == 'species-reverse':
        path = root / 'species_source/simulation_config.json'
        value = json.loads(path.read_text())
        value['temporal_functions'] = [dict(name='inlet_flow', kind='periodic_table',
                                          units='m3/s', period=1, file='flow.csv', interpolation='linear')]
        path.write_text(json.dumps(value, indent=2) + '\n')
        (path.parent / 'flow.csv').write_text('time,value\n0,0.001\n0.03,0.001\n0.04,0\n0.05,-0.001\n0.06,-0.001\n0.07,0.001\n')


class Acceptance:
    def __init__(self, options):
        self.options = options
        self.root = options.output_root.resolve()
        self.root.mkdir()
        self.environment = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                                PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps', IGA_PROFILE='1')
        self.report = dict(status='running', jobs=[], checks=[], binaries={
            str(options.binary): sha(options.binary), str(options.fixture_builder): sha(options.fixture_builder)})
        if not options.positive_only:
            self.report['binaries'][str(options.fault_binary)] = sha(options.fault_binary)

    def job(self, name, mode, extra=(), expected=0, fault=None, crash=False, case=None, binary=None, environment=None, ranks=None):
        log = self.root / name
        output = self.root / (name + '-output')
        log.mkdir()
        env = dict(self.environment)
        if fault:
            env['IGA_NATIVE_CHECKPOINT_FAULT'] = fault
        if environment:
            env.update(environment)
        argv = [str(binary or self.options.binary), '--graph-case', str(case or self.root / 'fixtures' / mode),
                '--output-dir', str(output)] + [str(x) for x in extra]
        job_ranks = self.options.ranks if ranks is None else ranks
        command = ['mpiexec', '--oversubscribe', '-np', str(job_ranks)]
        if expected == 0:
            command += ['python3', str(REPO / 'scripts/hpc_rank_run.py'), '--output-dir', str(log),
                        '--expected-ranks', str(job_ranks), '--timeout', str(self.options.timeout), '--']
        command += argv
        with (log / 'launcher.log').open('w') as stream:
            code = subprocess.call(['timeout', '--kill-after=5s', str(self.options.timeout + 15)] + command,
                                   cwd=REPO, env=env, stdout=stream, stderr=subprocess.STDOUT)
        self.report['jobs'].append(dict(name=name, command=command, returncode=code, expected=expected))
        print(name, 'exit', code, 'expected', expected, flush=True)
        assert code == expected, (name, code, expected)
        if expected == 0:
            for rank in range(job_ranks):
                directory = log / f'rank-{rank}'
                report = json.loads((directory / 'run.json').read_text())
                assert report['returncode'] == expected and not report['timed_out']
                assert report['resource'] is not None
                for file, digest in report['logs'].items():
                    assert sha(directory / file) == digest
                if not expected:
                    assert (directory / 'stderr.log').stat().st_size == 0
        if expected:
            assert not output.exists(), 'failed graph published an output directory'
        return output

    def copy_prefix(self, mode, name, steps=3):
        target = self.root / name
        target.mkdir()
        for count, epoch in generations(self.root / f'{mode}-full-bundle').items():
            if count <= steps:
                shutil.copytree(epoch, target / epoch.name)
        assert max(generations(target)) == steps
        return target

    def recover(self, name, mode, bundle, restored=3):
        output = self.job(name, mode, ['--restart-dir', bundle, '--checkpoint-dir', bundle])
        text = (self.root / name / 'rank-0/stdout.log').read_text()
        assert f'restored graph checkpoint step={restored} ' in text
        same_output(self.root / f'{mode}-full-output', output)
        expected = generations(self.root / f'{mode}-full-bundle')
        actual = generations(bundle)
        assert set(actual) == set(expected)
        for count, path in actual.items():
            assert payloads(path) == payloads(expected[count]), (name, count, 'payload differs')

    def positive(self, mode):
        fixture(self.options.fixture_builder, self.root / 'fixtures' / mode, mode, self.options.ranks)
        full = self.root / f'{mode}-full-bundle'
        split = self.root / f'{mode}-split-bundle'
        self.job(mode + '-full', mode, ['--checkpoint-dir', full])
        self.job(mode + '-save', mode, ['--checkpoint-dir', split, '--stop-after-step', '3'])
        self.recover(mode + '-resume', mode, split)
        assert images(full) == images(split), (mode, 'checkpoint envelope differs')
        if mode == 'species-reverse':
            with (self.root / f'{mode}-full-output/species_edge_amounts.csv').open() as source:
                rows = [row for row in csv.DictReader(source) if row['edge_id'] == 'edge_in' and row['species_id'] == 'red_logical']
            assert rows[2]['donor'] == rows[3]['donor'] != rows[4]['donor']
            assert abs(float(rows[3]['first_outward_amount'])) < 1e-12
        if mode in ('zero-fixed', 'flow-aitken'):
            with (self.root / f'{mode}-full-output/pressure_flow_steps.csv').open() as source:
                assert any(int(row['iterations']) > 1 for row in csv.DictReader(source))
            initial = {edge['id']: edge['initial_pressure_pa'] for edge in json.loads((self.root / 'fixtures' / mode / 'simulation_config.json').read_text())['couplings']}
            with (self.root / f'{mode}-full-output/pressure_flow_iterations.csv').open() as source:
                assert any(row['step'] == '4' and row['iteration'] == '1' and float(row['applied_pressure_pa']) != initial[row['edge_id']] for row in csv.DictReader(source))
        self.report['checks'].append(mode + ': complete histories and checkpoint bytes exact')

    def interruptions(self, mode):
        for fault, code in [('precommit', 1), ('stream', 86), ('manifest', 87), ('local-write', 1)]:
            name = mode + '-' + fault
            bundle = self.copy_prefix(mode, name + '-bundle')
            before = images(bundle, relative=False)
            env = {'TUBULARFLOWIGA_INJECT_BIFURCATION_FAILURE_STEP': '4'} if fault == 'precommit' else None
            self.job(name, mode, ['--restart-dir', bundle, '--checkpoint-dir', bundle, '--stop-after-step', '4'],
                     expected=code, fault=fault, crash=code in (86, 87), environment=env,
                     binary=self.options.binary if fault == 'precommit' else self.options.fault_binary)
            assert all(sha(Path(path)) == digest for path, digest in before.items())
            assert max(generations(bundle)) == 3
            if code in (86, 87):
                text = (self.root / name / 'launcher.log').read_text()
                assert ('stage=stream bytes=5' if code == 86 else 'stage=before-manifest') in text
                incomplete = [p for p in bundle.iterdir() if p.is_dir() and not (p / 'manifest').exists()]
                assert len(incomplete) == 1
                if code == 86:
                    assert any(p.name.endswith('.flow.rank-0.tmp') for p in incomplete[0].iterdir())
                else:
                    assert (incomplete[0] / 'manifest.tmp').is_file()
            self.recover(name + '-recovery', mode, bundle)
            self.report['checks'].append(name + ': previous generation preserved and fresh-job recovery exact')

    def compatibility(self, mode):
        for fault in ('corrupt', 'missing', 'mixed'):
            bundle = self.copy_prefix(mode, mode + '-' + fault + '-bundle')
            epochs = generations(bundle)
            victim = next(epochs[3].glob('*.flow.rank-0.shard'))
            if fault == 'corrupt':
                data = bytearray(victim.read_bytes()); data[-1] ^= 1; victim.write_bytes(data)
            elif fault == 'missing':
                victim.unlink()
            else:
                shutil.copy2(epochs[2] / victim.name, victim)
            # Invalid generation 3 remains on disk; recovery chooses 2 and writes
            # a different generation ID for 3 without overwriting the bad files.
            bad = images(epochs[3], relative=False)
            self.recover(mode + '-' + fault + '-recovery', mode, bundle, restored=2)
            assert all(sha(Path(path)) == digest for path, digest in bad.items())
            self.report['checks'].append(fault + ': rejected newest bundle and recovered previous complete state')
        source = self.root / f'{mode}-full-bundle'
        before = images(source)
        self.job('wrong-controls', mode, ['--restart-dir', source, '--three-d-max-newton', '31'], expected=1)
        self.job('wrong-source', mode, ['--restart-dir', source], expected=1, binary=self.options.fault_binary,
                 environment={'IGA_NATIVE_CHECKPOINT_TEST_SOURCE_SHA256': 'a' * 64})
        self.job('wrong-petsc', mode, ['--restart-dir', source], expected=1,
                 environment={'PETSC_OPTIONS': self.environment['PETSC_OPTIONS'] + ' -ksp_rtol 1e-9'})
        self.job('short-horizon', mode, ['--restart-dir', source, '--stop-after-step', '2'], expected=1)
        self.job('wrong-ranks', mode, ['--restart-dir', source], expected=1, ranks=self.options.ranks + 1)
        case = self.root / 'changed-case'
        shutil.copytree(self.root / 'fixtures' / mode, case)
        manifest = case / 'simulation_config.json'
        value = json.loads(manifest.read_text()); value['couplings'][0]['initial_pressure_pa'] += 0.001
        manifest.write_text(json.dumps(value, indent=2) + '\n')
        self.job('wrong-case', mode, ['--restart-dir', source], expected=1, case=case)
        self.job('interval-without-output', mode, ['--checkpoint-every', '1'], expected=1)
        empty = self.root / 'empty-bundle'; empty.mkdir()
        self.job('no-complete-checkpoint', mode, ['--restart-dir', empty], expected=1)
        if 'species' in self.options.modes:
            changed = self.root / 'changed-field-order'
            shutil.copytree(self.root / 'fixtures/species', changed)
            config = changed / 'junction/simulation_config.json'
            value = json.loads(config.read_text())
            next(system for system in value['equation_systems'] if system['kind'] == 'linear_transport')['unknowns'].reverse()
            config.write_text(json.dumps(value, indent=2) + '\n')
            self.job('wrong-field-order', 'species', ['--restart-dir', self.root / 'species-full-bundle'], expected=1, case=changed)
        if 'species-reverse' in self.options.modes:
            changed = self.root / 'changed-external-table'
            shutil.copytree(self.root / 'fixtures/species-reverse', changed)
            table = changed / 'species_source/flow.csv'
            table.write_text(table.read_text().replace('0,0.001\n', '0,0.0011\n', 1))
            self.job('wrong-external-table', 'species-reverse', ['--restart-dir', self.root / 'species-reverse-full-bundle'], expected=1, case=changed)
        assert images(source) == before
        relocated = self.root / 'relocated-case'; shutil.copytree(self.root / 'fixtures' / mode, relocated)
        output = self.job('relocated-case-restore', mode, ['--restart-dir', source], case=relocated)
        same_output(self.root / f'{mode}-full-output', output, relocated=True)
        self.report['checks'].append('compatibility rejection, read-only source preservation, relocated case restore')

    def signal(self, mode):
        bundle = self.root / 'signal-bundle'
        self.job('signal-stop', mode, ['--checkpoint-dir', bundle, '--checkpoint-every', '3'],
                 fault='signal', binary=self.options.fault_binary)
        assert set(generations(bundle)) == {3, 4}
        assert 'stopped after checkpoint request step=4' in (self.root / 'signal-stop/rank-0/stdout.log').read_text()
        output = self.job('signal-resume', mode, ['--restart-dir', bundle, '--checkpoint-dir', bundle])
        same_output(self.root / f'{mode}-full-output', output)
        self.report['checks'].append('single-rank SIGUSR1 saved at accepted step 4 outside interval and resumed exactly')
        early = self.root / 'signal-input-bundle'
        self.job('signal-input-stop', mode, ['--checkpoint-dir', early, '--checkpoint-every', '3'],
                 fault='signal-input', binary=self.options.fault_binary)
        assert set(generations(early)) == {1}
        assert 'stopped after checkpoint request step=1' in (self.root / 'signal-input-stop/rank-0/stdout.log').read_text()
        output = self.job('signal-input-resume', mode, ['--restart-dir', early, '--checkpoint-dir', early])
        same_output(self.root / f'{mode}-full-output', output)
        self.report['checks'].append('SIGUSR1 during asset input deferred until the first accepted checkpoint')

    def run(self):
        try:
            for mode in self.options.modes:
                self.positive(mode)
            if not self.options.positive_only:
                for mode in self.options.modes:
                    if mode in ('zero', 'flow', 'species'):
                        self.interruptions(mode)
                mode = 'flow' if 'flow' in self.options.modes else self.options.modes[0]
                self.compatibility(mode)
                self.signal('species' if 'species' in self.options.modes else mode)
            self.report['status'] = 'passed'
            self.report['scope'] = 'positive-only' if self.options.positive_only else 'native-checkpoint-integration'
        except Exception as error:
            self.report['status'] = 'failed'
            self.report['error'] = str(error)
            raise
        finally:
            (self.root / 'acceptance.json').write_text(json.dumps(self.report, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, default=REPO / 'solvers/coupling/iga_multidomain_flow')
    parser.add_argument('--fault-binary', type=Path, default=REPO / 'solvers/coupling/native_graph_checkpoint_fault_test')
    parser.add_argument('--fixture-builder', type=Path, default=REPO / 'solvers/coupling/native_graph_checkpoint_fixture_test')
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--ranks', type=int, default=3)
    parser.add_argument('--timeout', type=int, default=240)
    parser.add_argument('--modes', type=lambda text: text.split(','), default=list(MODES))
    parser.add_argument('--positive-only', action='store_true')
    options = parser.parse_args()
    if options.ranks < (1 if options.positive_only else 2) or options.timeout < 1 or not options.modes or len(set(options.modes)) != len(options.modes) or any(mode not in MODES for mode in options.modes):
        parser.error('require positive ranks (at least two for fault tests), a positive timeout and unique known modes')
    options.binary = options.binary.resolve(); options.fault_binary = options.fault_binary.resolve(); options.fixture_builder = options.fixture_builder.resolve()
    Acceptance(options).run()


if __name__ == '__main__':
    main()
