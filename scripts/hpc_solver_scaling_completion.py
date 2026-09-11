#!/usr/bin/env python3
"""Audit accepted HPC-04B/C candidate and rank-scaling evidence."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def digest_bytes(data):
	return hashlib.sha256(data).hexdigest()


def digest(path):
	return digest_bytes(path.read_bytes())


def load(path):
	return json.loads(path.read_text())


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('--evidence-root', type=Path, default=Path('outputs/hpc04'))
	parser.add_argument('--output-dir', type=Path, required=True)
	args = parser.parse_args()
	repo = Path(__file__).resolve().parents[1]
	root = (repo/args.evidence_root).resolve() if not args.evidence_root.is_absolute() else args.evidence_root.resolve()
	out = args.output_dir.resolve()
	out.mkdir(parents=True, exist_ok=False)

	manifests = {
		'current_small': root/'completion-current-small/acceptance.json',
		'current_small_audit': root/'completion-current-small/audit.json',
		'medium': root/'pde-candidates/refined-four-v1/acceptance.json',
		'medium_audit': root/'pde-candidates/refined-four-v1/audit.json',
		'large': root/'pde-candidates/large-two-v1/acceptance.json',
		'large_audit': root/'pde-candidates/large-two-v1/audit.json',
		'selfp': root/'pde-candidates/selfp-small-v2/acceptance.json',
		'selfp_audit': root/'pde-candidates/selfp-small-v2/audit.json',
		'immersed': root/'immersed-candidates/regular-acceptance.json',
		'immersed_small_cut': root/'immersed-candidates/small-cut-closed-v1/acceptance.json',
		'immersed_small_cut_audit': root/'immersed-candidates/small-cut-acceptance.json',
		'rank_scaling': root/'rank-scaling/fixed128-v1/acceptance.json',
		'rank_scaling_audit': root/'rank-scaling/fixed128-v1/audit.json',
	}
	data = {name: load(path) for name, path in manifests.items()}
	for name in ('current_small', 'medium', 'large', 'immersed', 'rank_scaling'):
		assert data[name]['status'] == 'passed', name
	for name in ('current_small_audit', 'medium_audit', 'large_audit', 'rank_scaling_audit'):
		assert data[name]['status'] == 'passed', name
	assert data['selfp']['status'] == 'completed_with_failures'
	assert data['selfp_audit']['status'] == 'evaluated'
	assert data['immersed_small_cut']['status'] == 'completed_with_failures'
	assert data['immersed_small_cut_audit']['status'] == 'evaluated'

	duct = []
	for name, expected in (('current_small', 4), ('medium', 4), ('large', 2)):
		row = data[name]
		assert len(row['candidates']) == expected
		assert all(candidate['status'] == 'passed' and candidate['returncode'] == 0
			for candidate in row['candidates'])
		duct.append({'suite': name, 'mesh': row['mesh'],
			'candidates': [candidate['name'] for candidate in row['candidates']]})

	selfp = data['selfp']['candidates']
	assert len(selfp) == 4 and sum(x['status'] == 'passed' for x in selfp) == 3
	assert sum(x['status'] == 'failed' for x in selfp) == 1
	regular = data['immersed']
	assert regular['jobs'] == 9 and regular['rank_reports'] == 18
	small_cut = data['immersed_small_cut']['jobs']
	assert len(small_cut) == 3 and sum(x['status'] == 'passed' for x in small_cut) == 1
	assert sum(x['status'] == 'failed' for x in small_cut) == 2

	scaling = data['rank_scaling']
	assert len(scaling['runs']) == 6 and len(scaling['summary']) == 6
	assert data['rank_scaling_audit']['rank_reports'] == 28
	assert len(data['rank_scaling_audit']['runs']) == 12
	assert all(len(row['iterations']) == 2 for row in scaling['summary'])

	launch = load(root/'pde-candidates/large-two-v1/launch-source.json')
	for path, expected in launch['sources'].items():
		content = subprocess.check_output(['git', 'show', launch['launch_commit']+':'+path], cwd=repo)
		assert digest_bytes(content) == expected, path
	assert digest(root/'immersed/native-v4-binary') == launch['native_binary_sha256']
	core = 'solvers/cpu/include/TransientFlowRuntime.hpp'
	assert subprocess.check_output(['git', 'show', launch['native_runtime_commit']+':'+core], cwd=repo) == (repo/core).read_bytes()

	report = {
		'status': 'passed',
		'base_revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip(),
		'evaluations': {'passed': 35, 'expected_nonconverged': 3, 'total': 38},
		'duct_meshes': duct,
		'immersed': {'regular_passed': 9, 'small_cut_passed': 1,
			'small_cut_expected_nonconverged': 2},
		'rank_scaling': scaling['summary'],
		'evidence_sha256': {name: digest(path) for name, path in manifests.items()},
		'historical_source': {'launch_commit': launch['launch_commit'],
			'native_runtime_commit': launch['native_runtime_commit'], 'verified_files': len(launch['sources'])},
		'limitations': [
			'Historical medium and large runs use authenticated frozen binaries; the current small matrix checks present solver behavior.',
			'Workstation timings are not an uncontended or cross-node benchmark.',
			'Negative convergence results remain reported and do not change solver defaults.',
		],
	}
	(out/'acceptance.json').write_text(json.dumps(report, indent=2)+'\n')
	print('HPC-04B/C completion audit passed: 38 evaluations, 12 rank-scaling solves')


if __name__ == '__main__':
	main()
