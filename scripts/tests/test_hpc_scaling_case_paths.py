"""Regression coverage for Bash braces in rank-dependent case defaults."""
import os
from pathlib import Path
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[2]


class ScalingCasePaths(unittest.TestCase):
    def expand(self, overrides):
        wrapper = Path(os.environ.get('IGA_TEST_SCALING_WRAPPER', str(ROOT / 'solvers/cpu/slurm/cross_node_scaling.sbatch')))
        assignments = '\n'.join(line for line in wrapper.read_text().splitlines()
                                if line.startswith(('strong_case=', 'weak_case=')))
        self.assertEqual(len(assignments.splitlines()), 2)
        env = {k: v for k, v in os.environ.items() if k not in ('IGA_STRONG_CASE', 'IGA_WEAK_CASE')}
        env.update(overrides)
        return subprocess.check_output(['bash', '-c', assignments + '\nprintf "%s\\n" "$strong_case" "$weak_case"'], env=env, text=True).splitlines()

    def test_default_rank_patterns(self):
        patterns = self.expand({})
        self.assertEqual(patterns, ['strong-{ranks}/root3d', 'weak-{ranks}/root3d'])
        self.assertEqual([x.format(ranks=256) for x in patterns], ['strong-256/root3d', 'weak-256/root3d'])

    def test_explicit_patterns(self):
        patterns = ['custom strong/{ranks}/domain', 'custom weak/{ranks}/domain']
        self.assertEqual(self.expand(dict(zip(['IGA_STRONG_CASE', 'IGA_WEAK_CASE'], patterns))), patterns)


if __name__ == '__main__':
    unittest.main()
