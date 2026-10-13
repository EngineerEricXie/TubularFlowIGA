#!/usr/bin/env python3
"""Run a historical `make -C <directory> <target>` through CMake and CTest.

The per-directory Makefiles call this script for every goal. Executables are
built with `cmake --build`; test targets build what they need and run the
matching CTest tests (label `make=<target>`) from that directory's build tree.
Make command-line variables arrive through the environment:

  PETSC_DIR, PETSC_ARCH   PETSc to configure with (kept once configured)
  CUDA_ARCHS, NVCC        enable the CUDA backend for these architectures/compiler
  EIGEN_DIR               Eigen 3 include directory for the spline preprocessor
  HDF5_CFLAGS, HDF5_LIBS  explicit HDF5 flags instead of find_package(HDF5)
  TFI_BUILD_DIR           build directory (default: <repository>/build)
"""
import fcntl
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
CUDA_GOALS_EXEMPT = {'execution-test', 'cuda_execution_test', 'clean'}


def fail(message):
	print(f'forward.py: {message}', file=sys.stderr)
	sys.exit(2)


def dry_run():
	"""True under `make -n`; the recursive recipe still runs, so only print."""
	flags = os.environ.get('MAKEFLAGS', '').split(' ', 1)[0]
	return not flags.startswith('-') and 'n' in flags


def is_open(fd):
	try:
		os.fstat(fd)
		return True
	except OSError:
		return False


def run(command):
	print('+ ' + ' '.join(command), flush=True)
	if dry_run():
		return 0
	# Keep make's jobserver pipe (make < 4.4) open for the generated build.
	match = re.search(r'--jobserver-(?:auth|fds)=(\d+),(\d+)', os.environ.get('MAKEFLAGS', ''))
	fds = [int(fd) for fd in match.groups()] if match else []
	return subprocess.run(command, pass_fds=[fd for fd in fds if is_open(fd)]).returncode


def parallel_level():
	match = re.search(r'(?:^|\s)-j\s*(\d+)', os.environ.get('MAKEFLAGS', ''))
	if match:
		return int(match.group(1))
	return int(os.environ.get('CMAKE_BUILD_PARALLEL_LEVEL', min(8, os.cpu_count() or 1)))


def read_cache(build):
	values = {}
	cache = build / 'CMakeCache.txt'
	if cache.exists():
		for line in cache.read_text().splitlines():
			match = re.match(r'^([A-Za-z0-9_]+):[A-Z]+=(.*)$', line)
			if match:
				values[match.group(1)] = match.group(2)
	return values


def desired_settings(directory, goal):
	settings = {}
	environment = os.environ
	if environment.get('PETSC_DIR'):
		settings['PETSC_DIR'] = environment['PETSC_DIR']
		settings['PETSC_ARCH'] = environment.get('PETSC_ARCH', '')
	if environment.get('EIGEN_DIR'):
		settings['EIGEN_DIR'] = environment['EIGEN_DIR']
	if environment.get('HDF5_CFLAGS'):
		settings['TFI_HDF5_CFLAGS'] = environment['HDF5_CFLAGS']
	if environment.get('HDF5_LIBS'):
		settings['TFI_HDF5_LIBS'] = environment['HDF5_LIBS']
	wants_cuda = directory == 'solvers/cuda' and goal not in CUDA_GOALS_EXEMPT
	if wants_cuda or environment.get('CUDA_ARCHS') or environment.get('NVCC'):
		settings['TFI_ENABLE_CUDA'] = 'ON'
		if environment.get('CUDA_ARCHS'):
			settings['CMAKE_CUDA_ARCHITECTURES'] = ';'.join(environment['CUDA_ARCHS'].split())
		nvcc = environment.get('NVCC')
		if nvcc:
			resolved = shutil.which(os.path.expanduser(nvcc))
			if not resolved:
				fail(f'NVCC={nvcc} was not found')
			settings['CMAKE_CUDA_COMPILER'] = resolved
	return settings


def configure(build, settings):
	cache = read_cache(build)
	if 'CMAKE_CUDA_COMPILER' in settings and cache.get('CMAKE_CUDA_COMPILER') \
			and cache['CMAKE_CUDA_COMPILER'] != settings['CMAKE_CUDA_COMPILER'] \
			and cache.get('TFI_ENABLE_CUDA') == 'ON':
		fail(f'{build} already uses CUDA compiler {cache["CMAKE_CUDA_COMPILER"]}; '
			'remove the build directory or set TFI_BUILD_DIR to change it')
	if not cache and 'TFI_ENABLE_CUDA' in settings and 'CMAKE_CUDA_COMPILER' not in settings \
			and not shutil.which('nvcc'):
		fail('nvcc is not on PATH; pass NVCC=/path/to/nvcc')
	stale = not (build / 'tfi-forward.json').exists() or any(
		cache.get(key) != value for key, value in settings.items())
	if not stale:
		return
	command = ['cmake', '-S', str(REPOSITORY), '-B', str(build)]
	command += [f'-D{key}={value}' for key, value in settings.items()]
	if run(command):
		fail('CMake configuration failed')


def resolve(entry, directory, goal):
	"""Return (targets to build, run tests?) for one goal."""
	if goal in entry['aliases']:
		return entry['aliases'][goal], False
	if goal in entry['tests']:
		test = entry['tests'][goal]
		if test['skipped']:
			fail(f'{directory}:{goal} cannot run here: ' + '; '.join(test['skipped']))
		return test['build'], True
	if goal in entry['targets']:
		return [goal], False
	if goal in entry['skipped']:
		fail(f'{directory}:{goal} is not configured: requires {entry["skipped"][goal]}')
	fail(f'{directory} has no target {goal!r}; see {directory}/CMakeLists.txt')


def main():
	if len(sys.argv) != 3:
		fail('usage: forward.py DIRECTORY GOAL')
	directory, goal = sys.argv[1], sys.argv[2]
	build = Path(os.environ.get('TFI_BUILD_DIR') or REPOSITORY / 'build').resolve()
	build.mkdir(parents=True, exist_ok=True)
	with open(build / '.forward.lock', 'w') as lock:
		fcntl.flock(lock, fcntl.LOCK_EX)
		configure(build, desired_settings(directory, goal))
		if not (build / 'tfi-forward.json').exists():
			return 0
		manifest = json.loads((build / 'tfi-forward.json').read_text())
		entry = manifest.get(directory)
		if entry is None:
			fail(f'{directory} is not part of the CMake project')
		if goal == 'clean':
			for name in entry['outputs']:
				if dry_run():
					print(f'+ rm -f {directory}/{name}')
				else:
					(REPOSITORY / directory / name).unlink(missing_ok=True)
			return 0
		targets, test = resolve(entry, directory, goal)
		missing = {name: entry['skipped'][name] for name in targets if name in entry['skipped']}
		if missing:
			fail(f'{directory}:{goal} needs unconfigured targets: '
				+ '; '.join(f'{name} requires {reason}' for name, reason in missing.items()))
		if goal in entry.get('goals', {}):
			targets = [entry['goals'][goal]]
		if targets:
			# Under `make -j`, the recursive (+) recipe passes make's jobserver to
			# the generated build, which then shares the caller's job slots.
			jobs = [] if '--jobserver' in os.environ.get('MAKEFLAGS', '') \
				else ['--parallel', str(parallel_level())]
			status = run(['cmake', '--build', str(build)] + jobs + ['--target'] + list(targets))
			if status:
				return status
	if not test:
		return 0
	label = '^' + re.sub(r'([][.^$|()*+?\\])', r'\\\1', f'make={goal}') + '$'
	command = ['ctest', '--test-dir', str(build / directory), '--output-on-failure', '-L', label]
	if re.search(r'(?:^|\s)-j\s*\d+', os.environ.get('MAKEFLAGS', '')):
		command += ['-j', str(parallel_level())]
	return run(command)


if __name__ == '__main__':
	sys.exit(main())
