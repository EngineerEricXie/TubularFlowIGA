#!/usr/bin/env bash

set -u

component=${1:-all}
failures=0

pass() { printf 'ok      %s\n' "$1"; }
fail() { printf 'missing %s\n' "$1"; failures=$((failures+1)); }
note() { printf 'note    %s\n' "$1"; }

have_command() {
	if command -v "$1" >/dev/null 2>&1; then
		pass "$1: $(command -v "$1")"
	else
		fail "$1"
	fi
}

check_base() {
	have_command make
	have_command "${CXX:-g++}"
}

check_preprocessing() {
	local eigen_dir=${EIGEN_DIR:-/usr/include/eigen3}
	if [[ -f "$eigen_dir/Eigen/Core" ]]; then
		pass "Eigen headers: $eigen_dir"
	else
		fail "Eigen headers (set EIGEN_DIR to the directory containing Eigen/)"
	fi
	have_command mpmetis
}

check_cpu() {
	have_command "${MPICXX:-mpicxx}"
	local petsc_dir=${PETSC_DIR:-}
	local petsc_arch=${PETSC_ARCH:-}
	if [[ -z "$petsc_dir" ]]; then
		fail "PETSC_DIR"
		return
	fi
	local prefix=$petsc_dir
	if [[ -n "$petsc_arch" ]]; then prefix=$prefix/$petsc_arch; fi
	if [[ -f "$prefix/lib/petsc/conf/petscvariables" ]]; then
		pass "PETSc configuration: $prefix"
	else
		fail "PETSc configuration: $prefix/lib/petsc/conf/petscvariables"
	fi
}

check_cuda() {
	have_command "${NVCC:-nvcc}"
}

case "$component" in
	preprocessing)
		check_base
		check_preprocessing
		;;
	cpu)
		check_base
		check_preprocessing
		check_cpu
		;;
	cuda)
		check_base
		check_preprocessing
		check_cuda
		;;
	all)
		check_base
		check_preprocessing
		check_cpu
		check_cuda
		;;
	*)
		printf 'usage: %s preprocessing|cpu|cuda|all\n' "$0" >&2
		exit 2
		;;
esac

if (( failures > 0 )); then
	note "$failures required item(s) missing for component '$component'"
	exit 1
fi
note "dependencies for component '$component' are available"
