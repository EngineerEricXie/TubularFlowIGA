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

check_hdf5() {
	if [[ -n ${HDF5_CFLAGS:-} && -n ${HDF5_LIBS:-} ]]; then
		pass "HDF5 compiler/link overrides"
		return
	fi
	for prefix in "${CONDA_PREFIX:-}" "${CONDA_ROOT:-}"; do
		if [[ -n $prefix && -f $prefix/include/hdf5.h \
			&& ( -f $prefix/lib/libhdf5.so || -f $prefix/lib/libhdf5.dylib ) ]]; then
			pass "HDF5 Conda prefix: $prefix"
			return
		fi
	done
	if ! command -v pkg-config >/dev/null 2>&1; then
		fail "pkg-config (needed to locate HDF5)"
		return
	fi
	if pkg-config --exists hdf5; then
		pass "HDF5: $(pkg-config --modversion hdf5)"
	else
		fail "HDF5 development package (pkg-config hdf5 or HDF5_CFLAGS/HDF5_LIBS)"
	fi
}

check_cuda() {
	have_command "${NVCC:-nvcc}"
}

check_one_d() {
	have_command "${MPICXX:-mpicxx}"
	if "${CXX:-g++}" -std=c++17 -dM -E -x c++ /dev/null 2>/dev/null \
		| grep -Eq '^#define __cplusplus 201703L$'; then
		pass "C++17 compiler mode"
	else
		fail "C++17 compiler support"
	fi
	if command -v pkg-config >/dev/null 2>&1 && pkg-config --exists PETSc; then
		pass "PETSc pkg-config: $(pkg-config --modversion PETSc)"
	else
		check_cpu
	fi
	if "${CXX:-g++}" -fopenmp -dM -E -x c++ /dev/null 2>/dev/null \
		| grep -q '^#define _OPENMP '; then
		pass "OpenMP compiler flag: -fopenmp"
	else
		fail "OpenMP compiler support"
	fi
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
		check_hdf5
		;;
	cuda)
		check_base
		check_preprocessing
		check_cuda
		check_hdf5
		;;
	one-d)
		check_base
		check_one_d
		;;
	all)
		check_base
		check_preprocessing
		check_cpu
		check_cuda
		check_hdf5
		;;
	*)
		printf 'usage: %s preprocessing|cpu|one-d|cuda|all\n' "$0" >&2
		exit 2
		;;
esac

if (( failures > 0 )); then
	note "$failures required item(s) missing for component '$component'"
	exit 1
fi
note "dependencies for component '$component' are available"
