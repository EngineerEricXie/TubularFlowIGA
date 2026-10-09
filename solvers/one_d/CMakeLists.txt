# Native 0D/1D solvers and their tests.
#
# Generated once from the historical Makefile; edit freely. Helpers are in
# cmake/TubularFlowHelpers.cmake. `make -C <this directory> <target>` forwards
# here, so every former target name and test name keeps working.
set(TFI_DIRECTORY_INCLUDES include ../../include ../cpu/include)
set(TFI_DIRECTORY_OPTIONS -Wall -Wextra -Wpedantic)

# Tools and runtimes.
tfi_executable(iga_0d SOURCES src/iga_1d.cpp PETSC HDF5 OPENMP)
tfi_executable(iga_1d SOURCES src/iga_1d.cpp PETSC HDF5 OPENMP)

# Test executables.
tfi_executable(one_d_accepted_checkpoint_test
	SOURCES tests/test_one_d_accepted_checkpoint.cpp
	MPI OPENMP)
tfi_executable(one_d_adapter_failure_test SOURCES tests/test_one_d_adapter_failure.cpp PETSC HDF5 OPENMP)
tfi_executable(one_d_checkpoint_format_test SOURCES tests/test_one_d_checkpoint_format.cpp PETSC)
tfi_executable(one_d_core_test SOURCES tests/test_one_d_core.cpp MPI OPENMP)
tfi_executable(one_d_coupling_test SOURCES tests/test_one_d_coupling.cpp MPI OPENMP)
tfi_executable(one_d_filename_failure_test
	SOURCES tests/test_one_d_filename_failure.cpp
	PETSC HDF5 OPENMP)
tfi_executable(one_d_implicit_failure_test
	SOURCES tests/test_one_d_implicit_failure.cpp
	PETSC HDF5 OPENMP)
tfi_executable(one_d_output_test SOURCES tests/test_one_d_output.cpp MPI HDF5 OPENMP)
tfi_executable(one_d_petsc_test SOURCES tests/test_one_d_petsc.cpp PETSC HDF5 OPENMP)
tfi_executable(one_d_runtime_test SOURCES tests/test_one_d_runtime.cpp MPI OPENMP)
tfi_executable(one_d_solver_options_test
	SOURCES tests/test_one_d_solver_options.cpp
	PETSC HDF5 OPENMP)
tfi_executable(one_d_stdout_failure_test
	SOURCES ../cpu/tests/test_tool_stdout_failure.cpp
	PETSC HDF5 OPENMP
	DEFINITIONS IGA_TEST_ONE_D_STDOUT)
tfi_executable(one_d_subcommunicator_test SOURCES tests/test_one_d_subcommunicator.cpp PETSC OPENMP)
tfi_executable(one_d_transport_junction_test SOURCES tests/test_one_d_transport_junction.cpp MPI)
tfi_executable(one_d_trial_failure_test SOURCES tests/test_one_d_trial_failure.cpp PETSC HDF5 OPENMP)

# Build-only aggregates of the historical Makefile.
tfi_make_alias(all TARGETS one_d_core_test)
tfi_make_alias(petsc TARGETS iga_1d)
tfi_make_alias(zero-d TARGETS iga_0d)

# Tests run from this directory with the original shell commands. MAKE
# names the historical make targets that ran them; LABELS select tiers.
tfi_test(one_d_core_test COMMAND [==[./one_d_core_test]==]
	REQUIRES one_d_core_test MAKE core-test test LABELS unit)
tfi_test(one_d_coupling_test COMMAND [==[./one_d_coupling_test]==]
	REQUIRES one_d_coupling_test MAKE core-test test LABELS unit)
tfi_test(one_d_runtime_test COMMAND [==[./one_d_runtime_test]==]
	REQUIRES one_d_runtime_test MAKE core-test test LABELS unit)
tfi_test(one_d_transport_junction_test COMMAND [==[./one_d_transport_junction_test]==]
	REQUIRES one_d_transport_junction_test MAKE core-test test LABELS unit)
tfi_test(one_d_implicit_failure_test COMMAND [==[OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps' timeout 180 @TFI_MPIEXEC@ -np 3 ./one_d_implicit_failure_test ../../examples/one_d/rigid_straight/skeleton_initial.swc]==]
	REQUIRES one_d_implicit_failure_test MAKE implicit-failure-test LABELS petsc mpi)
tfi_test(one_d_subcommunicator_test COMMAND [==[OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps' timeout 90 @TFI_MPIEXEC@ -np 4 ./one_d_subcommunicator_test ../../examples/one_d/rigid_straight/skeleton_initial.swc]==]
	REQUIRES one_d_subcommunicator_test MAKE subcommunicator-test LABELS petsc mpi)
tfi_test(one_d_output_test COMMAND [==[./one_d_output_test]==]
	REQUIRES one_d_output_test MAKE test LABELS unit)
tfi_test(one_d_petsc_test COMMAND [==[./one_d_petsc_test]==]
	REQUIRES one_d_petsc_test MAKE test LABELS petsc)
tfi_test(one_d_petsc_test.mpi COMMAND [==[@TFI_MPIEXEC@ -np "${MPI_TEST_RANKS:-2}" ./one_d_petsc_test]==]
	REQUIRES one_d_petsc_test MAKE test LABELS petsc mpi)
tfi_test(vtkhdf-paraview-test COMMAND [==[(TUBULARFLOWIGA_KEEP_TEST_OUTPUT=1 ./one_d_output_test) && ("${PVPYTHON:-pvpython}" tests/validate_one_d_vtkhdf.py) && (./one_d_output_test)]==]
	REQUIRES one_d_output_test MAKE vtkhdf-paraview-test LABELS paraview)

# Coordinated-failure, checkpoint, and MPI tests that no make target ran before.
# Each MPI test needs exactly the rank count shown; "$TMPDIR/out" is a fresh output path.
tfi_test(one_d_accepted_checkpoint_test COMMAND [==[./one_d_accepted_checkpoint_test "$TMPDIR/out"]==]
	REQUIRES one_d_accepted_checkpoint_test MAKE test one-d-accepted-checkpoint-test LABELS unit)
tfi_test(one_d_solver_options_test COMMAND [==[./one_d_solver_options_test ../../examples/one_d/rigid_straight/skeleton_initial.swc && @TFI_MPIEXEC@ -np 2 ./one_d_solver_options_test ../../examples/one_d/rigid_straight/skeleton_initial.swc]==]
	REQUIRES one_d_solver_options_test MAKE test one-d-solver-options-test LABELS petsc mpi)

# Multi-rank output-name regression driver from scripts/hpc.
tfi_test(hpc-one-d-filename-regression COMMAND [==[python3 ../../scripts/hpc/hpc_one_d_filename_regression.py --output-dir "$TMPDIR/out"]==]
	REQUIRES one_d_filename_failure_test MAKE hpc-regression-test LABELS petsc mpi)

# Coordinated 1D adapter and trial failures on the multispecies example network.
foreach(test one_d_adapter_failure_test one_d_trial_failure_test)
	tfi_test(${test} COMMAND "@TFI_MPIEXEC@ -np 3 ./${test} ../../examples/one_d/multispecies_physiology"
		REQUIRES ${test} MAKE test LABELS petsc mpi)
endforeach()
