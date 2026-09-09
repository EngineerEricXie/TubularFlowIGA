CXX ?= g++
PETSC_CXX ?= mpicxx
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS ?= -Iinclude -I../../include -I../cpu/include
OPENMP_FLAGS ?= -fopenmp
MPIEXEC ?= mpiexec
MPI_TEST_RANKS ?= 2
PETSC_DIR ?=
PETSC_ARCH ?=

ifeq ($(strip $(PETSC_DIR)),)
PETSC_CFLAGS ?= $(shell pkg-config --cflags PETSc 2>/dev/null)
PETSC_LIBS ?= $(shell pkg-config --libs PETSc 2>/dev/null)
else
PETSC_PREFIX := $(PETSC_DIR)$(if $(strip $(PETSC_ARCH)),/$(PETSC_ARCH))
PETSC_VARIABLES ?= $(PETSC_PREFIX)/lib/petsc/conf/petscvariables
include $(PETSC_VARIABLES)
PETSC_CFLAGS ?= $(PETSC_CC_INCLUDES)
PETSC_LIBS ?= $(PETSC_LIB)
endif

.PHONY: all petsc core-test test clean

all: one_d_core_test

one_d_filename_failure_test: tests/test_one_d_filename_failure.cpp src/iga_1d.cpp ../cpu/tests/StringStreamFailure.hpp $(wildcard include/*.hpp ../cpu/include/*.hpp ../../include/*.hpp)
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)

iga_1d one_d_core_test one_d_coupling_test one_d_runtime_test one_d_petsc_test: ../../include/CheckedText.hpp

petsc: iga_1d

iga_1d: ../cpu/include/CollectiveAssetInput.hpp ../cpu/include/Sha256.hpp

iga_1d one_d_coupling_test: ../../include/CouplingHistory.hpp

iga_1d: src/iga_1d.cpp ../cpu/include/CollectivePetscOptions.hpp ../cpu/include/CollectiveFailure.hpp ../../include/PhaseProfile.hpp include/OneDRuntime.hpp include/OneDCoupling.hpp include/OneDNetwork.hpp include/OneDFlow.hpp include/OneDTransport.hpp include/OneDImplicit.hpp include/OneDOutput.hpp include/OneDCheckpoint.hpp ../../include/CouplingPort.hpp ../../include/CouplingConfig.hpp ../../include/CouplingReplay.hpp ../../include/OxygenCapacity.hpp ../../include/VascularCoupling.hpp ../../include/OneDConfig.hpp ../../include/RadiusAnnotatedObj.hpp ../../include/SkeletonOutput.hpp ../../include/SimulationConfig.hpp ../../include/TemporalFunction.hpp ../../include/CaseConfig.hpp
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)

one_d_core_test: tests/test_one_d_core.cpp include/OneDNetwork.hpp include/OneDFlow.hpp include/OneDTransport.hpp ../../include/CouplingConfig.hpp ../../include/OxygenCapacity.hpp ../../include/VascularCoupling.hpp ../../include/OneDConfig.hpp ../../include/RadiusAnnotatedObj.hpp ../../include/SkeletonOutput.hpp ../../include/SimulationConfig.hpp ../../include/TemporalFunction.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $< -o $@

one_d_coupling_test: tests/test_one_d_coupling.cpp include/OneDCoupling.hpp include/OneDNetwork.hpp include/OneDFlow.hpp include/OneDTransport.hpp ../../include/CouplingConfig.hpp ../../include/CouplingReplay.hpp ../../include/OxygenCapacity.hpp ../../include/VascularCoupling.hpp ../../include/OneDConfig.hpp ../../include/RadiusAnnotatedObj.hpp ../../include/SkeletonOutput.hpp ../../include/SimulationConfig.hpp ../../include/TemporalFunction.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $< -o $@

one_d_runtime_test: tests/test_one_d_runtime.cpp include/OneDRuntime.hpp include/OneDCoupling.hpp include/OneDNetwork.hpp include/OneDFlow.hpp include/OneDTransport.hpp ../../include/OneDFlowDomainAdapter.hpp ../../include/FlowDomainPortMetadata.hpp ../../include/CoupledDomainRuntime.hpp ../../include/CouplingPort.hpp ../../include/CouplingConfig.hpp ../../include/OxygenCapacity.hpp ../../include/VascularCoupling.hpp ../../include/OneDConfig.hpp ../../include/RadiusAnnotatedObj.hpp ../../include/SkeletonOutput.hpp ../../include/SimulationConfig.hpp ../../include/TemporalFunction.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $< -o $@

one_d_petsc_test: tests/test_one_d_petsc.cpp ../cpu/include/CollectiveFailure.hpp ../../include/PhaseProfile.hpp include/OneDNetwork.hpp include/OneDFlow.hpp include/OneDTransport.hpp include/OneDImplicit.hpp include/OneDOutput.hpp include/OneDCheckpoint.hpp ../../include/OneDConfig.hpp ../../include/RadiusAnnotatedObj.hpp ../../include/SkeletonOutput.hpp ../../include/SimulationConfig.hpp ../../include/TemporalFunction.hpp ../../include/CaseConfig.hpp
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)

core-test: one_d_core_test one_d_coupling_test one_d_runtime_test
	./one_d_core_test
	./one_d_coupling_test
	./one_d_runtime_test

one_d_subcommunicator_test: tests/test_one_d_subcommunicator.cpp include/OneDImplicit.hpp include/OneDFlow.hpp include/OneDNetwork.hpp
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)

.PHONY: subcommunicator-test
one_d_checkpoint_format_test: tests/test_one_d_checkpoint_format.cpp
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)

one_d_trial_failure_test: tests/test_one_d_trial_failure.cpp include/OneDRuntime.hpp include/OneDTransport.hpp include/OneDFlow.hpp include/OneDCheckpoint.hpp include/OneDImplicit.hpp ../cpu/include/CollectiveFailure.hpp
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)

subcommunicator-test: one_d_subcommunicator_test
	OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps' timeout 90 $(MPIEXEC) -np 4 ./one_d_subcommunicator_test ../../examples/one_d/rigid_straight/skeleton_initial.swc

test: one_d_core_test one_d_coupling_test one_d_runtime_test one_d_petsc_test
	./one_d_core_test
	./one_d_coupling_test
	./one_d_runtime_test
	./one_d_petsc_test
	$(MPIEXEC) -np $(MPI_TEST_RANKS) ./one_d_petsc_test

clean:
	$(RM) one_d_adapter_failure_test
	$(RM) one_d_implicit_failure_test
	$(RM) one_d_trial_failure_test
	$(RM) one_d_checkpoint_format_test
	$(RM) one_d_filename_failure_test
	$(RM) one_d_subcommunicator_test
	$(RM) one_d_core_test one_d_coupling_test one_d_runtime_test one_d_petsc_test iga_1d

one_d_implicit_failure_test: tests/test_one_d_implicit_failure.cpp include/OneDImplicit.hpp include/OneDCheckpoint.hpp include/OneDFlow.hpp include/OneDNetwork.hpp
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)

iga_1d one_d_petsc_test one_d_subcommunicator_test one_d_trial_failure_test one_d_implicit_failure_test: include/OneDPetscSupport.hpp ../cpu/include/CollectiveFailure.hpp ../cpu/include/PetscReadArray.hpp
iga_1d one_d_petsc_test one_d_subcommunicator_test one_d_trial_failure_test one_d_implicit_failure_test one_d_adapter_failure_test: ../cpu/include/ExecutionResources.hpp

.PHONY: implicit-failure-test
implicit-failure-test: one_d_implicit_failure_test
	OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps' timeout 180 $(MPIEXEC) -np 3 ./one_d_implicit_failure_test ../../examples/one_d/rigid_straight/skeleton_initial.swc

one_d_adapter_failure_test: tests/test_one_d_adapter_failure.cpp ../../include/OneDFlowDomainAdapter.hpp include/OneDRuntime.hpp include/OneDImplicit.hpp include/OneDPetscSupport.hpp include/OneDCheckpoint.hpp include/OneDTransport.hpp ../cpu/include/CollectiveFailure.hpp ../cpu/include/PetscReadArray.hpp
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)
