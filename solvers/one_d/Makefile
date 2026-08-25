CXX ?= g++
PETSC_CXX ?= mpicxx
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS ?= -Iinclude -I../../include
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

petsc: iga_1d

iga_1d: src/iga_1d.cpp include/OneDNetwork.hpp include/OneDFlow.hpp include/OneDTransport.hpp include/OneDImplicit.hpp include/OneDOutput.hpp include/OneDCheckpoint.hpp ../../include/OneDConfig.hpp ../../include/SimulationConfig.hpp ../../include/TemporalFunction.hpp ../../include/CaseConfig.hpp
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)

one_d_core_test: tests/test_one_d_core.cpp include/OneDNetwork.hpp include/OneDFlow.hpp include/OneDTransport.hpp ../../include/OneDConfig.hpp ../../include/SimulationConfig.hpp ../../include/TemporalFunction.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $< -o $@

one_d_petsc_test: tests/test_one_d_petsc.cpp include/OneDNetwork.hpp include/OneDFlow.hpp include/OneDTransport.hpp include/OneDImplicit.hpp include/OneDOutput.hpp include/OneDCheckpoint.hpp ../../include/OneDConfig.hpp ../../include/SimulationConfig.hpp ../../include/TemporalFunction.hpp ../../include/CaseConfig.hpp
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)

core-test: one_d_core_test
	./one_d_core_test

test: one_d_core_test one_d_petsc_test
	./one_d_core_test
	./one_d_petsc_test
	$(MPIEXEC) -np $(MPI_TEST_RANKS) ./one_d_petsc_test

clean:
	$(RM) one_d_core_test one_d_petsc_test iga_1d
