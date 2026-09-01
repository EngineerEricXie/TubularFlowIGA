CXX ?= g++
PETSC_CXX ?= mpicxx
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS ?= -I../cpu/include -I../one_d/include -I../../include
OPENMP_FLAGS ?= -fopenmp
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

.PHONY: all test petsc clean

all: explicit_coupling_test

explicit_coupling_test: tests/test_explicit_coupling.cpp ../../include/ExplicitOneDThreeDCoupling.hpp ../../include/CouplingPort.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

test: explicit_coupling_test
	./explicit_coupling_test

petsc: iga_1d_3d_explicit

iga_1d_3d_explicit: src/iga_1d_3d_explicit.cpp ../../include/ExplicitOneDThreeDCoupling.hpp ../../include/ThreeDFlowCoupling.hpp ../../include/ThreeDVcaCoupling.hpp ../../include/SimulationConfig.hpp ../../include/CouplingConfig.hpp ../../include/CouplingPort.hpp ../cpu/include/BoundarySupport.hpp ../cpu/include/IgaDatabase.hpp ../cpu/include/TransientFlowRuntime.hpp ../cpu/include/TransportElement.hpp ../one_d/include/OneDImplicit.hpp ../one_d/include/OneDRuntime.hpp
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)

clean:
	$(RM) explicit_coupling_test iga_1d_3d_explicit
