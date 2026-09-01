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

.PHONY: all test petsc petsc-test convergence-test convergence-fixture-test convergence-axial-diagnostic convergence-isotropic-diagnostic convergence-length-diagnostic convergence-bulk-diagnostic temporal-convergence-test temporal-fixture-test clean

all: explicit_coupling_test simulation_graph_test

explicit_coupling_test: tests/test_explicit_coupling.cpp ../../include/ExplicitOneDThreeDCoupling.hpp ../../include/StrongOneDThreeDCoupling.hpp ../../include/AitkenRelaxation.hpp ../../include/CouplingPort.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

simulation_graph_test: tests/test_simulation_graph.cpp ../../include/SimulationGraph.hpp ../../include/CouplingEdge.hpp ../../include/CouplingPort.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

test: explicit_coupling_test simulation_graph_test
	./explicit_coupling_test
	./simulation_graph_test

petsc: iga_1d_3d_explicit

petsc-test: iga_1d_3d_explicit explicit_coupling_smoke_test
	PATH="$(PETSC_PREFIX)/bin:$$PATH" ./explicit_coupling_smoke_test

# Opt-in fixture/convergence prerequisite.  The full numerical ladder is kept
# out of the default smoke target because it is intentionally multi-minute.
convergence-test: iga_1d_3d_explicit straight_vessel_convergence_test
	PATH="$(PETSC_PREFIX)/bin:$$PATH" ./straight_vessel_convergence_test

convergence-fixture-test: straight_vessel_convergence_test
	TUBULARFLOWIGA_FIXTURE_ONLY=1 ./straight_vessel_convergence_test

convergence-axial-diagnostic: iga_1d_3d_explicit straight_vessel_convergence_test
	PATH="$(PETSC_PREFIX)/bin:$$PATH" TUBULARFLOWIGA_AXIAL_ONLY=1 ./straight_vessel_convergence_test

convergence-isotropic-diagnostic: iga_1d_3d_explicit straight_vessel_convergence_test
	PATH="$(PETSC_PREFIX)/bin:$$PATH" TUBULARFLOWIGA_ISOTROPIC_ONLY=1 ./straight_vessel_convergence_test

convergence-length-diagnostic: iga_1d_3d_explicit straight_vessel_convergence_test
	PATH="$(PETSC_PREFIX)/bin:$$PATH" TUBULARFLOWIGA_LENGTH_ONLY=1 ./straight_vessel_convergence_test

convergence-bulk-diagnostic: iga_1d_3d_explicit straight_vessel_convergence_test
	PATH="$(PETSC_PREFIX)/bin:$$PATH" TUBULARFLOWIGA_BULK_ONLY=1 ./straight_vessel_convergence_test

# Opt-in temporal verification.  This deliberately performs several pulsatile
# PETSc solves and remains outside both the default and smoke test targets.
temporal-convergence-test: iga_1d_3d_explicit pulsatile_straight_vessel_temporal_test
	PATH="$(PETSC_PREFIX)/bin:$$PATH" ./pulsatile_straight_vessel_temporal_test

temporal-fixture-test: pulsatile_straight_vessel_temporal_test
	TUBULARFLOWIGA_TEMPORAL_FIXTURE_ONLY=1 ./pulsatile_straight_vessel_temporal_test

iga_1d_3d_explicit: src/iga_1d_3d_explicit.cpp ../../include/ExplicitOneDThreeDCoupling.hpp ../../include/StrongOneDThreeDCoupling.hpp ../../include/AitkenRelaxation.hpp ../../include/SimulationGraph.hpp ../../include/CouplingEdge.hpp ../../include/ThreeDFlowCoupling.hpp ../../include/ThreeDVcaCoupling.hpp ../../include/SimulationConfig.hpp ../../include/CouplingConfig.hpp ../../include/CouplingPort.hpp ../cpu/include/BoundarySupport.hpp ../cpu/include/IgaDatabase.hpp ../cpu/include/TransientFlowRuntime.hpp ../cpu/include/TransportElement.hpp ../one_d/include/OneDImplicit.hpp ../one_d/include/OneDRuntime.hpp
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)

explicit_coupling_smoke_test: tests/test_explicit_coupling_smoke.cpp iga_1d_3d_explicit ../cpu/include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

straight_vessel_convergence_test: tests/test_straight_vessel_convergence.cpp tests/CouplingFixture.hpp ../cpu/include/BoundarySupport.hpp ../cpu/include/IgaDatabase.hpp ../one_d/include/OneDRuntime.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

pulsatile_straight_vessel_temporal_test: tests/test_pulsatile_straight_vessel_temporal.cpp tests/CouplingFixture.hpp ../cpu/include/BoundarySupport.hpp ../cpu/include/IgaDatabase.hpp ../one_d/include/OneDRuntime.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

clean:
	$(RM) explicit_coupling_test simulation_graph_test explicit_coupling_smoke_test straight_vessel_convergence_test pulsatile_straight_vessel_temporal_test iga_1d_3d_explicit
