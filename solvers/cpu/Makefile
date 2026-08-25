CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS ?= -Iinclude -I../../include
LDLIBS ?= -lstdc++fs

PETSC_DIR ?=
PETSC_ARCH ?=
PETSC_TARGETS := petsc petsc-test generic_transport_test navier_stokes_test iga_solve iga_assembly_smoke iga_transport iga_navier_stokes iga_mesh_check

ifneq ($(strip $(filter $(PETSC_TARGETS),$(MAKECMDGOALS))),)
ifeq ($(strip $(PETSC_DIR)),)
$(error PETSC_DIR must be set when building PETSc targets)
endif
PETSC_PREFIX := $(PETSC_DIR)$(if $(strip $(PETSC_ARCH)),/$(PETSC_ARCH))
PETSC_VARIABLES ?= $(PETSC_PREFIX)/lib/petsc/conf/petscvariables
include $(PETSC_VARIABLES)
endif

.PHONY: all petsc petsc-test test clean

all: iga_pack iga_inspect iga_case_check iga_config_check iga_flow_validate iga_womersley_reference

petsc: iga_assembly_smoke iga_solve iga_transport iga_navier_stokes iga_mesh_check

iga_pack: src/iga_pack.cpp include/IgaDatabase.hpp ../../include/IgaPreprocessCache.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

iga_inspect: src/iga_inspect.cpp include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

iga_case_check: src/iga_case_check.cpp include/IgaDatabase.hpp include/CaseInput.hpp include/GenericCaseInput.hpp ../../include/CaseConfig.hpp ../../include/SimulationConfig.hpp ../../include/TemporalFunction.hpp ../../include/OutletModel.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
iga_config_check: src/iga_config_check.cpp ../../include/SimulationConfig.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

iga_flow_validate: src/iga_flow_validate.cpp include/BoundaryFlow.hpp include/CaseInput.hpp include/IgaDatabase.hpp include/WomersleyReference.hpp ../../include/VelocitySeries.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

iga_womersley_reference: src/iga_womersley_reference.cpp include/WomersleyReference.hpp include/IgaDatabase.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)


case_config_test: tests/test_case_config.cpp include/CaseInput.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

simulation_config_test: tests/test_simulation_config.cpp include/GenericCaseInput.hpp ../../include/SimulationConfig.hpp ../../include/TemporalFunction.hpp ../../include/VelocitySeries.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

generic_transport_test: tests/test_generic_transport.cpp include/GenericTransportElement.hpp include/TransportElement.hpp include/CaseInput.hpp ../../include/SimulationConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

navier_stokes_test: tests/test_navier_stokes.cpp include/NavierStokesElement.hpp include/BoundaryFlow.hpp include/TransportElement.hpp include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

test: case_config_test simulation_config_test flow_checkpoint_test transport_checkpoint_test outlet_model_test boundary_flow_test womersley_reference_test
	./case_config_test
	./simulation_config_test
	./flow_checkpoint_test
	./transport_checkpoint_test
	./outlet_model_test
	./boundary_flow_test
	./womersley_reference_test

flow_checkpoint_test: tests/test_flow_checkpoint.cpp ../../include/FlowCheckpoint.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

transport_checkpoint_test: tests/test_transport_checkpoint.cpp ../../include/TransportCheckpoint.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

outlet_model_test: tests/test_outlet_model.cpp ../../include/OutletModel.hpp ../../include/OutletCheckpoint.hpp ../../include/FlowCheckpoint.hpp ../../include/SimulationConfig.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

boundary_flow_test: tests/test_boundary_flow.cpp include/BoundaryFlow.hpp include/OutletFlow.hpp include/PressureTractionFlow.hpp include/IgaDatabase.hpp ../../include/OutletModel.hpp ../../include/PressureTraction.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

womersley_reference_test: tests/test_womersley_reference.cpp include/WomersleyReference.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

petsc-test: generic_transport_test navier_stokes_test
	./generic_transport_test
	./navier_stokes_test

iga_assembly_smoke: src/iga_assembly_smoke.cpp include/IgaDatabase.hpp include/OwnedRowAssembler.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

iga_solve: src/iga_solve.cpp include/IgaDatabase.hpp include/OwnedRowAssembler.hpp include/CaseInput.hpp include/GenericCaseInput.hpp include/GenericTransportElement.hpp ../../include/SimulationConfig.hpp ../../include/TemporalFunction.hpp ../../include/VelocitySeries.hpp ../../include/FlowCheckpoint.hpp ../../include/TransportCheckpoint.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

iga_transport: src/iga_transport.cpp include/IgaDatabase.hpp include/OwnedRowAssembler.hpp include/CaseInput.hpp include/TransportElement.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

iga_navier_stokes: src/iga_navier_stokes.cpp include/IgaDatabase.hpp include/OwnedRowAssembler.hpp include/CaseInput.hpp include/GenericCaseInput.hpp include/TransportElement.hpp include/NavierStokesElement.hpp include/BoundaryFlow.hpp include/PressureTractionFlow.hpp ../../include/CaseConfig.hpp ../../include/SimulationConfig.hpp ../../include/TemporalFunction.hpp ../../include/FlowCheckpoint.hpp ../../include/OutletModel.hpp ../../include/OutletCheckpoint.hpp ../../include/PressureTraction.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

iga_mesh_check: src/iga_mesh_check.cpp include/IgaDatabase.hpp include/TransportElement.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

clean:
	$(RM) iga_pack iga_inspect iga_case_check iga_config_check iga_flow_validate iga_womersley_reference case_config_test simulation_config_test flow_checkpoint_test transport_checkpoint_test outlet_model_test boundary_flow_test womersley_reference_test generic_transport_test navier_stokes_test iga_solve iga_assembly_smoke iga_transport iga_navier_stokes iga_mesh_check
