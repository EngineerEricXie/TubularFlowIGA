CXX ?= g++
PETSC_CXX ?= mpicxx
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS ?= -I../cpu/include -I../one_d/include -I../../include
SPECIES_COUPLING_HEADER := ../../include/SpeciesCoupling.hpp
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

.PHONY: all test zero-d-flow-test simulation-graph-test surface-contracts-test fsi-runtime-contracts-test petsc petsc-test bifurcation-test multidomain-test immersed-case-factory-test phase6-aneurysm-closure-test phase9-multiscale-closure-test convergence-test convergence-fixture-test convergence-axial-diagnostic convergence-isotropic-diagnostic convergence-length-diagnostic convergence-bulk-diagnostic temporal-convergence-test temporal-fixture-test clean

all: explicit_coupling_test simulation_graph_test pressure_flow_executor_test species_pressure_flow_executor_test multidomain_config_test

iga_1d_3d_explicit iga_1d_3d_bifurcation iga_multidomain_flow multidomain_config_test: ../../include/CheckedText.hpp

collective_asset_input_test: tests/test_collective_asset_input.cpp ../cpu/include/CollectiveAssetInput.hpp ../cpu/include/CollectiveFailure.hpp ../cpu/include/Sha256.hpp ../../include/PhaseProfile.hpp
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG -DIGA_ASSET_INPUT_TESTING $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)

three_d_staged_failure_test: tests/test_three_d_staged_failure.cpp tests/test_multidomain_flow_smoke.cpp tests/test_bifurcation_coupling_smoke.cpp $(wildcard ../../include/*.hpp ../cpu/include/*.hpp)
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) -DIGA_TRANSPORT_RUNTIME_TESTING $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)

three_d_trial_failure_test: tests/test_three_d_trial_failure.cpp tests/test_bifurcation_coupling_smoke.cpp $(wildcard ../../include/*.hpp ../cpu/include/*.hpp)
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)

three_d_port_failure_test: tests/test_three_d_port_failure.cpp tests/test_bifurcation_coupling_smoke.cpp $(wildcard ../../include/*.hpp ../cpu/include/*.hpp)
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) -DIGA_FLOW_RUNTIME_TESTING $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)

iga_1d_3d_bifurcation iga_multidomain_flow: ../../include/MultidomainRunner.hpp
iga_1d_3d_explicit iga_1d_3d_bifurcation iga_multidomain_flow: ../cpu/include/PetscPhaseProfile.hpp ../../include/PhaseProfile.hpp
iga_1d_3d_explicit iga_1d_3d_bifurcation iga_multidomain_flow: ../cpu/include/OwnedRowAssembler.hpp ../cpu/include/ParallelOwnershipValidation.hpp ../../include/ParallelOwnership.hpp
iga_1d_3d_explicit iga_1d_3d_bifurcation iga_multidomain_flow: ../cpu/include/CollectiveFailure.hpp ../cpu/include/PetscReadArray.hpp
iga_1d_3d_bifurcation iga_multidomain_flow: ../cpu/include/CollectivePetscOptions.hpp
iga_1d_3d_bifurcation iga_multidomain_flow: ../cpu/include/CollectiveAssetInput.hpp ../cpu/include/Sha256.hpp

zero_d_flow_test: tests/test_zero_d_flow.cpp ../../include/ZeroDFlowDomain.hpp ../../include/DomainRuntimeRegistry.hpp ../../include/CoupledDomainRuntime.hpp ../../include/SimulationGraph.hpp ../../include/CouplingPort.hpp ../cpu/include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

zero-d-flow-test: zero_d_flow_test
	./zero_d_flow_test

distributed_surface_contracts_test: tests/test_distributed_surface_contracts.cpp ../../include/DistributedSurfaceInterface.hpp ../../include/DynamicWeightedAitkenRelaxation.hpp ../../include/AitkenRelaxation.hpp ../cpu/include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

surface-contracts-test: distributed_surface_contracts_test
	./distributed_surface_contracts_test

fsi_runtime_contracts_test: tests/test_fsi_runtime_contracts.cpp ../../include/FsiCouplingEdge.hpp ../../include/FsiDomainRuntime.hpp ../../include/DistributedSurfaceInterface.hpp ../../include/CoupledDomainRuntime.hpp ../../include/CouplingEdge.hpp ../../include/CouplingPort.hpp ../cpu/include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

fsi-runtime-contracts-test: fsi_runtime_contracts_test
	./fsi_runtime_contracts_test

explicit_coupling_test: tests/test_explicit_coupling.cpp ../../include/ExplicitOneDThreeDCoupling.hpp ../../include/StrongOneDThreeDCoupling.hpp ../../include/AitkenRelaxation.hpp ../../include/CouplingPort.hpp $(SPECIES_COUPLING_HEADER)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

simulation_graph_test: tests/test_simulation_graph.cpp ../../include/SimulationGraph.hpp ../../include/ZeroDFlowDomain.hpp ../../include/FsiCouplingEdge.hpp ../../include/DistributedSurfaceInterface.hpp ../../include/CouplingEdge.hpp ../../include/CouplingPort.hpp ../../include/SpeciesCoupling.hpp ../cpu/include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

simulation-graph-test: simulation_graph_test
	./simulation_graph_test

pressure_flow_executor_test: tests/test_pressure_flow_executor.cpp ../../include/PressureFlowComponentExecutor.hpp ../../include/DomainRuntimeRegistry.hpp ../../include/CoupledDomainRuntime.hpp ../../include/SimulationGraph.hpp ../../include/CouplingEdge.hpp ../../include/CouplingPort.hpp $(SPECIES_COUPLING_HEADER)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

species_pressure_flow_executor_test: tests/test_species_pressure_flow_executor.cpp ../../include/SpeciesPressureFlowComponentExecutor.hpp ../../include/PressureFlowComponentExecutor.hpp ../../include/DomainRuntimeRegistry.hpp ../../include/CoupledDomainRuntime.hpp ../../include/SimulationGraph.hpp ../../include/CouplingEdge.hpp ../../include/CouplingPort.hpp $(SPECIES_COUPLING_HEADER)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

multidomain_config_test: tests/test_multidomain_config.cpp ../../include/BifurcationMultidomainCase.hpp ../../include/SequentialMultidomainCase.hpp ../../include/MultidomainConfig.hpp ../../include/ZeroDFlowDomain.hpp ../../include/FlowDomainPortMetadata.hpp ../../include/PressureFlowComponentExecutor.hpp ../../include/DomainRuntimeRegistry.hpp ../../include/CoupledDomainRuntime.hpp ../../include/SimulationGraph.hpp ../../include/CouplingEdge.hpp ../../include/CouplingPort.hpp ../../include/SpeciesCoupling.hpp ../../include/CaseConfig.hpp ../cpu/include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

test: explicit_coupling_test simulation_graph_test pressure_flow_executor_test species_pressure_flow_executor_test multidomain_config_test zero_d_flow_test distributed_surface_contracts_test fsi_runtime_contracts_test
	./explicit_coupling_test
	./simulation_graph_test
	./pressure_flow_executor_test
	./species_pressure_flow_executor_test
	./multidomain_config_test
	./zero_d_flow_test
	./distributed_surface_contracts_test
	./fsi_runtime_contracts_test

petsc: iga_1d_3d_explicit iga_1d_3d_bifurcation iga_multidomain_flow

petsc-test: iga_1d_3d_explicit explicit_coupling_smoke_test
	PATH="$(PETSC_PREFIX)/bin:$$PATH" ./explicit_coupling_smoke_test

bifurcation-test: iga_1d_3d_bifurcation bifurcation_coupling_smoke_test
	PATH="$(PETSC_PREFIX)/bin:$$PATH" ./bifurcation_coupling_smoke_test

multidomain-test: iga_multidomain_flow multidomain_flow_smoke_test
	PATH="$(PETSC_PREFIX)/bin:$$PATH" ./multidomain_flow_smoke_test

multidomain_subcommunicator_test: tests/test_multidomain_subcommunicator.cpp tests/test_bifurcation_coupling_smoke.cpp tests/test_multidomain_flow_smoke.cpp src/iga_1d_3d_bifurcation.cpp ../../include/MultidomainRunner.hpp $(wildcard ../../include/*.hpp ../cpu/include/*.hpp ../one_d/include/*.hpp)
	$(PETSC_CXX) $(CPPFLAGS) -DIGA_MULTIDOMAIN_NO_MAIN $(CXXFLAGS) $(OPENMP_FLAGS) $(PETSC_CFLAGS) tests/test_multidomain_subcommunicator.cpp src/iga_1d_3d_bifurcation.cpp -o $@ $(PETSC_LIBS)

multidomain_failure_test: tests/test_multidomain_failure.cpp tests/test_bifurcation_coupling_smoke.cpp tests/test_multidomain_flow_smoke.cpp src/iga_1d_3d_bifurcation.cpp $(wildcard ../../include/*.hpp ../cpu/include/*.hpp ../one_d/include/*.hpp)
	$(PETSC_CXX) $(CPPFLAGS) -DIGA_MULTIDOMAIN_NO_MAIN $(CXXFLAGS) $(OPENMP_FLAGS) $(PETSC_CFLAGS) tests/test_multidomain_failure.cpp src/iga_1d_3d_bifurcation.cpp -o $@ $(PETSC_LIBS)

# Phase 6 closure combines direct production immersed solves (physical and
# transactional gates) with real explicit and strong-fixed 1D--3D--1D runs.
phase6_aneurysm_closure_test: tests/test_phase6_aneurysm_closure.cpp iga_multidomain_flow
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

phase6-aneurysm-closure-test: phase6_aneurysm_closure_test
	$(MAKE) -C ../cpu phase6-aneurysm-runtime-closure-test PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=$(PETSC_ARCH)
	PATH="$(PETSC_PREFIX)/bin:$$PATH" ./phase6_aneurysm_closure_test ./iga_multidomain_flow

phase9_multiscale_closure_test: tests/test_phase9_multiscale_closure.cpp tests/CouplingFixture.hpp iga_multidomain_flow ../cpu/include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

phase9-multiscale-closure-test: phase9_multiscale_closure_test
	PATH="$(PETSC_PREFIX)/bin:$$PATH" ./phase9_multiscale_closure_test ./iga_multidomain_flow

immersed-case-factory-test: immersed_case_factory_test
	PATH="$(PETSC_PREFIX)/bin:$$PATH" ./immersed_case_factory_test

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

iga_1d_3d_explicit: src/iga_1d_3d_explicit.cpp ../../include/ExplicitOneDThreeDCoupling.hpp ../../include/StrongOneDThreeDCoupling.hpp ../../include/AitkenRelaxation.hpp ../../include/SequentialMultidomainCase.hpp ../../include/MultidomainConfig.hpp ../../include/PressureFlowComponentExecutor.hpp ../../include/DomainRuntimeRegistry.hpp ../../include/CoupledDomainRuntime.hpp ../../include/OneDFlowDomainAdapter.hpp ../../include/ThreeDBodyFittedFlowDomainAdapter.hpp ../../include/FlowDomainPortMetadata.hpp ../../include/SimulationGraph.hpp ../../include/CouplingEdge.hpp ../../include/ThreeDFlowCoupling.hpp ../../include/ThreeDVcaCoupling.hpp ../../include/SimulationConfig.hpp ../../include/CouplingConfig.hpp ../../include/CouplingPort.hpp $(SPECIES_COUPLING_HEADER) ../cpu/include/BoundarySupport.hpp ../cpu/include/IgaDatabase.hpp ../cpu/include/TransientFlowRuntime.hpp ../cpu/include/TransportElement.hpp ../one_d/include/OneDImplicit.hpp ../one_d/include/OneDRuntime.hpp
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)

iga_1d_3d_bifurcation: src/iga_1d_3d_bifurcation.cpp ../../include/BifurcationMultidomainCase.hpp ../../include/SequentialMultidomainCase.hpp ../../include/ImmersedFlowCase.hpp ../../include/ThreeDImmersedFlowDomain.hpp ../../include/MultidomainConfig.hpp ../../include/PressureFlowComponentExecutor.hpp ../../include/SpeciesPressureFlowComponentExecutor.hpp ../../include/DomainRuntimeRegistry.hpp ../../include/CoupledDomainRuntime.hpp ../../include/OneDFlowDomainAdapter.hpp ../../include/ThreeDBodyFittedFlowDomainAdapter.hpp ../../include/ThreeDBodyFittedFlowTransportDomainAdapter.hpp ../../include/FlowDomainPortMetadata.hpp ../../include/SimulationGraph.hpp ../../include/CouplingEdge.hpp ../../include/ThreeDFlowCoupling.hpp ../../include/ThreeDVcaCoupling.hpp ../../include/SimulationConfig.hpp ../../include/CouplingConfig.hpp ../../include/CouplingPort.hpp $(SPECIES_COUPLING_HEADER) ../cpu/include/BoundarySupport.hpp ../cpu/include/IgaDatabase.hpp ../cpu/include/TransientFlowRuntime.hpp ../cpu/include/TransientTransportRuntime.hpp ../cpu/include/TransportElement.hpp ../cpu/include/ImmersedStaticFlowRuntime.hpp ../one_d/include/OneDImplicit.hpp ../one_d/include/OneDRuntime.hpp ../one_d/include/OneDTransport.hpp
	$(PETSC_CXX) $(CPPFLAGS) -DIGA_BIFURCATION_ENTRY=1 $(CXXFLAGS) $(OPENMP_FLAGS) $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)

iga_multidomain_flow: src/iga_1d_3d_bifurcation.cpp ../../include/SequentialMultidomainCase.hpp ../../include/ImmersedFlowCase.hpp ../../include/ThreeDImmersedFlowDomain.hpp ../../include/MultidomainConfig.hpp ../../include/PressureFlowComponentExecutor.hpp ../../include/SpeciesPressureFlowComponentExecutor.hpp ../../include/DomainRuntimeRegistry.hpp ../../include/CoupledDomainRuntime.hpp ../../include/OneDFlowDomainAdapter.hpp ../../include/ThreeDBodyFittedFlowDomainAdapter.hpp ../../include/ThreeDBodyFittedFlowTransportDomainAdapter.hpp ../../include/FlowDomainPortMetadata.hpp ../../include/SimulationGraph.hpp ../../include/CouplingEdge.hpp ../../include/ThreeDFlowCoupling.hpp ../../include/ThreeDVcaCoupling.hpp ../../include/SimulationConfig.hpp ../../include/CouplingConfig.hpp ../../include/CouplingPort.hpp $(SPECIES_COUPLING_HEADER) ../cpu/include/BoundarySupport.hpp ../cpu/include/IgaDatabase.hpp ../cpu/include/TransientFlowRuntime.hpp ../cpu/include/TransientTransportRuntime.hpp ../cpu/include/TransportElement.hpp ../cpu/include/ImmersedStaticFlowRuntime.hpp ../one_d/include/OneDImplicit.hpp ../one_d/include/OneDRuntime.hpp ../one_d/include/OneDTransport.hpp
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)

immersed_case_factory_test: tests/test_immersed_case_factory.cpp ../../include/ImmersedFlowCase.hpp ../../include/ThreeDImmersedFlowDomain.hpp ../../include/DomainRuntimeRegistry.hpp ../../include/SimulationGraph.hpp ../../include/SimulationConfig.hpp ../../include/CaseConfig.hpp ../cpu/include/ImmersedStaticFlowRuntime.hpp ../cpu/include/ImmersedFlowPort.hpp ../cpu/include/CutCellGhostPenalty.hpp ../cpu/include/ImmersedSurfaceQuadrature.hpp ../cpu/include/CutCellVolumeQuadrature.hpp ../cpu/include/CartesianDomainClassification.hpp ../cpu/include/SurfaceReaders.hpp
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)

explicit_coupling_smoke_test: tests/test_explicit_coupling_smoke.cpp iga_1d_3d_explicit ../cpu/include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

bifurcation_coupling_smoke_test: tests/test_bifurcation_coupling_smoke.cpp iga_1d_3d_bifurcation ../cpu/include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

multidomain_flow_smoke_test: tests/test_multidomain_flow_smoke.cpp tests/test_bifurcation_coupling_smoke.cpp iga_multidomain_flow iga_1d_3d_bifurcation ../cpu/include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

straight_vessel_convergence_test: tests/test_straight_vessel_convergence.cpp tests/CouplingFixture.hpp ../cpu/include/BoundarySupport.hpp ../cpu/include/IgaDatabase.hpp ../one_d/include/OneDRuntime.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

pulsatile_straight_vessel_temporal_test: tests/test_pulsatile_straight_vessel_temporal.cpp tests/CouplingFixture.hpp ../cpu/include/BoundarySupport.hpp ../cpu/include/IgaDatabase.hpp ../one_d/include/OneDRuntime.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

clean:
	$(RM) sequential_strong_failure_test
	$(RM) sequential_text_failure_test
	$(RM) sequential_initialization_failure_test
	$(RM) text_close_preload.so
	$(RM) registry_failure_test graph_registry_failure_test
	$(RM) species_collective_failure_test
	$(RM) pressure_flow_collective_failure_test
	$(RM) collective_asset_input_test
	$(RM) runtime_construction_failure_test
	$(RM) three_d_staged_failure_test
	$(RM) three_d_trial_failure_test
	$(RM) three_d_port_failure_test
	$(RM) multidomain_failure_test
	$(RM) multidomain_subcommunicator_test
	$(RM) explicit_coupling_test simulation_graph_test pressure_flow_executor_test species_pressure_flow_executor_test multidomain_config_test zero_d_flow_test distributed_surface_contracts_test fsi_runtime_contracts_test explicit_coupling_smoke_test bifurcation_coupling_smoke_test multidomain_flow_smoke_test immersed_case_factory_test phase6_aneurysm_closure_test phase9_multiscale_closure_test straight_vessel_convergence_test pulsatile_straight_vessel_temporal_test iga_1d_3d_explicit iga_1d_3d_bifurcation iga_multidomain_flow

iga_1d_3d_explicit iga_1d_3d_bifurcation iga_multidomain_flow: ../one_d/include/OneDPetscSupport.hpp ../cpu/include/CollectiveFailure.hpp ../cpu/include/PetscReadArray.hpp ../cpu/include/RuntimeConstruction.hpp

runtime_construction_failure_test: tests/test_runtime_construction_failure.cpp tests/test_bifurcation_coupling_smoke.cpp tests/CouplingFixture.hpp $(wildcard ../../include/*.hpp ../cpu/include/*.hpp)
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) -UNDEBUG -DIGA_RUNTIME_CONSTRUCTION_TESTING $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)

iga_1d_3d_explicit iga_1d_3d_bifurcation iga_multidomain_flow: ../cpu/include/CollectivePressureFlowExecution.hpp

pressure_flow_collective_failure_test: tests/test_pressure_flow_collective_failure.cpp tests/test_pressure_flow_executor.cpp $(wildcard ../../include/*.hpp ../cpu/include/*.hpp)
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)

iga_1d_3d_bifurcation iga_multidomain_flow: ../cpu/include/CollectiveSpeciesPressureFlowExecution.hpp
iga_1d_3d_explicit: ../../include/SpeciesPressureFlowComponentExecutor.hpp

species_collective_failure_test: tests/test_species_collective_failure.cpp tests/test_species_pressure_flow_executor.cpp $(wildcard ../../include/*.hpp ../cpu/include/*.hpp)
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)

.PHONY: multidomain-species-test
multidomain-species-test: iga_multidomain_flow iga_1d_3d_bifurcation multidomain_flow_smoke_test
	./multidomain_flow_smoke_test --species-only

iga_multidomain_flow iga_1d_3d_bifurcation iga_1d_3d_explicit: ../cpu/include/CollectiveDomainRuntimeRegistry.hpp

text_close_preload.so: tests/text_close_preload.cpp
	$(CXX) $(CXXFLAGS) -fPIC -shared $< -o $@ -ldl

iga_1d_3d_explicit: ../../include/MultidomainRunner.hpp ../cpu/include/CollectiveAssetInput.hpp ../cpu/include/CollectivePetscOptions.hpp
iga_1d_3d_explicit: ../cpu/include/SequentialExecution.hpp

sequential_strong_failure_test: tests/test_sequential_strong_failure.cpp tests/test_explicit_coupling_smoke.cpp src/iga_1d_3d_explicit.cpp $(wildcard ../../include/*.hpp ../cpu/include/*.hpp ../one_d/include/*.hpp)
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) -UNDEBUG -DIGA_RUNTIME_CONSTRUCTION_TESTING -DIGA_SEQUENTIAL_EXECUTION_TESTING -DIGA_SEQUENTIAL_NO_MAIN $(PETSC_CFLAGS) tests/test_sequential_strong_failure.cpp src/iga_1d_3d_explicit.cpp -o $@ $(PETSC_LIBS)

sequential_text_failure_test: tests/test_sequential_text_failure.cpp src/iga_1d_3d_explicit.cpp $(wildcard ../../include/*.hpp ../cpu/include/*.hpp ../one_d/include/*.hpp)
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) -UNDEBUG $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)

sequential_initialization_failure_test: tests/test_sequential_initialization_failure.cpp tests/test_explicit_coupling_smoke.cpp src/iga_1d_3d_explicit.cpp $(wildcard ../../include/*.hpp ../cpu/include/*.hpp ../one_d/include/*.hpp)
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) -UNDEBUG -DIGA_RUNTIME_CONSTRUCTION_TESTING -DIGA_SEQUENTIAL_NO_MAIN $(PETSC_CFLAGS) tests/test_sequential_initialization_failure.cpp src/iga_1d_3d_explicit.cpp -o $@ $(PETSC_LIBS)

registry_failure_test: tests/test_registry_failure.cpp tests/test_pressure_flow_executor.cpp $(wildcard ../../include/*.hpp ../cpu/include/*.hpp)
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG -DIGA_RUNTIME_CONSTRUCTION_TESTING $(PETSC_CFLAGS) $< -o $@ $(PETSC_LIBS)

graph_registry_failure_test: tests/test_graph_registry_failure.cpp tests/test_multidomain_flow_smoke.cpp tests/test_bifurcation_coupling_smoke.cpp src/iga_1d_3d_bifurcation.cpp $(wildcard ../../include/*.hpp ../cpu/include/*.hpp ../one_d/include/*.hpp)
	$(PETSC_CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) -UNDEBUG -DIGA_RUNTIME_CONSTRUCTION_TESTING -DIGA_MULTIDOMAIN_NO_MAIN $(PETSC_CFLAGS) tests/test_graph_registry_failure.cpp src/iga_1d_3d_bifurcation.cpp -o $@ $(PETSC_LIBS)
three_d_staged_failure_test three_d_trial_failure_test three_d_port_failure_test: ../cpu/tests/StringStreamFailure.hpp
pressure_flow_collective_failure_test species_collective_failure_test: tests/AbortReportAllocation.hpp
