CXX ?= g++
.DEFAULT_GOAL := all
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS ?= -Iinclude -I../../include
SPECIES_COUPLING_HEADER := ../../include/SpeciesCoupling.hpp
LDLIBS ?= -lstdc++fs
HDF5_CFLAGS ?= $(firstword $(shell pkg-config --cflags-only-I hdf5 2>/dev/null))
HDF5_ALL_CFLAGS ?= $(shell pkg-config --cflags hdf5 2>/dev/null)
HDF5_SERIAL_LIBRARY ?= $(firstword $(wildcard \
	/usr/lib/*/libhdf5_serial.so /usr/lib/*/libhdf5_serial.so.[0-9]* \
	/lib/*/libhdf5_serial.so /lib/*/libhdf5_serial.so.[0-9]*))
HDF5_LIBS ?= $(if $(HDF5_SERIAL_LIBRARY),$(HDF5_SERIAL_LIBRARY),$(shell pkg-config --libs hdf5 2>/dev/null))

PETSC_DIR ?=
PETSC_ARCH ?=
PETSC_TARGETS := body_fitted_subcommunicator_test petsc petsc-test phase5-closure-test phase5-closure-manufactured phase5-closure-manufactured-triplet phase5-closure-sliver phase5-closure-all phase5_closure_test phase6-aneurysm-depth2-regression phase6-aneurysm-runtime-closure-test phase6_aneurysm_runtime_closure_test phase7-focused-test phase7-lv-closure-test phase7_lv_closure_test phase8-compliant-channel-fsi-paraview phase8_compliant_channel_fsi_paraview quadrature_test generic_transport_test navier_stokes_test nitsche_wall_test ghost_penalty_test immersed_transient_state_test immersed-transient-state-test immersed_velocity_extension_test immersed-velocity-extension-test immersed_static_flow_test immersed_transient_flow_test immersed-transient-flow-test moving_immersed_transient_flow_test moving-immersed-transient-flow-test moving_immersed_flow_snapshot_test moving-immersed-flow-snapshot-test moving_immersed_flow_snapshot_publisher_test moving-immersed-flow-snapshot-publisher-test immersed_flow_port_test immersed_flow_domain_test immersed_aneurysm_jacobian_test transport_boundary_preflight_test vca_3d_runtime_test vca_3d_smoke_test iga_solve iga_assembly_smoke iga_transport iga_navier_stokes iga_mesh_check moving_cut_geometry_test moving-cut-geometry-test fluid_surface_traction_test fluid-surface-traction-test moving_immersed_transient_flow_fsi_runtime_test moving-immersed-transient-flow-fsi-runtime-test compliant_channel_fsi_test compliant-channel-fsi-test

PETSC_TARGETS += parallel_ownership_test parallel-ownership-test
PETSC_TARGETS += collective_failure_test
PETSC_TARGETS += collective_petsc_options_test petsc_solver_options_test petsc_factor_matrix_test moving_immersed_solver_options_test
PETSC_TARGETS += petsc_gather_test
PETSC_TARGETS += body_fitted_checkpoint_bundle_test body_fitted_accepted_checkpoint_test
PETSC_TARGETS += petsc_checkpoint_read_test
PETSC_TARGETS += petsc_checkpoint_write_test
PETSC_TARGETS += execution_resources_test execution_resources_serial_test
PETSC_TARGETS += memory_report_test
PETSC_TARGETS += immersed_distributed_assembly_test
PETSC_TARGETS += immersed_distributed_history_test
PETSC_TARGETS += immersed_distributed_transient_volume_test
PETSC_TARGETS += immersed_transient_distributed_operator_test
PETSC_TARGETS += immersed_conservation_roundoff_test
PETSC_TARGETS += immersed_transient_distributed_runtime_test immersed_transient_distributed_domain_test
PETSC_TARGETS += immersed_distributed_physics_test immersed_distributed_static_flow_test immersed_distributed_domain_test

immersed_distributed_history_test: tests/test_immersed_distributed_history.cpp include/ImmersedDistributedVelocityHistory.hpp include/ImmersedTransientState.hpp include/NavierStokesElement.hpp include/CollectiveFailure.hpp include/PetscReadArray.hpp include/RuntimeCleanup.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

immersed_distributed_transient_volume_test: tests/test_immersed_distributed_transient_volume.cpp include/ImmersedDistributedTransientVolume.hpp include/ImmersedDistributedVelocityHistory.hpp include/ImmersedDistributedAssembly.hpp include/ImmersedTransientState.hpp include/NavierStokesElement.hpp include/RuntimeConstruction.hpp include/CollectiveFailure.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

immersed_transient_distributed_operator_test: tests/test_immersed_transient_distributed_operator.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

immersed_conservation_roundoff_test: CPPFLAGS += -DIGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING
immersed_conservation_roundoff_test: tests/test_immersed_conservation_roundoff.cpp include/ImmersedTransientFlowRuntime.hpp include/MovingImmersedTransientFlowRuntime.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

immersed_transient_distributed_domain_test: tests/test_immersed_transient_distributed_domain.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

immersed_transient_distributed_runtime_test: tests/test_immersed_transient_distributed_runtime.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

immersed_distributed_static_flow_test immersed_distributed_domain_test: include/ImmersedDistributedNewtonRuntime.hpp

immersed_distributed_static_flow_test: tests/test_immersed_distributed_static_flow.cpp include/ImmersedStaticDistributedRuntime.hpp include/ImmersedStaticDistributedOperator.hpp include/ImmersedDistributedAssembly.hpp include/ImmersedStaticFlowRuntime.hpp include/CollectivePetscOptions.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

immersed_distributed_physics_test: tests/test_immersed_distributed_physics.cpp include/ImmersedStaticDistributedOperator.hpp include/ImmersedDistributedAssembly.hpp include/ImmersedStaticFlowRuntime.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

immersed_distributed_assembly_test: tests/test_immersed_distributed_assembly.cpp include/ImmersedDistributedAssembly.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

memory_report_test: tests/test_memory_report.cpp include/MemoryReport.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)
PETSC_TARGETS += parallel_immersed_volume_test parallel_immersed_volume_serial_test element_assembly_execution_test element_assembly_execution_openmp_test
PETSC_TARGETS += geometry_inspection_failure_test compliant_channel_fsi_openmp_test immersed_transient_flow_openmp_test
PETSC_TARGETS += parallel_body_fitted_flow_test iga_navier_stokes_openmp
PETSC_TARGETS += immersed_small_cut_solver_test
PETSC_TARGETS += flow_output_failure_test
PETSC_TARGETS += flow_step_failure_test legacy_transport_failure_test
PETSC_TARGETS += text_helper_failure_test text-helper-failure-test
PETSC_TARGETS += snapshot_publication_failure_test
PETSC_TARGETS += mesh_stdout_failure_test assembly_stdout_failure_test
PETSC_TARGETS += flow_stdout_failure_test transport_stdout_failure_test legacy_stdout_failure_test runtime_cleanup_cli_test

flow_stdout_failure_test: CPPFLAGS += -DIGA_TEST_FLOW_STDOUT
transport_stdout_failure_test: CPPFLAGS += -DIGA_TEST_TRANSPORT_STDOUT
legacy_stdout_failure_test: CPPFLAGS += -DIGA_TEST_LEGACY_STDOUT
flow_stdout_failure_test transport_stdout_failure_test legacy_stdout_failure_test: tests/test_tool_stdout_failure.cpp src/iga_navier_stokes.cpp src/iga_solve.cpp src/iga_transport.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $(HDF5_ALL_CFLAGS) $< -o $@ $(PETSC_LIB) $(HDF5_LIBS)
PETSC_TARGETS += womersley_stream_failure_test transport_budget_stream_failure_test
PETSC_TARGETS += fsi_export_stream_failure_test

fsi_export_stream_failure_test: tests/test_fsi_export_stream_failure.cpp tests/StringStreamFailure.hpp src/phase8_compliant_channel_fsi_paraview.cpp $(wildcard include/*.hpp ../../include/*.hpp tests/*.hpp)
	$(CXX) $(CPPFLAGS) -Itests $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

text_read_preload.so: tests/text_read_preload.cpp
	$(CXX) -shared -fPIC -Wall -Wextra -Wpedantic $< -o $@ -ldl

womersley_stream_failure_test: CPPFLAGS += -DIGA_TEST_WOMERSLEY_STREAM
womersley_stream_failure_test transport_budget_stream_failure_test: tests/test_serial_tool_stream_failure.cpp tests/StringStreamFailure.hpp src/iga_womersley_reference.cpp src/iga_transport_validate.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

mesh_stdout_failure_test: tests/test_tool_stdout_failure.cpp src/iga_mesh_check.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) -DIGA_TEST_MESH_STDOUT $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

assembly_stdout_failure_test: tests/test_tool_stdout_failure.cpp src/iga_assembly_smoke.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

snapshot_publication_failure_test: tests/test_snapshot_publication_failure.cpp tests/test_moving_immersed_flow_snapshot_publisher.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) -DIGA_MOVING_IMMERSED_FLOW_SNAPSHOT_PUBLISHER_TESTING $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

text_helper_failure_test: tests/test_text_helper_failure.cpp $(wildcard include/*.hpp ../one_d/include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) -I../one_d/include $(CXXFLAGS) $(PETSC_CC_INCLUDES) $(HDF5_ALL_CFLAGS) $< -o $@ $(PETSC_LIB) $(HDF5_LIBS)

.PHONY: text-helper-failure-test
text-helper-failure-test: text_helper_failure_test
	OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout --kill-after=5s 90s mpiexec -np 3 ./text_helper_failure_test

flow_step_failure_test: tests/test_flow_step_failure.cpp src/iga_navier_stokes.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $(HDF5_ALL_CFLAGS) $< -o $@ $(PETSC_LIB) $(HDF5_LIBS)

flow_output_failure_test: tests/test_flow_output_failure.cpp src/iga_navier_stokes.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $(HDF5_ALL_CFLAGS) $< -o $@ $(PETSC_LIB) $(HDF5_LIBS)

ifneq ($(strip $(filter $(PETSC_TARGETS),$(MAKECMDGOALS))),)
ifeq ($(strip $(PETSC_DIR)),)
$(error PETSC_DIR must be set when building PETSc targets)
endif
PETSC_PREFIX := $(PETSC_DIR)$(if $(strip $(PETSC_ARCH)),/$(PETSC_ARCH))
PETSC_VARIABLES ?= $(PETSC_PREFIX)/lib/petsc/conf/petscvariables
include $(PETSC_VARIABLES)
endif

.PHONY: all petsc petsc-test test config-check-dispatch-test vtkhdf-paraview-test clean phase5-closure-test phase5-closure-manufactured phase5-closure-manufactured-triplet phase5-closure-sliver phase5-closure-all phase6-aneurysm-depth2-regression phase6-aneurysm-runtime-closure-test phase7-focused-test phase7-lv-closure-test phase8-compliant-channel-fsi-paraview prescribed-surface-motion-test material-surface-patch-kinematics-test moving-cut-geometry-test immersed-transient-state-test immersed-transient-flow-test moving-immersed-transient-flow-test moving-immersed-transient-flow-fsi-runtime-test moving-immersed-flow-snapshot-test moving-immersed-flow-snapshot-publisher-test immersed-velocity-extension-test pretensioned-membrane-test pretensioned-membrane-fsi-runtime-test

# PETSc-backed runtimes share opt-in phase accounting and explicit KSP setup.
$(PETSC_TARGETS): include/PetscPhaseProfile.hpp ../../include/PhaseProfile.hpp
$(PETSC_TARGETS): include/ParallelOwnershipValidation.hpp ../../include/ParallelOwnership.hpp
$(PETSC_TARGETS): include/CollectiveFailure.hpp include/PetscReadArray.hpp include/RuntimeConstruction.hpp include/RuntimeCleanup.hpp
$(PETSC_TARGETS): include/ImmersedStaticFlowSetup.hpp
immersed_distributed_physics_test immersed_distributed_static_flow_test immersed_distributed_domain_test: ../../include/WeightedWorkPartition.hpp
$(PETSC_TARGETS): include/PetscGather.hpp
$(PETSC_TARGETS): include/PetscCheckpointRead.hpp
$(PETSC_TARGETS): include/PetscCheckpointWrite.hpp
$(PETSC_TARGETS) iga_pack iga_transport_validate iga_config_check iga_case_check iga_womersley_reference case_config_test simulation_config_test flow_checkpoint_test transport_checkpoint_test vca_checkpoint_test womersley_reference_test boundary_support_test temporal_vtkhdf_test: ../../include/CheckedText.hpp

checked_text_test: tests/test_checked_text.cpp ../../include/CheckedText.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

.PHONY: checked-text-test
checked-text-test: checked_text_test
	./checked_text_test

test: checked-text-test
$(PETSC_TARGETS): ../../include/ExecutionEnvironment.hpp include/ExecutionResources.hpp include/ElementAssemblyExecution.hpp include/ParallelElementBatch.hpp

parallel_body_fitted_flow_test: tests/test_parallel_body_fitted_flow.cpp include/TransientFlowRuntime.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -fopenmp -UNDEBUG -DIGA_FLOW_RUNTIME_TESTING $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

geometry_inspection_failure_test: tests/test_geometry_inspection_failure.cpp include/TransportElement.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

execution_resources_test: tests/test_execution_resources.cpp include/ExecutionResources.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -fopenmp $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

execution_resources_serial_test: tests/test_execution_resources.cpp include/ExecutionResources.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

petsc_checkpoint_write_test: tests/test_petsc_checkpoint_write.cpp include/PetscCheckpointWrite.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

petsc_checkpoint_read_test: tests/test_petsc_checkpoint_read.cpp include/PetscCheckpointRead.hpp include/CollectiveAssetInput.hpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

petsc_gather_test: tests/test_petsc_gather.cpp include/PetscGather.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

collective_failure_test: tests/test_collective_failure.cpp tests/test_parallel_ownership.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

collective_petsc_options_test: tests/test_collective_petsc_options.cpp include/CollectivePetscOptions.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

parallel_ownership_test: tests/test_parallel_ownership.cpp include/OwnedRowAssembler.hpp include/IgaDatabase.hpp include/SurfaceOwnershipValidation.hpp ../../include/DistributedSurfaceInterface.hpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

.PHONY: parallel-ownership-test
parallel-ownership-test: parallel_ownership_test
	OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout 30 $(MPIEXEC) -np 3 ./parallel_ownership_test
strong_fluid_structure_coupling_test pretensioned_membrane_test pretensioned_membrane_fsi_runtime_test: ../../include/PhaseProfile.hpp

# Shared Cartesian spline-face ordering is part of all immersed cut/ghost and
# wall paths; make header changes rebuild every affected executable.
immersed_face_targets := navier_stokes_test nitsche_wall_test ghost_penalty_test immersed_transient_state_test immersed_static_flow_test immersed_transient_flow_test moving_immersed_flow_snapshot_test immersed_flow_port_test immersed_flow_domain_test immersed_aneurysm_jacobian_test phase6_aneurysm_runtime_closure_test phase7_lv_closure_test phase5_closure_test moving_cut_geometry_test
$(immersed_face_targets): include/CubicCartesianSplineFace.hpp

runtime_cleanup_cli_test: tests/test_runtime_cleanup_cli.cpp src/iga_navier_stokes.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $(HDF5_ALL_CFLAGS) $< -o $@ $(PETSC_LIB) $(HDF5_LIBS)

all: iga_pack iga_inspect iga_case_check iga_config_check iga_flow_validate iga_transport_validate iga_womersley_reference

iga_inspect iga_flow_validate: ../../include/CheckedText.hpp

compliant_channel_fsi_test compliant_channel_fsi_openmp_test immersed_aneurysm_jacobian_test: tests/ReferenceStateOutput.hpp include/Sha256.hpp

reference_state_output_test: tests/test_reference_state_output.cpp tests/ReferenceStateOutput.hpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

element_assembly_execution_openmp_test: CXXFLAGS += -fopenmp
element_assembly_execution_test element_assembly_execution_openmp_test: tests/test_element_assembly_execution.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

parallel_immersed_volume_test: CXXFLAGS += -fopenmp
parallel_immersed_volume_test parallel_immersed_volume_serial_test: CPPFLAGS += -DIGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING
parallel_immersed_volume_test parallel_immersed_volume_serial_test: tests/test_parallel_immersed_volume.cpp tests/CompliantChannelFsiFixture.hpp include/ImmersedTransientFlowRuntime.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

parallel_element_batch_test: tests/test_parallel_element_batch.cpp include/ParallelElementBatch.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

parallel_element_batch_openmp_test: tests/test_parallel_element_batch.cpp include/ParallelElementBatch.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -fopenmp $< -o $@ $(LDLIBS)

.PHONY: parallel-element-batch-test
parallel-element-batch-test: parallel_element_batch_test parallel_element_batch_openmp_test
	./parallel_element_batch_test
	OMP_THREAD_LIMIT=4 ./parallel_element_batch_openmp_test

petsc: iga_assembly_smoke iga_solve iga_transport iga_navier_stokes iga_mesh_check

iga_pack: src/iga_pack.cpp include/IgaDatabase.hpp ../../include/IgaPreprocessCache.hpp ../../include/MeshConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

iga_inspect: src/iga_inspect.cpp include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

iga_case_check: src/iga_case_check.cpp include/IgaDatabase.hpp include/CaseInput.hpp include/GenericCaseInput.hpp ../../include/CaseConfig.hpp ../../include/SimulationConfig.hpp ../../include/TemporalFunction.hpp ../../include/OutletModel.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
iga_config_check: src/iga_config_check.cpp ../../include/MultidomainConfig.hpp ../../include/FlowDomainPortMetadata.hpp ../../include/SimulationGraph.hpp ../../include/CouplingEdge.hpp ../../include/CouplingPort.hpp ../../include/SpeciesCoupling.hpp ../../include/SimulationConfig.hpp ../../include/MeshConfig.hpp ../../include/OneDConfig.hpp ../../include/RadiusAnnotatedObj.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

config-check-dispatch-test: iga_config_check tests/test_iga_config_check.sh
	sh tests/test_iga_config_check.sh ./iga_config_check

iga_flow_validate: src/iga_flow_validate.cpp include/BoundaryFlow.hpp include/CaseInput.hpp include/IgaDatabase.hpp include/WomersleyReference.hpp ../../include/VelocitySeries.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

iga_transport_validate: src/iga_transport_validate.cpp include/TransportBudget.hpp include/BoundaryFlow.hpp include/GenericCaseInput.hpp include/CaseInput.hpp include/IgaDatabase.hpp include/Quadrature.hpp include/ElementGeometry.hpp ../../include/SimulationConfig.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

transport_budget_test: tests/test_transport_budget.cpp include/TransportBudget.hpp include/BoundaryFlow.hpp include/GenericCaseInput.hpp include/CaseInput.hpp include/IgaDatabase.hpp include/Quadrature.hpp include/ElementGeometry.hpp ../../include/SimulationConfig.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

.PHONY: transport-budget-test
transport-budget-test: transport_budget_test iga_transport_validate
	./transport_budget_test

test: transport-budget-test

iga_womersley_reference: src/iga_womersley_reference.cpp include/WomersleyReference.hpp include/IgaDatabase.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)


case_config_test: tests/test_case_config.cpp include/CaseInput.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

simulation_config_test: tests/test_simulation_config.cpp include/GenericCaseInput.hpp ../../include/CouplingPort.hpp $(SPECIES_COUPLING_HEADER) ../../include/FlowDomainPortMetadata.hpp ../../include/SimulationConfig.hpp ../../include/MeshConfig.hpp ../../include/ThreeDFlowCoupling.hpp ../../include/ThreeDVcaCoupling.hpp ../../include/TemporalFunction.hpp ../../include/VelocitySeries.hpp ../../include/CaseConfig.hpp ../../include/PhysiologyOutput.hpp ../../include/VtkOutput.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

boundary_support_test: tests/test_boundary_support.cpp include/BoundarySupport.hpp include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

coupling_port_test: tests/test_coupling_port.cpp ../../include/CouplingPort.hpp ../../include/SpeciesCoupling.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

quadrature_test: tests/test_quadrature.cpp include/Quadrature.hpp include/ElementGeometry.hpp include/TransportElement.hpp include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

generic_transport_test: tests/test_generic_transport.cpp include/GenericTransportElement.hpp include/FieldCouplingPattern.hpp include/TransportElement.hpp include/CaseInput.hpp ../../include/SimulationConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

navier_stokes_test: tests/test_navier_stokes.cpp include/ImmersedFlowPort.hpp include/ImmersedSurfaceQuadrature.hpp include/CartesianDomainClassification.hpp include/CartesianBackground.hpp include/SurfaceSpatialIndex.hpp include/SurfaceGeometry.hpp include/NavierStokesElement.hpp include/BoundaryFlow.hpp include/TransportElement.hpp include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

nitsche_wall_test: tests/test_nitsche_wall.cpp include/ImmersedNitscheWall.hpp include/PrescribedSurfaceMotion.hpp include/CutCellGhostPenalty.hpp include/ImmersedSurfaceQuadrature.hpp include/CutCellVolumeQuadrature.hpp include/NavierStokesElement.hpp include/CartesianDomainClassification.hpp include/CartesianBackground.hpp include/SurfaceSpatialIndex.hpp include/SurfaceGeometry.hpp include/Quadrature.hpp include/ElementGeometry.hpp include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

ghost_penalty_test: tests/test_ghost_penalty.cpp include/CutCellGhostPenalty.hpp include/ImmersedNitscheWall.hpp include/ImmersedSurfaceQuadrature.hpp include/CutCellVolumeQuadrature.hpp include/NavierStokesElement.hpp include/CartesianDomainClassification.hpp include/CartesianBackground.hpp include/SurfaceSpatialIndex.hpp include/SurfaceGeometry.hpp include/Quadrature.hpp include/ElementGeometry.hpp include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

immersed_transient_state_test: tests/test_immersed_transient_state.cpp include/ImmersedTransientState.hpp include/CutCellVolumeQuadrature.hpp include/CartesianDomainClassification.hpp include/CartesianBackground.hpp include/SurfaceSpatialIndex.hpp include/SurfaceGeometry.hpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

immersed-transient-state-test: immersed_transient_state_test
	./immersed_transient_state_test

immersed_velocity_extension_test: tests/test_immersed_velocity_extension.cpp include/ImmersedVelocityExtension.hpp include/ImmersedTransientState.hpp include/MovingCutGeometry.hpp include/CubicCartesianSplineFace.hpp include/CutCellVolumeQuadrature.hpp include/CartesianBackground.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

immersed_velocity_extension_test: CPPFLAGS += -DIGA_IMMERSED_VELOCITY_EXTENSION_TESTING

immersed-velocity-extension-test: immersed_velocity_extension_test
	./immersed_velocity_extension_test

immersed_static_flow_test: tests/test_immersed_static_flow.cpp include/ImmersedStaticFlowRuntime.hpp include/ImmersedFlowPort.hpp include/CutCellGhostPenalty.hpp include/ImmersedNitscheWall.hpp include/ImmersedSurfaceQuadrature.hpp include/CutCellVolumeQuadrature.hpp include/NavierStokesElement.hpp include/CartesianDomainClassification.hpp include/CartesianBackground.hpp include/SurfaceSpatialIndex.hpp include/SurfaceGeometry.hpp include/Quadrature.hpp include/ElementGeometry.hpp include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

immersed_transient_flow_openmp_test: CXXFLAGS += -fopenmp -UNDEBUG
immersed_transient_flow_test immersed_transient_flow_openmp_test: tests/test_immersed_transient_flow.cpp include/ImmersedTransientFlowRuntime.hpp include/ImmersedVelocityExtension.hpp include/MovingCutGeometry.hpp include/ImmersedTransientState.hpp include/ImmersedFlowPort.hpp include/CutCellGhostPenalty.hpp include/ImmersedNitscheWall.hpp include/NavierStokesElement.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

immersed-transient-flow-test: immersed_transient_flow_test
	./immersed_transient_flow_test

moving_immersed_transient_flow_test: tests/test_moving_immersed_transient_flow.cpp include/MaterialSurfaceKinematics.hpp include/MovingImmersedFlowSnapshotCapture.hpp include/MovingImmersedFlowSnapshotPublisher.hpp include/MovingImmersedFlowSnapshot.hpp include/MovingImmersedTransientFlowRuntime.hpp include/ImmersedTransientFlowRuntime.hpp include/ImmersedVelocityExtension.hpp include/MovingCutGeometry.hpp include/ImmersedTransientState.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

moving_immersed_transient_flow_test: CPPFLAGS += -DIGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING

moving-immersed-transient-flow-test: moving_immersed_transient_flow_test
	./moving_immersed_transient_flow_test

moving_immersed_flow_snapshot_test: tests/test_moving_immersed_flow_snapshot.cpp include/MovingImmersedFlowSnapshot.hpp include/MovingCutGeometry.hpp include/ImmersedTransientState.hpp include/TransportElement.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

moving-immersed-flow-snapshot-test: moving_immersed_flow_snapshot_test
	./moving_immersed_flow_snapshot_test

moving_immersed_flow_snapshot_publisher_test: tests/test_moving_immersed_flow_snapshot_publisher.cpp include/MovingImmersedFlowSnapshotPublisher.hpp include/MovingImmersedFlowSnapshot.hpp include/MovingImmersedTransientFlowRuntime.hpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

moving_immersed_flow_snapshot_publisher_test: CPPFLAGS += -DIGA_MOVING_IMMERSED_FLOW_SNAPSHOT_PUBLISHER_TESTING
moving_immersed_flow_snapshot_publisher_test: include/PrescribedSurfaceMotion.hpp ../../include/CheckedText.hpp

moving-immersed-flow-snapshot-publisher-test: moving_immersed_flow_snapshot_publisher_test
	./moving_immersed_flow_snapshot_publisher_test

immersed_flow_port_test: tests/test_immersed_flow_port.cpp include/ImmersedStaticFlowRuntime.hpp include/ImmersedFlowPort.hpp include/CutCellGhostPenalty.hpp include/ImmersedNitscheWall.hpp include/ImmersedSurfaceQuadrature.hpp include/CutCellVolumeQuadrature.hpp include/NavierStokesElement.hpp include/CartesianDomainClassification.hpp include/CartesianBackground.hpp include/SurfaceSpatialIndex.hpp include/SurfaceGeometry.hpp include/Quadrature.hpp include/ElementGeometry.hpp include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

immersed_flow_domain_test: tests/test_immersed_flow_domain.cpp ../../include/ThreeDImmersedFlowDomain.hpp ../../include/CoupledDomainRuntime.hpp ../../include/DomainRuntimeRegistry.hpp ../../include/SimulationGraph.hpp ../../include/FlowDomainPortMetadata.hpp ../../include/CouplingEdge.hpp ../../include/CouplingPort.hpp include/ImmersedStaticFlowRuntime.hpp include/ImmersedFlowPort.hpp include/CutCellGhostPenalty.hpp include/ImmersedNitscheWall.hpp include/ImmersedSurfaceQuadrature.hpp include/CutCellVolumeQuadrature.hpp include/NavierStokesElement.hpp include/CartesianDomainClassification.hpp include/CartesianBackground.hpp include/SurfaceSpatialIndex.hpp include/SurfaceGeometry.hpp include/Quadrature.hpp include/ElementGeometry.hpp include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

immersed_aneurysm_jacobian_test: tests/test_immersed_aneurysm_jacobian.cpp ../../include/ImmersedFlowCase.hpp ../../include/ThreeDImmersedFlowDomain.hpp ../../include/CoupledDomainRuntime.hpp ../../include/FlowDomainPortMetadata.hpp ../../include/CouplingPort.hpp include/ImmersedStaticFlowRuntime.hpp include/ImmersedFlowPort.hpp include/CutCellGhostPenalty.hpp include/ImmersedNitscheWall.hpp include/ImmersedSurfaceQuadrature.hpp include/CutCellVolumeQuadrature.hpp include/NavierStokesElement.hpp include/CartesianDomainClassification.hpp include/CartesianBackground.hpp include/SurfaceSpatialIndex.hpp include/SurfaceGeometry.hpp include/SurfaceReaders.hpp include/Quadrature.hpp include/ElementGeometry.hpp include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

# This is the focused depth-2 Jacobian/conservation regression, not a full
# coupled-chain closure target.  The checked-in production fixture stays at
# depth 3; the target makes an isolated depth-2 copy for this regression.
phase6-aneurysm-depth2-regression: immersed_aneurysm_jacobian_test
	@case_dir=$$(mktemp -d /tmp/immersed-aneurysm-depth2.XXXXXX); trap 'rm -rf "$$case_dir"' EXIT; \
	cp ../../examples/vascular_flow/immersed_aneurysm_chain/immersed/simulation_config.json "$$case_dir"/; \
	cp ../../examples/vascular_flow/immersed_aneurysm_chain/immersed/surface.vtp "$$case_dir"/; \
	sed 's/"max_depth":3/"max_depth":2/' ../../examples/vascular_flow/immersed_aneurysm_chain/immersed/immersed_geometry.json > "$$case_dir"/immersed_geometry.json; \
	./immersed_aneurysm_jacobian_test "$$case_dir"

# Real quasi-static immersed closure measurements and transactional replay.
# The coupled explicit/fixed runner evidence is composed by the coupling target.
phase6_aneurysm_runtime_closure_test: tests/test_phase6_aneurysm_runtime_closure.cpp ../../include/ImmersedFlowCase.hpp ../../include/ThreeDImmersedFlowDomain.hpp ../../include/CoupledDomainRuntime.hpp ../../include/FlowDomainPortMetadata.hpp ../../include/CouplingPort.hpp include/ImmersedStaticFlowRuntime.hpp include/ImmersedFlowPort.hpp include/CutCellGhostPenalty.hpp include/ImmersedNitscheWall.hpp include/ImmersedSurfaceQuadrature.hpp include/CutCellVolumeQuadrature.hpp include/NavierStokesElement.hpp include/CartesianDomainClassification.hpp include/CartesianBackground.hpp include/SurfaceSpatialIndex.hpp include/SurfaceGeometry.hpp include/SurfaceReaders.hpp include/Quadrature.hpp include/ElementGeometry.hpp include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

phase6-aneurysm-runtime-closure-test: phase6_aneurysm_runtime_closure_test
	./phase6_aneurysm_runtime_closure_test

# Focused real 16-step prescribed-motion closure.  It publishes only committed
# epochs into a temporary directory by default; pass an output stem directly to
# the executable when retaining VTU/JSON/PVD evidence is desired.
phase7_lv_closure_test: tests/test_phase7_lv_closure.cpp include/IdealizedLeftVentricleFixture.hpp include/MovingImmersedFlowSnapshotCapture.hpp include/MovingImmersedFlowSnapshotPublisher.hpp include/MovingImmersedTransientFlowRuntime.hpp include/ImmersedTransientFlowRuntime.hpp include/ImmersedVelocityExtension.hpp include/MovingCutGeometry.hpp include/PrescribedSurfaceMotion.hpp include/SurfaceGeometry.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

phase7_lv_closure_test: CPPFLAGS += -DIGA_MOVING_IMMERSED_FLOW_SNAPSHOT_PUBLISHER_TESTING

phase7-lv-closure-test: phase7_lv_closure_test
	./phase7_lv_closure_test

# Standard Phase 7 gate: focused moving-cut/runtime/snapshot coverage and its
# cut, ghost, Nitsche, and port prerequisites.  The expensive 16-step LV
# closure remains the separate phase7-lv-closure-test target above.
phase7-focused-test: cut_cell_volume_quadrature_test immersed_surface_quadrature_test nitsche_wall_test ghost_penalty_test immersed_flow_port_test moving_cut_geometry_test immersed_transient_state_test immersed_velocity_extension_test immersed_transient_flow_test moving_immersed_transient_flow_test moving_immersed_flow_snapshot_test moving_immersed_flow_snapshot_publisher_test
	./cut_cell_volume_quadrature_test
	./immersed_surface_quadrature_test
	./nitsche_wall_test
	./ghost_penalty_test
	./immersed_flow_port_test
	./moving_cut_geometry_test
	./immersed_transient_state_test
	./immersed_velocity_extension_test
	./immersed_transient_flow_test
	./moving_immersed_transient_flow_test
	./moving_immersed_flow_snapshot_test
	./moving_immersed_flow_snapshot_publisher_test

# Deliberately outside petsc-test: this is the bounded but substantive Phase 5
# real-PETSc closure evidence (three solves plus eight global sliver spectra).
phase5_closure_test: tests/test_phase5_closure.cpp include/ImmersedStaticFlowRuntime.hpp include/ImmersedFlowPort.hpp include/CutCellGhostPenalty.hpp include/ImmersedNitscheWall.hpp include/ImmersedSurfaceQuadrature.hpp include/CutCellVolumeQuadrature.hpp include/NavierStokesElement.hpp include/CartesianDomainClassification.hpp include/CartesianBackground.hpp include/SurfaceSpatialIndex.hpp include/SurfaceGeometry.hpp include/Quadrature.hpp include/ElementGeometry.hpp include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

phase5-closure-test: phase5_closure_test
	./phase5_closure_test --all

# Independent selector targets are suitable for scheduler-array tasks, e.g.
# make phase5-closure-sliver PHASE5_M=4 PETSC_DIR=... .
phase5-closure-manufactured: phase5_closure_test
	@test -n "$(PHASE5_LEVEL)"
	./phase5_closure_test --manufactured-level $(PHASE5_LEVEL)

# Complete manufactured-convergence evidence: runs 3, 6, and 12, and applies
# all pairwise/order/finest-cap gates.  Use scheduler resources for this target.
phase5-closure-manufactured-triplet: phase5_closure_test
	./phase5_closure_test --manufactured-triplet

phase5-closure-sliver: phase5_closure_test
	@test -n "$(PHASE5_M)"
	./phase5_closure_test --sliver-m $(PHASE5_M)

phase5-closure-all: phase5_closure_test
	./phase5_closure_test --all

test: case_config_test simulation_config_test boundary_support_test coupling_port_test sha256_test surface_geometry_test surface_spatial_index_test cartesian_domain_classification_test cut_cell_volume_quadrature_test immersed_surface_quadrature_test surface_stl_test surface_vtp_test cartesian_background_test vtk_output_test bezier_visualization_test temporal_vtkhdf_test flow_checkpoint_test transport_checkpoint_test vca_checkpoint_test outlet_model_test boundary_flow_test womersley_reference_test config-check-dispatch-test
	./case_config_test
	./simulation_config_test
	./boundary_support_test
	./coupling_port_test
	./sha256_test
	./surface_geometry_test
	./surface_spatial_index_test
	./cartesian_domain_classification_test
	./cut_cell_volume_quadrature_test
	./immersed_surface_quadrature_test
	./surface_stl_test
	./surface_vtp_test
	./cartesian_background_test
	./vtk_output_test
	./bezier_visualization_test
	./temporal_vtkhdf_test
	./flow_checkpoint_test
	./transport_checkpoint_test
	./vca_checkpoint_test
	./outlet_model_test
	./boundary_flow_test
	./womersley_reference_test

cartesian_background_test: tests/test_cartesian_background.cpp include/CartesianBackground.hpp include/ElementGeometry.hpp include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

sha256_test: tests/test_sha256.cpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

surface_geometry_test: tests/test_surface_geometry.cpp include/SurfaceGeometry.hpp include/ExactDyadicPredicates.hpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

prescribed_surface_motion_test: tests/test_prescribed_surface_motion.cpp include/MaterialSurfaceKinematics.hpp include/PrescribedSurfaceMotion.hpp include/SurfaceGeometry.hpp include/ExactDyadicPredicates.hpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

prescribed-surface-motion-test: prescribed_surface_motion_test
	./prescribed_surface_motion_test

material_surface_patch_kinematics_test: tests/test_material_surface_patch_kinematics.cpp include/MaterialSurfacePatchKinematics.hpp include/MaterialSurfacePatchMap.hpp include/MaterialSurfaceKinematics.hpp include/SurfaceGeometry.hpp ../../include/DistributedSurfaceInterface.hpp ../../include/FsiDomainRuntime.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

material-surface-patch-kinematics-test: material_surface_patch_kinematics_test
	./material_surface_patch_kinematics_test

moving_cut_geometry_test: tests/test_moving_cut_geometry.cpp include/MaterialSurfaceKinematics.hpp include/MovingCutGeometry.hpp include/PrescribedSurfaceMotion.hpp include/CutCellGhostPenalty.hpp include/ImmersedSurfaceQuadrature.hpp include/CutCellVolumeQuadrature.hpp include/CartesianDomainClassification.hpp include/CartesianBackground.hpp include/SurfaceSpatialIndex.hpp include/SurfaceGeometry.hpp include/Quadrature.hpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

moving-cut-geometry-test: moving_cut_geometry_test
	./moving_cut_geometry_test

fluid_surface_traction_test: tests/test_fluid_surface_traction.cpp include/FluidSurfaceTraction.hpp include/FluidCauchyStress.hpp include/ImmersedFlowPort.hpp include/ImmersedSurfaceQuadrature.hpp include/MaterialSurfaceKinematics.hpp include/MaterialSurfacePatchMap.hpp include/NavierStokesElement.hpp ../../include/DistributedSurfaceInterface.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

fluid-surface-traction-test: fluid_surface_traction_test
	./fluid_surface_traction_test

moving_immersed_transient_flow_fsi_runtime_test: CPPFLAGS += -DIGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING
moving_immersed_transient_flow_fsi_runtime_test: tests/test_moving_immersed_transient_flow_fsi_runtime.cpp include/MovingImmersedTransientFlowFsiRuntime.hpp include/MovingImmersedTransientFlowRuntime.hpp include/FluidSurfaceTraction.hpp include/MaterialSurfacePatchKinematics.hpp ../../include/FsiDomainRuntime.hpp ../../include/FsiCouplingEdge.hpp ../../include/DistributedSurfaceInterface.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

moving-immersed-transient-flow-fsi-runtime-test: moving_immersed_transient_flow_fsi_runtime_test
	./moving_immersed_transient_flow_fsi_runtime_test

compliant_channel_fsi_test compliant_channel_fsi_openmp_test: CPPFLAGS += -DIGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING
compliant_channel_fsi_openmp_test: CXXFLAGS += -fopenmp
compliant_channel_fsi_test compliant_channel_fsi_openmp_test: include/ImmersedTransientFlowRuntime.hpp
compliant_channel_fsi_test compliant_channel_fsi_openmp_test: tests/test_compliant_channel_fsi.cpp tests/CompliantChannelFsiFixture.hpp include/StrongFluidStructureCoupling.hpp include/MovingImmersedTransientFlowFsiRuntime.hpp include/PretensionedMembraneFsiRuntime.hpp include/PretensionedMembrane.hpp include/MovingImmersedTransientFlowRuntime.hpp include/FluidSurfaceTraction.hpp include/MaterialSurfacePatchKinematics.hpp ../../include/FsiDomainRuntime.hpp ../../include/FsiCouplingEdge.hpp ../../include/DistributedSurfaceInterface.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

compliant-channel-fsi-test: compliant_channel_fsi_test
	./compliant_channel_fsi_test

phase8_compliant_channel_fsi_paraview: CPPFLAGS += -Itests
phase8_compliant_channel_fsi_paraview: src/phase8_compliant_channel_fsi_paraview.cpp tests/CompliantChannelFsiFixture.hpp include/StrongFluidStructureCoupling.hpp include/MovingImmersedTransientFlowFsiRuntime.hpp include/MovingImmersedFlowSnapshotCapture.hpp include/MovingImmersedFlowSnapshot.hpp include/PretensionedMembraneFsiRuntime.hpp include/PretensionedMembrane.hpp include/MovingImmersedTransientFlowRuntime.hpp include/FluidSurfaceTraction.hpp include/MaterialSurfacePatchKinematics.hpp ../../include/FsiDomainRuntime.hpp ../../include/FsiCouplingEdge.hpp ../../include/DistributedSurfaceInterface.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

PHASE8_PARAVIEW_DIR ?= ../../results/phase8_compliant_channel_fsi_paraview
phase8-compliant-channel-fsi-paraview: phase8_compliant_channel_fsi_paraview
	./phase8_compliant_channel_fsi_paraview "$(PHASE8_PARAVIEW_DIR)"

pretensioned_membrane_test: tests/test_pretensioned_membrane.cpp include/PretensionedMembrane.hpp ../../include/DistributedSurfaceInterface.hpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

pretensioned-membrane-test: pretensioned_membrane_test
	./pretensioned_membrane_test

pretensioned_membrane_fsi_runtime_test: tests/test_pretensioned_membrane_fsi_runtime.cpp include/PretensionedMembraneFsiRuntime.hpp include/PretensionedMembrane.hpp ../../include/FsiDomainRuntime.hpp ../../include/FsiCouplingEdge.hpp ../../include/DistributedSurfaceInterface.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@ $(LDLIBS)

pretensioned-membrane-fsi-runtime-test: pretensioned_membrane_fsi_runtime_test
	./pretensioned_membrane_fsi_runtime_test

strong_fluid_structure_coupling_test: tests/test_strong_fluid_structure_coupling.cpp include/StrongFluidStructureCoupling.hpp ../../include/DynamicWeightedAitkenRelaxation.hpp ../../include/FsiCouplingEdge.hpp ../../include/DistributedSurfaceInterface.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@ $(LDLIBS)

strong-fluid-structure-coupling-test: strong_fluid_structure_coupling_test
	./strong_fluid_structure_coupling_test

surface_spatial_index_test: tests/test_surface_spatial_index.cpp include/SurfaceSpatialIndex.hpp include/SurfaceGeometry.hpp include/ExactDyadicPredicates.hpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

cartesian_domain_classification_test: tests/test_cartesian_domain_classification.cpp include/CartesianDomainClassification.hpp include/CartesianBackground.hpp include/SurfaceSpatialIndex.hpp include/SurfaceGeometry.hpp include/ExactDyadicPredicates.hpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

cut_cell_volume_quadrature_test: tests/test_cut_cell_volume_quadrature.cpp include/CutCellVolumeQuadrature.hpp include/CartesianDomainClassification.hpp include/CartesianBackground.hpp include/SurfaceSpatialIndex.hpp include/SurfaceGeometry.hpp include/ExactDyadicPredicates.hpp include/Quadrature.hpp include/ElementGeometry.hpp include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

immersed_surface_quadrature_test: tests/test_immersed_surface_quadrature.cpp include/ImmersedSurfaceQuadrature.hpp include/CartesianDomainClassification.hpp include/CartesianBackground.hpp include/SurfaceSpatialIndex.hpp include/SurfaceGeometry.hpp include/ExactDyadicPredicates.hpp include/Quadrature.hpp include/ElementGeometry.hpp include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

immersed_surface_quadrature_test: CPPFLAGS += -DIGA_EXACT_DYADIC_TESTING

cartesian_domain_classification_test: CPPFLAGS += -DIGA_EXACT_DYADIC_TESTING

cut_cell_volume_quadrature_test: CPPFLAGS += -DIGA_EXACT_DYADIC_TESTING

surface_spatial_index_test: CPPFLAGS += -DIGA_EXACT_DYADIC_TESTING

surface_stl_test: tests/test_surface_stl.cpp include/SurfaceReaders.hpp include/SurfaceGeometry.hpp include/ExactDyadicPredicates.hpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

surface_vtp_test: tests/test_surface_vtp.cpp include/SurfaceReaders.hpp include/SurfaceGeometry.hpp include/ExactDyadicPredicates.hpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

vtk_output_test: tests/test_vtk_output.cpp ../../include/VtkOutput.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

text_output_failure_test: tests/test_text_output_failure.cpp ../../include/VtkOutput.hpp ../../include/PhysiologyOutput.hpp ../../include/VelocitySeries.hpp ../../include/SimulationConfig.hpp ../../include/OxygenCapacity.hpp ../../include/BezierVisualization.hpp ../../include/CouplingHistory.hpp ../../include/VascularCoupling.hpp include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

.PHONY: text-output-failure-test
text-output-failure-test: text_output_failure_test
	@directory=$$(mktemp -d) || exit 1; trap 'rm -rf "$$directory"' EXIT; timeout 30 ./text_output_failure_test "$$directory/results"

test: text-output-failure-test

bezier_visualization_test: tests/test_bezier_visualization.cpp include/IgaDatabase.hpp ../../include/BezierVisualization.hpp ../../include/VtkOutput.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

temporal_vtkhdf_test: tests/test_temporal_vtkhdf.cpp include/IgaDatabase.hpp ../../include/BezierVisualization.hpp ../../include/TemporalVtkHdf.hpp ../../include/VtkOutput.hpp
	$(CXX) $(CPPFLAGS) $(HDF5_ALL_CFLAGS) $(CXXFLAGS) $< -o $@ $(HDF5_LIBS) $(LDLIBS)

# Error-return injection requires a linker supporting GNU --wrap.
vtkhdf_close_failure_test: tests/test_vtkhdf_close_failure.cpp include/IgaDatabase.hpp ../../include/BezierVisualization.hpp ../../include/TemporalVtkHdf.hpp ../../include/VtkOutput.hpp
	$(CXX) $(CPPFLAGS) $(HDF5_ALL_CFLAGS) $(CXXFLAGS) $< -o $@ -Wl,--wrap=H5Fflush,--wrap=H5Gclose,--wrap=H5Fclose,--wrap=H5Dclose,--wrap=H5Fget_obj_count $(HDF5_LIBS) $(LDLIBS)

.PHONY: vtkhdf-close-failure-test
vtkhdf-close-failure-test: vtkhdf_close_failure_test
	@directory=$$(mktemp -d) || exit 1; trap 'rm -rf "$$directory"' EXIT; timeout 30 ./vtkhdf_close_failure_test "$$directory/results"


PVPYTHON ?= pvpython
vtkhdf-paraview-test: temporal_vtkhdf_test
	TUBULARFLOWIGA_KEEP_TEST_OUTPUT=1 ./temporal_vtkhdf_test
	$(PVPYTHON) tests/validate_temporal_vtkhdf.py

flow_checkpoint_test: tests/test_flow_checkpoint.cpp ../../include/FlowCheckpoint.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

transport_checkpoint_test: tests/test_transport_checkpoint.cpp ../../include/TransportCheckpoint.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

vca_checkpoint_test: tests/test_vca_checkpoint.cpp ../../include/VcaCheckpoint.hpp ../../include/VascularCoupling.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

outlet_model_test: tests/test_outlet_model.cpp ../../include/OutletModel.hpp ../../include/OutletCheckpoint.hpp ../../include/FlowCheckpoint.hpp ../../include/SimulationConfig.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

boundary_flow_test: tests/test_boundary_flow.cpp include/BoundaryFlow.hpp include/OutletFlow.hpp include/PressureTractionFlow.hpp include/IgaDatabase.hpp ../../include/OutletModel.hpp ../../include/PressureTraction.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

womersley_reference_test: tests/test_womersley_reference.cpp include/WomersleyReference.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

petsc-test: quadrature_test generic_transport_test navier_stokes_test nitsche_wall_test ghost_penalty_test immersed_static_flow_test immersed_flow_port_test immersed_flow_domain_test transport_boundary_preflight_test vca_3d_runtime_test vca_3d_smoke_test
	./quadrature_test
	./generic_transport_test
	./navier_stokes_test
	./nitsche_wall_test
	./ghost_penalty_test
	./immersed_static_flow_test
	./immersed_flow_port_test
	./immersed_flow_domain_test
	PATH="$(PETSC_PREFIX)/bin:$$PATH" mpiexec -np 2 ./transport_boundary_preflight_test
	./vca_3d_runtime_test
	PATH="$(PETSC_PREFIX)/bin:$$PATH" ./vca_3d_smoke_test

vca_3d_runtime_test: tests/test_vca_3d_runtime.cpp include/IgaDatabase.hpp include/OwnedRowAssembler.hpp include/FieldCouplingPattern.hpp include/GenericCaseInput.hpp include/GenericTransportElement.hpp include/TransientFlowRuntime.hpp include/TransientTransportRuntime.hpp ../../include/ThreeDBodyFittedFlowDomainAdapter.hpp ../../include/ThreeDBodyFittedFlowTransportDomainAdapter.hpp ../../include/FlowDomainPortMetadata.hpp ../../include/ThreeDFlowCoupling.hpp ../../include/CoupledDomainRuntime.hpp ../../include/CouplingPort.hpp $(SPECIES_COUPLING_HEADER) ../../include/SimulationConfig.hpp ../../include/TemporalFunction.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

vca_3d_smoke_test: tests/test_vca_3d_smoke.cpp iga_navier_stokes include/IgaDatabase.hpp ../../include/SimulationConfig.hpp ../../include/VcaCheckpoint.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

transport_boundary_preflight_test: tests/test_transport_boundary_preflight.cpp include/IgaDatabase.hpp include/OwnedRowAssembler.hpp include/TransportBoundaryPreflight.hpp include/TransportElement.hpp ../../include/SimulationConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

iga_assembly_smoke: src/iga_assembly_smoke.cpp include/IgaDatabase.hpp include/OwnedRowAssembler.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

body_fitted_subcommunicator_test: tests/test_body_fitted_subcommunicator.cpp include/TransientFlowRuntime.hpp include/TransientTransportRuntime.hpp include/OwnedRowAssembler.hpp include/BoundarySupport.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

phase_profile_test: tests/test_phase_profile.cpp ../../include/PhaseProfile.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

.PHONY: phase-profile-test
phase-profile-test: phase_profile_test
	./phase_profile_test

iga_solve: src/iga_solve.cpp ../../include/PhaseProfile.hpp include/IgaDatabase.hpp include/MemoryReport.hpp include/OwnedRowAssembler.hpp include/TransportBoundaryPreflight.hpp include/FieldCouplingPattern.hpp include/CaseInput.hpp include/GenericCaseInput.hpp include/GenericTransportElement.hpp ../../include/SimulationConfig.hpp ../../include/TemporalFunction.hpp ../../include/VelocitySeries.hpp ../../include/FlowCheckpoint.hpp ../../include/TransportCheckpoint.hpp ../../include/VtkOutput.hpp ../../include/PhysiologyOutput.hpp ../../include/BezierVisualization.hpp ../../include/TemporalVtkHdf.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $(HDF5_CFLAGS) $< -o $@ $(PETSC_LIB) $(HDF5_LIBS)

iga_transport: src/iga_transport.cpp include/IgaDatabase.hpp include/OwnedRowAssembler.hpp include/FieldCouplingPattern.hpp include/CaseInput.hpp include/TransportElement.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

iga_navier_stokes_openmp: CXXFLAGS += -fopenmp
iga_navier_stokes iga_navier_stokes_openmp: include/CollectiveAssetInput.hpp include/Sha256.hpp include/PetscCheckpointRead.hpp include/PetscCheckpointWrite.hpp
iga_navier_stokes iga_navier_stokes_openmp: src/iga_navier_stokes.cpp include/BoundarySupport.hpp include/IgaDatabase.hpp include/OwnedRowAssembler.hpp include/CaseInput.hpp include/GenericCaseInput.hpp include/GenericTransportElement.hpp include/TransientFlowRuntime.hpp include/TransientTransportRuntime.hpp include/TransportElement.hpp include/NavierStokesElement.hpp include/BoundaryFlow.hpp include/PressureTractionFlow.hpp ../../include/CaseConfig.hpp ../../include/CouplingPort.hpp $(SPECIES_COUPLING_HEADER) ../../include/SimulationConfig.hpp ../../include/ThreeDFlowCoupling.hpp ../../include/ThreeDVcaCoupling.hpp ../../include/VascularCoupling.hpp ../../include/CouplingConfig.hpp ../../include/TemporalFunction.hpp ../../include/FlowCheckpoint.hpp ../../include/VcaCheckpoint.hpp ../../include/CouplingHistory.hpp ../../include/OutletModel.hpp ../../include/OutletCheckpoint.hpp ../../include/PressureTraction.hpp ../../include/VelocitySeries.hpp ../../include/VtkOutput.hpp ../../include/BezierVisualization.hpp ../../include/TemporalVtkHdf.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $(HDF5_CFLAGS) $< -o $@ $(PETSC_LIB) $(HDF5_LIBS)

iga_mesh_check: src/iga_mesh_check.cpp include/IgaDatabase.hpp include/TransportElement.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

iga_solve: include/CollectiveAssetInput.hpp include/CollectivePetscOptions.hpp include/MemoryReport.hpp include/PetscCheckpointRead.hpp include/PetscCheckpointWrite.hpp include/PetscGather.hpp include/PetscReadArray.hpp include/Sha256.hpp

clean:
	$(RM) petsc_solver_options_test
	$(RM) body_fitted_checkpoint_bundle_test
	$(RM) body_fitted_accepted_checkpoint_test
	$(RM) runtime_cleanup_cli_test
	$(RM) snapshot_publication_failure_test mesh_stdout_failure_test assembly_stdout_failure_test
	$(RM) flow_stdout_failure_test transport_stdout_failure_test legacy_stdout_failure_test
	$(RM) womersley_stream_failure_test transport_budget_stream_failure_test fsi_export_stream_failure_test text_read_preload.so
	$(RM) checked_text_test
	$(RM) parallel_body_fitted_flow_test iga_navier_stokes_openmp parallel_immersed_volume_serial_test element_assembly_execution_test element_assembly_execution_openmp_test parallel_immersed_volume_test compliant_channel_fsi_openmp_test immersed_transient_flow_openmp_test
	$(RM) iga_transport_validate transport_budget_test reference_state_output_test parallel_element_batch_test parallel_element_batch_openmp_test
	$(RM) petsc_checkpoint_read_test
	$(RM) petsc_checkpoint_write_test
	$(RM) execution_resources_test execution_resources_serial_test
	$(RM) memory_report_test
	$(RM) immersed_distributed_assembly_test weighted_work_partition_test
	$(RM) immersed_distributed_history_test
	$(RM) immersed_distributed_transient_volume_test
	$(RM) immersed_transient_distributed_operator_test
	$(RM) immersed_conservation_roundoff_test immersed_transient_distributed_runtime_test immersed_transient_distributed_domain_test
	$(RM) immersed_distributed_physics_test immersed_distributed_static_flow_test immersed_distributed_domain_test
	$(RM) text_helper_failure_test
	$(RM) geometry_inspection_failure_test
	$(RM) petsc_gather_test
	$(RM) collective_petsc_options_test
	$(RM) collective_failure_test
	$(RM) parallel_ownership_test
	$(RM) phase_profile_test body_fitted_subcommunicator_test
	$(RM) immersed_small_cut_solver_test legacy_transport_failure_test vtkhdf_close_failure_test flow_step_failure_test flow_output_failure_test text_output_failure_test iga_pack iga_inspect iga_case_check iga_config_check iga_flow_validate iga_womersley_reference case_config_test simulation_config_test boundary_support_test coupling_port_test sha256_test surface_geometry_test prescribed_surface_motion_test material_surface_patch_kinematics_test moving_cut_geometry_test fluid_surface_traction_test moving_immersed_transient_flow_fsi_runtime_test compliant_channel_fsi_test pretensioned_membrane_test pretensioned_membrane_fsi_runtime_test strong_fluid_structure_coupling_test surface_spatial_index_test cartesian_domain_classification_test cut_cell_volume_quadrature_test immersed_surface_quadrature_test surface_stl_test surface_vtp_test cartesian_background_test quadrature_test vtk_output_test bezier_visualization_test temporal_vtkhdf_test flow_checkpoint_test transport_checkpoint_test vca_checkpoint_test outlet_model_test boundary_flow_test womersley_reference_test generic_transport_test navier_stokes_test nitsche_wall_test ghost_penalty_test immersed_transient_state_test immersed_velocity_extension_test immersed_static_flow_test immersed_transient_flow_test moving_immersed_transient_flow_test moving_immersed_flow_snapshot_test moving_immersed_flow_snapshot_publisher_test immersed_flow_port_test immersed_flow_domain_test immersed_aneurysm_jacobian_test phase5_closure_test phase6_aneurysm_runtime_closure_test phase7_lv_closure_test transport_boundary_preflight_test vca_3d_runtime_test vca_3d_smoke_test iga_solve iga_assembly_smoke iga_transport iga_navier_stokes iga_mesh_check

iga_transport: include/PetscSolverOptions.hpp include/CollectiveAssetInput.hpp include/CollectivePetscOptions.hpp include/PetscGather.hpp include/PetscReadArray.hpp

iga_mesh_check iga_assembly_smoke: include/CollectiveAssetInput.hpp include/Sha256.hpp
iga_assembly_smoke: include/CollectivePetscOptions.hpp

legacy_transport_failure_test: tests/test_legacy_transport_failure.cpp include/PetscSolverOptions.hpp src/iga_transport.cpp include/ExecutionResources.hpp include/CollectiveAssetInput.hpp include/CollectivePetscOptions.hpp include/CollectiveFailure.hpp include/OwnedRowAssembler.hpp include/PetscGather.hpp include/PetscReadArray.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ -Wl,--wrap=MatCreateAIJ,--wrap=VecCreateMPI,--wrap=KSPCreate,--wrap=VecScatterCreateToZero,--wrap=MatAssemblyEnd,--wrap=MatZeroRows,--wrap=MatMult,--wrap=VecCopy,--wrap=KSPSetFromOptions,--wrap=KSPSolve,--wrap=VecNorm,--wrap=VecScatterEnd,--wrap=KSPDestroy,--wrap=MatSetValues,--wrap=VecSetValues,--wrap=VecGetArrayRead,--wrap=VecRestoreArrayRead $(PETSC_LIB)

weighted_work_partition_test: tests/test_weighted_work_partition.cpp ../../include/WeightedWorkPartition.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

.PHONY: weighted-work-partition-test
weighted-work-partition-test: weighted_work_partition_test
	./weighted_work_partition_test

test: weighted-work-partition-test

immersed_distributed_domain_test: ../../include/ThreeDImmersedDistributedDomain.hpp

immersed_distributed_domain_test: tests/test_immersed_distributed_domain.cpp ../../include/ThreeDImmersedDistributedFlowDomain.hpp ../../include/ThreeDImmersedFlowDomain.hpp include/ImmersedStaticDistributedRuntime.hpp include/ImmersedStaticDistributedOperator.hpp include/ImmersedDistributedAssembly.hpp include/ImmersedStaticFlowRuntime.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

immersed_aneurysm_jacobian_test phase6_aneurysm_runtime_closure_test: include/ImmersedDistributedNewtonRuntime.hpp ../../include/ThreeDImmersedDistributedFlowDomain.hpp ../cpu/include/ImmersedStaticDistributedRuntime.hpp ../cpu/include/ImmersedStaticDistributedOperator.hpp ../cpu/include/ImmersedDistributedAssembly.hpp ../../include/WeightedWorkPartition.hpp

immersed_aneurysm_jacobian_test phase6_aneurysm_runtime_closure_test: ../../include/ThreeDImmersedDistributedDomain.hpp

immersed_aneurysm_jacobian_test phase6_aneurysm_runtime_closure_test: ../../include/ThreeDImmersedTransientDistributedFlowDomain.hpp include/ImmersedTransientDistributedRuntime.hpp include/ImmersedTransientDistributedOperator.hpp include/ImmersedDistributedTransientVolume.hpp include/ImmersedDistributedVelocityHistory.hpp include/ImmersedTransientFlowRuntime.hpp include/ImmersedTransientState.hpp

body_fitted_accepted_checkpoint_test: tests/test_body_fitted_accepted_checkpoint.cpp ../coupling/tests/test_bifurcation_coupling_smoke.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

body_fitted_checkpoint_bundle_test: tests/test_body_fitted_checkpoint_bundle.cpp tests/test_body_fitted_accepted_checkpoint.cpp ../coupling/tests/test_bifurcation_coupling_smoke.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

petsc_solver_options_test: tests/test_petsc_solver_options.cpp include/PetscSolverOptions.hpp include/CollectivePetscOptions.hpp include/ExecutionResources.hpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

iga_navier_stokes: include/PetscSolverOptions.hpp include/CollectivePetscOptions.hpp

iga_solve: include/PetscSolverOptions.hpp include/RuntimeConstruction.hpp

# Immersed solver snapshots are shared by static, transient and moving runtimes.
$(PETSC_TARGETS): include/PetscSolverOptions.hpp

moving_immersed_solver_options_test: CPPFLAGS += -DIGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING
moving_immersed_solver_options_test: tests/test_moving_immersed_solver_options.cpp tests/test_moving_immersed_transient_flow.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

petsc_factor_matrix_test: tests/test_petsc_factor_matrix.cpp include/PetscSolverOptions.hpp include/CollectivePetscOptions.hpp include/ExecutionResources.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

immersed_small_cut_solver_test: tests/test_immersed_small_cut_solver.cpp tests/test_immersed_distributed_static_flow.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)
