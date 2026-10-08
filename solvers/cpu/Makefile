CXX ?= g++
.DEFAULT_GOAL := all
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS ?= -Iinclude -I../../include
SPECIES_COUPLING_HEADER := ../../include/SpeciesCoupling.hpp
LDLIBS ?= -lstdc++fs
# Keep headers and libraries from the same active Conda installation when
# pkg-config is unavailable; otherwise use the system HDF5 development files.
HDF5_CONDA_PREFIX ?= $(firstword $(foreach prefix,$(CONDA_PREFIX) $(CONDA_ROOT),\
	$(if $(and $(wildcard $(prefix)/include/hdf5.h),\
		$(wildcard $(prefix)/lib/libhdf5.so $(prefix)/lib/libhdf5.dylib)),$(prefix))))
HDF5_CFLAGS ?= $(if $(HDF5_CONDA_PREFIX),-I$(HDF5_CONDA_PREFIX)/include,\
	$(shell pkg-config --cflags hdf5 2>/dev/null))
HDF5_ALL_CFLAGS ?= $(HDF5_CFLAGS)
HDF5_SERIAL_LIBRARY ?= $(firstword $(wildcard \
	/usr/lib/*/libhdf5_serial.so /usr/lib/*/libhdf5_serial.so.[0-9]* \
	/lib/*/libhdf5_serial.so /lib/*/libhdf5_serial.so.[0-9]* \
	/usr/lib64/libhdf5.so /usr/lib/libhdf5.so /lib64/libhdf5.so /lib/libhdf5.so))
HDF5_CONDA_LIBS = -L$(HDF5_CONDA_PREFIX)/lib -Wl,-rpath,$(HDF5_CONDA_PREFIX)/lib -lhdf5
HDF5_LIBS ?= $(if $(HDF5_CONDA_PREFIX),$(HDF5_CONDA_LIBS),\
	$(if $(HDF5_SERIAL_LIBRARY),$(HDF5_SERIAL_LIBRARY),$(shell pkg-config --libs hdf5 2>/dev/null)))

PETSC_DIR ?=
PETSC_ARCH ?=
PETSC_TARGETS := body_fitted_subcommunicator_test petsc petsc-test quadrature_test generic_transport_test navier_stokes_test transport_boundary_preflight_test vca_3d_runtime_test vca_3d_smoke_test iga_solve iga_assembly_smoke iga_transport iga_navier_stokes iga_mesh_check

PETSC_TARGETS += parallel_ownership_test parallel-ownership-test
PETSC_TARGETS += native_tet_flow
PETSC_TARGETS += native_tet_ale_petsc_assembly_test native-tet-ale-petsc-assembly-test
PETSC_TARGETS += native_tet_ale_petsc_runtime_test native-tet-ale-petsc-runtime-test
PETSC_TARGETS += native_tet_ale_idealized_lv_flow_test native-tet-ale-idealized-lv-flow-test native_tet_ale_idealized_lv_radial_flow_test native-tet-ale-idealized-lv-radial-flow-test native_tet_ale_idealized_lv_backflow_flow_test native-tet-ale-idealized-lv-backflow-flow-test
PETSC_TARGETS += native_tet_ale_idealized_lv_radial_backflow_flow_test native-tet-ale-idealized-lv-radial-backflow-flow-test
PETSC_TARGETS += native_tet_ale_idealized_lv_two_cycle_backflow_flow_test native_tet_ale_idealized_lv_radial_two_cycle_backflow_flow_test
PETSC_TARGETS += native_tet_ale_idealized_lv_two_shell_single_step_backflow_flow_test
PETSC_TARGETS += native_tet_ale_idealized_lv_two_shell_backflow_flow_test
PETSC_TARGETS += native_tet_moving_species_petsc_runtime_test
PETSC_TARGETS += native_tet_wall_reservoir_exchange_test
PETSC_TARGETS += native_tet_darcy_petsc_test
PETSC_TARGETS += native_tet_darcy
PETSC_TARGETS += native_tet_matching_darcy_smoke
PETSC_TARGETS += native_tet_matching_solved_roi_smoke
PETSC_TARGETS += native_tet_cube_dual_tree_fixed_flow
PETSC_TARGETS += native_tet_cube_oxygen
PETSC_TARGETS += native_tet_flow_darcy_chain
PETSC_TARGETS += native_tet_fsi_steps
PETSC_TARGETS += native_tet_moving_species_spatial_convergence_test
PETSC_TARGETS += native_tet_moving_species_checkpoint_test
PETSC_TARGETS += native_tet_ale_graph_petsc_runtime_test
PETSC_TARGETS += native_tet_ale_flow_transport_domain_adapter_test
PETSC_TARGETS += native_tet_ale_species_graph_1d_test
PETSC_TARGETS += native_tet_ale_flow_0d_graph_test
PETSC_TARGETS += native_tet_hydraulic_graph
PETSC_TARGETS += native_tet_species_transport
PETSC_TARGETS += native-tet-ale-petsc-fieldsplit-test
PETSC_TARGETS += native_tet_ale_petsc_prescribed_runtime_test native-tet-ale-petsc-prescribed-runtime-test
PETSC_TARGETS += native_tet_ale_petsc_mesh_runtime_test native-tet-ale-petsc-mesh-runtime-test
PETSC_TARGETS += native_tet_bifurcation_runtime_test native-tet-bifurcation-runtime-test
PETSC_TARGETS += parallel_vtk_output_test
PETSC_TARGETS += distributed_point_identity_test
PETSC_TARGETS += parallel_bezier_visualization_test

parallel_bezier_visualization_test: tests/test_parallel_bezier_visualization.cpp tests/test_bezier_visualization.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

distributed_point_identity_test: tests/test_distributed_point_identity.cpp tests/test_bezier_visualization.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

parallel_vtk_output_test: tests/test_parallel_vtk_output.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

partitioned_vtk_output_test: tests/test_partitioned_vtk_output.cpp ../../include/PartitionedVtkOutput.hpp ../../include/VtkOutput.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

PETSC_TARGETS += collective_failure_test
PETSC_TARGETS += collective_petsc_options_test petsc_solver_options_test petsc_factor_matrix_test
PETSC_TARGETS += petsc_gather_test
PETSC_TARGETS += body_fitted_checkpoint_bundle_test body_fitted_accepted_checkpoint_test
PETSC_TARGETS += petsc_checkpoint_read_test
PETSC_TARGETS += petsc_checkpoint_write_test
PETSC_TARGETS += execution_resources_test execution_resources_serial_test
PETSC_TARGETS += memory_report_test

memory_report_test: tests/test_memory_report.cpp include/MemoryReport.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)
PETSC_TARGETS += element_assembly_execution_test element_assembly_execution_openmp_test
PETSC_TARGETS += geometry_inspection_failure_test
PETSC_TARGETS += parallel_body_fitted_flow_test iga_navier_stokes_openmp
PETSC_TARGETS += flow_output_failure_test
PETSC_TARGETS += flow_step_failure_test legacy_transport_failure_test
PETSC_TARGETS += text_helper_failure_test text-helper-failure-test
PETSC_TARGETS += mesh_stdout_failure_test assembly_stdout_failure_test
PETSC_TARGETS += flow_stdout_failure_test transport_stdout_failure_test legacy_stdout_failure_test runtime_cleanup_cli_test

flow_stdout_failure_test: CPPFLAGS += -DIGA_TEST_FLOW_STDOUT
transport_stdout_failure_test: CPPFLAGS += -DIGA_TEST_TRANSPORT_STDOUT
legacy_stdout_failure_test: CPPFLAGS += -DIGA_TEST_LEGACY_STDOUT
flow_stdout_failure_test transport_stdout_failure_test legacy_stdout_failure_test: tests/test_tool_stdout_failure.cpp src/iga_navier_stokes.cpp src/iga_solve.cpp src/iga_transport.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $(HDF5_ALL_CFLAGS) $< -o $@ $(PETSC_LIB) $(HDF5_LIBS)
PETSC_TARGETS += womersley_stream_failure_test transport_budget_stream_failure_test

text_read_preload.so: tests/text_read_preload.cpp
	$(CXX) -shared -fPIC -Wall -Wextra -Wpedantic $< -o $@ -ldl

womersley_stream_failure_test: CPPFLAGS += -DIGA_TEST_WOMERSLEY_STREAM
womersley_stream_failure_test transport_budget_stream_failure_test: tests/test_serial_tool_stream_failure.cpp tests/StringStreamFailure.hpp src/iga_womersley_reference.cpp src/iga_transport_validate.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

mesh_stdout_failure_test: tests/test_tool_stdout_failure.cpp src/iga_mesh_check.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) -DIGA_TEST_MESH_STDOUT $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

assembly_stdout_failure_test: tests/test_tool_stdout_failure.cpp src/iga_assembly_smoke.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

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

.PHONY: all petsc petsc-test test config-check-dispatch-test vtkhdf-paraview-test clean prescribed-surface-motion-test

# PETSc-backed runtimes share opt-in phase accounting and explicit KSP setup.
$(PETSC_TARGETS): include/PetscPhaseProfile.hpp ../../include/PhaseProfile.hpp
$(PETSC_TARGETS): include/ParallelOwnershipValidation.hpp ../../include/ParallelOwnership.hpp
$(PETSC_TARGETS): include/CollectiveFailure.hpp include/PetscReadArray.hpp include/RuntimeConstruction.hpp include/RuntimeCleanup.hpp
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

parallel_ownership_test: tests/test_parallel_ownership.cpp include/SurfaceGhostKinematics.hpp include/OwnedPointValues.hpp include/DistributedPointIdentity.hpp ../../include/DynamicWeightedAitkenRelaxation.hpp ../../include/AitkenRelaxation.hpp include/OwnedRowAssembler.hpp include/IgaDatabase.hpp include/SurfaceOwnershipValidation.hpp ../../include/DistributedSurfaceInterface.hpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

.PHONY: parallel-ownership-test
parallel-ownership-test: parallel_ownership_test
	OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout 30 $(MPIEXEC) -np 3 ./parallel_ownership_test
strong_fluid_structure_coupling_test: ../../include/PhaseProfile.hpp

runtime_cleanup_cli_test: tests/test_runtime_cleanup_cli.cpp src/iga_navier_stokes.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $(HDF5_ALL_CFLAGS) $< -o $@ $(PETSC_LIB) $(HDF5_LIBS)

all: iga_pack iga_inspect iga_bezier_export iga_case_check iga_config_check iga_flow_validate iga_transport_validate iga_womersley_reference surface_fem_preflight

native_tet_fem_test: tests/test_native_tet_fem.cpp include/NativeTetFem.hpp include/NativeTetBoundaryFlow.hpp include/NativeTetHydraulicVisualization.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

native_tet_mesh_components_test: tests/test_native_tet_mesh_components.cpp include/NativeTetMeshComponents.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

native_tet_mesh_inspect: src/native_tet_mesh_inspect.cpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

native_tet_matching_mesh_inspect: src/native_tet_matching_mesh_inspect.cpp include/NativeTetMatchingInterfaceSource.hpp include/NativeTetVesselTissueSourceMap.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

native_tet_matching_darcy_smoke: src/native_tet_matching_darcy_smoke.cpp include/NativeTetDarcyPetsc.hpp include/NativeTetMatchingInterfaceSource.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

native_tet_matching_solved_roi_smoke: src/native_tet_matching_solved_roi_smoke.cpp include/NativeTetAlePetscRuntime.hpp include/NativeTetAlePetscPreallocation.hpp include/NativeTetBoundaryFlow.hpp include/NativeTetDarcyPetsc.hpp include/NativeTetMatchingInterfaceSource.hpp include/NativeTetMeshComponents.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

native_tet_cube_dual_tree_fixed_flow: src/native_tet_cube_dual_tree_fixed_flow.cpp include/NativeTetAlePetscRuntime.hpp include/NativeTetAlePetscPreallocation.hpp include/NativeTetBoundaryFlow.hpp include/NativeTetDarcyPetsc.hpp include/NativeTetMatchingInterfaceSource.hpp include/NativeTetMatchingFsiInterface.hpp include/NativeTetLinearElasticPetsc.hpp include/NativeTetLinearElasticVisualization.hpp include/NativeTetMeshComponents.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

native_tet_cube_oxygen: src/native_tet_cube_oxygen.cpp $(wildcard include/NativeTet*.hpp include/ParallelVtkOutput.hpp include/Collective*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

native_tet_flow_darcy_chain: src/native_tet_flow_darcy_chain.cpp $(wildcard include/NativeTet*.hpp include/ParallelVtkOutput.hpp include/Collective*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

native_tet_fsi_steps: src/native_tet_fsi_steps.cpp $(wildcard include/NativeTet*.hpp include/ParallelVtkOutput.hpp include/Collective*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

native_tet_hyperelastic_solid_test: tests/test_native_tet_hyperelastic_solid.cpp include/NativeTetHyperelasticSolid.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

native_tet_solid_static_solver_test: tests/test_native_tet_solid_static_solver.cpp include/NativeTetSolidStaticSolver.hpp include/NativeTetHyperelasticSolid.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

native_tet_solid_prestress_test: tests/test_native_tet_solid_prestress.cpp include/NativeTetSolidPrestress.hpp include/NativeTetSolidStaticSolver.hpp include/NativeTetHyperelasticSolid.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

native_tet_mixed_hyperelastic_solid_test: tests/test_native_tet_mixed_hyperelastic_solid.cpp include/NativeTetMixedHyperelasticSolid.hpp include/NativeTetHyperelasticSolid.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

native_tet_mixed_solid_static_solver_test: tests/test_native_tet_mixed_solid_static_solver.cpp include/NativeTetMixedSolidStaticSolver.hpp include/NativeTetMixedHyperelasticSolid.hpp include/NativeTetSolidStaticSolver.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

native_tet_mixed_locking_test: tests/test_native_tet_mixed_locking.cpp include/NativeTetMixedSolidStaticSolver.hpp include/NativeTetMixedHyperelasticSolid.hpp include/NativeTetAleKinematics.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

native_tet_matching_fsi_interface_test: tests/test_native_tet_matching_fsi_interface.cpp include/NativeTetMatchingFsiInterface.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

native_tet_solid_fsi_runtime_test: tests/test_native_tet_solid_fsi_runtime.cpp include/NativeTetSolidFsiRuntime.hpp include/NativeTetSolidStaticSolver.hpp include/NativeTetHyperelasticSolid.hpp include/StrongFluidStructureCoupling.hpp ../../include/FsiDomainRuntime.hpp ../../include/FsiCouplingEdge.hpp ../../include/DistributedSurfaceInterface.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

native_tet_ale_solid_fsi_test: tests/test_native_tet_ale_solid_fsi.cpp include/NativeTetAleFsiRuntime.hpp include/NativeTetSolidFsiRuntime.hpp include/NativeTetAleDenseRuntime.hpp include/NativeTetMatchingFsiInterface.hpp include/StrongFluidStructureCoupling.hpp ../../include/FsiDomainRuntime.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

native_tet_compliant_channel_fsi_test: tests/test_native_tet_compliant_channel_fsi.cpp include/NativeTetFsiCheckpoint.hpp include/NativeTetAleFsiRuntime.hpp include/NativeTetSolidFsiRuntime.hpp include/NativeTetAleDenseRuntime.hpp include/NativeTetMatchingFsiInterface.hpp include/StrongFluidStructureCoupling.hpp ../../include/FsiDomainRuntime.hpp ../../include/CheckpointMetadataCodec.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

native_tet_fixed_spatial_convergence_test: tests/test_native_tet_fixed_spatial_convergence.cpp include/NativeTetAleDenseRuntime.hpp include/NativeTetAleKinematics.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

native_tet_ale_immersed_chamber_motion_test: tests/test_native_tet_ale_immersed_chamber_motion.cpp include/IdealizedLeftVentricleFixture.hpp include/NativeTetAleConservation.hpp include/NativeTetAleKinematics.hpp include/NativeTetAleMeshMotion.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

native_tet_star_radial_refinement_test: tests/test_native_tet_star_radial_refinement.cpp include/NativeTetStarRadialRefinement.hpp include/IdealizedLeftVentricleFixture.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

native_tet_moving_species_transport_test: tests/test_native_tet_moving_species_transport.cpp include/NativeTetMovingSpeciesTransport.hpp include/NativeTetSpeciesVisualization.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

native_tet_moving_species_ports_test: tests/test_native_tet_moving_species_ports.cpp include/NativeTetMovingSpeciesPorts.hpp include/NativeTetMovingSpeciesTransport.hpp include/NativeTetFem.hpp ../../include/CoupledDomainRuntime.hpp ../../include/FlowDomainPortMetadata.hpp ../../include/SimulationGraph.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

native_tet_ale_graph_ports_test: tests/test_native_tet_ale_graph_ports.cpp include/NativeTetAleGraphPorts.hpp include/NativeTetAleBoundaryControl.hpp include/NativeTetAleBackflow.hpp ../../include/CoupledDomainRuntime.hpp ../../include/FlowDomainPortMetadata.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

native_tet_ale_graph_petsc_runtime_test: tests/test_native_tet_ale_graph_petsc_runtime.cpp include/NativeTetAleGraphPorts.hpp include/NativeTetAlePetscRuntime.hpp include/NativeTetAlePetscPreallocation.hpp include/NativeTetAleBoundaryControl.hpp include/NativeTetMovingSpeciesPetscRuntime.hpp include/NativeTetMovingSpeciesPorts.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) -UNDEBUG $< -o $@ $(PETSC_LIB)

native_tet_ale_flow_transport_domain_adapter_test: tests/test_native_tet_ale_flow_transport_domain_adapter.cpp include/NativeTetAleFlowTransportDomainAdapter.hpp include/NativeTetAleGraphPorts.hpp include/NativeTetAlePetscRuntime.hpp include/NativeTetAlePetscPreallocation.hpp include/NativeTetMovingSpeciesCheckpoint.hpp include/NativeTetMovingSpeciesPorts.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) -UNDEBUG $< -o $@ $(PETSC_LIB)

native_tet_ale_species_graph_1d_test: tests/test_native_tet_ale_species_graph_1d.cpp $(wildcard include/NativeTet*.hpp include/Collective*.hpp ../one_d/include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) -I../one_d/include $(CXXFLAGS) $(PETSC_CC_INCLUDES) -UNDEBUG $< -o $@ $(PETSC_LIB)

native_tet_ale_flow_0d_graph_test: tests/test_native_tet_ale_flow_0d_graph.cpp $(wildcard include/NativeTet*.hpp include/Collective*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) -UNDEBUG $< -o $@ $(PETSC_LIB)

native_tet_hydraulic_graph: src/native_tet_hydraulic_graph.cpp include/ParallelVtkOutput.hpp $(wildcard include/NativeTet*.hpp include/Collective*.hpp ../one_d/include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) -I../one_d/include $(CXXFLAGS) $(PETSC_CC_INCLUDES) $(HDF5_ALL_CFLAGS) $< -o $@ $(PETSC_LIB) $(HDF5_LIBS)

native_tet_species_transport: src/native_tet_species_transport.cpp include/ParallelVtkOutput.hpp $(wildcard include/NativeTet*.hpp include/Collective*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $(HDF5_ALL_CFLAGS) $< -o $@ $(PETSC_LIB) $(HDF5_LIBS)

native_tet_moving_species_petsc_runtime_test: tests/test_native_tet_moving_species_petsc_runtime.cpp include/NativeTetMovingSpeciesPetscRuntime.hpp include/NativeTetMovingSpeciesTransport.hpp include/NativeTetRt0Flux.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) -UNDEBUG $< -o $@ $(PETSC_LIB)

native_tet_wall_reservoir_exchange_test: tests/test_native_tet_wall_reservoir_exchange.cpp include/NativeTetWallReservoirExchange.hpp include/NativeTetMovingSpeciesPetscRuntime.hpp include/NativeTetMovingSpeciesTransport.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) -UNDEBUG $< -o $@ $(PETSC_LIB)

native_tet_darcy_petsc_test: tests/test_native_tet_darcy_petsc.cpp include/NativeTetDarcyPetsc.hpp include/NativeTetMovingSpeciesPetscRuntime.hpp include/NativeTetMovingSpeciesTransport.hpp include/NativeTetVesselTissueSourceMap.hpp include/NativeTetMatchingInterfaceSource.hpp include/NativeTetBoundaryFlow.hpp include/NativeTetRt0Flux.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) -UNDEBUG $< -o $@ $(PETSC_LIB)

native_tet_darcy: src/native_tet_darcy.cpp include/NativeTetDarcyVisualization.hpp include/NativeTetDarcyPetsc.hpp include/NativeTetRt0Flux.hpp include/NativeTetFem.hpp include/ParallelVtkOutput.hpp ../../include/PartitionedVtkOutput.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

native_tet_moving_species_spatial_convergence_test: tests/test_native_tet_moving_species_spatial_convergence.cpp include/NativeTetMovingSpeciesPetscRuntime.hpp include/NativeTetMovingSpeciesTransport.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) -UNDEBUG $< -o $@ $(PETSC_LIB)

native_tet_moving_species_checkpoint_test: tests/test_native_tet_moving_species_checkpoint.cpp include/NativeTetMovingSpeciesCheckpoint.hpp include/NativeTetMovingSpeciesPetscRuntime.hpp include/NativeTetMovingSpeciesTransport.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) -UNDEBUG $< -o $@ $(PETSC_LIB)

.PHONY: native-tet-moving-species-transport-test
native-tet-moving-species-transport-test: native_tet_moving_species_transport_test
	./native_tet_moving_species_transport_test

.PHONY: native-tet-moving-species-ports-test
native-tet-moving-species-ports-test: native_tet_moving_species_ports_test
	./native_tet_moving_species_ports_test

.PHONY: native-tet-ale-graph-ports-test
native-tet-ale-graph-ports-test: native_tet_ale_graph_ports_test
	./native_tet_ale_graph_ports_test

native_tet_ale_kinematics_test: tests/test_native_tet_ale_kinematics.cpp include/NativeTetAleKinematics.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

native_tet_velocity_field_test: tests/test_native_tet_velocity_field.cpp include/NativeTetVelocityField.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

.PHONY: native-tet-velocity-field-test
native-tet-velocity-field-test: native_tet_velocity_field_test
	./native_tet_velocity_field_test

native_tet_ale_mesh_motion_test: tests/test_native_tet_ale_mesh_motion.cpp include/NativeTetAleMeshMotion.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

native_tet_ale_boundary_test: tests/test_native_tet_ale_boundary.cpp include/NativeTetAleBoundary.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

native_tet_ale_boundary_control_test: tests/test_native_tet_ale_boundary_control.cpp include/NativeTetAleBoundaryControl.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

native_tet_ale_backflow_test: tests/test_native_tet_ale_backflow.cpp include/NativeTetAleBackflow.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@

native_tet_ale_step_control_test: tests/test_native_tet_ale_step_control.cpp include/NativeTetAleStepControl.hpp include/NativeTetAleKinematics.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

native_tet_ale_conservation_test: tests/test_native_tet_ale_conservation.cpp include/NativeTetAleConservation.hpp include/NativeTetAleKinematics.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

native_tet_ale_transient_test: tests/test_native_tet_ale_transient.cpp include/NativeTetAleTransient.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

native_tet_ale_dense_runtime_test: tests/test_native_tet_ale_dense_runtime.cpp include/NativeTetAleDenseRuntime.hpp include/NativeTetAleTransient.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

native_tet_ale_prescribed_runtime_test: tests/test_native_tet_ale_prescribed_runtime.cpp include/NativeTetAleDenseRuntime.hpp include/NativeTetAleTransient.hpp include/NativeTetAleConservation.hpp include/NativeTetAleKinematics.hpp include/NativeTetAleMeshMotion.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

native_tet_ale_temporal_convergence_test: tests/test_native_tet_ale_temporal_convergence.cpp include/NativeTetAleDenseRuntime.hpp include/NativeTetAleTransient.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

native_tet_fixed_temporal_convergence_test: tests/test_native_tet_fixed_temporal_convergence.cpp include/NativeTetAleDenseRuntime.hpp include/NativeTetAleTransient.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

native_tet_flow: src/native_tet_flow.cpp include/NativeTetFem.hpp ../../include/CaseConfig.hpp ../../include/CheckedText.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

native_tet_ale_petsc_assembly_test: tests/test_native_tet_ale_petsc_assembly.cpp include/NativeTetAlePetscAssembly.hpp include/NativeTetAleTransient.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

native_tet_ale_petsc_runtime_test: tests/test_native_tet_ale_petsc_runtime.cpp include/NativeTetAlePetscRuntime.hpp include/NativeTetAlePetscPreallocation.hpp include/NativeTetAlePetscAssembly.hpp include/NativeTetAleTransient.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

native_tet_ale_idealized_lv_flow_test: tests/test_native_tet_ale_idealized_lv_flow.cpp include/IdealizedLeftVentricleFixture.hpp include/NativeTetAlePetscRuntime.hpp include/NativeTetAlePetscPreallocation.hpp include/NativeTetBoundaryFlow.hpp include/NativeTetStarRadialRefinement.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

native_tet_ale_idealized_lv_radial_flow_test: tests/test_native_tet_ale_idealized_lv_flow.cpp include/IdealizedLeftVentricleFixture.hpp include/NativeTetAlePetscRuntime.hpp include/NativeTetAlePetscPreallocation.hpp include/NativeTetBoundaryFlow.hpp include/NativeTetStarRadialRefinement.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) -DIGA_NATIVE_LV_RADIAL_REFINEMENT $< -o $@ $(PETSC_LIB)

native_tet_ale_idealized_lv_backflow_flow_test: tests/test_native_tet_ale_idealized_lv_flow.cpp include/IdealizedLeftVentricleFixture.hpp include/NativeTetAlePetscRuntime.hpp include/NativeTetAlePetscPreallocation.hpp include/NativeTetAleBackflow.hpp include/NativeTetBoundaryFlow.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) -DIGA_NATIVE_LV_BACKFLOW_STABILIZATION $< -o $@ $(PETSC_LIB)

native_tet_ale_idealized_lv_radial_backflow_flow_test: tests/test_native_tet_ale_idealized_lv_flow.cpp include/IdealizedLeftVentricleFixture.hpp include/NativeTetAlePetscRuntime.hpp include/NativeTetAlePetscPreallocation.hpp include/NativeTetAleBackflow.hpp include/NativeTetBoundaryFlow.hpp include/NativeTetStarRadialRefinement.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) -DIGA_NATIVE_LV_RADIAL_REFINEMENT -DIGA_NATIVE_LV_BACKFLOW_STABILIZATION $< -o $@ $(PETSC_LIB)

native_tet_ale_idealized_lv_two_cycle_backflow_flow_test: tests/test_native_tet_ale_idealized_lv_flow.cpp include/IdealizedLeftVentricleFixture.hpp include/NativeTetAlePetscRuntime.hpp include/NativeTetAlePetscPreallocation.hpp include/NativeTetAleBackflow.hpp include/NativeTetBoundaryFlow.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) -DIGA_NATIVE_LV_TWO_CYCLES -DIGA_NATIVE_LV_BACKFLOW_STABILIZATION $< -o $@ $(PETSC_LIB)

native_tet_ale_idealized_lv_radial_two_cycle_backflow_flow_test: tests/test_native_tet_ale_idealized_lv_flow.cpp include/IdealizedLeftVentricleFixture.hpp include/NativeTetAlePetscRuntime.hpp include/NativeTetAlePetscPreallocation.hpp include/NativeTetAleBackflow.hpp include/NativeTetBoundaryFlow.hpp include/NativeTetStarRadialRefinement.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) -DIGA_NATIVE_LV_RADIAL_REFINEMENT -DIGA_NATIVE_LV_TWO_CYCLES -DIGA_NATIVE_LV_BACKFLOW_STABILIZATION $< -o $@ $(PETSC_LIB)

native_tet_ale_idealized_lv_two_shell_single_step_backflow_flow_test: tests/test_native_tet_ale_idealized_lv_flow.cpp include/IdealizedLeftVentricleFixture.hpp include/NativeTetAlePetscRuntime.hpp include/NativeTetAlePetscPreallocation.hpp include/NativeTetAleBackflow.hpp include/NativeTetBoundaryFlow.hpp include/NativeTetStarRadialRefinement.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) -DIGA_NATIVE_LV_TWO_RADIAL_SHELLS -DIGA_NATIVE_LV_SINGLE_STEP -DIGA_NATIVE_LV_BACKFLOW_STABILIZATION $< -o $@ $(PETSC_LIB)

native_tet_ale_idealized_lv_two_shell_backflow_flow_test: tests/test_native_tet_ale_idealized_lv_flow.cpp include/IdealizedLeftVentricleFixture.hpp include/NativeTetAlePetscRuntime.hpp include/NativeTetAlePetscPreallocation.hpp include/NativeTetAleBackflow.hpp include/NativeTetBoundaryFlow.hpp include/NativeTetStarRadialRefinement.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) -DIGA_NATIVE_LV_TWO_RADIAL_SHELLS -DIGA_NATIVE_LV_BACKFLOW_STABILIZATION $< -o $@ $(PETSC_LIB)

.PHONY: native-tet-ale-idealized-lv-radial-backflow-flow-test
native-tet-ale-idealized-lv-radial-backflow-flow-test: native_tet_ale_idealized_lv_radial_backflow_flow_test
	timeout --kill-after=10s 600s ./native_tet_ale_idealized_lv_radial_backflow_flow_test

native_tet_ale_petsc_runtime_test native_tet_ale_petsc_prescribed_runtime_test native_tet_ale_petsc_mesh_runtime_test native_tet_ale_idealized_lv_flow_test native_tet_ale_idealized_lv_radial_flow_test native_tet_ale_idealized_lv_backflow_flow_test native_tet_ale_idealized_lv_radial_backflow_flow_test native_tet_ale_idealized_lv_two_cycle_backflow_flow_test native_tet_ale_idealized_lv_radial_two_cycle_backflow_flow_test native_tet_ale_idealized_lv_two_shell_single_step_backflow_flow_test native_tet_ale_idealized_lv_two_shell_backflow_flow_test: include/NativeTetAleBackflow.hpp include/NativeTetAleBoundaryControl.hpp

native_tet_ale_petsc_prescribed_runtime_test: tests/test_native_tet_ale_petsc_prescribed_runtime.cpp include/NativeTetAlePetscRuntime.hpp include/NativeTetAlePetscPreallocation.hpp include/NativeTetAlePetscAssembly.hpp include/NativeTetAleConservation.hpp include/NativeTetAleKinematics.hpp include/NativeTetAleMeshMotion.hpp include/NativeTetAleTransient.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

native_tet_ale_petsc_mesh_runtime_test: tests/test_native_tet_ale_petsc_mesh_runtime.cpp include/NativeTetAlePetscRuntime.hpp include/NativeTetAlePetscPreallocation.hpp include/NativeTetAlePetscAssembly.hpp include/NativeTetAleConservation.hpp include/NativeTetAleKinematics.hpp include/NativeTetAleTransient.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

native_tet_bifurcation_runtime_test: tests/test_native_tet_bifurcation_runtime.cpp include/NativeTetBoundaryFlow.hpp include/NativeTetDarcyPetsc.hpp include/NativeTetRt0Flux.hpp include/NativeTetVesselTissueSourceMap.hpp include/NativeTetMatchingInterfaceSource.hpp include/NativeTetAlePetscRuntime.hpp include/NativeTetAlePetscPreallocation.hpp include/NativeTetAlePetscAssembly.hpp include/NativeTetAleTransient.hpp include/NativeTetFem.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

.PHONY: native-tet-ale-petsc-assembly-test
native-tet-ale-petsc-assembly-test: native_tet_ale_petsc_assembly_test
	timeout --kill-after=5s 90s mpiexec -np 2 ./native_tet_ale_petsc_assembly_test

.PHONY: native-tet-ale-petsc-runtime-test
native-tet-ale-petsc-runtime-test: native_tet_ale_petsc_runtime_test
	timeout --kill-after=5s 90s mpiexec -np 2 ./native_tet_ale_petsc_runtime_test

.PHONY: native-tet-ale-idealized-lv-flow-test
native-tet-ale-idealized-lv-flow-test: native_tet_ale_idealized_lv_flow_test
	timeout --kill-after=5s 180s mpiexec -np 2 ./native_tet_ale_idealized_lv_flow_test

.PHONY: native-tet-ale-idealized-lv-radial-flow-test
native-tet-ale-idealized-lv-radial-flow-test: native_tet_ale_idealized_lv_radial_flow_test
	timeout --kill-after=5s 180s mpiexec -np 2 ./native_tet_ale_idealized_lv_radial_flow_test

.PHONY: native-tet-ale-idealized-lv-backflow-flow-test
native-tet-ale-idealized-lv-backflow-flow-test: native_tet_ale_idealized_lv_backflow_flow_test
	timeout --kill-after=5s 180s mpiexec -np 2 ./native_tet_ale_idealized_lv_backflow_flow_test

.PHONY: native-tet-ale-petsc-fieldsplit-test
native-tet-ale-petsc-fieldsplit-test: native_tet_ale_petsc_runtime_test
	timeout --kill-after=5s 90s mpiexec -np 2 ./native_tet_ale_petsc_runtime_test -native_ale_scale true -native_ale_user_schur false -ksp_type fgmres -ksp_rtol 1e-10 -ksp_max_it 300 -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type full -pc_fieldsplit_schur_precondition selfp -fieldsplit_velocity_ksp_type preonly -fieldsplit_velocity_pc_type hypre -fieldsplit_pressure_ksp_type preonly -fieldsplit_pressure_pc_type jacobi

.PHONY: native-tet-ale-petsc-prescribed-runtime-test
native-tet-ale-petsc-prescribed-runtime-test: native_tet_ale_petsc_prescribed_runtime_test
	timeout --kill-after=5s 90s mpiexec -np 2 ./native_tet_ale_petsc_prescribed_runtime_test

.PHONY: native-tet-ale-petsc-mesh-runtime-test
native-tet-ale-petsc-mesh-runtime-test: native_tet_ale_petsc_mesh_runtime_test
	test -n "$(T4_ALE_MESH)" || (echo "set T4_ALE_MESH=/path/to/volume.msh" >&2; exit 2)
	timeout --kill-after=5s 300s mpiexec -np 2 ./native_tet_ale_petsc_mesh_runtime_test "$(T4_ALE_MESH)" $(if $(strip $(T4_ALE_RESULT)),"$(T4_ALE_RESULT)","/tmp/native-tet-ale-mesh-result.json") "$(if $(strip $(T4_ALE_MODE)),$(T4_ALE_MODE),translation)"

.PHONY: native-tet-bifurcation-runtime-test
native-tet-bifurcation-runtime-test: native_tet_bifurcation_runtime_test
	test -n "$(T2_BIFURCATION_MESH)" || (echo "set T2_BIFURCATION_MESH=/path/to/y.msh" >&2; exit 2)
	timeout --kill-after=5s 300s mpiexec -np $(if $(strip $(T2_BIFURCATION_RANKS)),$(T2_BIFURCATION_RANKS),2) ./native_tet_bifurcation_runtime_test "$(T2_BIFURCATION_MESH)" $(if $(strip $(T2_BIFURCATION_RESULT)),"$(T2_BIFURCATION_RESULT)") -bifurcation_steps $(if $(strip $(T2_BIFURCATION_STEPS)),$(T2_BIFURCATION_STEPS),5) -bifurcation_dt $(if $(strip $(T2_BIFURCATION_DT)),$(T2_BIFURCATION_DT),0.01)

.PHONY: native-tet-fem-test
native-tet-fem-test: native_tet_fem_test
	./native_tet_fem_test

test: native-tet-fem-test

.PHONY: native-tet-hyperelastic-solid-test
native-tet-hyperelastic-solid-test: native_tet_hyperelastic_solid_test
	./native_tet_hyperelastic_solid_test

test: native-tet-hyperelastic-solid-test

.PHONY: native-tet-solid-static-solver-test
native-tet-solid-static-solver-test: native_tet_solid_static_solver_test
	./native_tet_solid_static_solver_test

test: native-tet-solid-static-solver-test

.PHONY: native-tet-matching-fsi-interface-test
native-tet-matching-fsi-interface-test: native_tet_matching_fsi_interface_test
	./native_tet_matching_fsi_interface_test

test: native-tet-matching-fsi-interface-test

.PHONY: native-tet-solid-fsi-runtime-test
native-tet-solid-fsi-runtime-test: native_tet_solid_fsi_runtime_test
	./native_tet_solid_fsi_runtime_test

test: native-tet-solid-fsi-runtime-test

.PHONY: native-tet-solid-prestress-test
native-tet-solid-prestress-test: native_tet_solid_prestress_test
	./native_tet_solid_prestress_test

test: native-tet-solid-prestress-test

.PHONY: native-tet-mixed-hyperelastic-solid-test
native-tet-mixed-hyperelastic-solid-test: native_tet_mixed_hyperelastic_solid_test
	./native_tet_mixed_hyperelastic_solid_test

test: native-tet-mixed-hyperelastic-solid-test

.PHONY: native-tet-mixed-solid-static-solver-test
native-tet-mixed-solid-static-solver-test: native_tet_mixed_solid_static_solver_test
	./native_tet_mixed_solid_static_solver_test

test: native-tet-mixed-solid-static-solver-test

.PHONY: native-tet-mixed-locking-test
native-tet-mixed-locking-test: native_tet_mixed_locking_test
	./native_tet_mixed_locking_test

test: native-tet-mixed-locking-test

.PHONY: native-tet-ale-solid-fsi-test
native-tet-ale-solid-fsi-test: native_tet_ale_solid_fsi_test
	./native_tet_ale_solid_fsi_test

test: native-tet-ale-solid-fsi-test

.PHONY: native-tet-compliant-channel-fsi-test
native-tet-compliant-channel-fsi-test: native_tet_compliant_channel_fsi_test
	./native_tet_compliant_channel_fsi_test

test: native-tet-compliant-channel-fsi-test

.PHONY: native-tet-fixed-spatial-convergence-test
native-tet-fixed-spatial-convergence-test: native_tet_fixed_spatial_convergence_test
	./native_tet_fixed_spatial_convergence_test

test: native-tet-fixed-spatial-convergence-test

.PHONY: native-tet-ale-immersed-chamber-motion-test
native-tet-ale-immersed-chamber-motion-test: native_tet_ale_immersed_chamber_motion_test
	./native_tet_ale_immersed_chamber_motion_test

.PHONY: native-tet-star-radial-refinement-test
native-tet-star-radial-refinement-test: native_tet_star_radial_refinement_test
	./native_tet_star_radial_refinement_test

test: native-tet-star-radial-refinement-test

test: native-tet-ale-immersed-chamber-motion-test

.PHONY: native-tet-ale-kinematics-test
native-tet-ale-kinematics-test: native_tet_ale_kinematics_test
	./native_tet_ale_kinematics_test

test: native-tet-ale-kinematics-test

.PHONY: native-tet-ale-mesh-motion-test
native-tet-ale-mesh-motion-test: native_tet_ale_mesh_motion_test
	./native_tet_ale_mesh_motion_test

test: native-tet-ale-mesh-motion-test

.PHONY: native-tet-ale-boundary-test
native-tet-ale-boundary-test: native_tet_ale_boundary_test
	./native_tet_ale_boundary_test

test: native-tet-ale-boundary-test

.PHONY: native-tet-ale-boundary-control-test
native-tet-ale-boundary-control-test: native_tet_ale_boundary_control_test
	./native_tet_ale_boundary_control_test

test: native-tet-ale-boundary-control-test

.PHONY: native-tet-ale-backflow-test
native-tet-ale-backflow-test: native_tet_ale_backflow_test
	./native_tet_ale_backflow_test

test: native-tet-ale-backflow-test

.PHONY: native-tet-ale-step-control-test
native-tet-ale-step-control-test: native_tet_ale_step_control_test
	./native_tet_ale_step_control_test

test: native-tet-ale-step-control-test

.PHONY: native-tet-ale-conservation-test
native-tet-ale-conservation-test: native_tet_ale_conservation_test
	./native_tet_ale_conservation_test

test: native-tet-ale-conservation-test

.PHONY: native-tet-ale-transient-test
native-tet-ale-transient-test: native_tet_ale_transient_test
	./native_tet_ale_transient_test

test: native-tet-ale-transient-test

.PHONY: native-tet-ale-dense-runtime-test
native-tet-ale-dense-runtime-test: native_tet_ale_dense_runtime_test
	./native_tet_ale_dense_runtime_test

test: native-tet-ale-dense-runtime-test

.PHONY: native-tet-ale-prescribed-runtime-test
native-tet-ale-prescribed-runtime-test: native_tet_ale_prescribed_runtime_test
	./native_tet_ale_prescribed_runtime_test

test: native-tet-ale-prescribed-runtime-test

.PHONY: native-tet-ale-temporal-convergence-test
native-tet-ale-temporal-convergence-test: native_tet_ale_temporal_convergence_test
	./native_tet_ale_temporal_convergence_test

test: native-tet-ale-temporal-convergence-test

.PHONY: native-tet-fixed-temporal-convergence-test
native-tet-fixed-temporal-convergence-test: native_tet_fixed_temporal_convergence_test
	./native_tet_fixed_temporal_convergence_test

test: native-tet-fixed-temporal-convergence-test

iga_inspect iga_flow_validate: ../../include/CheckedText.hpp
iga_inspect: include/Sha256.hpp

surface_fem_preflight: src/surface_fem_preflight.cpp include/SurfaceReaders.hpp include/SurfaceGeometry.hpp include/ExactDyadicPredicates.hpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

reference_state_output_test: tests/test_reference_state_output.cpp tests/ReferenceStateOutput.hpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

element_assembly_execution_openmp_test: CXXFLAGS += -fopenmp
element_assembly_execution_test element_assembly_execution_openmp_test: tests/test_element_assembly_execution.cpp
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

iga_bezier_export: src/iga_bezier_export.cpp include/IgaDatabase.hpp ../../include/BezierVisualization.hpp ../../include/TemporalVtkHdf.hpp ../../include/VtkOutput.hpp
	$(CXX) $(CPPFLAGS) $(HDF5_ALL_CFLAGS) $(CXXFLAGS) $< -o $@ $(HDF5_LIBS) $(LDLIBS)

iga_case_check: src/iga_case_check.cpp include/IgaDatabase.hpp include/CaseInput.hpp include/GenericCaseInput.hpp ../../include/CaseConfig.hpp ../../include/SimulationConfig.hpp ../../include/MeshConfig.hpp ../../include/TemporalFunction.hpp ../../include/OutletModel.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
iga_config_check: src/iga_config_check.cpp ../../include/MultidomainConfig.hpp ../../include/FlowDomainPortMetadata.hpp ../../include/SimulationGraph.hpp ../../include/CouplingEdge.hpp ../../include/CouplingPort.hpp ../../include/SpeciesCoupling.hpp ../../include/SimulationConfig.hpp ../../include/MeshConfig.hpp ../../include/OneDConfig.hpp ../../include/RadiusAnnotatedObj.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

config-check-dispatch-test: iga_config_check tests/test_iga_config_check.sh
	sh tests/test_iga_config_check.sh ./iga_config_check

iga_flow_validate: src/iga_flow_validate.cpp include/BoundaryFlow.hpp include/CaseInput.hpp include/IgaDatabase.hpp include/WomersleyReference.hpp ../../include/VelocitySeries.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

iga_transport_validate: src/iga_transport_validate.cpp include/TransportBudget.hpp include/BoundaryFlow.hpp include/GenericCaseInput.hpp include/CaseInput.hpp include/IgaDatabase.hpp include/Quadrature.hpp include/ElementGeometry.hpp ../../include/SimulationConfig.hpp ../../include/MeshConfig.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

transport_budget_test: tests/test_transport_budget.cpp include/TransportBudget.hpp include/BoundaryFlow.hpp include/GenericCaseInput.hpp include/CaseInput.hpp include/IgaDatabase.hpp include/Quadrature.hpp include/ElementGeometry.hpp ../../include/SimulationConfig.hpp ../../include/MeshConfig.hpp ../../include/CaseConfig.hpp
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

generic_transport_test: tests/test_generic_transport.cpp include/GenericTransportElement.hpp include/FieldCouplingPattern.hpp include/TransportElement.hpp include/CaseInput.hpp ../../include/SimulationConfig.hpp ../../include/MeshConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

navier_stokes_test: tests/test_navier_stokes.cpp include/SurfaceSpatialIndex.hpp include/SurfaceGeometry.hpp include/NavierStokesElement.hpp include/BoundaryFlow.hpp include/TransportElement.hpp include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

test: case_config_test simulation_config_test boundary_support_test coupling_port_test sha256_test surface_geometry_test surface_spatial_index_test surface_stl_test surface_vtp_test vtk_output_test bezier_visualization_test temporal_vtkhdf_test flow_checkpoint_test transport_checkpoint_test vca_checkpoint_test outlet_model_test boundary_flow_test womersley_reference_test config-check-dispatch-test
	./case_config_test
	./simulation_config_test
	./boundary_support_test
	./coupling_port_test
	./sha256_test
	./surface_geometry_test
	./surface_spatial_index_test
	./surface_stl_test
	./surface_vtp_test
	./vtk_output_test
	./bezier_visualization_test
	./temporal_vtkhdf_test
	./flow_checkpoint_test
	./transport_checkpoint_test
	./vca_checkpoint_test
	./outlet_model_test
	./boundary_flow_test
	./womersley_reference_test

sha256_test: tests/test_sha256.cpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

nonnegative_least_squares_test: tests/test_nonnegative_least_squares.cpp include/NonnegativeLeastSquares.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

surface_geometry_test: tests/test_surface_geometry.cpp include/SurfaceGeometry.hpp include/ExactDyadicPredicates.hpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

prescribed_surface_motion_test: tests/test_prescribed_surface_motion.cpp include/MaterialSurfaceKinematics.hpp include/PrescribedSurfaceMotion.hpp include/SurfaceGeometry.hpp include/ExactDyadicPredicates.hpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

prescribed-surface-motion-test: prescribed_surface_motion_test
	./prescribed_surface_motion_test

strong_fluid_structure_coupling_test: tests/test_strong_fluid_structure_coupling.cpp include/StrongFluidStructureCoupling.hpp ../../include/DynamicWeightedAitkenRelaxation.hpp ../../include/FsiCouplingEdge.hpp ../../include/DistributedSurfaceInterface.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $< -o $@ $(LDLIBS)

strong-fluid-structure-coupling-test: strong_fluid_structure_coupling_test
	./strong_fluid_structure_coupling_test

surface_spatial_index_test: tests/test_surface_spatial_index.cpp include/SurfaceSpatialIndex.hpp include/SurfaceGeometry.hpp include/ExactDyadicPredicates.hpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

surface_spatial_index_test: CPPFLAGS += -DIGA_EXACT_DYADIC_TESTING

surface_stl_test: tests/test_surface_stl.cpp include/SurfaceReaders.hpp include/SurfaceGeometry.hpp include/ExactDyadicPredicates.hpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

surface_vtp_test: tests/test_surface_vtp.cpp include/SurfaceReaders.hpp include/SurfaceGeometry.hpp include/ExactDyadicPredicates.hpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

vtk_output_test: tests/test_vtk_output.cpp ../../include/VtkOutput.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

text_output_failure_test: tests/test_text_output_failure.cpp ../../include/VtkOutput.hpp ../../include/PhysiologyOutput.hpp ../../include/VelocitySeries.hpp ../../include/SimulationConfig.hpp ../../include/MeshConfig.hpp ../../include/OxygenCapacity.hpp ../../include/BezierVisualization.hpp ../../include/CouplingHistory.hpp ../../include/VascularCoupling.hpp include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

.PHONY: text-output-failure-test
text-output-failure-test: text_output_failure_test
	@directory=$$(mktemp -d) || exit 1; trap 'rm -rf "$$directory"' EXIT; timeout 30 ./text_output_failure_test "$$directory/results"

test: text-output-failure-test

bezier_point_signature_test: tests/test_bezier_point_signature.cpp tests/test_bezier_visualization.cpp include/IgaDatabase.hpp ../../include/BezierVisualization.hpp ../../include/VtkOutput.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

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
	$(PVPYTHON) tests/validate_temporal_vtkhdf.py /tmp/tubularflowiga-temporal-vtkhdf-test/geometry.vtkhdf --times none --arrays none --value none
	$(PVPYTHON) tests/validate_temporal_vtkhdf.py --unit-cube-interpolation

flow_checkpoint_test: tests/test_flow_checkpoint.cpp ../../include/FlowCheckpoint.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

transport_checkpoint_test: tests/test_transport_checkpoint.cpp ../../include/TransportCheckpoint.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

vca_checkpoint_test: tests/test_vca_checkpoint.cpp ../../include/VcaCheckpoint.hpp ../../include/VascularCoupling.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

outlet_model_test: tests/test_outlet_model.cpp ../../include/OutletModel.hpp ../../include/OutletCheckpoint.hpp ../../include/FlowCheckpoint.hpp ../../include/SimulationConfig.hpp ../../include/MeshConfig.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

boundary_flow_test: tests/test_boundary_flow.cpp include/BoundaryFlow.hpp include/OutletFlow.hpp include/PressureTractionFlow.hpp include/IgaDatabase.hpp ../../include/OutletModel.hpp ../../include/PressureTraction.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

womersley_reference_test: tests/test_womersley_reference.cpp include/WomersleyReference.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

petsc-test: quadrature_test generic_transport_test navier_stokes_test transport_boundary_preflight_test vca_3d_runtime_test vca_3d_smoke_test
	./quadrature_test
	./generic_transport_test
	./navier_stokes_test
	PATH="$(PETSC_PREFIX)/bin:$$PATH" mpiexec -np 2 ./transport_boundary_preflight_test
	./vca_3d_runtime_test
	PATH="$(PETSC_PREFIX)/bin:$$PATH" ./vca_3d_smoke_test

vca_3d_runtime_test: tests/test_vca_3d_runtime.cpp include/IgaDatabase.hpp include/OwnedRowAssembler.hpp include/FieldCouplingPattern.hpp include/GenericCaseInput.hpp include/GenericTransportElement.hpp include/TransientFlowRuntime.hpp include/TransientTransportRuntime.hpp ../../include/ThreeDBodyFittedFlowDomainAdapter.hpp ../../include/ThreeDBodyFittedFlowTransportDomainAdapter.hpp ../../include/FlowDomainPortMetadata.hpp ../../include/ThreeDFlowCoupling.hpp ../../include/CoupledDomainRuntime.hpp ../../include/CouplingPort.hpp $(SPECIES_COUPLING_HEADER) ../../include/SimulationConfig.hpp ../../include/MeshConfig.hpp ../../include/TemporalFunction.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

vca_3d_smoke_test: tests/test_vca_3d_smoke.cpp iga_navier_stokes include/IgaDatabase.hpp ../../include/SimulationConfig.hpp ../../include/MeshConfig.hpp ../../include/VcaCheckpoint.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

transport_boundary_preflight_test: tests/test_transport_boundary_preflight.cpp include/IgaDatabase.hpp include/OwnedRowAssembler.hpp include/TransportBoundaryPreflight.hpp include/TransportElement.hpp ../../include/SimulationConfig.hpp ../../include/MeshConfig.hpp
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

iga_solve: src/iga_solve.cpp ../../include/PhaseProfile.hpp include/IgaDatabase.hpp include/MemoryReport.hpp include/OwnedRowAssembler.hpp include/TransportBoundaryPreflight.hpp include/FieldCouplingPattern.hpp include/CaseInput.hpp include/GenericCaseInput.hpp include/GenericTransportElement.hpp include/PetscBezierVisualization.hpp include/ParallelVtkOutput.hpp ../../include/SimulationConfig.hpp ../../include/MeshConfig.hpp ../../include/TemporalFunction.hpp ../../include/VelocitySeries.hpp ../../include/FlowCheckpoint.hpp ../../include/TransportCheckpoint.hpp ../../include/VtkOutput.hpp ../../include/PhysiologyOutput.hpp ../../include/BezierVisualization.hpp ../../include/PartitionedVtkOutput.hpp ../../include/TemporalVtkHdf.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $(HDF5_CFLAGS) $< -o $@ $(PETSC_LIB) $(HDF5_LIBS)

iga_transport: src/iga_transport.cpp include/IgaDatabase.hpp include/OwnedRowAssembler.hpp include/FieldCouplingPattern.hpp include/CaseInput.hpp include/TransportElement.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

iga_navier_stokes_openmp: CXXFLAGS += -fopenmp
iga_navier_stokes iga_navier_stokes_openmp: include/CollectiveAssetInput.hpp include/Sha256.hpp include/PetscCheckpointRead.hpp include/PetscCheckpointWrite.hpp
iga_navier_stokes iga_navier_stokes_openmp: src/iga_navier_stokes.cpp include/BoundarySupport.hpp include/IgaDatabase.hpp include/OwnedRowAssembler.hpp include/CaseInput.hpp include/GenericCaseInput.hpp include/GenericTransportElement.hpp include/TransientFlowRuntime.hpp include/TransientTransportRuntime.hpp include/TransportElement.hpp include/NavierStokesElement.hpp include/BoundaryFlow.hpp include/PressureTractionFlow.hpp ../../include/CaseConfig.hpp ../../include/CouplingPort.hpp $(SPECIES_COUPLING_HEADER) ../../include/SimulationConfig.hpp ../../include/MeshConfig.hpp ../../include/ThreeDFlowCoupling.hpp ../../include/ThreeDVcaCoupling.hpp ../../include/VascularCoupling.hpp ../../include/CouplingConfig.hpp ../../include/TemporalFunction.hpp ../../include/FlowCheckpoint.hpp ../../include/VcaCheckpoint.hpp ../../include/CouplingHistory.hpp ../../include/OutletModel.hpp ../../include/OutletCheckpoint.hpp ../../include/PressureTraction.hpp ../../include/VelocitySeries.hpp ../../include/VtkOutput.hpp ../../include/BezierVisualization.hpp ../../include/TemporalVtkHdf.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $(HDF5_CFLAGS) $< -o $@ $(PETSC_LIB) $(HDF5_LIBS)

iga_navier_stokes iga_navier_stokes_openmp iga_solve iga_config_check simulation_config_test: ../../include/MeshConfig.hpp

iga_mesh_check: src/iga_mesh_check.cpp include/IgaDatabase.hpp include/TransportElement.hpp include/CollectiveAssetInput.hpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

iga_solve: include/CollectiveAssetInput.hpp include/CollectivePetscOptions.hpp include/MemoryReport.hpp include/PetscCheckpointRead.hpp include/PetscCheckpointWrite.hpp include/PetscGather.hpp include/PetscReadArray.hpp include/Sha256.hpp

clean:
	$(RM) native_tet_flow native_tet_fem_test native_tet_mesh_components_test native_tet_mesh_inspect native_tet_matching_mesh_inspect native_tet_matching_darcy_smoke native_tet_matching_solved_roi_smoke
	$(RM) parallel_bezier_visualization_test
	$(RM) distributed_point_identity_test
	$(RM) bezier_point_signature_test
	$(RM) parallel_vtk_output_test
	$(RM) partitioned_vtk_output_test
	$(RM) nonnegative_least_squares_test
	$(RM) petsc_solver_options_test
	$(RM) body_fitted_checkpoint_bundle_test
	$(RM) body_fitted_accepted_checkpoint_test
	$(RM) runtime_cleanup_cli_test
	$(RM) mesh_stdout_failure_test assembly_stdout_failure_test
	$(RM) flow_stdout_failure_test transport_stdout_failure_test legacy_stdout_failure_test
	$(RM) womersley_stream_failure_test transport_budget_stream_failure_test text_read_preload.so
	$(RM) checked_text_test
	$(RM) parallel_body_fitted_flow_test iga_navier_stokes_openmp element_assembly_execution_test element_assembly_execution_openmp_test
	$(RM) iga_transport_validate transport_budget_test reference_state_output_test parallel_element_batch_test parallel_element_batch_openmp_test
	$(RM) petsc_checkpoint_read_test
	$(RM) petsc_checkpoint_write_test
	$(RM) execution_resources_test execution_resources_serial_test
	$(RM) memory_report_test
	$(RM) weighted_work_partition_test
	$(RM) text_helper_failure_test
	$(RM) geometry_inspection_failure_test
	$(RM) petsc_gather_test
	$(RM) collective_petsc_options_test
	$(RM) collective_failure_test
	$(RM) parallel_ownership_test
	$(RM) phase_profile_test body_fitted_subcommunicator_test
	$(RM) legacy_transport_failure_test vtkhdf_close_failure_test flow_step_failure_test flow_output_failure_test text_output_failure_test iga_pack iga_inspect iga_bezier_export iga_case_check iga_config_check iga_flow_validate iga_womersley_reference surface_fem_preflight case_config_test simulation_config_test boundary_support_test coupling_port_test sha256_test surface_geometry_test prescribed_surface_motion_test strong_fluid_structure_coupling_test surface_spatial_index_test surface_stl_test surface_vtp_test quadrature_test vtk_output_test bezier_visualization_test temporal_vtkhdf_test flow_checkpoint_test transport_checkpoint_test vca_checkpoint_test outlet_model_test boundary_flow_test womersley_reference_test generic_transport_test navier_stokes_test transport_boundary_preflight_test vca_3d_runtime_test vca_3d_smoke_test native_tet_hyperelastic_solid_test native_tet_solid_static_solver_test native_tet_solid_prestress_test native_tet_mixed_hyperelastic_solid_test native_tet_mixed_locking_test native_tet_matching_fsi_interface_test native_tet_solid_fsi_runtime_test native_tet_ale_solid_fsi_test native_tet_compliant_channel_fsi_test native_tet_ale_idealized_lv_flow_test native_tet_star_radial_refinement_test native_tet_ale_boundary_control_test native_tet_ale_immersed_chamber_motion_test native_tet_fixed_spatial_convergence_test native_tet_fixed_temporal_convergence_test iga_solve iga_assembly_smoke iga_transport iga_navier_stokes iga_mesh_check

	$(RM) native_tet_mixed_solid_static_solver_test native_tet_ale_idealized_lv_radial_flow_test native_tet_ale_idealized_lv_backflow_flow_test native_tet_ale_backflow_test
	$(RM) native_tet_ale_idealized_lv_radial_backflow_flow_test native_tet_ale_idealized_lv_two_cycle_backflow_flow_test native_tet_ale_idealized_lv_radial_two_cycle_backflow_flow_test native_tet_ale_idealized_lv_two_shell_single_step_backflow_flow_test native_tet_ale_idealized_lv_two_shell_backflow_flow_test
	$(RM) native_tet_moving_species_transport_test
	$(RM) native_tet_moving_species_ports_test
	$(RM) native_tet_ale_graph_ports_test
	$(RM) native_tet_ale_graph_petsc_runtime_test
	$(RM) native_tet_ale_flow_transport_domain_adapter_test
	$(RM) native_tet_ale_species_graph_1d_test
	$(RM) native_tet_ale_flow_0d_graph_test
	$(RM) native_tet_hydraulic_graph native_tet_flow_export temporal_unstructured_vtkhdf_test
	$(RM) native_tet_flow_darcy_chain
	$(RM) native_tet_fsi_steps
	$(RM) native_tet_species_transport
	$(RM) native_tet_moving_species_petsc_runtime_test
	$(RM) native_tet_wall_reservoir_exchange_test
	$(RM) native_tet_moving_species_spatial_convergence_test
	$(RM) native_tet_moving_species_checkpoint_test

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

body_fitted_accepted_checkpoint_test: tests/test_body_fitted_accepted_checkpoint.cpp ../coupling/tests/test_bifurcation_coupling_smoke.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

body_fitted_checkpoint_bundle_test: tests/test_body_fitted_checkpoint_bundle.cpp tests/test_body_fitted_accepted_checkpoint.cpp ../coupling/tests/test_bifurcation_coupling_smoke.cpp $(wildcard include/*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

petsc_solver_options_test: tests/test_petsc_solver_options.cpp include/PetscSolverOptions.hpp include/CollectivePetscOptions.hpp include/ExecutionResources.hpp include/Sha256.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

iga_navier_stokes: include/PetscSolverOptions.hpp include/CollectivePetscOptions.hpp

iga_solve: include/PetscSolverOptions.hpp include/RuntimeConstruction.hpp

$(PETSC_TARGETS): include/PetscSolverOptions.hpp

petsc_factor_matrix_test: tests/test_petsc_factor_matrix.cpp include/PetscSolverOptions.hpp include/CollectivePetscOptions.hpp include/ExecutionResources.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

native_tet_flow_export: src/native_tet_flow_export.cpp $(wildcard include/NativeTet*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(HDF5_ALL_CFLAGS) $(CXXFLAGS) $< -o $@ $(HDF5_LIBS)

temporal_unstructured_vtkhdf_test: tests/test_temporal_unstructured_vtkhdf.cpp $(wildcard include/NativeTet*.hpp ../../include/*.hpp)
	$(CXX) $(CPPFLAGS) $(HDF5_ALL_CFLAGS) $(CXXFLAGS) $< -o $@ -Wl,--wrap=H5Fflush $(HDF5_LIBS)
