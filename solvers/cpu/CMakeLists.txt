# CPU backend: IGA tools and MPI/PETSc runtimes, native tetrahedral FEM
# runtimes, and their tests.
#
# Generated once from the historical Makefile; edit freely. Helpers are in
# cmake/TubularFlowHelpers.cmake. `make -C <this directory> <target>` forwards
# here, so every former target name and test name keeps working.
set(TFI_DIRECTORY_INCLUDES include ../../include)
set(TFI_DIRECTORY_OPTIONS -Wall -Wextra -Wpedantic)

# Failure-injection libraries loaded with LD_PRELOAD.
tfi_preload_library(text_read_preload.so SOURCES tests/text_read_preload.cpp LIBRARIES dl)

# Tools and runtimes.
tfi_executable(iga_assembly_smoke SOURCES src/iga_assembly_smoke.cpp PETSC)
tfi_executable(iga_bezier_export SOURCES src/iga_bezier_export.cpp HDF5)
tfi_executable(iga_case_check SOURCES src/iga_case_check.cpp)
tfi_executable(iga_config_check SOURCES src/iga_config_check.cpp)
tfi_executable(iga_flow_validate SOURCES src/iga_flow_validate.cpp)
tfi_executable(iga_inspect SOURCES src/iga_inspect.cpp)
tfi_executable(iga_mesh_check SOURCES src/iga_mesh_check.cpp PETSC)
tfi_executable(iga_navier_stokes SOURCES src/iga_navier_stokes.cpp PETSC HDF5)
tfi_executable(iga_navier_stokes_openmp SOURCES src/iga_navier_stokes.cpp PETSC HDF5 OPENMP)
tfi_executable(iga_pack SOURCES src/iga_pack.cpp)
tfi_executable(iga_solve SOURCES src/iga_solve.cpp PETSC HDF5)
tfi_executable(iga_transport SOURCES src/iga_transport.cpp PETSC)
tfi_executable(iga_transport_validate SOURCES src/iga_transport_validate.cpp)
tfi_executable(iga_womersley_reference SOURCES src/iga_womersley_reference.cpp)
tfi_executable(native_tet_cube_dual_tree_fixed_flow
	SOURCES src/native_tet_cube_dual_tree_fixed_flow.cpp
	PETSC)
tfi_executable(native_tet_cube_oxygen SOURCES src/native_tet_cube_oxygen.cpp PETSC)
tfi_executable(native_tet_darcy SOURCES src/native_tet_darcy.cpp PETSC)
tfi_executable(native_tet_flow SOURCES src/native_tet_flow.cpp PETSC)
tfi_executable(native_tet_flow_darcy_chain SOURCES src/native_tet_flow_darcy_chain.cpp PETSC)
tfi_executable(native_tet_flow_export SOURCES src/native_tet_flow_export.cpp HDF5)
tfi_executable(native_tet_fsi_steps SOURCES src/native_tet_fsi_steps.cpp PETSC)
tfi_executable(native_tet_hydraulic_graph
	SOURCES src/native_tet_hydraulic_graph.cpp
	PETSC HDF5
	INCLUDES ../one_d/include)
tfi_executable(native_tet_matching_darcy_smoke
	SOURCES src/native_tet_matching_darcy_smoke.cpp
	PETSC)
tfi_executable(native_tet_matching_mesh_inspect SOURCES src/native_tet_matching_mesh_inspect.cpp)
tfi_executable(native_tet_matching_solved_roi_smoke
	SOURCES src/native_tet_matching_solved_roi_smoke.cpp
	PETSC)
tfi_executable(native_tet_mesh_inspect SOURCES src/native_tet_mesh_inspect.cpp)
tfi_executable(native_tet_species_transport SOURCES src/native_tet_species_transport.cpp PETSC HDF5)
tfi_executable(surface_fem_preflight SOURCES src/surface_fem_preflight.cpp)

# Test executables.
tfi_executable(assembly_stdout_failure_test SOURCES tests/test_tool_stdout_failure.cpp PETSC)
tfi_executable(bezier_point_signature_test SOURCES tests/test_bezier_point_signature.cpp)
tfi_executable(bezier_visualization_test SOURCES tests/test_bezier_visualization.cpp)
tfi_executable(body_fitted_accepted_checkpoint_test
	SOURCES tests/test_body_fitted_accepted_checkpoint.cpp
	PETSC)
tfi_executable(body_fitted_checkpoint_bundle_test
	SOURCES tests/test_body_fitted_checkpoint_bundle.cpp
	PETSC)
tfi_executable(body_fitted_subcommunicator_test
	SOURCES tests/test_body_fitted_subcommunicator.cpp
	PETSC)
tfi_executable(boundary_flow_test SOURCES tests/test_boundary_flow.cpp)
tfi_executable(boundary_support_test SOURCES tests/test_boundary_support.cpp)
tfi_executable(case_config_test SOURCES tests/test_case_config.cpp)
tfi_executable(checked_text_test SOURCES tests/test_checked_text.cpp)
tfi_executable(collective_failure_test SOURCES tests/test_collective_failure.cpp PETSC)
tfi_executable(collective_petsc_options_test SOURCES tests/test_collective_petsc_options.cpp PETSC)
tfi_executable(coupling_port_test SOURCES tests/test_coupling_port.cpp)
tfi_executable(distributed_point_identity_test
	SOURCES tests/test_distributed_point_identity.cpp
	PETSC)
tfi_executable(element_assembly_execution_openmp_test
	SOURCES tests/test_element_assembly_execution.cpp
	PETSC OPENMP
	OPTIONS -UNDEBUG)
tfi_executable(element_assembly_execution_test
	SOURCES tests/test_element_assembly_execution.cpp
	PETSC
	OPTIONS -UNDEBUG)
tfi_executable(execution_resources_serial_test SOURCES tests/test_execution_resources.cpp PETSC)
tfi_executable(execution_resources_test SOURCES tests/test_execution_resources.cpp PETSC OPENMP)
tfi_executable(flow_checkpoint_test SOURCES tests/test_flow_checkpoint.cpp)
tfi_executable(flow_output_failure_test SOURCES tests/test_flow_output_failure.cpp PETSC HDF5)
tfi_executable(flow_stdout_failure_test
	SOURCES tests/test_tool_stdout_failure.cpp
	PETSC HDF5
	DEFINITIONS IGA_TEST_FLOW_STDOUT)
tfi_executable(flow_step_failure_test SOURCES tests/test_flow_step_failure.cpp PETSC HDF5)
tfi_executable(generic_transport_test SOURCES tests/test_generic_transport.cpp PETSC)
tfi_executable(geometry_inspection_failure_test
	SOURCES tests/test_geometry_inspection_failure.cpp
	PETSC)
tfi_executable(legacy_stdout_failure_test
	SOURCES tests/test_tool_stdout_failure.cpp
	PETSC HDF5
	DEFINITIONS IGA_TEST_LEGACY_STDOUT)
tfi_executable(legacy_transport_failure_test
	SOURCES tests/test_legacy_transport_failure.cpp
	PETSC
	LINK_OPTIONS
		-Wl,--wrap=MatCreateAIJ,--wrap=VecCreateMPI,--wrap=KSPCreate,--wrap=VecScatterCreateToZero,--wrap=MatAssemblyEnd,--wrap=MatZeroRows,--wrap=MatMult,--wrap=VecCopy,--wrap=KSPSetFromOptions,--wrap=KSPSolve,--wrap=VecNorm,--wrap=VecScatterEnd,--wrap=KSPDestroy,--wrap=MatSetValues,--wrap=VecSetValues,--wrap=VecGetArrayRead,--wrap=VecRestoreArrayRead)
tfi_executable(memory_report_test SOURCES tests/test_memory_report.cpp PETSC)
tfi_executable(mesh_stdout_failure_test
	SOURCES tests/test_tool_stdout_failure.cpp
	PETSC
	DEFINITIONS IGA_TEST_MESH_STDOUT)
tfi_executable(native_tet_ale_backflow_test
	SOURCES tests/test_native_tet_ale_backflow.cpp
	OPTIONS -UNDEBUG)
tfi_executable(native_tet_ale_boundary_control_test
	SOURCES tests/test_native_tet_ale_boundary_control.cpp)
tfi_executable(native_tet_ale_boundary_test SOURCES tests/test_native_tet_ale_boundary.cpp)
tfi_executable(native_tet_ale_conservation_test SOURCES tests/test_native_tet_ale_conservation.cpp)
tfi_executable(native_tet_ale_dense_runtime_test
	SOURCES tests/test_native_tet_ale_dense_runtime.cpp)
tfi_executable(native_tet_ale_flow_0d_graph_test
	SOURCES tests/test_native_tet_ale_flow_0d_graph.cpp
	PETSC
	OPTIONS -UNDEBUG)
tfi_executable(native_tet_ale_flow_transport_domain_adapter_test
	SOURCES tests/test_native_tet_ale_flow_transport_domain_adapter.cpp
	PETSC
	OPTIONS -UNDEBUG)
tfi_executable(native_tet_ale_graph_petsc_runtime_test
	SOURCES tests/test_native_tet_ale_graph_petsc_runtime.cpp
	PETSC
	OPTIONS -UNDEBUG)
tfi_executable(native_tet_ale_graph_ports_test
	SOURCES tests/test_native_tet_ale_graph_ports.cpp
	OPTIONS -UNDEBUG)
tfi_executable(native_tet_ale_idealized_lv_backflow_flow_test
	SOURCES tests/test_native_tet_ale_idealized_lv_flow.cpp
	PETSC
	DEFINITIONS IGA_NATIVE_LV_BACKFLOW_STABILIZATION)
tfi_executable(native_tet_ale_idealized_lv_flow_test
	SOURCES tests/test_native_tet_ale_idealized_lv_flow.cpp
	PETSC)
tfi_executable(native_tet_ale_idealized_lv_radial_backflow_flow_test
	SOURCES tests/test_native_tet_ale_idealized_lv_flow.cpp
	PETSC
	DEFINITIONS IGA_NATIVE_LV_RADIAL_REFINEMENT IGA_NATIVE_LV_BACKFLOW_STABILIZATION)
tfi_executable(native_tet_ale_idealized_lv_radial_flow_test
	SOURCES tests/test_native_tet_ale_idealized_lv_flow.cpp
	PETSC
	DEFINITIONS IGA_NATIVE_LV_RADIAL_REFINEMENT)
tfi_executable(native_tet_ale_idealized_lv_radial_two_cycle_backflow_flow_test
	SOURCES tests/test_native_tet_ale_idealized_lv_flow.cpp
	PETSC
	DEFINITIONS IGA_NATIVE_LV_RADIAL_REFINEMENT IGA_NATIVE_LV_TWO_CYCLES
		IGA_NATIVE_LV_BACKFLOW_STABILIZATION)
tfi_executable(native_tet_ale_idealized_lv_two_cycle_backflow_flow_test
	SOURCES tests/test_native_tet_ale_idealized_lv_flow.cpp
	PETSC
	DEFINITIONS IGA_NATIVE_LV_TWO_CYCLES IGA_NATIVE_LV_BACKFLOW_STABILIZATION)
tfi_executable(native_tet_ale_idealized_lv_two_shell_backflow_flow_test
	SOURCES tests/test_native_tet_ale_idealized_lv_flow.cpp
	PETSC
	DEFINITIONS IGA_NATIVE_LV_TWO_RADIAL_SHELLS IGA_NATIVE_LV_BACKFLOW_STABILIZATION)
tfi_executable(native_tet_ale_idealized_lv_two_shell_single_step_backflow_flow_test
	SOURCES tests/test_native_tet_ale_idealized_lv_flow.cpp
	PETSC
	DEFINITIONS IGA_NATIVE_LV_TWO_RADIAL_SHELLS IGA_NATIVE_LV_SINGLE_STEP
		IGA_NATIVE_LV_BACKFLOW_STABILIZATION)
tfi_executable(native_tet_ale_immersed_chamber_motion_test
	SOURCES tests/test_native_tet_ale_immersed_chamber_motion.cpp
	OPTIONS -UNDEBUG)
tfi_executable(native_tet_ale_kinematics_test SOURCES tests/test_native_tet_ale_kinematics.cpp)
tfi_executable(native_tet_ale_mesh_motion_test SOURCES tests/test_native_tet_ale_mesh_motion.cpp)
tfi_executable(native_tet_ale_petsc_assembly_test
	SOURCES tests/test_native_tet_ale_petsc_assembly.cpp
	PETSC)
tfi_executable(native_tet_ale_petsc_mesh_runtime_test
	SOURCES tests/test_native_tet_ale_petsc_mesh_runtime.cpp
	PETSC)
tfi_executable(native_tet_ale_petsc_prescribed_runtime_test
	SOURCES tests/test_native_tet_ale_petsc_prescribed_runtime.cpp
	PETSC)
tfi_executable(native_tet_ale_petsc_runtime_test
	SOURCES tests/test_native_tet_ale_petsc_runtime.cpp
	PETSC)
tfi_executable(native_tet_ale_prescribed_runtime_test
	SOURCES tests/test_native_tet_ale_prescribed_runtime.cpp)
tfi_executable(native_tet_ale_solid_fsi_test
	SOURCES tests/test_native_tet_ale_solid_fsi.cpp
	OPTIONS -UNDEBUG)
tfi_executable(native_tet_ale_species_graph_1d_test
	SOURCES tests/test_native_tet_ale_species_graph_1d.cpp
	PETSC
	OPTIONS -UNDEBUG
	INCLUDES ../one_d/include)
tfi_executable(native_tet_ale_step_control_test SOURCES tests/test_native_tet_ale_step_control.cpp)
tfi_executable(native_tet_ale_temporal_convergence_test
	SOURCES tests/test_native_tet_ale_temporal_convergence.cpp)
tfi_executable(native_tet_ale_transient_test SOURCES tests/test_native_tet_ale_transient.cpp)
tfi_executable(native_tet_bifurcation_runtime_test
	SOURCES tests/test_native_tet_bifurcation_runtime.cpp
	PETSC)
tfi_executable(native_tet_compliant_channel_fsi_test
	SOURCES tests/test_native_tet_compliant_channel_fsi.cpp
	OPTIONS -UNDEBUG)
tfi_executable(native_tet_darcy_petsc_test
	SOURCES tests/test_native_tet_darcy_petsc.cpp
	PETSC
	OPTIONS -UNDEBUG)
tfi_executable(native_tet_fem_test SOURCES tests/test_native_tet_fem.cpp)
tfi_executable(native_tet_fixed_spatial_convergence_test
	SOURCES tests/test_native_tet_fixed_spatial_convergence.cpp
	OPTIONS -UNDEBUG)
tfi_executable(native_tet_fixed_temporal_convergence_test
	SOURCES tests/test_native_tet_fixed_temporal_convergence.cpp)
tfi_executable(native_tet_hyperelastic_solid_test
	SOURCES tests/test_native_tet_hyperelastic_solid.cpp)
tfi_executable(native_tet_matching_fsi_interface_test
	SOURCES tests/test_native_tet_matching_fsi_interface.cpp)
tfi_executable(native_tet_mesh_components_test SOURCES tests/test_native_tet_mesh_components.cpp)
tfi_executable(native_tet_mixed_hyperelastic_solid_test
	SOURCES tests/test_native_tet_mixed_hyperelastic_solid.cpp
	OPTIONS -UNDEBUG)
tfi_executable(native_tet_mixed_locking_test
	SOURCES tests/test_native_tet_mixed_locking.cpp
	OPTIONS -UNDEBUG)
tfi_executable(native_tet_mixed_solid_static_solver_test
	SOURCES tests/test_native_tet_mixed_solid_static_solver.cpp
	OPTIONS -UNDEBUG)
tfi_executable(native_tet_moving_species_checkpoint_test
	SOURCES tests/test_native_tet_moving_species_checkpoint.cpp
	PETSC
	OPTIONS -UNDEBUG)
tfi_executable(native_tet_moving_species_petsc_runtime_test
	SOURCES tests/test_native_tet_moving_species_petsc_runtime.cpp
	PETSC
	OPTIONS -UNDEBUG)
tfi_executable(native_tet_moving_species_ports_test
	SOURCES tests/test_native_tet_moving_species_ports.cpp
	OPTIONS -UNDEBUG)
tfi_executable(native_tet_moving_species_spatial_convergence_test
	SOURCES tests/test_native_tet_moving_species_spatial_convergence.cpp
	PETSC
	OPTIONS -UNDEBUG)
tfi_executable(native_tet_moving_species_transport_test
	SOURCES tests/test_native_tet_moving_species_transport.cpp
	OPTIONS -UNDEBUG)
tfi_executable(native_tet_solid_fsi_runtime_test
	SOURCES tests/test_native_tet_solid_fsi_runtime.cpp
	OPTIONS -UNDEBUG)
tfi_executable(native_tet_solid_prestress_test
	SOURCES tests/test_native_tet_solid_prestress.cpp
	OPTIONS -UNDEBUG)
tfi_executable(native_tet_solid_static_solver_test
	SOURCES tests/test_native_tet_solid_static_solver.cpp)
tfi_executable(native_tet_star_radial_refinement_test
	SOURCES tests/test_native_tet_star_radial_refinement.cpp
	OPTIONS -UNDEBUG)
tfi_executable(native_tet_velocity_field_test SOURCES tests/test_native_tet_velocity_field.cpp)
tfi_executable(native_tet_wall_reservoir_exchange_test
	SOURCES tests/test_native_tet_wall_reservoir_exchange.cpp
	PETSC
	OPTIONS -UNDEBUG)
tfi_executable(navier_stokes_test SOURCES tests/test_navier_stokes.cpp PETSC)
tfi_executable(nonnegative_least_squares_test SOURCES tests/test_nonnegative_least_squares.cpp)
tfi_executable(outlet_model_test SOURCES tests/test_outlet_model.cpp)
tfi_executable(parallel_bezier_visualization_test
	SOURCES tests/test_parallel_bezier_visualization.cpp
	PETSC)
tfi_executable(parallel_body_fitted_flow_test
	SOURCES tests/test_parallel_body_fitted_flow.cpp
	PETSC OPENMP
	DEFINITIONS IGA_FLOW_RUNTIME_TESTING
	OPTIONS -UNDEBUG)
tfi_executable(parallel_element_batch_openmp_test
	SOURCES tests/test_parallel_element_batch.cpp
	OPENMP)
tfi_executable(parallel_element_batch_test SOURCES tests/test_parallel_element_batch.cpp)
tfi_executable(parallel_ownership_test SOURCES tests/test_parallel_ownership.cpp PETSC)
tfi_executable(parallel_vtk_output_test SOURCES tests/test_parallel_vtk_output.cpp PETSC)
tfi_executable(partitioned_vtk_output_test SOURCES tests/test_partitioned_vtk_output.cpp)
tfi_executable(petsc_checkpoint_read_test SOURCES tests/test_petsc_checkpoint_read.cpp PETSC)
tfi_executable(petsc_checkpoint_write_test SOURCES tests/test_petsc_checkpoint_write.cpp PETSC)
tfi_executable(petsc_factor_matrix_test SOURCES tests/test_petsc_factor_matrix.cpp PETSC)
tfi_executable(petsc_gather_test SOURCES tests/test_petsc_gather.cpp PETSC)
tfi_executable(petsc_solver_options_test SOURCES tests/test_petsc_solver_options.cpp PETSC)
tfi_executable(phase_profile_test SOURCES tests/test_phase_profile.cpp)
tfi_executable(prescribed_surface_motion_test SOURCES tests/test_prescribed_surface_motion.cpp)
tfi_executable(quadrature_test SOURCES tests/test_quadrature.cpp PETSC)
tfi_executable(reference_state_output_test SOURCES tests/test_reference_state_output.cpp)
tfi_executable(runtime_cleanup_cli_test SOURCES tests/test_runtime_cleanup_cli.cpp PETSC HDF5)
tfi_executable(sha256_test SOURCES tests/test_sha256.cpp)
tfi_executable(simulation_config_test SOURCES tests/test_simulation_config.cpp)
tfi_executable(strong_fluid_structure_coupling_test
	SOURCES tests/test_strong_fluid_structure_coupling.cpp
	OPTIONS -UNDEBUG)
tfi_executable(surface_geometry_test SOURCES tests/test_surface_geometry.cpp)
tfi_executable(surface_spatial_index_test
	SOURCES tests/test_surface_spatial_index.cpp
	DEFINITIONS IGA_EXACT_DYADIC_TESTING)
tfi_executable(surface_stl_test SOURCES tests/test_surface_stl.cpp)
tfi_executable(surface_vtp_test SOURCES tests/test_surface_vtp.cpp)
tfi_executable(temporal_unstructured_vtkhdf_test
	SOURCES tests/test_temporal_unstructured_vtkhdf.cpp
	HDF5
	LINK_OPTIONS -Wl,--wrap=H5Fflush)
tfi_executable(temporal_vtkhdf_test SOURCES tests/test_temporal_vtkhdf.cpp HDF5)
tfi_executable(text_helper_failure_test
	SOURCES tests/test_text_helper_failure.cpp
	PETSC HDF5
	INCLUDES ../one_d/include)
tfi_executable(text_output_failure_test SOURCES tests/test_text_output_failure.cpp)
tfi_executable(transport_boundary_preflight_test
	SOURCES tests/test_transport_boundary_preflight.cpp
	PETSC)
tfi_executable(transport_budget_stream_failure_test
	SOURCES tests/test_serial_tool_stream_failure.cpp
	PETSC)
tfi_executable(transport_budget_test SOURCES tests/test_transport_budget.cpp)
tfi_executable(transport_checkpoint_test SOURCES tests/test_transport_checkpoint.cpp)
tfi_executable(transport_stdout_failure_test
	SOURCES tests/test_tool_stdout_failure.cpp
	PETSC HDF5
	DEFINITIONS IGA_TEST_TRANSPORT_STDOUT)
tfi_executable(vca_3d_runtime_test SOURCES tests/test_vca_3d_runtime.cpp PETSC)
tfi_executable(vca_3d_smoke_test SOURCES tests/test_vca_3d_smoke.cpp MPI)
tfi_executable(vca_checkpoint_test SOURCES tests/test_vca_checkpoint.cpp)
tfi_executable(vtk_output_test SOURCES tests/test_vtk_output.cpp)
tfi_executable(vtkhdf_close_failure_test
	SOURCES tests/test_vtkhdf_close_failure.cpp
	HDF5
	LINK_OPTIONS
		-Wl,--wrap=H5Fflush,--wrap=H5Gclose,--wrap=H5Fclose,--wrap=H5Dclose,--wrap=H5Fget_obj_count)
tfi_executable(weighted_work_partition_test SOURCES tests/test_weighted_work_partition.cpp)
tfi_executable(womersley_reference_test SOURCES tests/test_womersley_reference.cpp)
tfi_executable(womersley_stream_failure_test
	SOURCES tests/test_serial_tool_stream_failure.cpp
	PETSC
	DEFINITIONS IGA_TEST_WOMERSLEY_STREAM)

# Build-only aggregates of the historical Makefile.
tfi_make_alias(all TARGETS iga_pack iga_inspect iga_bezier_export iga_case_check iga_config_check iga_flow_validate iga_transport_validate iga_womersley_reference surface_fem_preflight)
tfi_make_alias(petsc TARGETS iga_assembly_smoke iga_solve iga_transport iga_navier_stokes iga_mesh_check)

# Tests run from this directory with the original shell commands. MAKE
# names the historical make targets that ran them; LABELS select tiers.
tfi_test(checked_text_test COMMAND [==[./checked_text_test]==]
	REQUIRES checked_text_test MAKE checked-text-test test LABELS unit)
tfi_test(iga_config_check COMMAND [==[sh tests/test_iga_config_check.sh ./iga_config_check]==]
	REQUIRES iga_config_check MAKE config-check-dispatch-test test LABELS unit)
tfi_test(native_tet_ale_backflow_test COMMAND [==[./native_tet_ale_backflow_test]==]
	REQUIRES native_tet_ale_backflow_test MAKE native-tet-ale-backflow-test test LABELS unit)
tfi_test(native_tet_ale_boundary_control_test COMMAND [==[./native_tet_ale_boundary_control_test]==]
	REQUIRES native_tet_ale_boundary_control_test MAKE native-tet-ale-boundary-control-test test LABELS unit)
tfi_test(native_tet_ale_boundary_test COMMAND [==[./native_tet_ale_boundary_test]==]
	REQUIRES native_tet_ale_boundary_test MAKE native-tet-ale-boundary-test test LABELS unit)
tfi_test(native_tet_ale_conservation_test COMMAND [==[./native_tet_ale_conservation_test]==]
	REQUIRES native_tet_ale_conservation_test MAKE native-tet-ale-conservation-test test LABELS unit)
tfi_test(native_tet_ale_dense_runtime_test COMMAND [==[./native_tet_ale_dense_runtime_test]==]
	REQUIRES native_tet_ale_dense_runtime_test MAKE native-tet-ale-dense-runtime-test test LABELS unit)
tfi_test(native_tet_ale_graph_ports_test COMMAND [==[./native_tet_ale_graph_ports_test]==]
	REQUIRES native_tet_ale_graph_ports_test MAKE native-tet-ale-graph-ports-test LABELS unit)
tfi_test(native_tet_ale_idealized_lv_backflow_flow_test COMMAND [==[timeout --kill-after=5s 180s mpiexec -np 2 ./native_tet_ale_idealized_lv_backflow_flow_test]==]
	REQUIRES native_tet_ale_idealized_lv_backflow_flow_test MAKE native-tet-ale-idealized-lv-backflow-flow-test LABELS petsc mpi)
tfi_test(native_tet_ale_idealized_lv_flow_test COMMAND [==[timeout --kill-after=5s 180s mpiexec -np 2 ./native_tet_ale_idealized_lv_flow_test]==]
	REQUIRES native_tet_ale_idealized_lv_flow_test MAKE native-tet-ale-idealized-lv-flow-test LABELS petsc mpi)
tfi_test(native_tet_ale_idealized_lv_radial_backflow_flow_test COMMAND [==[timeout --kill-after=10s 600s ./native_tet_ale_idealized_lv_radial_backflow_flow_test]==]
	REQUIRES native_tet_ale_idealized_lv_radial_backflow_flow_test MAKE native-tet-ale-idealized-lv-radial-backflow-flow-test LABELS petsc)
tfi_test(native_tet_ale_idealized_lv_radial_flow_test COMMAND [==[timeout --kill-after=5s 180s mpiexec -np 2 ./native_tet_ale_idealized_lv_radial_flow_test]==]
	REQUIRES native_tet_ale_idealized_lv_radial_flow_test MAKE native-tet-ale-idealized-lv-radial-flow-test LABELS petsc mpi)
tfi_test(native_tet_ale_immersed_chamber_motion_test COMMAND [==[./native_tet_ale_immersed_chamber_motion_test]==]
	REQUIRES native_tet_ale_immersed_chamber_motion_test MAKE native-tet-ale-immersed-chamber-motion-test test LABELS unit)
tfi_test(native_tet_ale_kinematics_test COMMAND [==[./native_tet_ale_kinematics_test]==]
	REQUIRES native_tet_ale_kinematics_test MAKE native-tet-ale-kinematics-test test LABELS unit)
tfi_test(native_tet_ale_mesh_motion_test COMMAND [==[./native_tet_ale_mesh_motion_test]==]
	REQUIRES native_tet_ale_mesh_motion_test MAKE native-tet-ale-mesh-motion-test test LABELS unit)
tfi_test(native_tet_ale_petsc_assembly_test COMMAND [==[timeout --kill-after=5s 90s mpiexec -np 2 ./native_tet_ale_petsc_assembly_test]==]
	REQUIRES native_tet_ale_petsc_assembly_test MAKE native-tet-ale-petsc-assembly-test LABELS petsc mpi)
tfi_test(native_tet_ale_petsc_runtime_test COMMAND [==[timeout --kill-after=5s 90s mpiexec -np 2 ./native_tet_ale_petsc_runtime_test -native_ale_scale true -native_ale_user_schur false -ksp_type fgmres -ksp_rtol 1e-10 -ksp_max_it 300 -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type full -pc_fieldsplit_schur_precondition selfp -fieldsplit_velocity_ksp_type preonly -fieldsplit_velocity_pc_type hypre -fieldsplit_pressure_ksp_type preonly -fieldsplit_pressure_pc_type jacobi]==]
	REQUIRES native_tet_ale_petsc_runtime_test MAKE native-tet-ale-petsc-fieldsplit-test LABELS petsc mpi)
tfi_test(native-tet-ale-petsc-mesh-runtime-test COMMAND [==[test -n "${T4_ALE_MESH}" || (echo "set T4_ALE_MESH=/path/to/volume.msh" >&2; exit 2); timeout --kill-after=5s 300s mpiexec -np 2 ./native_tet_ale_petsc_mesh_runtime_test "${T4_ALE_MESH}" "${T4_ALE_RESULT:-/tmp/native-tet-ale-mesh-result.json}" "${T4_ALE_MODE:-translation}"]==]
	REQUIRES native_tet_ale_petsc_mesh_runtime_test MAKE native-tet-ale-petsc-mesh-runtime-test LABELS petsc mpi manual)
tfi_test(native_tet_ale_petsc_prescribed_runtime_test COMMAND [==[timeout --kill-after=5s 90s mpiexec -np 2 ./native_tet_ale_petsc_prescribed_runtime_test]==]
	REQUIRES native_tet_ale_petsc_prescribed_runtime_test MAKE native-tet-ale-petsc-prescribed-runtime-test LABELS petsc mpi)
tfi_test(native_tet_ale_petsc_runtime_test.mpi COMMAND [==[timeout --kill-after=5s 90s mpiexec -np 2 ./native_tet_ale_petsc_runtime_test]==]
	REQUIRES native_tet_ale_petsc_runtime_test MAKE native-tet-ale-petsc-runtime-test LABELS petsc mpi)
tfi_test(native_tet_ale_prescribed_runtime_test COMMAND [==[./native_tet_ale_prescribed_runtime_test]==]
	REQUIRES native_tet_ale_prescribed_runtime_test MAKE native-tet-ale-prescribed-runtime-test test LABELS unit)
tfi_test(native_tet_ale_solid_fsi_test COMMAND [==[./native_tet_ale_solid_fsi_test]==]
	REQUIRES native_tet_ale_solid_fsi_test MAKE native-tet-ale-solid-fsi-test test LABELS unit)
tfi_test(native_tet_ale_step_control_test COMMAND [==[./native_tet_ale_step_control_test]==]
	REQUIRES native_tet_ale_step_control_test MAKE native-tet-ale-step-control-test test LABELS unit)
tfi_test(native_tet_ale_temporal_convergence_test COMMAND [==[./native_tet_ale_temporal_convergence_test]==]
	REQUIRES native_tet_ale_temporal_convergence_test MAKE native-tet-ale-temporal-convergence-test test LABELS unit)
tfi_test(native_tet_ale_transient_test COMMAND [==[./native_tet_ale_transient_test]==]
	REQUIRES native_tet_ale_transient_test MAKE native-tet-ale-transient-test test LABELS unit)
tfi_test(native-tet-bifurcation-runtime-test COMMAND [==[test -n "${T2_BIFURCATION_MESH}" || (echo "set T2_BIFURCATION_MESH=/path/to/y.msh" >&2; exit 2); timeout --kill-after=5s 300s mpiexec -np "${T2_BIFURCATION_RANKS:-2}" ./native_tet_bifurcation_runtime_test "${T2_BIFURCATION_MESH}" ${T2_BIFURCATION_RESULT:+"$T2_BIFURCATION_RESULT"} -bifurcation_steps "${T2_BIFURCATION_STEPS:-5}" -bifurcation_dt "${T2_BIFURCATION_DT:-0.01}"]==]
	REQUIRES native_tet_bifurcation_runtime_test MAKE native-tet-bifurcation-runtime-test LABELS petsc mpi manual)
tfi_test(native_tet_compliant_channel_fsi_test COMMAND [==[./native_tet_compliant_channel_fsi_test]==]
	REQUIRES native_tet_compliant_channel_fsi_test MAKE native-tet-compliant-channel-fsi-test test LABELS unit)
tfi_test(native_tet_fem_test COMMAND [==[./native_tet_fem_test]==]
	REQUIRES native_tet_fem_test MAKE native-tet-fem-test test LABELS unit)
tfi_test(native_tet_fixed_spatial_convergence_test COMMAND [==[./native_tet_fixed_spatial_convergence_test]==]
	REQUIRES native_tet_fixed_spatial_convergence_test MAKE native-tet-fixed-spatial-convergence-test test LABELS unit)
tfi_test(native_tet_fixed_temporal_convergence_test COMMAND [==[./native_tet_fixed_temporal_convergence_test]==]
	REQUIRES native_tet_fixed_temporal_convergence_test MAKE native-tet-fixed-temporal-convergence-test test LABELS unit)
tfi_test(native_tet_hyperelastic_solid_test COMMAND [==[./native_tet_hyperelastic_solid_test]==]
	REQUIRES native_tet_hyperelastic_solid_test MAKE native-tet-hyperelastic-solid-test test LABELS unit)
tfi_test(native_tet_matching_fsi_interface_test COMMAND [==[./native_tet_matching_fsi_interface_test]==]
	REQUIRES native_tet_matching_fsi_interface_test MAKE native-tet-matching-fsi-interface-test test LABELS unit)
tfi_test(native_tet_mixed_hyperelastic_solid_test COMMAND [==[./native_tet_mixed_hyperelastic_solid_test]==]
	REQUIRES native_tet_mixed_hyperelastic_solid_test MAKE native-tet-mixed-hyperelastic-solid-test test LABELS unit)
tfi_test(native_tet_mixed_locking_test COMMAND [==[./native_tet_mixed_locking_test]==]
	REQUIRES native_tet_mixed_locking_test MAKE native-tet-mixed-locking-test test LABELS unit)
tfi_test(native_tet_mixed_solid_static_solver_test COMMAND [==[./native_tet_mixed_solid_static_solver_test]==]
	REQUIRES native_tet_mixed_solid_static_solver_test MAKE native-tet-mixed-solid-static-solver-test test LABELS unit)
tfi_test(native_tet_moving_species_ports_test COMMAND [==[./native_tet_moving_species_ports_test]==]
	REQUIRES native_tet_moving_species_ports_test MAKE native-tet-moving-species-ports-test LABELS unit)
tfi_test(native_tet_moving_species_transport_test COMMAND [==[./native_tet_moving_species_transport_test]==]
	REQUIRES native_tet_moving_species_transport_test MAKE native-tet-moving-species-transport-test LABELS unit)
tfi_test(native_tet_solid_fsi_runtime_test COMMAND [==[./native_tet_solid_fsi_runtime_test]==]
	REQUIRES native_tet_solid_fsi_runtime_test MAKE native-tet-solid-fsi-runtime-test test LABELS unit)
tfi_test(native_tet_solid_prestress_test COMMAND [==[./native_tet_solid_prestress_test]==]
	REQUIRES native_tet_solid_prestress_test MAKE native-tet-solid-prestress-test test LABELS unit)
tfi_test(native_tet_solid_static_solver_test COMMAND [==[./native_tet_solid_static_solver_test]==]
	REQUIRES native_tet_solid_static_solver_test MAKE native-tet-solid-static-solver-test test LABELS unit)
tfi_test(native_tet_star_radial_refinement_test COMMAND [==[./native_tet_star_radial_refinement_test]==]
	REQUIRES native_tet_star_radial_refinement_test MAKE native-tet-star-radial-refinement-test test LABELS unit)
tfi_test(native_tet_velocity_field_test COMMAND [==[./native_tet_velocity_field_test]==]
	REQUIRES native_tet_velocity_field_test MAKE native-tet-velocity-field-test LABELS unit)
tfi_test(parallel_element_batch_test COMMAND [==[./parallel_element_batch_test]==]
	REQUIRES parallel_element_batch_test MAKE parallel-element-batch-test LABELS unit)
tfi_test(parallel_element_batch_openmp_test COMMAND [==[OMP_THREAD_LIMIT=4 ./parallel_element_batch_openmp_test]==]
	REQUIRES parallel_element_batch_openmp_test MAKE parallel-element-batch-test LABELS unit)
tfi_test(parallel_ownership_test COMMAND [==[OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout 30 @TFI_MPIEXEC@ -np 3 ./parallel_ownership_test]==]
	REQUIRES parallel_ownership_test MAKE parallel-ownership-test LABELS petsc mpi)
tfi_test(quadrature_test COMMAND [==[./quadrature_test]==]
	REQUIRES quadrature_test MAKE petsc-test LABELS petsc)
tfi_test(generic_transport_test COMMAND [==[./generic_transport_test]==]
	REQUIRES generic_transport_test MAKE petsc-test LABELS petsc)
tfi_test(navier_stokes_test COMMAND [==[./navier_stokes_test]==]
	REQUIRES navier_stokes_test MAKE petsc-test LABELS petsc)
tfi_test(transport_boundary_preflight_test COMMAND [==[PATH="@TFI_PETSC_PREFIX@/bin:$PATH" mpiexec -np 2 ./transport_boundary_preflight_test]==]
	REQUIRES transport_boundary_preflight_test MAKE petsc-test LABELS petsc mpi)
tfi_test(vca_3d_runtime_test COMMAND [==[./vca_3d_runtime_test]==]
	REQUIRES vca_3d_runtime_test MAKE petsc-test LABELS petsc)
tfi_test(vca_3d_smoke_test COMMAND [==[PATH="@TFI_PETSC_PREFIX@/bin:$PATH" ./vca_3d_smoke_test]==]
	REQUIRES iga_navier_stokes vca_3d_smoke_test MAKE petsc-test LABELS petsc)
tfi_test(phase_profile_test COMMAND [==[./phase_profile_test]==]
	REQUIRES phase_profile_test MAKE phase-profile-test LABELS unit)
tfi_test(prescribed_surface_motion_test COMMAND [==[./prescribed_surface_motion_test]==]
	REQUIRES prescribed_surface_motion_test MAKE prescribed-surface-motion-test LABELS unit)
tfi_test(strong_fluid_structure_coupling_test COMMAND [==[./strong_fluid_structure_coupling_test]==]
	REQUIRES strong_fluid_structure_coupling_test MAKE strong-fluid-structure-coupling-test LABELS unit)
tfi_test(case_config_test COMMAND [==[./case_config_test]==]
	REQUIRES case_config_test MAKE test LABELS unit)
tfi_test(simulation_config_test COMMAND [==[./simulation_config_test]==]
	REQUIRES simulation_config_test MAKE test LABELS unit)
tfi_test(boundary_support_test COMMAND [==[./boundary_support_test]==]
	REQUIRES boundary_support_test MAKE test LABELS unit)
tfi_test(coupling_port_test COMMAND [==[./coupling_port_test]==]
	REQUIRES coupling_port_test MAKE test LABELS unit)
tfi_test(sha256_test COMMAND [==[./sha256_test]==]
	REQUIRES sha256_test MAKE test LABELS unit)
tfi_test(surface_geometry_test COMMAND [==[./surface_geometry_test]==]
	REQUIRES surface_geometry_test MAKE test LABELS unit)
tfi_test(surface_spatial_index_test COMMAND [==[./surface_spatial_index_test]==]
	REQUIRES surface_spatial_index_test MAKE test LABELS unit)
tfi_test(surface_stl_test COMMAND [==[./surface_stl_test]==]
	REQUIRES surface_stl_test MAKE test LABELS unit)
tfi_test(surface_vtp_test COMMAND [==[./surface_vtp_test]==]
	REQUIRES surface_vtp_test MAKE test LABELS unit)
tfi_test(vtk_output_test COMMAND [==[./vtk_output_test]==]
	REQUIRES vtk_output_test MAKE test LABELS unit)
tfi_test(bezier_visualization_test COMMAND [==[./bezier_visualization_test]==]
	REQUIRES bezier_visualization_test MAKE test LABELS unit)
tfi_test(temporal_vtkhdf_test COMMAND [==[./temporal_vtkhdf_test]==]
	REQUIRES temporal_vtkhdf_test MAKE test LABELS unit)
tfi_test(flow_checkpoint_test COMMAND [==[./flow_checkpoint_test]==]
	REQUIRES flow_checkpoint_test MAKE test LABELS unit)
tfi_test(transport_checkpoint_test COMMAND [==[./transport_checkpoint_test]==]
	REQUIRES transport_checkpoint_test MAKE test LABELS unit)
tfi_test(vca_checkpoint_test COMMAND [==[./vca_checkpoint_test]==]
	REQUIRES vca_checkpoint_test MAKE test LABELS unit)
tfi_test(outlet_model_test COMMAND [==[./outlet_model_test]==]
	REQUIRES outlet_model_test MAKE test LABELS unit)
tfi_test(boundary_flow_test COMMAND [==[./boundary_flow_test]==]
	REQUIRES boundary_flow_test MAKE test LABELS unit)
tfi_test(womersley_reference_test COMMAND [==[./womersley_reference_test]==]
	REQUIRES womersley_reference_test MAKE test LABELS unit)
tfi_test(text_helper_failure_test COMMAND [==[OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout --kill-after=5s 90s mpiexec -np 3 ./text_helper_failure_test]==]
	REQUIRES text_helper_failure_test MAKE text-helper-failure-test LABELS petsc mpi)
tfi_test(text_output_failure_test COMMAND [==[directory=$(mktemp -d) || exit 1; trap 'rm -rf "$directory"' EXIT; timeout 30 ./text_output_failure_test "$directory/results"]==]
	REQUIRES text_output_failure_test MAKE test text-output-failure-test LABELS unit)
tfi_test(transport_budget_test COMMAND [==[./transport_budget_test]==]
	REQUIRES iga_transport_validate transport_budget_test MAKE test transport-budget-test LABELS unit)
tfi_test(vtkhdf_close_failure_test COMMAND [==[directory=$(mktemp -d) || exit 1; trap 'rm -rf "$directory"' EXIT; timeout 30 ./vtkhdf_close_failure_test "$directory/results"]==]
	REQUIRES vtkhdf_close_failure_test MAKE vtkhdf-close-failure-test LABELS unit)
tfi_test(vtkhdf-paraview-test COMMAND [==[(TUBULARFLOWIGA_KEEP_TEST_OUTPUT=1 ./temporal_vtkhdf_test) && ("${PVPYTHON:-pvpython}" tests/validate_temporal_vtkhdf.py) && ("${PVPYTHON:-pvpython}" tests/validate_temporal_vtkhdf.py /tmp/tubularflowiga-temporal-vtkhdf-test/geometry.vtkhdf --times none --arrays none --value none) && ("${PVPYTHON:-pvpython}" tests/validate_temporal_vtkhdf.py --unit-cube-interpolation)]==]
	REQUIRES temporal_vtkhdf_test MAKE vtkhdf-paraview-test LABELS paraview)
tfi_test(weighted_work_partition_test COMMAND [==[./weighted_work_partition_test]==]
	REQUIRES weighted_work_partition_test MAKE test weighted-work-partition-test LABELS unit)
