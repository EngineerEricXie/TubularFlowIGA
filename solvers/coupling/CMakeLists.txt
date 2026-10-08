# Multidomain coupling runners and their tests.
#
# Generated once from the historical Makefile; edit freely. Helpers are in
# cmake/TubularFlowHelpers.cmake. `make -C <this directory> <target>` forwards
# here, so every former target name and test name keeps working.
set(TFI_DIRECTORY_INCLUDES ../cpu/include ../one_d/include ../../include)
set(TFI_DIRECTORY_OPTIONS -Wall -Wextra -Wpedantic)

# Failure-injection libraries loaded with LD_PRELOAD.
tfi_preload_library(text_close_preload.so SOURCES tests/text_close_preload.cpp MPI LIBRARIES dl)

# Tools and runtimes.
tfi_executable(hpc_duct_solver_fixture
	SOURCES tests/test_hpc_duct_solver_fixture.cpp
	MPI
	OPTIONS -UNDEBUG)
tfi_executable(iga_1d_3d_bifurcation
	SOURCES src/iga_1d_3d_bifurcation.cpp
	PETSC OPENMP NATIVE_CHECKPOINT_IDENTITY
	DEFINITIONS IGA_BIFURCATION_ENTRY=1)
tfi_executable(iga_1d_3d_explicit SOURCES src/iga_1d_3d_explicit.cpp PETSC OPENMP)
tfi_executable(iga_multidomain_flow
	SOURCES src/iga_1d_3d_bifurcation.cpp
	PETSC OPENMP NATIVE_CHECKPOINT_IDENTITY)

# Test executables.
tfi_executable(bifurcation_cleanup_cli_test
	SOURCES ../cpu/tests/test_runtime_cleanup_cli.cpp
	PETSC OPENMP
	DEFINITIONS IGA_CLEANUP_GRAPH IGA_BIFURCATION_ENTRY=1)
tfi_executable(bifurcation_coupling_smoke_test
	SOURCES tests/test_bifurcation_coupling_smoke.cpp
	MPI)
tfi_executable(bifurcation_stdout_failure_test
	SOURCES ../cpu/tests/test_tool_stdout_failure.cpp
	PETSC OPENMP
	DEFINITIONS IGA_TEST_GRAPH_STDOUT IGA_BIFURCATION_ENTRY=1)
tfi_executable(collective_asset_input_test
	SOURCES tests/test_collective_asset_input.cpp
	PETSC
	DEFINITIONS IGA_ASSET_INPUT_TESTING
	OPTIONS -UNDEBUG)
tfi_executable(coupled_checkpoint_bundle_test
	SOURCES tests/test_coupled_checkpoint_bundle.cpp
	MPI
	DEFINITIONS IGA_COUPLED_CHECKPOINT_TESTING)
tfi_executable(distributed_surface_contracts_test
	SOURCES tests/test_distributed_surface_contracts.cpp
	MPI
	OPTIONS -UNDEBUG)
tfi_executable(domain_group_runtime_proxy_test
	SOURCES tests/test_domain_group_runtime_proxy.cpp
	PETSC)
tfi_executable(explicit_coupling_smoke_test SOURCES tests/test_explicit_coupling_smoke.cpp MPI)
tfi_executable(explicit_coupling_test SOURCES tests/test_explicit_coupling.cpp MPI)
tfi_executable(fsi_runtime_contracts_test
	SOURCES tests/test_fsi_runtime_contracts.cpp
	MPI
	OPTIONS -UNDEBUG)
tfi_executable(graph_checkpoint_history_test SOURCES tests/test_graph_checkpoint_history.cpp MPI)
tfi_executable(graph_cleanup_cli_test
	SOURCES ../cpu/tests/test_runtime_cleanup_cli.cpp
	PETSC OPENMP
	DEFINITIONS IGA_CLEANUP_GRAPH)
tfi_executable(graph_registry_failure_test
	SOURCES tests/test_graph_registry_failure.cpp src/iga_1d_3d_bifurcation.cpp
	PETSC OPENMP
	DEFINITIONS IGA_RUNTIME_CONSTRUCTION_TESTING IGA_MULTIDOMAIN_NO_MAIN
	OPTIONS -UNDEBUG)
tfi_executable(graph_stdout_failure_test
	SOURCES ../cpu/tests/test_tool_stdout_failure.cpp
	PETSC OPENMP
	DEFINITIONS IGA_TEST_GRAPH_STDOUT)
tfi_executable(hpc_duct_partition_test
	SOURCES tests/test_hpc_duct_partition.cpp
	MPI
	OPTIONS -UNDEBUG)
tfi_executable(multidomain_config_test
	SOURCES tests/test_multidomain_config.cpp
	MPI
	OPTIONS -UNDEBUG)
tfi_executable(multidomain_domain_groups_test
	SOURCES tests/test_multidomain_domain_groups.cpp src/iga_1d_3d_bifurcation.cpp
	PETSC OPENMP
	DEFINITIONS IGA_MULTIDOMAIN_NO_MAIN)
tfi_executable(multidomain_failure_test
	SOURCES tests/test_multidomain_failure.cpp src/iga_1d_3d_bifurcation.cpp
	PETSC OPENMP
	DEFINITIONS IGA_MULTIDOMAIN_NO_MAIN)
tfi_executable(multidomain_flow_smoke_test SOURCES tests/test_multidomain_flow_smoke.cpp MPI)
tfi_executable(multidomain_subcommunicator_test
	SOURCES tests/test_multidomain_subcommunicator.cpp src/iga_1d_3d_bifurcation.cpp
	PETSC OPENMP
	DEFINITIONS IGA_MULTIDOMAIN_NO_MAIN)
tfi_executable(native_graph_checkpoint_fault_test
	SOURCES tests/test_native_graph_checkpoint_fault.cpp
	PETSC OPENMP NATIVE_CHECKPOINT_IDENTITY)
tfi_executable(native_graph_checkpoint_fixture_test
	SOURCES tests/test_native_graph_checkpoint_fixture.cpp
	MPI)
tfi_executable(native_graph_checkpoint_groups_test
	SOURCES tests/test_native_graph_checkpoint_groups.cpp src/iga_1d_3d_bifurcation.cpp
	PETSC OPENMP NATIVE_CHECKPOINT_IDENTITY
	DEFINITIONS IGA_MULTIDOMAIN_NO_MAIN)
tfi_executable(phase9_multiscale_closure_test
	SOURCES tests/test_phase9_multiscale_closure.cpp
	MPI
	OPTIONS -UNDEBUG)
tfi_executable(pressure_flow_collective_failure_test
	SOURCES tests/test_pressure_flow_collective_failure.cpp
	PETSC
	OPTIONS -UNDEBUG)
tfi_executable(pressure_flow_executor_test SOURCES tests/test_pressure_flow_executor.cpp MPI)
tfi_executable(pulsatile_straight_vessel_temporal_test
	SOURCES tests/test_pulsatile_straight_vessel_temporal.cpp
	MPI
	OPTIONS -UNDEBUG)
tfi_executable(registry_failure_test
	SOURCES tests/test_registry_failure.cpp
	PETSC
	DEFINITIONS IGA_RUNTIME_CONSTRUCTION_TESTING
	OPTIONS -UNDEBUG)
tfi_executable(runtime_construction_failure_test
	SOURCES tests/test_runtime_construction_failure.cpp
	PETSC OPENMP
	DEFINITIONS IGA_RUNTIME_CONSTRUCTION_TESTING IGA_RUNTIME_CLEANUP_TESTING
	OPTIONS -UNDEBUG)
tfi_executable(sequential_cleanup_cli_test
	SOURCES ../cpu/tests/test_runtime_cleanup_cli.cpp
	PETSC OPENMP
	DEFINITIONS IGA_CLEANUP_SEQUENTIAL)
tfi_executable(sequential_initialization_failure_test
	SOURCES tests/test_sequential_initialization_failure.cpp src/iga_1d_3d_explicit.cpp
	PETSC OPENMP
	DEFINITIONS IGA_RUNTIME_CONSTRUCTION_TESTING IGA_SEQUENTIAL_NO_MAIN
	OPTIONS -UNDEBUG)
tfi_executable(sequential_stdout_failure_test
	SOURCES ../cpu/tests/test_tool_stdout_failure.cpp
	PETSC OPENMP
	DEFINITIONS IGA_TEST_SEQUENTIAL_STDOUT)
tfi_executable(sequential_strong_failure_test
	SOURCES tests/test_sequential_strong_failure.cpp src/iga_1d_3d_explicit.cpp
	PETSC OPENMP
	DEFINITIONS IGA_RUNTIME_CONSTRUCTION_TESTING IGA_SEQUENTIAL_EXECUTION_TESTING
		IGA_SEQUENTIAL_NO_MAIN
	OPTIONS -UNDEBUG)
tfi_executable(sequential_text_failure_test
	SOURCES tests/test_sequential_text_failure.cpp
	PETSC OPENMP
	OPTIONS -UNDEBUG)
tfi_executable(simulation_graph_test SOURCES tests/test_simulation_graph.cpp MPI OPTIONS -UNDEBUG)
tfi_executable(species_collective_failure_test
	SOURCES tests/test_species_collective_failure.cpp
	PETSC
	OPTIONS -UNDEBUG)
tfi_executable(species_pressure_flow_executor_test
	SOURCES tests/test_species_pressure_flow_executor.cpp
	MPI
	OPTIONS -UNDEBUG)
tfi_executable(straight_vessel_convergence_test
	SOURCES tests/test_straight_vessel_convergence.cpp
	MPI
	OPTIONS -UNDEBUG)
tfi_executable(three_d_port_failure_test
	SOURCES tests/test_three_d_port_failure.cpp
	PETSC OPENMP
	DEFINITIONS IGA_FLOW_RUNTIME_TESTING)
tfi_executable(three_d_staged_failure_test
	SOURCES tests/test_three_d_staged_failure.cpp
	PETSC OPENMP
	DEFINITIONS IGA_TRANSPORT_RUNTIME_TESTING)
tfi_executable(three_d_trial_failure_test SOURCES tests/test_three_d_trial_failure.cpp PETSC OPENMP)
tfi_executable(zero_d_checkpoint_test SOURCES tests/test_zero_d_checkpoint.cpp MPI)
tfi_executable(zero_d_flow_species_checkpoint_test
	SOURCES tests/test_zero_d_flow_species_checkpoint.cpp
	MPI
	OPTIONS -UNDEBUG)
tfi_executable(zero_d_flow_test SOURCES tests/test_zero_d_flow.cpp MPI OPTIONS -UNDEBUG)
tfi_executable(zero_d_source_reservoir_species_runtime_test
	SOURCES tests/test_zero_d_source_reservoir_species_runtime.cpp
	MPI
	OPTIONS -UNDEBUG)
tfi_executable(zero_d_source_reservoir_species_test
	SOURCES tests/test_zero_d_source_reservoir_species.cpp
	MPI
	OPTIONS -UNDEBUG)
tfi_executable(zero_d_species_reservoir_test
	SOURCES tests/test_zero_d_species_reservoir.cpp
	MPI
	OPTIONS -UNDEBUG)
tfi_executable(zero_d_terminal_rcr_species_graph_test
	SOURCES tests/test_zero_d_terminal_rcr_species_graph.cpp
	MPI
	OPTIONS -UNDEBUG)
tfi_executable(zero_d_terminal_rcr_species_runtime_test
	SOURCES tests/test_zero_d_terminal_rcr_species_runtime.cpp
	MPI
	OPTIONS -UNDEBUG)
tfi_executable(zero_d_terminal_rcr_species_test
	SOURCES tests/test_zero_d_terminal_rcr_species.cpp
	MPI
	OPTIONS -UNDEBUG)

# Build-only aggregates of the historical Makefile.
tfi_make_alias(all TARGETS explicit_coupling_test simulation_graph_test pressure_flow_executor_test species_pressure_flow_executor_test multidomain_config_test)
tfi_make_alias(petsc TARGETS iga_1d_3d_explicit iga_1d_3d_bifurcation iga_multidomain_flow)

# Tests run from this directory with the original shell commands. MAKE
# names the historical make targets that ran them; LABELS select tiers.
tfi_test(bifurcation_coupling_smoke_test COMMAND [==[PATH="@TFI_PETSC_PREFIX@/bin:$PATH" ./bifurcation_coupling_smoke_test]==]
	REQUIRES bifurcation_coupling_smoke_test iga_1d_3d_bifurcation MAKE bifurcation-test LABELS petsc)
tfi_test(coupled_checkpoint_bundle_test COMMAND [==[./coupled_checkpoint_bundle_test "$(mktemp -d)/checkpoint"]==]
	REQUIRES coupled_checkpoint_bundle_test MAKE checkpoint-bundle-test LABELS unit)
tfi_test(straight_vessel_convergence_test COMMAND [==[PATH="@TFI_PETSC_PREFIX@/bin:$PATH" TUBULARFLOWIGA_AXIAL_ONLY=1 ./straight_vessel_convergence_test]==]
	REQUIRES iga_1d_3d_explicit straight_vessel_convergence_test MAKE convergence-axial-diagnostic LABELS petsc diagnostic)
tfi_test(straight_vessel_convergence_test.2 COMMAND [==[PATH="@TFI_PETSC_PREFIX@/bin:$PATH" TUBULARFLOWIGA_BULK_ONLY=1 ./straight_vessel_convergence_test]==]
	REQUIRES iga_1d_3d_explicit straight_vessel_convergence_test MAKE convergence-bulk-diagnostic LABELS petsc diagnostic)
tfi_test(straight_vessel_convergence_test.3 COMMAND [==[TUBULARFLOWIGA_FIXTURE_ONLY=1 ./straight_vessel_convergence_test]==]
	REQUIRES straight_vessel_convergence_test MAKE convergence-fixture-test LABELS unit)
tfi_test(straight_vessel_convergence_test.4 COMMAND [==[PATH="@TFI_PETSC_PREFIX@/bin:$PATH" TUBULARFLOWIGA_ISOTROPIC_ONLY=1 ./straight_vessel_convergence_test]==]
	REQUIRES iga_1d_3d_explicit straight_vessel_convergence_test MAKE convergence-isotropic-diagnostic LABELS petsc diagnostic)
tfi_test(straight_vessel_convergence_test.5 COMMAND [==[PATH="@TFI_PETSC_PREFIX@/bin:$PATH" TUBULARFLOWIGA_LENGTH_ONLY=1 ./straight_vessel_convergence_test]==]
	REQUIRES iga_1d_3d_explicit straight_vessel_convergence_test MAKE convergence-length-diagnostic LABELS petsc diagnostic)
tfi_test(straight_vessel_convergence_test.6 COMMAND [==[PATH="@TFI_PETSC_PREFIX@/bin:$PATH" ./straight_vessel_convergence_test]==]
	REQUIRES iga_1d_3d_explicit straight_vessel_convergence_test MAKE convergence-test LABELS petsc slow)
tfi_test(fsi_runtime_contracts_test COMMAND [==[./fsi_runtime_contracts_test]==]
	REQUIRES fsi_runtime_contracts_test MAKE fsi-runtime-contracts-test test LABELS unit)
tfi_test(multidomain_flow_smoke_test COMMAND [==[./multidomain_flow_smoke_test --species-only]==]
	REQUIRES iga_1d_3d_bifurcation iga_multidomain_flow multidomain_flow_smoke_test MAKE multidomain-species-test LABELS petsc)
tfi_test(multidomain_flow_smoke_test.2 COMMAND [==[PATH="@TFI_PETSC_PREFIX@/bin:$PATH" ./multidomain_flow_smoke_test]==]
	REQUIRES iga_1d_3d_bifurcation iga_multidomain_flow multidomain_flow_smoke_test MAKE multidomain-test LABELS petsc)
tfi_test(explicit_coupling_smoke_test COMMAND [==[PATH="@TFI_PETSC_PREFIX@/bin:$PATH" ./explicit_coupling_smoke_test]==]
	REQUIRES explicit_coupling_smoke_test iga_1d_3d_explicit MAKE petsc-test LABELS petsc)
tfi_test(phase9_multiscale_closure_test COMMAND [==[PATH="@TFI_PETSC_PREFIX@/bin:$PATH" ./phase9_multiscale_closure_test ./iga_multidomain_flow]==]
	REQUIRES iga_multidomain_flow phase9_multiscale_closure_test MAKE phase9-multiscale-closure-test LABELS petsc)
tfi_test(simulation_graph_test COMMAND [==[./simulation_graph_test]==]
	REQUIRES simulation_graph_test MAKE simulation-graph-test test LABELS unit)
tfi_test(distributed_surface_contracts_test COMMAND [==[./distributed_surface_contracts_test]==]
	REQUIRES distributed_surface_contracts_test MAKE surface-contracts-test test LABELS unit)
tfi_test(pulsatile_straight_vessel_temporal_test COMMAND [==[PATH="@TFI_PETSC_PREFIX@/bin:$PATH" ./pulsatile_straight_vessel_temporal_test]==]
	REQUIRES iga_1d_3d_explicit pulsatile_straight_vessel_temporal_test MAKE temporal-convergence-test LABELS petsc slow)
tfi_test(pulsatile_straight_vessel_temporal_test.2 COMMAND [==[TUBULARFLOWIGA_TEMPORAL_FIXTURE_ONLY=1 ./pulsatile_straight_vessel_temporal_test]==]
	REQUIRES pulsatile_straight_vessel_temporal_test MAKE temporal-fixture-test LABELS unit)
tfi_test(explicit_coupling_test COMMAND [==[./explicit_coupling_test]==]
	REQUIRES explicit_coupling_test MAKE test LABELS unit)
tfi_test(pressure_flow_executor_test COMMAND [==[./pressure_flow_executor_test]==]
	REQUIRES pressure_flow_executor_test MAKE test LABELS unit)
tfi_test(species_pressure_flow_executor_test COMMAND [==[./species_pressure_flow_executor_test]==]
	REQUIRES species_pressure_flow_executor_test MAKE test LABELS unit)
tfi_test(multidomain_config_test COMMAND [==[./multidomain_config_test]==]
	REQUIRES multidomain_config_test MAKE test LABELS unit)
tfi_test(zero_d_flow_test COMMAND [==[./zero_d_flow_test]==]
	REQUIRES zero_d_flow_test MAKE test zero-d-flow-test LABELS unit)
tfi_test(zero_d_species_reservoir_test COMMAND [==[./zero_d_species_reservoir_test]==]
	REQUIRES zero_d_species_reservoir_test MAKE test zero-d-species-reservoir-test LABELS unit)
tfi_test(zero_d_source_reservoir_species_test COMMAND [==[./zero_d_source_reservoir_species_test]==]
	REQUIRES zero_d_source_reservoir_species_test MAKE test zero-d-source-reservoir-species-test LABELS unit)
tfi_test(zero_d_source_reservoir_species_runtime_test COMMAND [==[./zero_d_source_reservoir_species_runtime_test]==]
	REQUIRES zero_d_source_reservoir_species_runtime_test MAKE test zero-d-source-reservoir-species-runtime-test LABELS unit)
tfi_test(zero_d_flow_species_checkpoint_test COMMAND [==[checkpoint_dir=$(mktemp -d); ./zero_d_flow_species_checkpoint_test save "$checkpoint_dir/paired" && ./zero_d_flow_species_checkpoint_test resume "$checkpoint_dir/paired"]==]
	REQUIRES zero_d_flow_species_checkpoint_test MAKE test zero-d-flow-species-checkpoint-test LABELS unit)
tfi_test(zero_d_terminal_rcr_species_test COMMAND [==[./zero_d_terminal_rcr_species_test]==]
	REQUIRES zero_d_terminal_rcr_species_test MAKE test zero-d-terminal-rcr-species-test LABELS unit)
tfi_test(zero_d_terminal_rcr_species_runtime_test COMMAND [==[./zero_d_terminal_rcr_species_runtime_test]==]
	REQUIRES zero_d_terminal_rcr_species_runtime_test MAKE test zero-d-terminal-rcr-species-runtime-test LABELS unit)
tfi_test(zero_d_terminal_rcr_species_graph_test COMMAND [==[./zero_d_terminal_rcr_species_graph_test]==]
	REQUIRES zero_d_terminal_rcr_species_graph_test MAKE test zero-d-terminal-rcr-species-graph-test LABELS unit)

# Coordinated-failure, checkpoint, and MPI tests that no make target ran before.
# Each MPI test needs exactly the rank count shown; "$TMPDIR/out" is a fresh output path.
tfi_test(collective_asset_input_test COMMAND [==[@TFI_MPIEXEC@ -np 3 ./collective_asset_input_test "$TMPDIR/out"]==]
	REQUIRES collective_asset_input_test MAKE petsc-test collective-asset-input-test LABELS petsc mpi)
tfi_test(domain_group_runtime_proxy_test COMMAND [==[@TFI_MPIEXEC@ -np 3 ./domain_group_runtime_proxy_test]==]
	REQUIRES domain_group_runtime_proxy_test MAKE petsc-test domain-group-runtime-proxy-test LABELS petsc mpi)
tfi_test(graph_checkpoint_history_test COMMAND [==[./graph_checkpoint_history_test]==]
	REQUIRES graph_checkpoint_history_test MAKE test graph-checkpoint-history-test LABELS unit)
tfi_test(hpc_duct_partition_test COMMAND [==[./hpc_duct_partition_test "$TMPDIR/out"]==]
	REQUIRES hpc_duct_partition_test MAKE test hpc-duct-partition-test LABELS unit)
tfi_test(multidomain_domain_groups_test COMMAND [==[@TFI_MPIEXEC@ -np 5 ./multidomain_domain_groups_test "$TMPDIR/out"]==]
	REQUIRES multidomain_domain_groups_test MAKE petsc-test multidomain-domain-groups-test LABELS petsc mpi)
tfi_test(pressure_flow_collective_failure_test COMMAND [==[@TFI_MPIEXEC@ -np 3 ./pressure_flow_collective_failure_test]==]
	REQUIRES pressure_flow_collective_failure_test MAKE petsc-test pressure-flow-collective-failure-test LABELS petsc mpi)
tfi_test(registry_failure_test COMMAND [==[@TFI_MPIEXEC@ -np 3 ./registry_failure_test]==]
	REQUIRES registry_failure_test MAKE petsc-test registry-failure-test LABELS petsc mpi)
tfi_test(sequential_text_failure_test COMMAND [==[@TFI_MPIEXEC@ -np 3 ./sequential_text_failure_test "$TMPDIR/out"]==]
	REQUIRES sequential_text_failure_test MAKE petsc-test sequential-text-failure-test LABELS petsc mpi)
tfi_test(species_collective_failure_test COMMAND [==[@TFI_MPIEXEC@ -np 3 ./species_collective_failure_test]==]
	REQUIRES species_collective_failure_test MAKE petsc-test species-collective-failure-test LABELS petsc mpi)
tfi_test(three_d_port_failure_test COMMAND [==[@TFI_MPIEXEC@ -np 3 ./three_d_port_failure_test "$TMPDIR/out"]==]
	REQUIRES three_d_port_failure_test MAKE petsc-test three-d-port-failure-test LABELS petsc mpi)
tfi_test(zero_d_checkpoint_test COMMAND [==[./zero_d_checkpoint_test "$TMPDIR/out"]==]
	REQUIRES zero_d_checkpoint_test MAKE test zero-d-checkpoint-test LABELS unit)

# Multi-rank regression drivers from scripts/hpc that need no retained fixture.
tfi_test(hpc-failure-regression COMMAND [==[python3 ../../scripts/hpc/hpc_failure_regression.py --output-dir "$TMPDIR/out"]==]
	REQUIRES collective_failure_test collective_petsc_options_test multidomain_failure_test
	MAKE hpc-regression-test LABELS petsc mpi)
tfi_test(hpc-solver-views COMMAND [==[python3 ../../scripts/hpc/hpc_solver_views.py --output-dir "$TMPDIR/out"]==]
	REQUIRES iga_multidomain_flow native_graph_checkpoint_fixture_test
	MAKE hpc-regression-test LABELS petsc mpi)
tfi_test(hpc-native-graph-checkpoint COMMAND [==[python3 ../../scripts/hpc/hpc_native_graph_checkpoint.py --output-root "$TMPDIR/out"]==]
	REQUIRES iga_multidomain_flow native_graph_checkpoint_fault_test native_graph_checkpoint_fixture_test
	MAKE hpc-regression-test LABELS petsc mpi slow)
