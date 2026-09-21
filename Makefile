.PHONY: all mesh mesh-test solver-test cpu cpu-test cpu-petsc t1-audit t1-fem-mesh-test t1-multiregion-audit-test t1-iga-manifest-test t1-immersed-manifest-test t2-contract-audit t3-contract-audit t3-solid-contract-audit t4-contract-audit t4-petsc-assembly-test t5-matching-interface-audit hpc-build-manifest hpc-prepare-scaling hpc-cross-node-binaries hpc-test-unit hpc-test-mpi hpc-test-gpu fluid-surface-traction-test moving-immersed-transient-flow-fsi-runtime-test compliant-channel-fsi-test phase8-compliant-channel-fsi-paraview pretensioned-membrane-test pretensioned-membrane-fsi-runtime-test material-surface-patch-kinematics-test phase7-focused-test phase7-lv-closure-test one-d-petsc one-d-test zero-d-petsc zero-d-test coupling coupling-test coupling-zero-d-flow-test coupling-simulation-graph-test coupling-surface-contracts-test coupling-fsi-runtime-contracts-test coupling-petsc-test coupling-convergence-test coupling-convergence-axial-diagnostic coupling-convergence-isotropic-diagnostic coupling-convergence-length-diagnostic coupling-convergence-bulk-diagnostic coupling-temporal-convergence-test coupling-temporal-fixture-test cuda spline workflow-test clean

all: cpu

mesh:
	$(MAKE) -C preprocessing/mesh

mesh-test:
	$(MAKE) -C preprocessing/mesh test

solver-test:
	$(MAKE) -C solvers/cpu surface_fem_preflight
	python3 scripts/tests/test_solver_entry.py
	python3 scripts/tests/test_label_surface_bc_gui.py

cpu:
	$(MAKE) -C solvers/cpu

cpu-test:
	$(MAKE) -C solvers/cpu test

t1-audit: t1-fem-mesh-test t1-ftetwild-mesh-test t1-multiregion-audit-test t1-iga-manifest-test
	python3 scripts/validate_t1_geometry_inventory.py
	python3 scripts/validate_t1_mesher_comparison.py
	python3 scripts/tests/test_t1_mesher_comparison.py
	python3 -m unittest scripts.tests.test_dicom_seg_region_overlap

t1-fem-mesh-test:
	$(MAKE) -C solvers/cpu surface_fem_preflight
	python3 scripts/tests/test_surface_to_fem_volume.py

t1-multiregion-audit-test:
	$(MAKE) -C solvers/cpu native_tet_mesh_inspect native_tet_matching_mesh_inspect
	python3 scripts/tests/test_multiregion_tet_mesh_audit.py
	python3 scripts/tests/test_dicom_seg_multiregion_tet.py
	python3 scripts/tests/test_split_multiregion_tet_for_native.py

.PHONY: t1-ftetwild-mesh-test
t1-ftetwild-mesh-test:
	$(MAKE) -C solvers/cpu surface_fem_preflight native-tet-fem-test
	python3 scripts/tests/test_bent_pipe_surface.py
	python3 scripts/tests/test_y_pipe_surface.py
	python3 scripts/tests/test_ftetwild_to_fem_volume.py

t1-iga-manifest-test:
	$(MAKE) -C solvers/cpu iga_inspect
	python3 scripts/tests/test_iga_geometry_manifest.py

t1-immersed-manifest-test:
	$(MAKE) -C solvers/cpu moving-cut-geometry-test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))

t2-contract-audit:
	$(MAKE) -C solvers/cpu surface_fem_preflight native-tet-fem-test
	python3 scripts/tests/test_circular_pipe_surface.py
	python3 scripts/tests/test_prepare_t2_bifurcation_meshes.py
	python3 scripts/validate_t2_fixed_flow_contract.py
	python3 scripts/tests/test_t2_fixed_flow_contract.py
	python3 scripts/tests/test_t2_fixed_flow_result.py
	python3 scripts/tests/test_native_t2_result_pipeline.py
	python3 scripts/validate_t2_native_functional_evidence.py
	python3 scripts/tests/test_t2_native_functional_evidence.py
	python3 scripts/validate_t2_native_bifurcation_evidence.py
	python3 scripts/validate_t2_native_temporal_evidence.py
	python3 scripts/validate_t2_native_spatial_evidence.py
	$(MAKE) -C solvers/cpu native-tet-fixed-temporal-convergence-test
	$(MAKE) -C solvers/cpu native-tet-fixed-spatial-convergence-test

t3-contract-audit:
	python3 scripts/validate_t3_native_iga_shell_contract.py
	$(MAKE) -C solvers/cpu native-iga-shell-surface-test
	$(MAKE) -C solvers/cpu native-iga-kirchhoff-love-test
	$(MAKE) -C solvers/cpu native-iga-shell-patch-system-test
	$(MAKE) -C solvers/cpu native-iga-shell-static-solver-test
	$(MAKE) -C solvers/cpu native-iga-shell-spatial-convergence-test
	$(MAKE) -C solvers/cpu native-iga-shell-flat-plate-test
	$(MAKE) -C solvers/cpu native-iga-shell-patch-interface-test

t3-solid-contract-audit:
	python3 scripts/validate_t3_native_tet_solid_contract.py
	$(MAKE) -C solvers/cpu native-tet-hyperelastic-solid-test
	$(MAKE) -C solvers/cpu native-tet-solid-static-solver-test
	$(MAKE) -C solvers/cpu native-tet-solid-prestress-test
	$(MAKE) -C solvers/cpu native-tet-mixed-hyperelastic-solid-test
	$(MAKE) -C solvers/cpu native-tet-mixed-solid-static-solver-test
	$(MAKE) -C solvers/cpu native-tet-mixed-locking-test

t5-matching-interface-audit:
	python3 scripts/validate_t5_native_matching_interface_contract.py
	$(MAKE) -C solvers/cpu native-tet-matching-fsi-interface-test
	$(MAKE) -C solvers/cpu native-tet-solid-fsi-runtime-test
	$(MAKE) -C solvers/cpu native-tet-ale-solid-fsi-test
	$(MAKE) -C solvers/cpu native-tet-compliant-channel-fsi-test

t4-contract-audit:
	python3 scripts/validate_t4_native_ale_contract.py
	python3 scripts/validate_t4_native_ale_evidence.py
	python3 scripts/validate_t4_matched_lv_qoi_evidence.py
	$(MAKE) -C solvers/cpu native-tet-ale-kinematics-test
	$(MAKE) -C solvers/cpu native-tet-velocity-field-test
	$(MAKE) -C solvers/cpu native-tet-ale-mesh-motion-test
	$(MAKE) -C solvers/cpu native-tet-ale-boundary-test
	$(MAKE) -C solvers/cpu native-tet-ale-boundary-control-test
	$(MAKE) -C solvers/cpu native-tet-ale-backflow-test
	$(MAKE) -C solvers/cpu native-tet-ale-step-control-test
	$(MAKE) -C solvers/cpu native-tet-ale-conservation-test
	$(MAKE) -C solvers/cpu native-tet-ale-transient-test
	$(MAKE) -C solvers/cpu native-tet-ale-dense-runtime-test
	$(MAKE) -C solvers/cpu native-tet-ale-prescribed-runtime-test
	$(MAKE) -C solvers/cpu native-tet-ale-temporal-convergence-test
	$(MAKE) -C solvers/cpu native-tet-ale-immersed-chamber-motion-test
	$(MAKE) -C solvers/cpu native-tet-star-radial-refinement-test

.PHONY: t6-native-lv-qoi-audit
t6-native-lv-qoi-audit:
	python3 scripts/validate_t6_native_lv_qoi_evidence.py

.PHONY: t7-native-moving-species-test
t7-native-moving-species-test:
	$(MAKE) -C solvers/cpu native-tet-moving-species-transport-test
	$(MAKE) -C solvers/cpu native-tet-moving-species-ports-test
	$(MAKE) -C solvers/cpu native-tet-ale-graph-ports-test
	python3 scripts/validate_t7_native_moving_species_evidence.py

.PHONY: t7-native-moving-species-petsc-build
t7-native-moving-species-petsc-build:
	$(MAKE) -C solvers/cpu native_tet_moving_species_petsc_runtime_test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))

T7_SPECIES_MPI_RANKS ?= 2
.PHONY: t7-native-moving-species-petsc-test
t7-native-moving-species-petsc-test: t7-native-moving-species-petsc-build
	timeout --kill-after=10s 120s mpiexec -np $(T7_SPECIES_MPI_RANKS) ./solvers/cpu/native_tet_moving_species_petsc_runtime_test

.PHONY: t7-native-moving-species-spatial-convergence-test
t7-native-moving-species-spatial-convergence-test:
	$(MAKE) -C solvers/cpu native_tet_moving_species_spatial_convergence_test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	timeout --kill-after=10s 120s mpiexec -np $(T7_SPECIES_MPI_RANKS) ./solvers/cpu/native_tet_moving_species_spatial_convergence_test

.PHONY: t7-native-moving-species-checkpoint-test
t7-native-moving-species-checkpoint-test:
	$(MAKE) -C solvers/cpu native_tet_moving_species_checkpoint_test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	timeout --kill-after=10s 120s mpiexec -np $(T7_SPECIES_MPI_RANKS) ./solvers/cpu/native_tet_moving_species_checkpoint_test

.PHONY: t7-native-ale-graph-petsc-test
t7-native-ale-graph-petsc-test:
	$(MAKE) -C solvers/cpu native_tet_ale_graph_petsc_runtime_test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	timeout --kill-after=10s 120s mpiexec -np $(T7_SPECIES_MPI_RANKS) ./solvers/cpu/native_tet_ale_graph_petsc_runtime_test

.PHONY: t7-native-staged-ale-test
t7-native-staged-ale-test:
	$(MAKE) -C solvers/cpu native_tet_ale_flow_transport_domain_adapter_test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	timeout --kill-after=10s 120s mpiexec -np $(T7_SPECIES_MPI_RANKS) ./solvers/cpu/native_tet_ale_flow_transport_domain_adapter_test

.PHONY: t7-native-ale-species-graph-1d-test
t7-native-ale-species-graph-1d-test:
	$(MAKE) -C solvers/cpu native_tet_ale_species_graph_1d_test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	timeout --kill-after=10s 120s mpiexec -np $(T7_SPECIES_MPI_RANKS) ./solvers/cpu/native_tet_ale_species_graph_1d_test
	timeout --kill-after=10s 120s mpiexec -np $(T7_SPECIES_MPI_RANKS) ./solvers/cpu/native_tet_ale_species_graph_1d_test --channel
	timeout --kill-after=10s 120s mpiexec -np $(T7_SPECIES_MPI_RANKS) ./solvers/cpu/native_tet_ale_species_graph_1d_test --channel-strong

.PHONY: t7-native-ale-species-graph-rcr-test
t7-native-ale-species-graph-rcr-test:
	$(MAKE) -C solvers/cpu native_tet_ale_species_graph_1d_test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	timeout --kill-after=10s 120s mpiexec -np $(T7_SPECIES_MPI_RANKS) ./solvers/cpu/native_tet_ale_species_graph_1d_test --channel-rcr
	timeout --kill-after=10s 120s mpiexec -np $(T7_SPECIES_MPI_RANKS) ./solvers/cpu/native_tet_ale_species_graph_1d_test --channel-rcr-strong
	timeout --kill-after=10s 120s mpiexec -np $(T7_SPECIES_MPI_RANKS) ./solvers/cpu/native_tet_ale_species_graph_1d_test --channel-rcr-distal-backflow
	timeout --kill-after=10s 120s mpiexec -np $(T7_SPECIES_MPI_RANKS) ./solvers/cpu/native_tet_ale_species_graph_1d_test --channel-rcr-distal-backflow-strong

.PHONY: t7-native-ale-species-graph-zero-d-source-test
t7-native-ale-species-graph-zero-d-source-test:
	$(MAKE) -C solvers/cpu native_tet_ale_species_graph_1d_test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	timeout --kill-after=10s 120s mpiexec -np $(T7_SPECIES_MPI_RANKS) ./solvers/cpu/native_tet_ale_species_graph_1d_test --channel-zero-d-source-rcr
	timeout --kill-after=10s 120s mpiexec -np $(T7_SPECIES_MPI_RANKS) ./solvers/cpu/native_tet_ale_species_graph_1d_test --channel-zero-d-source-rcr-strong

.PHONY: t7-native-ale-species-graph-checkpoint-test
t7-native-ale-species-graph-checkpoint-test:
	$(MAKE) -C solvers/cpu native_tet_ale_species_graph_1d_test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	@checkpoint_root=$$(mktemp -d); \
		timeout --kill-after=10s 120s mpiexec -np $(T7_SPECIES_MPI_RANKS) ./solvers/cpu/native_tet_ale_species_graph_1d_test --channel-zero-d-source-rcr-checkpoint-save "$$checkpoint_root" && \
		timeout --kill-after=10s 120s mpiexec -np $(T7_SPECIES_MPI_RANKS) ./solvers/cpu/native_tet_ale_species_graph_1d_test --channel-zero-d-source-rcr-checkpoint-resume "$$checkpoint_root"

.PHONY: t7-native-ale-species-monotone-front-checkpoint-test
t7-native-ale-species-monotone-front-checkpoint-test:
	$(MAKE) -C solvers/cpu native_tet_ale_species_graph_1d_test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	@checkpoint_root=$$(mktemp -d); \
		timeout --kill-after=10s 120s mpiexec -np $(T7_SPECIES_MPI_RANKS) ./solvers/cpu/native_tet_ale_species_graph_1d_test --channel-zero-d-source-rcr-monotone-front-checkpoint-save "$$checkpoint_root" && \
		timeout --kill-after=10s 120s mpiexec -np $(T7_SPECIES_MPI_RANKS) ./solvers/cpu/native_tet_ale_species_graph_1d_test --channel-zero-d-source-rcr-monotone-front-checkpoint-resume "$$checkpoint_root"

.PHONY: t7-native-ale-flow-0d-graph-test
t7-native-ale-flow-0d-graph-test:
	$(MAKE) -C solvers/cpu native_tet_ale_flow_0d_graph_test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	timeout --kill-after=10s 120s mpiexec -np $(T7_SPECIES_MPI_RANKS) ./solvers/cpu/native_tet_ale_flow_0d_graph_test
	timeout --kill-after=10s 120s mpiexec -np $(T7_SPECIES_MPI_RANKS) ./solvers/cpu/native_tet_ale_flow_0d_graph_test --strong
	timeout --kill-after=10s 120s mpiexec -np $(T7_SPECIES_MPI_RANKS) ./solvers/cpu/native_tet_ale_flow_0d_graph_test --balanced-explicit
	timeout --kill-after=10s 120s mpiexec -np $(T7_SPECIES_MPI_RANKS) ./solvers/cpu/native_tet_ale_flow_0d_graph_test --balanced-strong
	python3 scripts/validate_t7_native_0d_graph_evidence.py

.PHONY: t7-native-ale-flow-0d-cross-process-test
t7-native-ale-flow-0d-cross-process-test:
	$(MAKE) -C solvers/cpu native_tet_ale_flow_0d_graph_test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	python3 scripts/run_t7_native_0d_cross_process.py --ranks $(T7_SPECIES_MPI_RANKS) --binary solvers/cpu/native_tet_ale_flow_0d_graph_test
	python3 scripts/validate_t7_native_0d_graph_evidence.py

.PHONY: t7-native-tet-hydraulic-cli-test
t7-native-tet-hydraulic-cli-test:
	$(MAKE) -C solvers/cpu native_tet_hydraulic_graph PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	python3 scripts/test_native_tet_hydraulic_cli.py --ranks $(T7_SPECIES_MPI_RANKS) --binary solvers/cpu/native_tet_hydraulic_graph --case solvers/cpu/tests/data/native_tet_hydraulic_star.json
	python3 scripts/validate_t7_native_hydraulic_cli_evidence.py

.PHONY: t7-native-tet-ftetwild-smoke
t7-native-tet-ftetwild-smoke:
	$(MAKE) -C solvers/cpu native_tet_hydraulic_graph surface_fem_preflight PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	python3 scripts/test_native_tet_hydraulic_ftetwild.py --ranks $(T7_SPECIES_MPI_RANKS) --binary solvers/cpu/native_tet_hydraulic_graph --ftetwild $(FTETWILD_BIN)

.PHONY: t7-native-tet-species-ftetwild-test
t7-native-tet-species-ftetwild-test:
	$(MAKE) -C solvers/cpu native_tet_ale_species_graph_1d_test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	python3 scripts/test_native_tet_species_ftetwild.py --ranks $(T7_SPECIES_MPI_RANKS) --binary solvers/cpu/native_tet_ale_species_graph_1d_test --ftetwild $(FTETWILD_BIN)

.PHONY: t9-native-tet-species-cli-test
t9-native-tet-species-cli-test:
	$(MAKE) -C solvers/cpu native_tet_species_transport native_tet_moving_species_transport_test native_tet_moving_species_petsc_runtime_test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	./solvers/cpu/native_tet_moving_species_transport_test
	mpiexec -np 1 solvers/cpu/native_tet_moving_species_petsc_runtime_test
	mpiexec -np 2 solvers/cpu/native_tet_moving_species_petsc_runtime_test
	mpiexec -np 4 solvers/cpu/native_tet_moving_species_petsc_runtime_test
	python3 scripts/test_native_tet_species_cli.py --binary solvers/cpu/native_tet_species_transport $(if $(strip $(FTETWILD_BIN)),--ftetwild $(FTETWILD_BIN),)
	python3 scripts/validate_t9_native_species_cli_evidence.py

.PHONY: t6-native-wall-reservoir-exchange-test
t6-native-wall-reservoir-exchange-test:
	$(MAKE) -C solvers/cpu native_tet_wall_reservoir_exchange_test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	mpiexec -np 1 solvers/cpu/native_tet_wall_reservoir_exchange_test
	mpiexec -np 2 solvers/cpu/native_tet_wall_reservoir_exchange_test
	mpiexec -np 4 solvers/cpu/native_tet_wall_reservoir_exchange_test
	$(if $(strip $(FTETWILD_BIN)),python3 scripts/test_native_tet_wall_reservoir_ftetwild.py --binary solvers/cpu/native_tet_wall_reservoir_exchange_test --ftetwild $(FTETWILD_BIN) --ranks 2,)
	python3 scripts/validate_t6_native_wall_reservoir_evidence.py

.PHONY: t6-native-tet-darcy-test
t6-native-tet-darcy-test:
	$(MAKE) -C solvers/cpu native_tet_mesh_components_test
	solvers/cpu/native_tet_mesh_components_test
	$(MAKE) -C solvers/cpu native_tet_darcy_petsc_test native_tet_darcy native_tet_matching_darcy_smoke native_tet_matching_solved_roi_smoke surface_fem_preflight PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	mpiexec -np 1 solvers/cpu/native_tet_darcy_petsc_test
	mpiexec -np 2 solvers/cpu/native_tet_darcy_petsc_test
	mpiexec -np 4 solvers/cpu/native_tet_darcy_petsc_test
	$(if $(strip $(FTETWILD_BIN)),python3 scripts/test_native_tet_darcy_ftetwild.py --binary solvers/cpu/native_tet_darcy_petsc_test --ftetwild $(FTETWILD_BIN) --ranks 2,)
	python3 scripts/test_native_tet_darcy_cli.py --binary solvers/cpu/native_tet_darcy $(if $(strip $(FTETWILD_BIN)),--ftetwild $(FTETWILD_BIN),)
	python3 scripts/test_native_tet_darcy_workflow.py --solver solvers/cpu/native_tet_darcy $(if $(strip $(FTETWILD_BIN)),--ftetwild $(FTETWILD_BIN) --gmsh,)
	python3 -m unittest scripts.tests.test_run_native_matching_roi_functional
	python3 -m unittest scripts.tests.test_run_liver_roi_functional_case
	python3 -m unittest scripts.tests.test_liver_roi_facet_ledger
	python3 -m unittest scripts.tests.test_run_liver_roi_from_raw
	python3 -m unittest scripts.tests.test_time_native_mpi_rank
	python3 -m unittest scripts.tests.test_audit_seg_exposed_vessel_faces
	python3 -m unittest scripts.tests.test_liver_roi_summary_bindings
	python3 scripts/validate_t6_native_tet_darcy_evidence.py
	python3 scripts/validate_liver_roi_functional_evidence.py

.PHONY: t7-native-tet-species-ftetwild-same-mesh-mpi-test
t7-native-tet-species-ftetwild-same-mesh-mpi-test:
	$(MAKE) -C solvers/cpu native_tet_ale_species_graph_1d_test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	python3 scripts/test_native_tet_species_ftetwild.py --compare-ranks --binary solvers/cpu/native_tet_ale_species_graph_1d_test --ftetwild $(FTETWILD_BIN)

.PHONY: t9-native-tet-workflow-test
t9-native-tet-workflow-test:
	$(MAKE) -C solvers/cpu native_tet_hydraulic_graph surface_fem_preflight PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	python3 scripts/test_native_tet_workflow.py --solver solvers/cpu/native_tet_hydraulic_graph $(if $(strip $(FTETWILD_BIN)),--ftetwild $(FTETWILD_BIN) --gmsh,)

.PHONY: t9-native-tet-species-workflow-test
t9-native-tet-species-workflow-test:
	$(MAKE) -C solvers/cpu native_tet_species_transport surface_fem_preflight PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	python3 scripts/test_native_tet_species_workflow.py --solver solvers/cpu/native_tet_species_transport $(if $(strip $(FTETWILD_BIN)),--ftetwild $(FTETWILD_BIN) --gmsh,)

.PHONY: t9-native-y-rcr-smoke
t9-native-y-rcr-smoke:
	$(MAKE) -C solvers/cpu native_tet_hydraulic_graph surface_fem_preflight PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	python3 scripts/test_native_tet_hydraulic_y.py --solver solvers/cpu/native_tet_hydraulic_graph --ranks $(T7_SPECIES_MPI_RANKS) $(if $(strip $(FTETWILD_BIN)),--ftetwild $(FTETWILD_BIN),)

.PHONY: t9-native-y-rcr-strong-smoke
t9-native-y-rcr-strong-smoke:
	$(MAKE) -C solvers/cpu native_tet_hydraulic_graph surface_fem_preflight PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	python3 scripts/test_native_tet_hydraulic_y.py --solver solvers/cpu/native_tet_hydraulic_graph --ranks $(T7_SPECIES_MPI_RANKS) --coupling fixed $(if $(strip $(FTETWILD_BIN)),--ftetwild $(FTETWILD_BIN),)

t4-petsc-assembly-test:
	$(MAKE) -C solvers/cpu native-tet-ale-petsc-assembly-test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	$(MAKE) -C solvers/cpu native-tet-ale-petsc-runtime-test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	$(MAKE) -C solvers/cpu native-tet-ale-petsc-prescribed-runtime-test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))

.PHONY: t4-petsc-mesh-runtime-test
t4-petsc-mesh-runtime-test:
	$(MAKE) -C solvers/cpu native-tet-ale-petsc-mesh-runtime-test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH)) T4_ALE_MESH="$(T4_ALE_MESH)" T4_ALE_RESULT="$(T4_ALE_RESULT)" T4_ALE_MODE="$(T4_ALE_MODE)"

.PHONY: t4-petsc-fieldsplit-test
t4-petsc-fieldsplit-test:
	$(MAKE) -C solvers/cpu native-tet-ale-petsc-fieldsplit-test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))

.PHONY: t2-bifurcation-runtime-test
t2-bifurcation-runtime-test:
	$(MAKE) -C solvers/cpu native-tet-bifurcation-runtime-test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH)) T2_BIFURCATION_MESH="$(T2_BIFURCATION_MESH)" T2_BIFURCATION_RESULT="$(T2_BIFURCATION_RESULT)" T2_BIFURCATION_RANKS="$(T2_BIFURCATION_RANKS)" T2_BIFURCATION_STEPS="$(T2_BIFURCATION_STEPS)" T2_BIFURCATION_DT="$(T2_BIFURCATION_DT)"

.PHONY: t2-native-functional-coarse t2-native-functional-level
T2_MESH_SERIES ?=
T2_FUNCTIONAL_RESULT ?=$(CURDIR)/native-t2-functional-level-1.json
T2_FUNCTIONAL_LEVEL ?=1
T2_MPI_RANKS ?=1
t2-native-functional-coarse:
	test -n "$(T2_MESH_SERIES)" || (echo "set T2_MESH_SERIES=/path/to/mesh-series.json" >&2; exit 2)
	$(MAKE) -C solvers/cpu native_tet_flow PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	python3 scripts/run_native_t2_level.py "$(T2_MESH_SERIES)" 1 "$(T2_FUNCTIONAL_RESULT)" --mode functional --mpi-ranks $(T2_MPI_RANKS)

t2-native-functional-level:
	test -n "$(T2_MESH_SERIES)" || (echo "set T2_MESH_SERIES=/path/to/mesh-series.json" >&2; exit 2)
	$(MAKE) -C solvers/cpu native_tet_flow PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	python3 scripts/run_native_t2_level.py "$(T2_MESH_SERIES)" $(T2_FUNCTIONAL_LEVEL) "$(T2_FUNCTIONAL_RESULT)" --mode functional --mpi-ranks $(T2_MPI_RANKS)

cpu-petsc:
	$(MAKE) -C solvers/cpu petsc
	$(MAKE) -C solvers/one_d petsc

HPC_BUILD_MANIFEST ?= $(CURDIR)/hpc-build-manifest.json
hpc-build-manifest:
	./scripts/check_dependencies.sh cpu --report $(HPC_BUILD_MANIFEST)

HPC_SCALING_CASES ?= $(CURDIR)/hpc-scaling-cases
HPC_SCALING_RANKS ?= 1 64 128 256
hpc-prepare-scaling:
	$(MAKE) -C solvers/coupling hpc_duct_solver_fixture
	python3 scripts/hpc_prepare_scaling_cases.py --output-dir $(HPC_SCALING_CASES) --ranks $(HPC_SCALING_RANKS)

hpc-cross-node-binaries:
	$(MAKE) -C solvers/cpu petsc iga_flow_validate immersed_moving_distributed_fsi_runtime_test distributed_immersed_extension_test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
	$(MAKE) -C solvers/coupling petsc hpc_duct_solver_fixture PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))

HPC_TEST_OUTPUT ?= $(CURDIR)/hpc-test-results
HPC_TEST_TIMEOUT ?= 300
hpc-test-unit:
	python3 scripts/hpc_test_tiers.py --tier unit --timeout $(HPC_TEST_TIMEOUT) --output-dir $(HPC_TEST_OUTPUT)/unit

hpc-test-mpi:
	python3 scripts/hpc_test_tiers.py --tier mpi --timeout $(HPC_TEST_TIMEOUT) --output-dir $(HPC_TEST_OUTPUT)/mpi

hpc-test-gpu:
	python3 scripts/hpc_test_tiers.py --tier gpu --timeout $(HPC_TEST_TIMEOUT) --output-dir $(HPC_TEST_OUTPUT)/gpu

fluid-surface-traction-test:
	$(MAKE) -C solvers/cpu fluid-surface-traction-test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))

moving-immersed-transient-flow-fsi-runtime-test:
	$(MAKE) -C solvers/cpu moving-immersed-transient-flow-fsi-runtime-test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))

compliant-channel-fsi-test:
	$(MAKE) -C solvers/cpu compliant-channel-fsi-test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))

PHASE8_PARAVIEW_DIR ?= $(CURDIR)/artifacts/visualization/compliant-channel-fsi
phase8-compliant-channel-fsi-paraview:
	$(MAKE) -C solvers/cpu phase8-compliant-channel-fsi-paraview PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH)) PHASE8_PARAVIEW_DIR=$(abspath $(PHASE8_PARAVIEW_DIR))

pretensioned-membrane-test:
	$(MAKE) -C solvers/cpu pretensioned-membrane-test

pretensioned-membrane-fsi-runtime-test:
	$(MAKE) -C solvers/cpu pretensioned-membrane-fsi-runtime-test

material-surface-patch-kinematics-test:
	$(MAKE) -C solvers/cpu material-surface-patch-kinematics-test

phase7-focused-test:
	$(MAKE) -C solvers/cpu phase7-focused-test PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=$(PETSC_ARCH)

phase7-lv-closure-test:
	$(MAKE) -C solvers/cpu phase7-lv-closure-test PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=$(PETSC_ARCH)

one-d-petsc:
	$(MAKE) -C solvers/one_d petsc

one-d-test:
	$(MAKE) -C solvers/one_d test

zero-d-petsc:
	$(MAKE) -C solvers/one_d zero-d

zero-d-test:
	$(MAKE) -C solvers/one_d test

coupling:
	$(MAKE) -C solvers/coupling petsc

coupling-test:
	$(MAKE) -C solvers/coupling test

coupling-zero-d-flow-test:
	$(MAKE) -C solvers/coupling zero-d-flow-test

.PHONY: coupling-zero-d-species-reservoir-test
coupling-zero-d-species-reservoir-test:
	$(MAKE) -C solvers/coupling zero-d-species-reservoir-test

coupling-zero-d-source-reservoir-species-test:
	$(MAKE) -C solvers/coupling zero-d-source-reservoir-species-test

.PHONY: coupling-zero-d-source-reservoir-species-test

coupling-zero-d-source-reservoir-species-runtime-test:
	$(MAKE) -C solvers/coupling zero-d-source-reservoir-species-runtime-test

.PHONY: coupling-zero-d-source-reservoir-species-runtime-test

coupling-zero-d-flow-species-checkpoint-test:
	$(MAKE) -C solvers/coupling zero-d-flow-species-checkpoint-test

.PHONY: coupling-zero-d-flow-species-checkpoint-test

.PHONY: coupling-zero-d-terminal-rcr-species-test
coupling-zero-d-terminal-rcr-species-test:
	$(MAKE) -C solvers/coupling zero-d-terminal-rcr-species-test

.PHONY: coupling-zero-d-terminal-rcr-species-runtime-test
coupling-zero-d-terminal-rcr-species-runtime-test:
	$(MAKE) -C solvers/coupling zero-d-terminal-rcr-species-runtime-test

.PHONY: coupling-zero-d-terminal-rcr-species-graph-test
coupling-zero-d-terminal-rcr-species-graph-test:
	$(MAKE) -C solvers/coupling zero-d-terminal-rcr-species-graph-test

coupling-simulation-graph-test:
	$(MAKE) -C solvers/coupling simulation-graph-test

coupling-surface-contracts-test:
	$(MAKE) -C solvers/coupling surface-contracts-test

coupling-fsi-runtime-contracts-test:
	$(MAKE) -C solvers/coupling fsi-runtime-contracts-test

coupling-petsc-test:
	$(MAKE) -C solvers/coupling petsc-test

coupling-convergence-test:
	$(MAKE) -C solvers/coupling convergence-test

coupling-convergence-axial-diagnostic:
	$(MAKE) -C solvers/coupling convergence-axial-diagnostic

coupling-convergence-isotropic-diagnostic:
	$(MAKE) -C solvers/coupling convergence-isotropic-diagnostic

coupling-convergence-length-diagnostic:
	$(MAKE) -C solvers/coupling convergence-length-diagnostic

coupling-convergence-bulk-diagnostic:
	$(MAKE) -C solvers/coupling convergence-bulk-diagnostic

coupling-temporal-convergence-test:
	$(MAKE) -C solvers/coupling temporal-convergence-test

coupling-temporal-fixture-test:
	$(MAKE) -C solvers/coupling temporal-fixture-test

cuda:
	$(MAKE) -C solvers/cuda

spline:
	$(MAKE) -C preprocessing/spline

workflow-test: cpu
	./scripts/tests/test_run_cases.sh
	python3 ./scripts/tests/test_vtk_centerline_to_obj.py

clean:
	$(MAKE) -C preprocessing/mesh clean
	$(MAKE) -C preprocessing/spline clean
	$(MAKE) -C solvers/cpu clean
	$(MAKE) -C solvers/one_d clean
	$(MAKE) -C solvers/coupling clean
	$(MAKE) -C solvers/cuda clean
