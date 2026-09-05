CXX ?= g++
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
PETSC_TARGETS := petsc petsc-test phase5-closure-test phase5-closure-manufactured phase5-closure-manufactured-triplet phase5-closure-sliver phase5-closure-all phase5_closure_test phase6-aneurysm-depth2-regression phase6-aneurysm-runtime-closure-test phase6_aneurysm_runtime_closure_test phase7-focused-test phase7-lv-closure-test phase7_lv_closure_test quadrature_test generic_transport_test navier_stokes_test nitsche_wall_test ghost_penalty_test immersed_transient_state_test immersed-transient-state-test immersed_velocity_extension_test immersed-velocity-extension-test immersed_static_flow_test immersed_transient_flow_test immersed-transient-flow-test moving_immersed_transient_flow_test moving-immersed-transient-flow-test moving_immersed_flow_snapshot_test moving-immersed-flow-snapshot-test moving_immersed_flow_snapshot_publisher_test moving-immersed-flow-snapshot-publisher-test immersed_flow_port_test immersed_flow_domain_test immersed_aneurysm_jacobian_test transport_boundary_preflight_test vca_3d_runtime_test vca_3d_smoke_test iga_solve iga_assembly_smoke iga_transport iga_navier_stokes iga_mesh_check moving_cut_geometry_test moving-cut-geometry-test fluid_surface_traction_test fluid-surface-traction-test moving_immersed_transient_flow_fsi_runtime_test moving-immersed-transient-flow-fsi-runtime-test compliant_channel_fsi_test compliant-channel-fsi-test

ifneq ($(strip $(filter $(PETSC_TARGETS),$(MAKECMDGOALS))),)
ifeq ($(strip $(PETSC_DIR)),)
$(error PETSC_DIR must be set when building PETSc targets)
endif
PETSC_PREFIX := $(PETSC_DIR)$(if $(strip $(PETSC_ARCH)),/$(PETSC_ARCH))
PETSC_VARIABLES ?= $(PETSC_PREFIX)/lib/petsc/conf/petscvariables
include $(PETSC_VARIABLES)
endif

.PHONY: all petsc petsc-test test config-check-dispatch-test vtkhdf-paraview-test clean phase5-closure-test phase5-closure-manufactured phase5-closure-manufactured-triplet phase5-closure-sliver phase5-closure-all phase6-aneurysm-depth2-regression phase6-aneurysm-runtime-closure-test phase7-focused-test phase7-lv-closure-test prescribed-surface-motion-test material-surface-patch-kinematics-test moving-cut-geometry-test immersed-transient-state-test immersed-transient-flow-test moving-immersed-transient-flow-test moving-immersed-transient-flow-fsi-runtime-test moving-immersed-flow-snapshot-test moving-immersed-flow-snapshot-publisher-test immersed-velocity-extension-test pretensioned-membrane-test pretensioned-membrane-fsi-runtime-test

# Shared Cartesian spline-face ordering is part of all immersed cut/ghost and
# wall paths; make header changes rebuild every affected executable.
immersed_face_targets := navier_stokes_test nitsche_wall_test ghost_penalty_test immersed_transient_state_test immersed_static_flow_test immersed_transient_flow_test moving_immersed_flow_snapshot_test immersed_flow_port_test immersed_flow_domain_test immersed_aneurysm_jacobian_test phase6_aneurysm_runtime_closure_test phase7_lv_closure_test phase5_closure_test moving_cut_geometry_test
$(immersed_face_targets): include/CubicCartesianSplineFace.hpp

all: iga_pack iga_inspect iga_case_check iga_config_check iga_flow_validate iga_womersley_reference

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

immersed_transient_flow_test: tests/test_immersed_transient_flow.cpp include/ImmersedTransientFlowRuntime.hpp include/ImmersedVelocityExtension.hpp include/MovingCutGeometry.hpp include/ImmersedTransientState.hpp include/ImmersedFlowPort.hpp include/CutCellGhostPenalty.hpp include/ImmersedNitscheWall.hpp include/NavierStokesElement.hpp
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

compliant_channel_fsi_test: CPPFLAGS += -DIGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING
compliant_channel_fsi_test: tests/test_compliant_channel_fsi.cpp tests/CompliantChannelFsiFixture.hpp include/StrongFluidStructureCoupling.hpp include/MovingImmersedTransientFlowFsiRuntime.hpp include/PretensionedMembraneFsiRuntime.hpp include/PretensionedMembrane.hpp include/MovingImmersedTransientFlowRuntime.hpp include/FluidSurfaceTraction.hpp include/MaterialSurfacePatchKinematics.hpp ../../include/FsiDomainRuntime.hpp ../../include/FsiCouplingEdge.hpp ../../include/DistributedSurfaceInterface.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -UNDEBUG $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

compliant-channel-fsi-test: compliant_channel_fsi_test
	./compliant_channel_fsi_test

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

bezier_visualization_test: tests/test_bezier_visualization.cpp include/IgaDatabase.hpp ../../include/BezierVisualization.hpp ../../include/VtkOutput.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

temporal_vtkhdf_test: tests/test_temporal_vtkhdf.cpp include/IgaDatabase.hpp ../../include/BezierVisualization.hpp ../../include/TemporalVtkHdf.hpp ../../include/VtkOutput.hpp
	$(CXX) $(CPPFLAGS) $(HDF5_ALL_CFLAGS) $(CXXFLAGS) $< -o $@ $(HDF5_LIBS) $(LDLIBS)

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

iga_solve: src/iga_solve.cpp include/IgaDatabase.hpp include/MemoryReport.hpp include/OwnedRowAssembler.hpp include/TransportBoundaryPreflight.hpp include/FieldCouplingPattern.hpp include/CaseInput.hpp include/GenericCaseInput.hpp include/GenericTransportElement.hpp ../../include/SimulationConfig.hpp ../../include/TemporalFunction.hpp ../../include/VelocitySeries.hpp ../../include/FlowCheckpoint.hpp ../../include/TransportCheckpoint.hpp ../../include/VtkOutput.hpp ../../include/PhysiologyOutput.hpp ../../include/BezierVisualization.hpp ../../include/TemporalVtkHdf.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $(HDF5_CFLAGS) $< -o $@ $(PETSC_LIB) $(HDF5_LIBS)

iga_transport: src/iga_transport.cpp include/IgaDatabase.hpp include/OwnedRowAssembler.hpp include/FieldCouplingPattern.hpp include/CaseInput.hpp include/TransportElement.hpp ../../include/CaseConfig.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

iga_navier_stokes: src/iga_navier_stokes.cpp include/BoundarySupport.hpp include/IgaDatabase.hpp include/OwnedRowAssembler.hpp include/CaseInput.hpp include/GenericCaseInput.hpp include/GenericTransportElement.hpp include/TransientFlowRuntime.hpp include/TransientTransportRuntime.hpp include/TransportElement.hpp include/NavierStokesElement.hpp include/BoundaryFlow.hpp include/PressureTractionFlow.hpp ../../include/CaseConfig.hpp ../../include/CouplingPort.hpp $(SPECIES_COUPLING_HEADER) ../../include/SimulationConfig.hpp ../../include/ThreeDFlowCoupling.hpp ../../include/ThreeDVcaCoupling.hpp ../../include/VascularCoupling.hpp ../../include/CouplingConfig.hpp ../../include/TemporalFunction.hpp ../../include/FlowCheckpoint.hpp ../../include/VcaCheckpoint.hpp ../../include/CouplingHistory.hpp ../../include/OutletModel.hpp ../../include/OutletCheckpoint.hpp ../../include/PressureTraction.hpp ../../include/VelocitySeries.hpp ../../include/VtkOutput.hpp ../../include/BezierVisualization.hpp ../../include/TemporalVtkHdf.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $(HDF5_CFLAGS) $< -o $@ $(PETSC_LIB) $(HDF5_LIBS)

iga_mesh_check: src/iga_mesh_check.cpp include/IgaDatabase.hpp include/TransportElement.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

clean:
	$(RM) iga_pack iga_inspect iga_case_check iga_config_check iga_flow_validate iga_womersley_reference case_config_test simulation_config_test boundary_support_test coupling_port_test sha256_test surface_geometry_test prescribed_surface_motion_test material_surface_patch_kinematics_test moving_cut_geometry_test fluid_surface_traction_test moving_immersed_transient_flow_fsi_runtime_test compliant_channel_fsi_test pretensioned_membrane_test pretensioned_membrane_fsi_runtime_test strong_fluid_structure_coupling_test surface_spatial_index_test cartesian_domain_classification_test cut_cell_volume_quadrature_test immersed_surface_quadrature_test surface_stl_test surface_vtp_test cartesian_background_test quadrature_test vtk_output_test bezier_visualization_test temporal_vtkhdf_test flow_checkpoint_test transport_checkpoint_test vca_checkpoint_test outlet_model_test boundary_flow_test womersley_reference_test generic_transport_test navier_stokes_test nitsche_wall_test ghost_penalty_test immersed_transient_state_test immersed_velocity_extension_test immersed_static_flow_test immersed_transient_flow_test moving_immersed_transient_flow_test moving_immersed_flow_snapshot_test moving_immersed_flow_snapshot_publisher_test immersed_flow_port_test immersed_flow_domain_test immersed_aneurysm_jacobian_test phase5_closure_test phase6_aneurysm_runtime_closure_test phase7_lv_closure_test transport_boundary_preflight_test vca_3d_runtime_test vca_3d_smoke_test iga_solve iga_assembly_smoke iga_transport iga_navier_stokes iga_mesh_check
