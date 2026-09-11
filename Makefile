.PHONY: all mesh mesh-test cpu cpu-test cpu-petsc hpc-build-manifest hpc-prepare-scaling hpc-cross-node-binaries hpc-test-unit hpc-test-mpi hpc-test-gpu fluid-surface-traction-test moving-immersed-transient-flow-fsi-runtime-test compliant-channel-fsi-test phase8-compliant-channel-fsi-paraview pretensioned-membrane-test pretensioned-membrane-fsi-runtime-test material-surface-patch-kinematics-test phase7-focused-test phase7-lv-closure-test one-d-petsc one-d-test coupling coupling-test coupling-zero-d-flow-test coupling-simulation-graph-test coupling-surface-contracts-test coupling-fsi-runtime-contracts-test coupling-petsc-test coupling-convergence-test coupling-convergence-axial-diagnostic coupling-convergence-isotropic-diagnostic coupling-convergence-length-diagnostic coupling-convergence-bulk-diagnostic coupling-temporal-convergence-test coupling-temporal-fixture-test cuda spline clean

all: cpu

mesh:
	$(MAKE) -C preprocessing/mesh

mesh-test:
	$(MAKE) -C preprocessing/mesh test

cpu:
	$(MAKE) -C solvers/cpu

cpu-test:
	$(MAKE) -C solvers/cpu test

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
	$(MAKE) -C solvers/cpu petsc immersed_moving_distributed_fsi_runtime_test distributed_immersed_extension_test PETSC_DIR=$(PETSC_DIR) $(if $(strip $(PETSC_ARCH)),PETSC_ARCH=$(PETSC_ARCH))
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

PHASE8_PARAVIEW_DIR ?= $(CURDIR)/results/phase8_compliant_channel_fsi_paraview
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

coupling:
	$(MAKE) -C solvers/coupling petsc

coupling-test:
	$(MAKE) -C solvers/coupling test

coupling-zero-d-flow-test:
	$(MAKE) -C solvers/coupling zero-d-flow-test

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

clean:
	$(MAKE) -C preprocessing/mesh clean
	$(MAKE) -C preprocessing/spline clean
	$(MAKE) -C solvers/cpu clean
	$(MAKE) -C solvers/one_d clean
	$(MAKE) -C solvers/coupling clean
	$(MAKE) -C solvers/cuda clean
