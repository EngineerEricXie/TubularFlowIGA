.PHONY: all mesh mesh-test cpu cpu-test cpu-petsc phase7-focused-test phase7-lv-closure-test one-d-petsc one-d-test coupling coupling-test coupling-surface-contracts-test coupling-fsi-runtime-contracts-test coupling-petsc-test coupling-convergence-test coupling-convergence-axial-diagnostic coupling-convergence-isotropic-diagnostic coupling-convergence-length-diagnostic coupling-convergence-bulk-diagnostic coupling-temporal-convergence-test coupling-temporal-fixture-test cuda spline clean

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
