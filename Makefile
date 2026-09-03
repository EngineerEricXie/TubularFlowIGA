.PHONY: all mesh mesh-test cpu cpu-test cpu-petsc one-d-petsc one-d-test cuda spline workflow-test clean

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

one-d-petsc:
	$(MAKE) -C solvers/one_d petsc

one-d-test:
	$(MAKE) -C solvers/one_d test

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
	$(MAKE) -C solvers/cuda clean
