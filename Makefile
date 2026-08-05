.PHONY: all cpu cpu-petsc cuda spline clean

all: cpu

cpu:
	$(MAKE) -C solvers/cpu

cpu-petsc:
	$(MAKE) -C solvers/cpu petsc

cuda:
	$(MAKE) -C solvers/cuda

spline:
	$(MAKE) -C preprocessing/spline

clean:
	$(MAKE) -C preprocessing/spline clean
	$(MAKE) -C solvers/cpu clean
	$(MAKE) -C solvers/cuda clean
