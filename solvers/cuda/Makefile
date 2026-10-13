# Targets and tests for this directory are defined in CMakeLists.txt. This file
# keeps every historical `make -C solvers/cuda <target>` command working: executables
# are built with CMake and test targets run through CTest (see docs/BUILD.md).
.DEFAULT_GOAL := all
.PHONY: FORCE
Makefile: ;
FORCE: ;
%: FORCE
	+@python3 ../../cmake/forward.py solvers/cuda $@
