CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS ?= -Iinclude
LDLIBS ?= -lstdc++fs

PETSC_DIR ?=
PETSC_ARCH ?=
PETSC_TARGETS := petsc iga_assembly_smoke iga_transport iga_navier_stokes iga_mesh_check

ifneq ($(strip $(filter $(PETSC_TARGETS),$(MAKECMDGOALS))),)
ifeq ($(strip $(PETSC_DIR)),)
$(error PETSC_DIR must be set when building PETSc targets)
endif
PETSC_PREFIX := $(PETSC_DIR)$(if $(strip $(PETSC_ARCH)),/$(PETSC_ARCH))
PETSC_VARIABLES ?= $(PETSC_PREFIX)/lib/petsc/conf/petscvariables
include $(PETSC_VARIABLES)
endif

.PHONY: all petsc clean

all: iga_pack iga_inspect

petsc: iga_assembly_smoke iga_transport iga_navier_stokes iga_mesh_check

iga_pack: src/iga_pack.cpp include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

iga_inspect: src/iga_inspect.cpp include/IgaDatabase.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

iga_assembly_smoke: src/iga_assembly_smoke.cpp include/IgaDatabase.hpp include/OwnedRowAssembler.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

iga_transport: src/iga_transport.cpp include/IgaDatabase.hpp include/OwnedRowAssembler.hpp include/CaseInput.hpp include/TransportElement.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

iga_navier_stokes: src/iga_navier_stokes.cpp include/IgaDatabase.hpp include/OwnedRowAssembler.hpp include/CaseInput.hpp include/TransportElement.hpp include/NavierStokesElement.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

iga_mesh_check: src/iga_mesh_check.cpp include/IgaDatabase.hpp include/TransportElement.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PETSC_CC_INCLUDES) $< -o $@ $(PETSC_LIB)

clean:
	$(RM) iga_pack iga_inspect iga_assembly_smoke iga_transport iga_navier_stokes iga_mesh_check
