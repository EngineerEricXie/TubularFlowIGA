NVCC ?= nvcc
CUDA_ARCHS ?= 70 80 89 90
CPPFLAGS ?= -Iinclude -I../cpu/include -I../../include
NVCCFLAGS ?= -O3 -std=c++17 -lineinfo
HOSTWARN ?= -Wall,-Wextra
GENCODE := $(foreach arch,$(CUDA_ARCHS),-gencode arch=compute_$(arch),code=sm_$(arch))
LDLIBS ?= -lcublas -lstdc++fs

TARGET := iga_cuda
SOURCES := src/iga_cuda.cu
HEADERS := include/BlockCsr.hpp include/CudaRuntime.hpp include/DeviceMesh.hpp \
	include/Gmres.hpp include/IgaCudaKernels.cuh include/GenericTransportKernels.cuh \
	include/ConfiguredTransportHost.inc ../cpu/include/CaseInput.hpp ../cpu/include/GenericCaseInput.hpp \
	../cpu/include/IgaDatabase.hpp ../../include/CaseConfig.hpp ../../include/SimulationConfig.hpp
HEADERS += ../../include/TemporalFunction.hpp

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SOURCES) Makefile $(HEADERS)
	$(NVCC) $(CPPFLAGS) $(NVCCFLAGS) $(GENCODE) -Xcompiler=$(HOSTWARN) $< -o $@ $(LDLIBS)

clean:
	$(RM) $(TARGET)
