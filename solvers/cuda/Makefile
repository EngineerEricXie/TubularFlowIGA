NVCC ?= nvcc
CUDA_ARCHS ?= 70 80 89 90
CPPFLAGS += -Iinclude -I../cpu/include -I../../include
CUDA_CPPFLAGS := $(filter-out -O%,$(CPPFLAGS))
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
HEADERS += ../../include/TemporalFunction.hpp ../../include/FlowCheckpoint.hpp ../../include/VelocitySeries.hpp
HEADERS += ../../include/OutletModel.hpp ../cpu/include/BoundaryFlow.hpp
HEADERS += ../../include/OutletCheckpoint.hpp
HEADERS += ../cpu/include/OutletFlow.hpp
HEADERS += ../../include/PressureTraction.hpp ../cpu/include/PressureTractionFlow.hpp
HEADERS += ../../include/TransportCheckpoint.hpp

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SOURCES) Makefile $(HEADERS)
	$(NVCC) $(CUDA_CPPFLAGS) $(NVCCFLAGS) $(GENCODE) -Xcompiler=$(HOSTWARN) $< -o $@ $(LDLIBS)

clean:
	$(RM) $(TARGET)
