NVCC ?= nvcc
CUDA_ARCHS ?= 70 80 89 90
CPPFLAGS += -Iinclude -I../cpu/include -I../../include
CUDA_CPPFLAGS := $(filter-out -O%,$(CPPFLAGS))
HDF5_CONDA_PREFIX ?= $(firstword $(foreach prefix,$(CONDA_PREFIX) $(CONDA_ROOT),\
	$(if $(wildcard $(prefix)/include/hdf5.h),$(prefix))))
HDF5_CFLAGS ?= $(if $(HDF5_CONDA_PREFIX),-I$(HDF5_CONDA_PREFIX)/include,\
	$(shell pkg-config --cflags hdf5 2>/dev/null))
HDF5_SERIAL_LIBRARY ?= $(firstword $(wildcard \
	/usr/lib/*/libhdf5_serial.so /usr/lib/*/libhdf5_serial.so.[0-9]* \
	/lib/*/libhdf5_serial.so /lib/*/libhdf5_serial.so.[0-9]*))
HDF5_LIBS ?= $(if $(HDF5_CONDA_PREFIX),\
	-L$(HDF5_CONDA_PREFIX)/lib -Xlinker=-rpath -Xlinker=$(HDF5_CONDA_PREFIX)/lib -lhdf5,\
	$(if $(HDF5_SERIAL_LIBRARY),\
		-L$(dir $(HDF5_SERIAL_LIBRARY)) -l:$(notdir $(HDF5_SERIAL_LIBRARY)),\
		$(shell pkg-config --libs hdf5 2>/dev/null)))
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
HEADERS += ../../include/VtkOutput.hpp ../../include/PhysiologyOutput.hpp
HEADERS += ../../include/BezierVisualization.hpp ../../include/TemporalVtkHdf.hpp

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SOURCES) Makefile $(HEADERS)
	$(NVCC) $(CUDA_CPPFLAGS) $(HDF5_CFLAGS) $(NVCCFLAGS) $(GENCODE) -Xcompiler=$(HOSTWARN) $< -o $@ $(LDLIBS) $(HDF5_LIBS)

clean:
	$(RM) $(TARGET)
