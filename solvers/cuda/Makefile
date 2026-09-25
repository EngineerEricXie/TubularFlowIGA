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
	../cpu/include/IgaDatabase.hpp ../../include/CaseConfig.hpp ../../include/SimulationConfig.hpp ../../include/MeshConfig.hpp
HEADERS += ../../include/TemporalFunction.hpp ../../include/FlowCheckpoint.hpp ../../include/VelocitySeries.hpp
HEADERS += ../../include/OutletModel.hpp ../cpu/include/BoundaryFlow.hpp
HEADERS += ../../include/OutletCheckpoint.hpp
HEADERS += ../cpu/include/OutletFlow.hpp
HEADERS += ../../include/PressureTraction.hpp ../cpu/include/PressureTractionFlow.hpp
HEADERS += ../../include/TransportCheckpoint.hpp
HEADERS += ../../include/VtkOutput.hpp ../../include/PhysiologyOutput.hpp
HEADERS += ../../include/BezierVisualization.hpp ../../include/TemporalVtkHdf.hpp
HEADERS += include/DeviceAllocation.hpp
HEADERS += ../cpu/include/BoundarySupport.hpp
HEADERS += ../../include/PhaseProfile.hpp
HEADERS += ../../include/CheckedText.hpp ../../include/ExecutionEnvironment.hpp include/CudaExecution.hpp

.PHONY: all clean allocation-test execution-test

all: $(TARGET) native_tet_flow_cuda native_tet_darcy_cuda native_tet_species_cuda native_tet_fsi_cuda

native_tet_flow_cuda: src/native_tet_flow_cuda.cu include/NativeTetFlowKernels.cuh include/CudaRuntime.hpp ../cpu/include/NativeTetFem.hpp ../cpu/include/NativeTetAleTransient.hpp
	$(NVCC) $(CUDA_CPPFLAGS) $(NVCCFLAGS) $(GENCODE) -Xcompiler=$(HOSTWARN) $< -o $@ -lcusolver -lcublas -lstdc++fs

native_tet_darcy_cuda: src/native_tet_darcy_cuda.cu include/NativeTetDarcyKernels.cuh include/CudaRuntime.hpp ../cpu/include/NativeTetDarcyPetsc.hpp ../cpu/include/NativeTetDarcyVisualization.hpp ../cpu/include/NativeTetVesselTissueSourceMap.hpp
	$(NVCC) $(CUDA_CPPFLAGS) $(filter -I%,$(shell pkg-config --cflags PETSc)) $(shell mpicxx --showme:compile 2>/dev/null) $(NVCCFLAGS) $(GENCODE) -Xcompiler=$(HOSTWARN) $< -o $@ -lcusolver -lcublas -lstdc++fs

native_tet_species_cuda: src/native_tet_species_cuda.cu include/NativeTetSpeciesKernels.cuh include/CudaRuntime.hpp ../cpu/include/NativeTetSpeciesVisualization.hpp
	$(NVCC) $(CUDA_CPPFLAGS) $(NVCCFLAGS) $(GENCODE) -Xcompiler=$(HOSTWARN) $< -o $@ -lcusolver -lcublas -lstdc++fs

native_tet_fsi_cuda: src/native_tet_fsi_cuda.cu src/native_tet_flow_cuda.cu include/NativeTetFlowKernels.cuh include/CudaRuntime.hpp ../cpu/include/NativeTetSolidStaticSolver.hpp ../cpu/include/NativeTetAleMeshMotion.hpp
	$(NVCC) $(CUDA_CPPFLAGS) $(filter -I%,$(shell pkg-config --cflags PETSc)) $(shell mpicxx --showme:compile 2>/dev/null) $(NVCCFLAGS) $(GENCODE) -Xcompiler=$(HOSTWARN) $< -o $@ -lcusolver -lcublas -lstdc++fs

cuda_execution_test: tests/test_execution.cpp include/CudaExecution.hpp ../../include/ExecutionEnvironment.hpp ../../include/CheckedText.hpp
	$(CXX) -Iinclude -I../../include -std=c++17 -Wall -Wextra -Wpedantic $< -o $@

execution-test: cuda_execution_test
	./cuda_execution_test

$(TARGET): $(SOURCES) Makefile $(HEADERS)
	$(NVCC) $(CUDA_CPPFLAGS) $(HDF5_CFLAGS) $(NVCCFLAGS) $(GENCODE) -Xcompiler=$(HOSTWARN) $< -o $@ $(LDLIBS) $(HDF5_LIBS)

device_allocation_test: tests/test_device_allocation.cu include/CudaRuntime.hpp include/DeviceAllocation.hpp
	$(NVCC) -Iinclude $(NVCCFLAGS) $(GENCODE) -Xcompiler=$(HOSTWARN) $< -o $@ -lcublas

cuda_stdout_failure_test: tests/test_stdout_failure.cu ../cpu/tests/test_tool_stdout_failure.cpp $(SOURCES) Makefile $(HEADERS)
	$(NVCC) $(CUDA_CPPFLAGS) $(HDF5_CFLAGS) $(NVCCFLAGS) $(GENCODE) -Xcompiler=$(HOSTWARN) $< -o $@ $(LDLIBS) $(HDF5_LIBS)

allocation-test: device_allocation_test
	./device_allocation_test

clean:
	$(RM) $(TARGET) native_tet_flow_cuda native_tet_darcy_cuda native_tet_species_cuda native_tet_fsi_cuda device_allocation_test cuda_stdout_failure_test cuda_execution_test
