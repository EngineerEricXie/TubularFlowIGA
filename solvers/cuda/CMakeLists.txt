# FP64 single-GPU backend. Configure with TFI_ENABLE_CUDA=ON (`make cuda` does
# this) and CMAKE_CUDA_ARCHITECTURES, for example 89 for an RTX 4080 SUPER.
set(TFI_DIRECTORY_INCLUDES include ../cpu/include ../../include)
set(TFI_DIRECTORY_OPTIONS -Wall -Wextra -Wpedantic)

# The host-side execution contract needs no GPU or CUDA toolkit.
tfi_executable(cuda_execution_test SOURCES tests/test_execution.cpp)
tfi_test(cuda_execution_test COMMAND [==[./cuda_execution_test]==]
	REQUIRES cuda_execution_test MAKE execution-test LABELS unit)

set(cuda_options $<$<COMPILE_LANGUAGE:CUDA>:-lineinfo -Xcompiler=-Wall,-Wextra>)
if(TFI_ENABLE_CUDA)
	set(cuda_libraries CUDA::cublas CUDA::cusolver)
endif()
tfi_executable(iga_cuda SOURCES src/iga_cuda.cu CUDA HDF5 OPTIONS ${cuda_options}
	LIBRARIES ${cuda_libraries})
tfi_executable(native_tet_flow_cuda SOURCES src/native_tet_flow_cuda.cu CUDA
	OPTIONS ${cuda_options} LIBRARIES ${cuda_libraries})
tfi_executable(native_tet_species_cuda SOURCES src/native_tet_species_cuda.cu CUDA
	OPTIONS ${cuda_options} LIBRARIES ${cuda_libraries})
tfi_executable(cuda_stdout_failure_test SOURCES tests/test_stdout_failure.cu CUDA HDF5
	OPTIONS ${cuda_options} LIBRARIES ${cuda_libraries})
tfi_executable(device_allocation_test SOURCES tests/test_device_allocation.cu CUDA
	OPTIONS ${cuda_options} LIBRARIES ${cuda_libraries})

# Darcy and FSI share CPU headers that declare PETSc and MPI types; they use
# those headers only and do not link PETSc.
foreach(name native_tet_darcy_cuda native_tet_fsi_cuda)
	string(REPLACE "_cuda" "_cuda.cu" source "src/${name}")
	if(TFI_ENABLE_CUDA AND NOT TARGET TubularFlow::PETSc)
		tfi_skip_target(${name} "PETSc headers (pass PETSC_DIR and PETSC_ARCH)")
		continue()
	endif()
	tfi_executable(${name} SOURCES ${source} CUDA OPTIONS ${cuda_options}
		INCLUDES ${PETSC_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_DIRS} LIBRARIES ${cuda_libraries})
endforeach()

tfi_make_alias(all TARGETS iga_cuda native_tet_flow_cuda native_tet_darcy_cuda
	native_tet_species_cuda native_tet_fsi_cuda)
tfi_test(device_allocation_test COMMAND [==[./device_allocation_test]==]
	REQUIRES device_allocation_test MAKE allocation-test LABELS gpu)
