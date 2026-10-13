cmake_minimum_required(VERSION 3.20)
project(TubularFlowIGA LANGUAGES C CXX)

# Every target is optimized but keeps assertions: the historical Makefiles
# never defined NDEBUG, and the regression tests rely on assert().
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3" CACHE STRING "Release C++ flags (assertions stay enabled)" FORCE)
set(CMAKE_CUDA_FLAGS_RELEASE "-O3" CACHE STRING "Release CUDA flags (assertions stay enabled)" FORCE)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(TubularFlowDependencies)
include(TubularFlowHelpers)
enable_testing()

# Native graph checkpoints embed a hash of the graph sources and their build
# definition; the header is regenerated whenever any hashed file changes.
set(TFI_NATIVE_CHECKPOINT_IDENTITY_HEADER
	"${CMAKE_BINARY_DIR}/generated/NativeCheckpointSourceIdentity.hpp")
file(GLOB native_checkpoint_sources CONFIGURE_DEPENDS
	"${PROJECT_SOURCE_DIR}/include/*.hpp"
	"${PROJECT_SOURCE_DIR}/solvers/cpu/include/*.hpp"
	"${PROJECT_SOURCE_DIR}/solvers/one_d/include/*.hpp")
add_custom_command(OUTPUT "${TFI_NATIVE_CHECKPOINT_IDENTITY_HEADER}"
	COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/scripts/native_checkpoint_source_identity.py"
		--header "${TFI_NATIVE_CHECKPOINT_IDENTITY_HEADER}"
	DEPENDS ${native_checkpoint_sources}
		"${PROJECT_SOURCE_DIR}/solvers/coupling/src/iga_1d_3d_bifurcation.cpp"
		"${PROJECT_SOURCE_DIR}/solvers/coupling/CMakeLists.txt"
		"${PROJECT_SOURCE_DIR}/cmake/TubularFlowHelpers.cmake"
		"${PROJECT_SOURCE_DIR}/scripts/native_checkpoint_source_identity.py"
	COMMENT "Hashing native graph checkpoint sources")
add_custom_target(tfi_native_checkpoint_identity DEPENDS "${TFI_NATIVE_CHECKPOINT_IDENTITY_HEADER}")

add_subdirectory(preprocessing/mesh)
add_subdirectory(preprocessing/spline)
add_subdirectory(solvers/cpu)
add_subdirectory(solvers/one_d)
add_subdirectory(solvers/coupling)
add_subdirectory(solvers/cuda)

tfi_write_forward_manifest()
