# Dependency-free centerline smoothing and hexahedral control-mesh generator.
set(TFI_DIRECTORY_INCLUDES include ../../include)
set(TFI_DIRECTORY_OPTIONS -Wall -Wextra -Wpedantic)

# The mesher core is compiled once and linked into both the tool and its test.
add_library(tfi_mesh_core OBJECT
	src/SwcGraph.cpp src/BSpline.cpp src/GeometryDiagnostics.cpp
	src/HexMesh.cpp src/CrossSectionTemplate.cpp src/MeshGenerator.cpp)
target_include_directories(tfi_mesh_core PRIVATE ${TFI_DIRECTORY_INCLUDES})
target_compile_options(tfi_mesh_core PRIVATE ${TFI_DIRECTORY_OPTIONS})

tfi_executable(tubular_mesh SOURCES src/main.cpp $<TARGET_OBJECTS:tfi_mesh_core>)
tfi_executable(mesh_core_test SOURCES tests/test_core.cpp $<TARGET_OBJECTS:tfi_mesh_core>)

tfi_make_alias(all TARGETS tubular_mesh)
tfi_test(mesh_core_test COMMAND [==[./mesh_core_test ../../meshgeneration/template]==]
	REQUIRES mesh_core_test MAKE test LABELS unit)
