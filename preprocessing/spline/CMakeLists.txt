# Spline construction and Bezier extraction (C++11, Eigen, OpenMP).
if(NOT TFI_EIGEN_INCLUDE_DIR)
	tfi_skip_target(spline "Eigen 3 (pass EIGEN_DIR)")
	tfi_make_alias(all TARGETS spline)
	return()
endif()
set(TFI_DIRECTORY_INCLUDES "${TFI_EIGEN_INCLUDE_DIR}" ../../include)
set(TFI_DIRECTORY_OPTIONS -Wall -Wextra)

tfi_executable(spline SOURCES BasicDataStructure.cpp kernel.cpp KnotInsertion.cpp main.cpp
	OPENMP LIBRARIES pthread m dl)
set_target_properties(spline PROPERTIES CXX_STANDARD 11)

tfi_make_alias(all TARGETS spline)
