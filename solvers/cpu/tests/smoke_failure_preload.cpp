#include <petscmat.h>
#include <dlfcn.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace {
Mat observed_matrix = nullptr;
Mat insertion_matrix = nullptr;
bool injected = false;

bool Enabled(const char* name)
{
	const char* rank = std::getenv("OMPI_COMM_WORLD_RANK");
	const char* mode = std::getenv("TUBULARFLOWIGA_TEST_SMOKE_FAULT");
	return !injected && rank && std::strcmp(rank, "0") == 0
		&& mode && std::strcmp(mode, name) == 0;
}

void Record(const char* name)
{
	injected = true;
	std::fprintf(stderr, "[smoke-failure-injection] %s\n", name);
}
}

// Test-only Linux interposition. MPI collectives complete their real calls
// before a returned error is injected on rank zero. Local insertion also tests
// PETSc's actual mixed-insertion-mode error under the returning handler.
extern "C" PetscErrorCode MatSetValues(Mat matrix, PetscInt rows, const PetscInt* row,
	PetscInt columns, const PetscInt* column, const PetscScalar* values, InsertMode mode)
{
	using Function = PetscErrorCode (*)(Mat, PetscInt, const PetscInt*, PetscInt,
		const PetscInt*, const PetscScalar*, InsertMode);
	static const auto function = reinterpret_cast<Function>(dlsym(RTLD_NEXT, "MatSetValues"));
	if (!function) return PETSC_ERR_LIB;
	const auto status = function(matrix, rows, row, columns, column, values, mode);
	if (!status && rows && columns) insertion_matrix = matrix;
	if (!status && mode == ADD_VALUES && rows && columns && Enabled("mixed-insertion")) {
		Record("mixed-insertion");
		return function(matrix, rows, row, columns, column, values, INSERT_VALUES);
	}
	return status;
}

extern "C" PetscErrorCode MatAssemblyEnd(Mat matrix, MatAssemblyType type)
{
	using Function = PetscErrorCode (*)(Mat, MatAssemblyType);
	static const auto function = reinterpret_cast<Function>(dlsym(RTLD_NEXT, "MatAssemblyEnd"));
	if (!function) return PETSC_ERR_LIB;
	const auto status = function(matrix, type);
	// Select the matrix populated by the caller, excluding PETSc submatrices.
	if (!status && matrix == insertion_matrix && Enabled("assembly")) {
		Record("assembly");
		return PETSC_ERR_USER;
	}
	return status;
}

extern "C" PetscErrorCode MatGetInfo(Mat matrix, MatInfoType type, MatInfo* info)
{
	using Function = PetscErrorCode (*)(Mat, MatInfoType, MatInfo*);
	static const auto function = reinterpret_cast<Function>(dlsym(RTLD_NEXT, "MatGetInfo"));
	if (!function) return PETSC_ERR_LIB;
	const auto status = function(matrix, type, info);
	if (!status && type == MAT_GLOBAL_SUM) {
		observed_matrix = matrix;
		if (Enabled("info")) {
			Record("info");
			return PETSC_ERR_USER;
		}
	}
	return status;
}

extern "C" PetscErrorCode MatDestroy(Mat* matrix)
{
	using Function = PetscErrorCode (*)(Mat*);
	static const auto function = reinterpret_cast<Function>(dlsym(RTLD_NEXT, "MatDestroy"));
	if (!function) return PETSC_ERR_LIB;
	const bool target = matrix && *matrix && *matrix == observed_matrix;
	const auto status = function(matrix);
	if (!status && target && Enabled("destroy")) {
		Record("destroy");
		return PETSC_ERR_USER;
	}
	return status;
}
