#include "IgaDatabase.hpp"
#include "OwnedRowAssembler.hpp"

#include <petscksp.h>

#include <iostream>
#include <vector>

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, "Owned-row IGA assembly smoke test\n");
	int rank = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	int status = 0;
	try {
		if (argc != 3) throw std::runtime_error("usage: iga_assembly_smoke DATABASE.ntiga FIELDS");
		iga::Database database(argv[1]);
		const auto fields = static_cast<PetscInt>(std::stol(argv[2]));
		iga::OwnedRowAssembler assembler(database, PETSC_COMM_WORLD, fields);
		Mat matrix = assembler.CreateMatrix();
		for (const auto& element : assembler.elements()) {
			const auto n = element.connectivity.size() * static_cast<std::size_t>(fields);
			std::vector<PetscScalar> values(n * n, 1.0);
			assembler.AddElementMatrix(matrix, element, values);
		}
		iga::OwnedRowAssembler::Assemble(matrix);
		PetscBool missing = PETSC_FALSE;
	#if PETSC_VERSION_LT(3, 25, 0)
		PetscInt row = -1;
		MatMissingDiagonal(matrix, &missing, &row);
	#else
		IS zero_diagonals = nullptr;
		MatFindZeroDiagonals(matrix, &zero_diagonals);
		PetscInt zero_diagonal_count = 0;
		ISGetSize(zero_diagonals, &zero_diagonal_count);
		missing = zero_diagonal_count > 0 ? PETSC_TRUE : PETSC_FALSE;
		ISDestroy(&zero_diagonals);
	#endif
		MatInfo info{};
		MatGetInfo(matrix, MAT_GLOBAL_SUM, &info);
		if (rank == 0) {
			std::cout << "global_rows=" << assembler.global_rows()
				<< " nz_used=" << static_cast<long long>(info.nz_used)
				<< " nz_allocated=" << static_cast<long long>(info.nz_allocated)
				<< " mallocs=" << static_cast<long long>(info.mallocs)
				<< " missing_diagonal=" << static_cast<int>(missing) << '\n';
		}
		if (missing || info.mallocs != 0.0) status = 1;
		MatDestroy(&matrix);
	} catch (const std::exception& e) {
		std::cerr << "rank " << rank << ": " << e.what() << '\n';
		status = 1;
	}
	int global_status = 0;
	MPI_Allreduce(&status, &global_status, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
	PetscFinalize();
	return global_status;
}
