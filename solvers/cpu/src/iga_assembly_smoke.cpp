#include "ExecutionResources.hpp"
#include "CollectiveAssetInput.hpp"
#include "CollectivePetscOptions.hpp"
#include <memory>
#include "IgaDatabase.hpp"
#include "OwnedRowAssembler.hpp"

#include <petscksp.h>

#include <iostream>
#include <vector>

namespace {
class ReturnErrors {
public:
	void Enable(MPI_Comm communicator)
	{
		iga::CollectiveLocalStage(communicator, "assembly smoke error handler", [&] {
			if (PetscPushErrorHandler(PetscReturnErrorHandler, nullptr))
				throw std::runtime_error("cannot install returning PETSc error handler");
			active_ = true;
		});
	}
	~ReturnErrors() { if (active_) PetscPopErrorHandler(); }
private:
	bool active_ = false;
};

struct MatrixOwner {
	MatrixOwner() = default;
	MatrixOwner(const MatrixOwner&) = delete;
	MatrixOwner& operator=(const MatrixOwner&) = delete;
	Mat value = nullptr;
	~MatrixOwner() { if (value) MatDestroy(&value); }
};
}

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, "Owned-row IGA assembly smoke test\n");
	int rank = 0, ranks = 1;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	int status = 0;
	try {
		iga::RequireExecutionResources(PETSC_COMM_WORLD, &std::cout);
		ReturnErrors errors;
		errors.Enable(PETSC_COMM_WORLD);
		std::unique_ptr<iga::Database> database_owner;
		PetscInt fields = 0;
		std::string field_description;
		std::string database_fingerprint;
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "assembly smoke input", [&] {
			if (argc != 3) throw std::runtime_error("usage: iga_assembly_smoke DATABASE.ntiga FIELDS");
			const std::string text(argv[2]); std::size_t used = 0;
			const auto parsed = std::stoll(text, &used);
			if (used != text.size() || parsed <= 0 || static_cast<unsigned long long>(parsed)
				> static_cast<unsigned long long>(std::numeric_limits<PetscInt>::max()))
				throw std::runtime_error("FIELDS must be a positive integer within PetscInt capacity");
			fields = static_cast<PetscInt>(parsed);
			field_description = std::to_string(fields);
			database_fingerprint = iga::ReadAssetFingerprint(argv[1]);
			database_owner = std::make_unique<iga::Database>(argv[1]);
			iga::ValidatePackedExecution(database_owner->header().ranks, database_owner->header().nodes, fields, ranks);
		});
		iga::RequireCollectiveSameText(PETSC_COMM_WORLD, "assembly smoke field agreement", field_description);
		iga::RequireCollectiveSameText(PETSC_COMM_WORLD, "assembly smoke asset database", database_fingerprint);
		iga::RequireCollectivePetscOptions(PETSC_COMM_WORLD);
		auto& database = *database_owner;
		iga::OwnedRowAssembler assembler(database, PETSC_COMM_WORLD, fields);
		MatrixOwner owner;
		owner.value = assembler.CreateMatrix();
		const auto matrix = owner.value;
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "assembly smoke element insertion", [&] {
			for (const auto& element : assembler.elements()) {
				const auto n = element.connectivity.size() * static_cast<std::size_t>(fields);
				std::vector<PetscScalar> values(n * n, 1.0);
				assembler.AddElementMatrix(matrix, element, values);
			}
		});
		iga::OwnedRowAssembler::Assemble(matrix, PETSC_COMM_WORLD);
		PetscBool missing = PETSC_FALSE;
	#if PETSC_VERSION_LT(3, 25, 0)
		PetscInt row = -1;
		iga::RequireCollectivePetscSuccess(PETSC_COMM_WORLD, "assembly smoke diagonal", MatMissingDiagonal(matrix, &missing, &row));
	#else
		IS zero_diagonals = nullptr;
		iga::RequireCollectivePetscSuccess(PETSC_COMM_WORLD, "assembly smoke diagonal", MatFindZeroDiagonals(matrix, &zero_diagonals));
		PetscInt zero_diagonal_count = 0;
		iga::RequireCollectivePetscSuccess(PETSC_COMM_WORLD, "assembly smoke diagonal count", ISGetSize(zero_diagonals, &zero_diagonal_count));
		missing = zero_diagonal_count > 0 ? PETSC_TRUE : PETSC_FALSE;
		iga::RequireCollectivePetscSuccess(PETSC_COMM_WORLD, "assembly smoke diagonal cleanup", ISDestroy(&zero_diagonals));
	#endif
		MatInfo info{};
		iga::RequireCollectivePetscSuccess(PETSC_COMM_WORLD, "assembly smoke matrix info", MatGetInfo(matrix, MAT_GLOBAL_SUM, &info));
		iga::RequireCollectivePetscSuccess(PETSC_COMM_WORLD, "assembly smoke destroy", MatDestroy(&owner.value));
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "assembly smoke result logging", [&] {
			if (rank != 0) return;
			std::cout << "global_rows=" << assembler.global_rows()
				<< " nz_used=" << static_cast<long long>(info.nz_used)
				<< " nz_allocated=" << static_cast<long long>(info.nz_allocated)
				<< " mallocs=" << static_cast<long long>(info.mallocs)
				<< " missing_diagonal=" << static_cast<int>(missing) << '\n';
			iga::FlushCheckedText(std::cout);
		});
		if (missing || info.mallocs != 0.0) status = 1;
	} catch (const std::exception& e) {
		std::cerr << "rank " << rank << ": " << e.what() << '\n';
		status = 1;
	}
	int global_status = 0;
	MPI_Allreduce(&status, &global_status, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
	PetscFinalize();
	return global_status;
}
