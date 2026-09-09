#include "ExecutionResources.hpp"
#include "CollectiveAssetInput.hpp"
#include <memory>
#include "IgaDatabase.hpp"
#include "TransportElement.hpp"

#include <petscsys.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, "Check IGA element Jacobians at quadrature points\n");
	int rank = 0, size = 1;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	int status = 0;
	try {
		iga::RequireExecutionResources(PETSC_COMM_WORLD, &std::cout);
		std::unique_ptr<iga::Database> database_owner;
		std::vector<iga::Element> elements;
		std::string database_fingerprint;
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "mesh check database input", [&] {
			if (argc != 2) throw std::runtime_error("usage: iga_mesh_check DATABASE.ntiga");
			database_fingerprint = iga::ReadAssetFingerprint(argv[1]);
			database_owner = std::make_unique<iga::Database>(argv[1]);
			iga::ValidatePackedExecution(database_owner->header().ranks, database_owner->header().nodes, 1, size);
			for (const auto owner : database_owner->owners())
				if (owner < 0 || owner >= size)
					throw std::runtime_error("mesh check element owner is outside communicator");
		});
		iga::RequireCollectiveSameText(PETSC_COMM_WORLD, "mesh check asset database", database_fingerprint);
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "mesh check element input", [&] {
			elements = database_owner->LoadOwned(rank);
			for (const auto& element : elements)
				if (element.owner != rank)
					throw std::runtime_error("mesh check element owner differs from ownership index");
		});
		auto& database = *database_owner;
		const auto quality = iga::InspectGeometry(elements, rank, PETSC_COMM_WORLD);
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "mesh check result logging", [&] {
			if (rank != 0) return;
			std::cout << "elements=" << database.header().elements << " minimum_detJ=" << quality.minimum_determinant
				<< " bad_elements=" << quality.bad_elements << " bad_samples=" << quality.bad_samples;
			if (quality.bad_elements) std::cout << " first_bad_element=" << quality.first_bad_element;
			std::cout << '\n';
			iga::FlushCheckedText(std::cout);
		});
		status = quality.bad_elements ? 2 : 0;
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		status = 1;
	}
	int global_status = 0;
	MPI_Allreduce(&status, &global_status, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
	PetscFinalize();
	return global_status;
}
