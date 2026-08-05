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
		if (argc != 2) throw std::runtime_error("usage: iga_mesh_check DATABASE.ntiga");
		iga::Database database(argv[1]);
		if (database.header().ranks != static_cast<std::uint32_t>(size))
			throw std::runtime_error("MPI size must match database rank count");
		const auto quality = iga::InspectGeometry(database.LoadOwned(rank), rank, PETSC_COMM_WORLD);
		if (rank == 0) std::cout << "elements=" << database.header().elements << " minimum_detJ=" << quality.minimum_determinant
			<< " bad_elements=" << quality.bad_elements << " bad_samples=" << quality.bad_samples;
		if (rank == 0 && quality.bad_elements) std::cout << " first_bad_element=" << quality.first_bad_element;
		if (rank == 0) std::cout << '\n';
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
