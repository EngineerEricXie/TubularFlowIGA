#include "TransportElement.hpp"
#include <iostream>

namespace {
void Check(bool value, const char* message)
{
	if (!value) throw std::runtime_error(message);
}

iga::Element Cube(int rank)
{
	iga::Element element; element.id = rank; element.owner = rank;
	for (int k = 0; k < 4; ++k) for (int j = 0; j < 4; ++j) for (int i = 0; i < 4; ++i)
		element.bezier_points[16*k+4*j+i] = {{rank+i/3.0,j/3.0,k/3.0}};
	return element;
}

void Run(MPI_Comm comm)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&ranks);
	const auto cube = Cube(rank);
	std::vector<iga::Element> elements{cube};
	std::string message;
	try {
		iga::InspectGeometry(elements,[&](const iga::Element&) {
			if (rank == ranks-1) throw std::runtime_error("local ownership predicate failure");
			return true;
		},comm);
	} catch (const std::exception& error) { message = error.what(); }
	iga::CollectiveLocalStage(comm,"geometry failure diagnostic",[&] {
		Check(message.find("geometry inspection") != std::string::npos,"missing geometry failure diagnostic");
	});
	iga::RequireCollectiveSameText(comm,"geometry common diagnostic",message);
	const auto healthy = iga::InspectGeometry(elements,rank,comm);
	Check(healthy.bad_elements == 0 && std::abs(healthy.minimum_determinant-0.125)<1e-14,"healthy geometry changed");
	if (rank == 0) for (auto& point : elements[0].bezier_points) point[2] = 0.0;
	const auto flat = iga::InspectGeometry(elements,rank,comm);
	Check(flat.bad_elements == 1 && flat.bad_samples == 64 && flat.first_bad_element == 0
		&& flat.minimum_determinant == 0.0,"nonpositive determinant classification changed");
	elements = {cube};
	if (rank != 0) elements.clear();
	const auto sparse = iga::InspectGeometry(elements,rank,comm);
	Check(sparse.bad_elements == 0 && std::abs(sparse.minimum_determinant-0.125)<1e-14,"empty rank changed geometry reduction");
	if (rank == 0) std::cout << "geometry_inspection ranks=" << ranks << " failure_retry=passed bad_geometry=passed empty_rank=passed\n";
}
}

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int rank = 0, ranks = 1;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank); MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	try {
		Check(ranks == 3,"geometry inspection test requires 3 ranks");
		Run(PETSC_COMM_WORLD);
		MPI_Comm group = MPI_COMM_NULL;
		MPI_Comm_split(PETSC_COMM_WORLD,rank == 0 ? 0 : 1,rank,&group);
		Run(group); MPI_Comm_free(&group);
		if (rank == 0) std::cout << "geometry_inspection all comparisons passed\n";
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		MPI_Abort(PETSC_COMM_WORLD,1); return 1;
	}
	PetscFinalize(); return 0;
}
