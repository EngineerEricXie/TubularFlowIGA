#include "PetscGather.hpp"

#include <cmath>
#include <iostream>

namespace {

void Check(bool condition, const char* message)
{
	if (!condition) throw std::runtime_error(message);
}

template <class Work>
void Failure(MPI_Comm communicator, Work&& work)
{
	std::string message;
	try { work(); } catch (const std::exception& error) { message = error.what(); }
	iga::CollectiveLocalStage(communicator,"gather test diagnostic",[&] {
		Check(message.find("global gather preparation") != std::string::npos,"missing gather failure");
	});
	iga::RequireCollectiveSameText(communicator,"gather test common diagnostic",message);
}

void Run(MPI_Comm communicator, PetscInt rows, int seed)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(communicator,&rank); MPI_Comm_size(communicator,&ranks);
	Vec source = nullptr;
	iga::RequireCollectivePetscSuccess(communicator,"test vector create",VecCreateMPI(communicator,PETSC_DECIDE,rows,&source));
	PetscInt begin = 0, end = 0;
	VecGetOwnershipRange(source,&begin,&end);
	PetscScalar* values = nullptr;
	VecGetArray(source,&values);
	for (PetscInt row = begin; row < end; ++row) values[row-begin] = seed+row+1.0;
	VecRestoreArray(source,&values);
	const auto check = [&] {
		const auto result = iga::GatherAllPetscReal(source,communicator,rows);
		iga::CollectiveLocalStage(communicator,"gather test values",[&] {
			Check(result.size() == static_cast<std::size_t>(rows),"wrong global row count");
			for (PetscInt row = 0; row < rows; ++row) Check(result[row] == seed+row+1.0,"wrong global row order/value");
		});
	};
	check();
	for (int mode = 0; mode < 2; ++mode) {
		std::uint64_t expected = rows;
		if (rank == ranks-1) expected = mode == 0 ? expected+1 : std::numeric_limits<std::uint64_t>::max();
		Failure(communicator,[&] { iga::GatherAllPetscReal(source,communicator,expected); });
		check();
	}
	if (rows > 0) {
		for (double raw : {std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity(),-0.0}) {
			VecGetArray(source,&values);
			if (begin == 0 && end > begin) values[0] = raw;
			VecRestoreArray(source,&values);
			const auto result = iga::GatherAllPetscReal(source,communicator,rows);
			iga::CollectiveLocalStage(communicator,"gather test raw values",[&] {
				Check(std::isnan(raw) ? std::isnan(result[0]) : result[0] == raw && std::signbit(result[0]) == std::signbit(raw),
					"gather changed raw nonfinite/signed-zero value");
			});
		}
		VecGetArray(source,&values);
		if (begin == 0 && end > begin) values[0] = seed+1.0;
		VecRestoreArray(source,&values);
		check();
	}
	iga::RequireCollectivePetscSuccess(communicator,"test vector destroy",VecDestroy(&source));
	if (rank == 0) std::cout << "petsc_gather ranks=" << ranks << " rows=" << rows << " cases=2 passed\n";
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int rank = 0, ranks = 1;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank); MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	try {
		Check(ranks == 3,"test_petsc_gather requires 3 ranks");
		Run(PETSC_COMM_WORLD,2,10);
		Run(PETSC_COMM_WORLD,0,10);
		MPI_Comm group = MPI_COMM_NULL;
		MPI_Comm_split(PETSC_COMM_WORLD,rank == 0 ? 0 : 1,rank,&group);
		Run(group,rank == 0 ? 3 : 1,rank == 0 ? 100 : 200);
		Run(group,0,rank == 0 ? 100 : 200);
		MPI_Comm_free(&group);
		if (rank == 0) std::cout << "petsc_gather all comparisons passed\n";
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		MPI_Abort(PETSC_COMM_WORLD,1); return 1;
	}
	PetscFinalize(); return 0;
}
