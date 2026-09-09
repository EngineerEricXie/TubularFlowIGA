#include <petscvec.h>
#include <iostream>
#include <stdexcept>

// Interoperability oracle: the original distributed VecLoad/VecView path.
// Production checkpoint I/O now uses COMM_SELF, but must keep this format.
int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	auto check = [](PetscErrorCode code) {
		if (code) throw std::runtime_error("legacy distributed checkpoint I/O failed");
	};
	try {
		if (argc != 3) throw std::runtime_error("require INPUT.state OUTPUT.state");
		Vec state = nullptr;
		PetscViewer viewer = nullptr;
		check(VecCreate(PETSC_COMM_WORLD, &state));
		check(PetscViewerBinaryOpen(PETSC_COMM_WORLD, argv[1], FILE_MODE_READ, &viewer));
		check(VecLoad(state, viewer));
		check(PetscViewerDestroy(&viewer));
		check(PetscViewerBinaryOpen(PETSC_COMM_WORLD, argv[2], FILE_MODE_WRITE, &viewer));
		check(VecView(state, viewer));
		check(PetscViewerDestroy(&viewer));
		check(VecDestroy(&state));
		int rank = 0;
		MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
		std::cout << "checkpoint format rank=" << rank << " passed\n";
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n';
		MPI_Abort(PETSC_COMM_WORLD, 2);
	}
	PetscFinalize();
	return 0;
}
