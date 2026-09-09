#ifndef IGA_ONE_D_PETSC_SUPPORT_HPP
#define IGA_ONE_D_PETSC_SUPPORT_HPP

#include "CollectiveFailure.hpp"
#include "ExecutionResources.hpp"
#include "PetscReadArray.hpp"

#include <petscksp.h>
#include <petscsnes.h>
#include <stdexcept>
#include <string>
#include <string_view>

namespace iga {

inline void OneDPetscCheck(PetscErrorCode error, std::string_view operation)
{
	if (error) throw std::runtime_error(std::string(operation) + " failed with PETSc error " + std::to_string(error));
}

// Call after a collective has returned on every member. This does not recover
// from a backend failure which leaves a process inside that collective.
inline void OneDCollectivePetscCheck(MPI_Comm communicator, PetscErrorCode error,
	const char* operation)
{
	CollectiveLocalStage(communicator, operation, [&] { OneDPetscCheck(error, operation); });
}

inline void OneDRequireSameInteger(MPI_Comm communicator, int value, const char* stage)
{
	int minimum = 0, maximum = 0;
	MPI_Allreduce(&value, &minimum, 1, MPI_INT, MPI_MIN, communicator);
	MPI_Allreduce(&value, &maximum, 1, MPI_INT, MPI_MAX, communicator);
	CollectiveLocalStage(communicator, stage, [&] {
		if (minimum != maximum) throw std::runtime_error("inconsistent implicit solver layout or method");
	});
}

// Owners are stack-bound to a collective solve. Local work is agreed before
// unwinding, so established distributed objects are destroyed in common order.
// Partially completed PETSc object creation/process loss is outside this contract.
struct OneDPetscObjects {
	Mat matrix = nullptr;
	Vec rhs = nullptr, solution = nullptr;
	KSP ksp = nullptr;
	SNES snes = nullptr;
	VecScatter scatter = nullptr;
	OneDPetscObjects() = default;
	OneDPetscObjects(const OneDPetscObjects&) = delete;
	OneDPetscObjects& operator=(const OneDPetscObjects&) = delete;
	~OneDPetscObjects()
	{
		if (snes) SNESDestroy(&snes);
		if (ksp) KSPDestroy(&ksp);
		if (scatter) VecScatterDestroy(&scatter);
		if (matrix) MatDestroy(&matrix);
		if (solution) VecDestroy(&solution);
		if (rhs) VecDestroy(&rhs);
	}
	void Close(MPI_Comm communicator)
	{
		OneDCollectivePetscCheck(communicator, SNESDestroy(&snes), "SNESDestroy");
		OneDCollectivePetscCheck(communicator, KSPDestroy(&ksp), "KSPDestroy");
		OneDCollectivePetscCheck(communicator, VecScatterDestroy(&scatter), "VecScatterDestroy");
		OneDCollectivePetscCheck(communicator, MatDestroy(&matrix), "MatDestroy");
		OneDCollectivePetscCheck(communicator, VecDestroy(&solution), "VecDestroy solution");
		OneDCollectivePetscCheck(communicator, VecDestroy(&rhs), "VecDestroy rhs");
	}
};

} // namespace iga

#endif
