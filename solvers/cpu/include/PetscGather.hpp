#ifndef IGA_PETSC_GATHER_HPP
#define IGA_PETSC_GATHER_HPP

#include "CollectiveFailure.hpp"
#include "PetscReadArray.hpp"

#include <cstdint>
#include <limits>
#include <vector>

namespace iga {

// Temporary handles belong to a collectively entered gather. Local failures
// are agreed before unwinding, so established objects are destroyed in common
// order. Partially failed PETSc object creation/process loss is not recoverable
// by this owner. Close checks normal cleanup; the destructor is an unwind path.
struct PetscGatherObjects {
	Vec all = nullptr;
	VecScatter scatter = nullptr;
	PetscGatherObjects() = default;
	PetscGatherObjects(const PetscGatherObjects&) = delete;
	PetscGatherObjects& operator=(const PetscGatherObjects&) = delete;
	~PetscGatherObjects()
	{
		if (scatter) VecScatterDestroy(&scatter);
		if (all) VecDestroy(&all);
	}
	void Close(MPI_Comm communicator)
	{
		RequireCollectivePetscSuccess(communicator, "global gather destroy scatter", VecScatterDestroy(&scatter));
		RequireCollectivePetscSuccess(communicator, "global gather destroy vector", VecDestroy(&all));
	}
};

// Preserve global row order and raw values, including nonfinite values needed
// for diagnostics. This deliberately replicates the full vector; it is not a
// replacement for the owned/ghost path used by distributed element assembly.
inline std::vector<double> GatherAllPetscReal(Vec source, MPI_Comm communicator,
	std::uint64_t expected_rows)
{
	std::vector<double> result;
	CollectiveLocalStage(communicator, "global gather preparation", [&] {
		PetscInt rows = 0;
		if (VecGetSize(source, &rows)) throw std::runtime_error("VecGetSize failed");
		if (rows < 0 || static_cast<std::uint64_t>(rows) != expected_rows
			|| expected_rows > result.max_size())
			throw std::runtime_error("global gather vector size does not match expected rows");
		result.resize(static_cast<std::size_t>(expected_rows));
	});
	PetscGatherObjects objects;
	RequireCollectivePetscSuccess(communicator, "global gather create", VecScatterCreateToAll(source, &objects.scatter, &objects.all));
	RequireCollectivePetscSuccess(communicator, "global gather begin", VecScatterBegin(objects.scatter, source, objects.all, INSERT_VALUES, SCATTER_FORWARD));
	RequireCollectivePetscSuccess(communicator, "global gather end", VecScatterEnd(objects.scatter, source, objects.all, INSERT_VALUES, SCATTER_FORWARD));
	CollectiveLocalStage(communicator, "global gather result", [&] {
		PetscInt count = 0;
		if (VecGetLocalSize(objects.all, &count)) throw std::runtime_error("gathered VecGetLocalSize failed");
		if (count < 0 || static_cast<std::uint64_t>(count) != expected_rows)
			throw std::runtime_error("gathered vector has an unexpected layout");
		PetscReadArray view;
		view.Acquire(objects.all);
		for (std::size_t i = 0; i < result.size(); ++i) result[i] = PetscRealPart(view.Data()[i]);
		view.Restore();
	});
	objects.Close(communicator);
	return result;
}

} // namespace iga

#endif
