#ifndef IGA_ONE_D_PETSC_SUPPORT_HPP
#define IGA_ONE_D_PETSC_SUPPORT_HPP

#include "CollectiveFailure.hpp"
#include "ExecutionResources.hpp"
#include "PetscReadArray.hpp"
#include "PetscSolverOptions.hpp"

#include <petscksp.h>
#include <petscsnes.h>
#include <stdexcept>
#include <string>
#include <string_view>
#include <iomanip>
#include <ostream>

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

// Persistent runtime snapshot; PETSc objects created by an advance borrow it.
// Construct outside local-only preparation stages because options agreement
// uses this communicator. Diagnostics describe the last successful solve.
struct OneDPetscSolverContext {
	PetscSolverOptions options;
	PetscKspConfiguration linear;
	std::string nonlinear_type;
	PetscInt nonlinear_iterations = 0;
	SNESConvergedReason nonlinear_reason = SNES_CONVERGED_ITERATING;
	bool has_result = false, nonlinear = false;
	OneDPetscSolverContext(MPI_Comm communicator, const std::string& prefix = {},
		const std::set<std::string>& excluded = {}) : options(communicator, prefix, nullptr, {}, excluded) {}
	void Record(MPI_Comm communicator, KSP solver, SNES snes = nullptr)
	{
		CollectiveLocalStage(communicator, "1d effective solver configuration", [&] {
			linear = CaptureKspConfiguration(solver);
			nonlinear = snes != nullptr;
			if (snes) {
				const char* name = nullptr;
				OneDPetscCheck(SNESGetType(snes, &name), "SNESGetType");
				nonlinear_type = name ? name : "unset";
				OneDPetscCheck(SNESGetIterationNumber(snes, &nonlinear_iterations), "SNESGetIterationNumber");
				OneDPetscCheck(SNESGetConvergedReason(snes, &nonlinear_reason), "SNESGetConvergedReason");
			}
			has_result = true;
		});
		options.RecordUsed();
	}
};

// Local output: callers coordinate stream failures on the runtime communicator.
inline void WriteOneDPetscSolverConfiguration(std::ostream& output, const OneDPetscSolverContext& context, int step)
{
	if (!context.has_result) return;
	const auto precision = output.precision();
	const auto& linear = context.linear;
	output << std::setprecision(17) << "one_d_solver_configuration {\"prefix\":" << std::quoted(linear.prefix)
		<< ",\"step\":" << step << ",\"ksp\":" << std::quoted(linear.ksp)
		<< ",\"pc\":" << std::quoted(linear.pc) << ",\"factor_backend\":" << std::quoted(linear.factor_backend)
		<< ",\"rtol\":" << linear.relative_tolerance << ",\"atol\":" << linear.absolute_tolerance
		<< ",\"last_iterations\":" << linear.last_iterations << ",\"last_reason\":" << static_cast<int>(linear.last_reason);
	if (context.nonlinear) output << ",\"snes\":" << std::quoted(context.nonlinear_type)
		<< ",\"snes_iterations\":" << context.nonlinear_iterations << ",\"snes_reason\":" << static_cast<int>(context.nonlinear_reason);
	output << "}\n";
	output.precision(precision);
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
