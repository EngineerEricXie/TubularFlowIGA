#ifndef IGA_PETSC_PHASE_PROFILE_HPP
#define IGA_PETSC_PHASE_PROFILE_HPP

#include "PhaseProfile.hpp"
#include <petscksp.h>

namespace iga {

// Make the setup normally deferred by KSPSolve explicit, including block-PC
// setup. PETSc recommends KSPSetUpOnBlocks for separating these timings.
// Preserve PETSc error codes and the existing KSP/PC choice and tolerances.
inline PetscErrorCode SetupProfiledKsp(KSP solver)
{
	PhaseScope phase(ProfilePhase::SolverSetup);
	PetscErrorCode error = KSPSetUp(solver);
	if (error) return error;
	return KSPSetUpOnBlocks(solver);
}

inline PetscErrorCode SolveProfiledKsp(KSP solver, Vec rhs, Vec solution)
{
	const PetscErrorCode error = SetupProfiledKsp(solver);
	if (error) return error;
	PhaseScope phase(ProfilePhase::LinearSolve);
	return KSPSolve(solver, rhs, solution);
}

} // namespace iga

#endif
