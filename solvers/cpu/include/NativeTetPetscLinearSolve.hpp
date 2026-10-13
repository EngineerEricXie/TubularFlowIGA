#ifndef IGA_NATIVE_TET_PETSC_LINEAR_SOLVE_HPP
#define IGA_NATIVE_TET_PETSC_LINEAR_SOLVE_HPP

#include "NativeTetCsrMatrix.hpp"

#include <petscksp.h>

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

struct NativeTetPetscLinearSolveDefaults
{
	std::string ksp_type = KSPPREONLY;
	std::string pc_type = PCLU;
	// Factorization package for PCLU/PCCHOLESKY; empty keeps PETSc's choice.
	std::string factor_package =
#ifdef PETSC_HAVE_MUMPS
		MATSOLVERMUMPS;
#else
		"";
#endif
	double relative_tolerance = 1e-12;
	double absolute_tolerance = 1e-14;
	PetscInt maximum_iterations = 10000;
};

namespace native_tet_petsc_linear_solve_detail {

inline void Check(PetscErrorCode code, const std::string& context)
{
	if (code) throw std::runtime_error(context+" failed with PETSc error "+std::to_string(static_cast<int>(code)));
}

struct Objects
{
	Mat matrix = nullptr;
	Vec solution = nullptr, rhs = nullptr, gathered = nullptr;
	KSP solver = nullptr;
	VecScatter scatter = nullptr;

	~Objects()
	{
		if (scatter) VecScatterDestroy(&scatter);
		if (gathered) VecDestroy(&gathered);
		if (solver) KSPDestroy(&solver);
		if (rhs) VecDestroy(&rhs);
		if (solution) VecDestroy(&solution);
		if (matrix) MatDestroy(&matrix);
	}
};

} // namespace native_tet_petsc_linear_solve_detail

// A NativeTetLinearSolve that runs on PETSC_COMM_WORLD. Every rank holds the
// same reduced CSR system; each inserts only the rows it owns, so assembly
// needs no communication, and every rank receives the whole solution.
// Options carrying `prefix` (for example -fsi_solid_ksp_type) override the
// defaults; the solve fails unless PETSc reports convergence.
inline NativeTetLinearSolve MakeNativeTetPetscLinearSolve(std::string prefix,
	NativeTetPetscLinearSolveDefaults defaults = {})
{
	return [prefix = std::move(prefix), defaults = std::move(defaults)](
		const NativeTetCsrMatrix& matrix, const std::vector<double>& right_hand_side) {
		using native_tet_petsc_linear_solve_detail::Check;
		const std::string name = "native tetrahedral PETSc solve ("+prefix+")";
		if (right_hand_side.size() != matrix.size || matrix.row_offsets.size() != matrix.size+1)
			throw std::invalid_argument(name+": system size mismatch");
		native_tet_petsc_linear_solve_detail::Objects petsc;
		const auto size = static_cast<PetscInt>(matrix.size);
		Check(VecCreate(PETSC_COMM_WORLD, &petsc.rhs), name+" VecCreate");
		Check(VecSetSizes(petsc.rhs, PETSC_DECIDE, size), name+" VecSetSizes");
		Check(VecSetFromOptions(petsc.rhs), name+" VecSetFromOptions");
		Check(VecDuplicate(petsc.rhs, &petsc.solution), name+" VecDuplicate");
		PetscInt begin = 0, end = 0;
		Check(VecGetOwnershipRange(petsc.rhs, &begin, &end), name+" VecGetOwnershipRange");
		std::vector<PetscInt> diagonal(static_cast<std::size_t>(end-begin)), off_diagonal(diagonal.size());
		for (PetscInt row = begin; row < end; ++row)
			for (auto entry = matrix.row_offsets[row]; entry < matrix.row_offsets[row+1]; ++entry) {
				const PetscInt column = static_cast<PetscInt>(matrix.columns[entry]);
				++(column >= begin && column < end ? diagonal : off_diagonal)[static_cast<std::size_t>(row-begin)];
			}
		Check(MatCreate(PETSC_COMM_WORLD, &petsc.matrix), name+" MatCreate");
		Check(MatSetSizes(petsc.matrix, end-begin, end-begin, size, size), name+" MatSetSizes");
		Check(MatSetType(petsc.matrix, MATAIJ), name+" MatSetType");
		Check(MatSeqAIJSetPreallocation(petsc.matrix, 0, diagonal.data()), name+" MatSeqAIJSetPreallocation");
		Check(MatMPIAIJSetPreallocation(petsc.matrix, 0, diagonal.data(), 0, off_diagonal.data()),
			name+" MatMPIAIJSetPreallocation");
		std::vector<PetscInt> columns;
		for (PetscInt row = begin; row < end; ++row) {
			const auto first = matrix.row_offsets[row], last = matrix.row_offsets[row+1];
			columns.assign(matrix.columns.begin()+static_cast<std::ptrdiff_t>(first),
				matrix.columns.begin()+static_cast<std::ptrdiff_t>(last));
			Check(MatSetValues(petsc.matrix, 1, &row, static_cast<PetscInt>(columns.size()), columns.data(),
				matrix.values.data()+first, INSERT_VALUES), name+" MatSetValues");
			Check(VecSetValue(petsc.rhs, row, right_hand_side[static_cast<std::size_t>(row)], INSERT_VALUES),
				name+" VecSetValue");
		}
		Check(MatAssemblyBegin(petsc.matrix, MAT_FINAL_ASSEMBLY), name+" MatAssemblyBegin");
		Check(MatAssemblyEnd(petsc.matrix, MAT_FINAL_ASSEMBLY), name+" MatAssemblyEnd");
		Check(VecAssemblyBegin(petsc.rhs), name+" VecAssemblyBegin");
		Check(VecAssemblyEnd(petsc.rhs), name+" VecAssemblyEnd");
		Check(KSPCreate(PETSC_COMM_WORLD, &petsc.solver), name+" KSPCreate");
		Check(KSPSetOptionsPrefix(petsc.solver, prefix.c_str()), name+" KSPSetOptionsPrefix");
		Check(KSPSetOperators(petsc.solver, petsc.matrix, petsc.matrix), name+" KSPSetOperators");
		Check(KSPSetType(petsc.solver, defaults.ksp_type.c_str()), name+" KSPSetType");
		PC preconditioner = nullptr;
		Check(KSPGetPC(petsc.solver, &preconditioner), name+" KSPGetPC");
		Check(PCSetType(preconditioner, defaults.pc_type.c_str()), name+" PCSetType");
		if (!defaults.factor_package.empty())
			Check(PCFactorSetMatSolverType(preconditioner, defaults.factor_package.c_str()),
				name+" PCFactorSetMatSolverType");
		Check(KSPSetTolerances(petsc.solver, defaults.relative_tolerance, defaults.absolute_tolerance,
			PETSC_DEFAULT, defaults.maximum_iterations), name+" KSPSetTolerances");
		Check(KSPSetFromOptions(petsc.solver), name+" KSPSetFromOptions");
		Check(KSPSolve(petsc.solver, petsc.rhs, petsc.solution), name+" KSPSolve");
		KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
		Check(KSPGetConvergedReason(petsc.solver, &reason), name+" KSPGetConvergedReason");
		if (reason <= 0)
			throw std::runtime_error(name+" did not converge (reason "+std::to_string(static_cast<int>(reason))+")");
		Check(VecScatterCreateToAll(petsc.solution, &petsc.scatter, &petsc.gathered), name+" VecScatterCreateToAll");
		Check(VecScatterBegin(petsc.scatter, petsc.solution, petsc.gathered, INSERT_VALUES, SCATTER_FORWARD),
			name+" VecScatterBegin");
		Check(VecScatterEnd(petsc.scatter, petsc.solution, petsc.gathered, INSERT_VALUES, SCATTER_FORWARD),
			name+" VecScatterEnd");
		const PetscScalar* values = nullptr;
		Check(VecGetArrayRead(petsc.gathered, &values), name+" VecGetArrayRead");
		std::vector<double> solution(matrix.size);
		for (std::size_t row = 0; row < matrix.size; ++row) solution[row] = PetscRealPart(values[row]);
		Check(VecRestoreArrayRead(petsc.gathered, &values), name+" VecRestoreArrayRead");
		for (const double value : solution)
			if (!std::isfinite(value)) throw std::runtime_error(name+" returned a nonfinite solution");
		return solution;
	};
}

} // namespace iga

#endif
