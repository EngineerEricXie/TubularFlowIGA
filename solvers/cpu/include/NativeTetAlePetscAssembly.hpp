#ifndef IGA_NATIVE_TET_ALE_PETSC_ASSEMBLY_HPP
#define IGA_NATIVE_TET_ALE_PETSC_ASSEMBLY_HPP

#include "NativeTetFem.hpp"

#include <petscmat.h>
#include <petscvec.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace iga {

inline void NativeTetAlePetscCheck(PetscErrorCode error,const char* context)
{
	if(error) throw std::runtime_error(std::string("PETSc ALE assembly failure in ")+context);
}

inline std::array<PetscInt,NativeTaylorHoodElementSystem::dofs>
NativeTetAlePetscCellRows(const NativeTetMesh& mesh,
	const NativeTaylorHoodTopology& topology,std::size_t cell)
{
	const std::uint64_t velocity_nodes=mesh.points.size()+topology.edges.size();
	std::array<PetscInt,NativeTaylorHoodElementSystem::dofs> rows{};
	std::size_t entry=0;
	for(const auto node:topology.cell_velocity_nodes.at(cell))
		for(int component=0;component<3;++component)
			rows[entry++]=static_cast<PetscInt>(3*static_cast<std::uint64_t>(node)+component);
	for(const auto node:mesh.cells.at(cell).nodes)
		rows[entry++]=static_cast<PetscInt>(3*velocity_nodes+node);
	return rows;
}

// All element physics stays in the caller-provided native builder. This layer
// only partitions cells and performs collective PETSc ADD_VALUES assembly.
template<class ElementBuilder>
void AssembleNativeTetAlePetscElements(const NativeTetMesh& mesh,
	const NativeTaylorHoodTopology& topology,int rank,int ranks,
	ElementBuilder&& build_element,Mat jacobian,Vec residual)
{
	if(rank<0||ranks<=0||rank>=ranks)
		throw std::invalid_argument("native ALE PETSc rank partition is invalid");
	NativeTetAlePetscCheck(MatZeroEntries(jacobian),"MatZeroEntries");
	NativeTetAlePetscCheck(VecSet(residual,0.0),"VecSet residual");
	for(std::size_t cell=static_cast<std::size_t>(rank);cell<mesh.cells.size();
			cell+=static_cast<std::size_t>(ranks)) {
		const auto rows=NativeTetAlePetscCellRows(mesh,topology,cell);
		const auto system=build_element(cell,rows);
		NativeTetAlePetscCheck(MatSetValues(jacobian,static_cast<PetscInt>(rows.size()),
			rows.data(),static_cast<PetscInt>(rows.size()),rows.data(),system.jacobian.data(),
			ADD_VALUES),"MatSetValues element");
		NativeTetAlePetscCheck(VecSetValues(residual,static_cast<PetscInt>(rows.size()),
			rows.data(),system.residual.data(),ADD_VALUES),"VecSetValues element");
	}
	NativeTetAlePetscCheck(MatAssemblyBegin(jacobian,MAT_FINAL_ASSEMBLY),
		"MatAssemblyBegin");
	NativeTetAlePetscCheck(MatAssemblyEnd(jacobian,MAT_FINAL_ASSEMBLY),"MatAssemblyEnd");
	NativeTetAlePetscCheck(VecAssemblyBegin(residual),"VecAssemblyBegin");
	NativeTetAlePetscCheck(VecAssemblyEnd(residual),"VecAssemblyEnd");
}

} // namespace iga

#endif
