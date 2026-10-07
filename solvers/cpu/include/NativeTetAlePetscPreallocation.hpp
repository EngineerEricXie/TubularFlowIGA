#ifndef IGA_NATIVE_TET_ALE_PETSC_PREALLOCATION_HPP
#define IGA_NATIVE_TET_ALE_PETSC_PREALLOCATION_HPP

#include "NativeTetAlePetscAssembly.hpp"
#include "NativeTetAleBoundaryControl.hpp"

namespace iga {

// Count the exact union of native 34-DOF element stencils for each owned row.
// All cells must be inspected: remote element owners also contribute to a row.
// Pressure/velocity controls and backflow keep their original numerical terms.
inline void CreateNativeTetAlePetscMatrix(const NativeTetMesh& mesh,
	const NativeTaylorHoodTopology& topology,const NativeTetAleBoundaryConditions& conditions,
	const std::map<int,NativeTetAleBoundaryFluxOperator>& boundary_operators,Mat* matrix)
{
	const std::size_t physical=3*(mesh.points.size()+topology.edges.size())+mesh.points.size();
	const auto total=physical+conditions.flow_rate_controls.size();
	if(total>static_cast<std::size_t>(std::numeric_limits<PetscInt>::max()))
		throw std::invalid_argument("native tetra matrix exceeds PETSc index range");
	PetscInt global=static_cast<PetscInt>(total),local=PETSC_DECIDE,prefix=0;
	NativeTetAlePetscCheck(PetscSplitOwnership(PETSC_COMM_WORLD,&local,&global),"symbolic ownership");
	if(MPI_Scan(&local,&prefix,1,MPIU_INT,MPI_SUM,PETSC_COMM_WORLD)!=MPI_SUCCESS)
		throw std::runtime_error("native tetra symbolic ownership scan failed");
	const PetscInt begin=prefix-local,end=prefix;
	std::vector<std::vector<PetscInt>> columns(static_cast<std::size_t>(local));
	const auto add=[&](PetscInt row,PetscInt column){
		if(row>=begin&&row<end)columns[row-begin].push_back(column);
	};
	for(PetscInt row=begin;row<end;++row)add(row,row);
	for(std::size_t cell=0;cell<mesh.cells.size();++cell){
		const auto rows=NativeTetAlePetscCellRows(mesh,topology,cell);
		for(const auto row:rows)if(row>=begin&&row<end){
			auto& entries=columns[row-begin];
			entries.insert(entries.end(),rows.begin(),rows.end());
		}
	}
	for(std::size_t control=0;control<conditions.flow_rate_controls.size();++control){
		const PetscInt row=static_cast<PetscInt>(physical+control);
		const auto& coefficients=boundary_operators.at(
			conditions.flow_rate_controls[control].boundary_label).velocity_coefficients;
		for(std::size_t dof=0;dof<coefficients.size();++dof)if(coefficients[dof]!=0.){
			add(row,static_cast<PetscInt>(dof));add(static_cast<PetscInt>(dof),row);
		}
	}
	std::vector<PetscInt> diagonal(local),off_diagonal(local);
	for(std::size_t row=0;row<columns.size();++row){
		auto& entries=columns[row];
		std::sort(entries.begin(),entries.end());
		entries.erase(std::unique(entries.begin(),entries.end()),entries.end());
		for(const auto column:entries)
			++(column>=begin&&column<end?diagonal[row]:off_diagonal[row]);
	}
	NativeTetAlePetscCheck(MatCreateAIJ(PETSC_COMM_WORLD,local,local,global,global,
		0,diagonal.data(),0,off_diagonal.data(),matrix),"native tetra exact preallocation");
	NativeTetAlePetscCheck(MatSetOption(*matrix,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE),
		"native tetra preallocation guard");
	// Establish controller slots before the element-only MatAssemblyEnd, which
	// otherwise compresses capacity unused by that first assembly phase.
	NativeTetAlePetscCheck(MatSetOption(*matrix,MAT_IGNORE_ZERO_ENTRIES,PETSC_FALSE),
		"native tetra structural zeros");
	for(PetscInt row=begin;row<end;++row){
		const auto& entries=columns[row-begin];
		const std::vector<PetscScalar> zeros(entries.size(),0.);
		NativeTetAlePetscCheck(MatSetValues(*matrix,1,&row,entries.size(),entries.data(),
			zeros.data(),INSERT_VALUES),"native tetra establish stencil");
	}
	NativeTetAlePetscCheck(MatAssemblyBegin(*matrix,MAT_FINAL_ASSEMBLY),"native stencil begin");
	NativeTetAlePetscCheck(MatAssemblyEnd(*matrix,MAT_FINAL_ASSEMBLY),"native stencil end");
}

} // namespace iga
#endif
