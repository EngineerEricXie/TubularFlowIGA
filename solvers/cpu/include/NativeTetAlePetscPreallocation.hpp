#ifndef IGA_NATIVE_TET_ALE_PETSC_PREALLOCATION_HPP
#define IGA_NATIVE_TET_ALE_PETSC_PREALLOCATION_HPP

#include "NativeTetAlePetscAssembly.hpp"
#include "NativeTetAleBoundaryControl.hpp"

namespace iga {

// Count the exact union of native 34-DOF element stencils for each owned row.
// All cells must be inspected: remote element owners also contribute to a row.
// Rows of one velocity node (three components) or one pressure node share a
// stencil, so it is formed once per node from that node's incident cells
// instead of storing every duplicated column. Pressure/velocity controls and
// backflow keep their original numerical terms; flow and resistance controls
// are bordered rows/columns after the physical unknowns, in that order.
inline void CreateNativeTetAlePetscMatrix(const NativeTetMesh& mesh,
	const NativeTaylorHoodTopology& topology,const NativeTetAleBoundaryConditions& conditions,
	const std::map<int,NativeTetAleBoundaryFluxOperator>& boundary_operators,Mat* matrix)
{
	const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
	const std::size_t physical=3*velocity_nodes+mesh.points.size();
	const std::size_t flow_controls=conditions.flow_rate_controls.size();
	const auto total=physical+flow_controls+conditions.resistance_controls.size();
	if(total>static_cast<std::size_t>(std::numeric_limits<PetscInt>::max()))
		throw std::invalid_argument("native tetra matrix exceeds PETSc index range");
	PetscInt global=static_cast<PetscInt>(total),local=PETSC_DECIDE,prefix=0;
	NativeTetAlePetscCheck(PetscSplitOwnership(PETSC_COMM_WORLD,&local,&global),"symbolic ownership");
	if(MPI_Scan(&local,&prefix,1,MPIU_INT,MPI_SUM,PETSC_COMM_WORLD)!=MPI_SUCCESS)
		throw std::runtime_error("native tetra symbolic ownership scan failed");
	const PetscInt begin=prefix-local,end=prefix;
	const auto owned=[&](std::size_t row){
		return static_cast<PetscInt>(row)>=begin&&static_cast<PetscInt>(row)<end;
	};
	// Node n<velocity_nodes is a velocity node (rows 3n..3n+2); node
	// velocity_nodes+v is the pressure node of vertex v (row 3*velocity_nodes+v).
	const std::size_t stencil_nodes=velocity_nodes+mesh.points.size();
	std::vector<std::size_t> first(stencil_nodes+1,0);
	const auto touches=[&](std::size_t node){
		const std::size_t row=node<velocity_nodes?3*node:3*velocity_nodes+(node-velocity_nodes);
		const std::size_t last=node<velocity_nodes?row+2:row;
		return static_cast<PetscInt>(last)>=begin&&static_cast<PetscInt>(row)<end;
	};
	for(std::size_t cell=0;cell<mesh.cells.size();++cell){
		for(const auto node:topology.cell_velocity_nodes[cell])if(touches(node))++first[node+1];
		for(const auto vertex:mesh.cells[cell].nodes)
			if(touches(velocity_nodes+vertex))++first[velocity_nodes+vertex+1];
	}
	for(std::size_t node=0;node<stencil_nodes;++node)first[node+1]+=first[node];
	std::vector<std::uint32_t> incident(first.back());
	{
		auto cursor=first;
		for(std::size_t cell=0;cell<mesh.cells.size();++cell){
			for(const auto node:topology.cell_velocity_nodes[cell])
				if(touches(node))incident[cursor[node]++]=static_cast<std::uint32_t>(cell);
			for(const auto vertex:mesh.cells[cell].nodes)
				if(touches(velocity_nodes+vertex))
					incident[cursor[velocity_nodes+vertex]++]=static_cast<std::uint32_t>(cell);
		}
	}
	std::vector<PetscInt> columns;
	const auto stencil=[&](std::size_t node){
		columns.clear();
		for(std::size_t entry=first[node];entry<first[node+1];++entry){
			const auto rows=NativeTetAlePetscCellRows(mesh,topology,incident[entry]);
			columns.insert(columns.end(),rows.begin(),rows.end());
		}
		std::sort(columns.begin(),columns.end());
		columns.erase(std::unique(columns.begin(),columns.end()),columns.end());
	};
	// Bordered control couplings: row -> extra columns, plus the transposed
	// velocity entries of every controlled label.
	std::map<std::size_t,std::vector<PetscInt>> control_columns;
	const auto border=[&](std::size_t row,const NativeTetAleBoundaryFluxOperator& operator_value,
		bool pressure){
		auto& entries=control_columns[row];
		for(const auto dof:operator_value.velocity_dofs){
			entries.push_back(static_cast<PetscInt>(dof));
			if(owned(dof))control_columns[dof].push_back(static_cast<PetscInt>(row));
		}
		if(pressure)for(const auto node:operator_value.pressure_nodes)
			entries.push_back(static_cast<PetscInt>(3*velocity_nodes+node));
	};
	for(std::size_t control=0;control<flow_controls;++control){
		const auto& specification=conditions.flow_rate_controls[control];
		border(physical+control,boundary_operators.at(specification.boundary_label),
			specification.target_pressure_slope_m3_s_pa!=0.0);
	}
	for(std::size_t control=0;control<conditions.resistance_controls.size();++control)
		border(physical+flow_controls+control,boundary_operators.at(
			conditions.resistance_controls[control].boundary_label),false);
	const auto row_columns=[&](std::size_t row,const std::vector<PetscInt>& base,
		std::vector<PetscInt>& output){
		output=base;
		output.push_back(static_cast<PetscInt>(row));
		const auto extra=control_columns.find(row);
		if(extra!=control_columns.end())output.insert(output.end(),extra->second.begin(),extra->second.end());
		std::sort(output.begin(),output.end());
		output.erase(std::unique(output.begin(),output.end()),output.end());
	};
	const auto row_node=[&](std::size_t row){
		return row<3*velocity_nodes?row/3:velocity_nodes+(row-3*velocity_nodes);
	};
	std::vector<PetscInt> diagonal(local),off_diagonal(local);
	std::vector<PetscInt> entries;
	const std::vector<PetscInt> no_columns;
	std::size_t cached=std::numeric_limits<std::size_t>::max();
	for(PetscInt row=begin;row<end;++row){
		const auto unsigned_row=static_cast<std::size_t>(row);
		if(unsigned_row<physical){
			const auto node=row_node(unsigned_row);
			if(node!=cached){stencil(node);cached=node;}
			row_columns(unsigned_row,columns,entries);
		}else row_columns(unsigned_row,no_columns,entries);
		for(const auto column:entries)
			++(column>=begin&&column<end?diagonal[row-begin]:off_diagonal[row-begin]);
	}
	NativeTetAlePetscCheck(MatCreateAIJ(PETSC_COMM_WORLD,local,local,global,global,
		0,diagonal.data(),0,off_diagonal.data(),matrix),"native tetra exact preallocation");
	NativeTetAlePetscCheck(MatSetOption(*matrix,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE),
		"native tetra preallocation guard");
	// Establish controller slots before the element-only MatAssemblyEnd, which
	// otherwise compresses capacity unused by that first assembly phase.
	NativeTetAlePetscCheck(MatSetOption(*matrix,MAT_IGNORE_ZERO_ENTRIES,PETSC_FALSE),
		"native tetra structural zeros");
	cached=std::numeric_limits<std::size_t>::max();
	std::vector<PetscScalar> zeros;
	for(PetscInt row=begin;row<end;++row){
		const auto unsigned_row=static_cast<std::size_t>(row);
		if(unsigned_row<physical){
			const auto node=row_node(unsigned_row);
			if(node!=cached){stencil(node);cached=node;}
			row_columns(unsigned_row,columns,entries);
		}else row_columns(unsigned_row,no_columns,entries);
		zeros.assign(entries.size(),0.);
		NativeTetAlePetscCheck(MatSetValues(*matrix,1,&row,static_cast<PetscInt>(entries.size()),
			entries.data(),zeros.data(),INSERT_VALUES),"native tetra establish stencil");
	}
	NativeTetAlePetscCheck(MatAssemblyBegin(*matrix,MAT_FINAL_ASSEMBLY),"native stencil begin");
	NativeTetAlePetscCheck(MatAssemblyEnd(*matrix,MAT_FINAL_ASSEMBLY),"native stencil end");
}

} // namespace iga
#endif
