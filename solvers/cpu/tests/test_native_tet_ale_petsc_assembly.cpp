#include "NativeTetAlePetscAssembly.hpp"
#include "NativeTetAleTransient.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <vector>

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int exit_code=0;
	try {
		int rank=0,ranks=1;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
		iga::NativeTetMesh mesh;
		mesh.points={{{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}},{{0.25,0.25,0.25}}}};
		mesh.cells={iga::NativeTetCell{1,{{4,1,2,3}}},
			iga::NativeTetCell{2,{{0,4,2,3}}},iga::NativeTetCell{3,{{0,1,4,3}}},
			iga::NativeTetCell{4,{{0,1,2,4}}}};
		const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
		const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
		const PetscInt dofs=static_cast<PetscInt>(3*velocity_nodes+mesh.points.size());
		std::vector<double> current(static_cast<std::size_t>(dofs));
		std::vector<double> committed(static_cast<std::size_t>(dofs));
		for(std::size_t dof=0;dof<current.size();++dof) {
			current[dof]=0.01*std::sin(static_cast<double>(dof+1));
			committed[dof]=0.008*std::cos(static_cast<double>(2*dof+1));
		}
		const std::vector<std::array<double,3>> grid(mesh.points.size(),{{0.03,-0.01,0.02}});
		Mat matrix=nullptr;Vec residual=nullptr;
		iga::NativeTetAlePetscCheck(MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,
			dofs,dofs,50,nullptr,50,nullptr,&matrix),"MatCreateAIJ test");
		iga::NativeTetAlePetscCheck(MatCreateVecs(matrix,nullptr,&residual),"MatCreateVecs test");
		auto builder=[&](std::size_t cell,const auto& rows) {
			std::array<double,iga::NativeTaylorHoodElementSystem::dofs> local{},old{};
			for(std::size_t entry=0;entry<rows.size();++entry) {
				local[entry]=current[static_cast<std::size_t>(rows[entry])];
				old[entry]=committed[static_cast<std::size_t>(rows[entry])];
			}
			std::array<std::array<double,3>,4> local_grid{};
			for(std::size_t node=0;node<4;++node)
				local_grid[node]=grid[mesh.cells[cell].nodes[node]];
			return iga::BuildNativeTaylorHoodAleTransientElement(mesh,mesh.cells[cell],
				local,old,local_grid,{997.0,0.0035},0.02);
		};
		iga::AssembleNativeTetAlePetscElements(mesh,topology,rank,ranks,builder,matrix,residual);

		std::vector<double> expected_residual(static_cast<std::size_t>(dofs),0.0);
		std::vector<double> expected_matrix(static_cast<std::size_t>(dofs)*dofs,0.0);
		for(std::size_t cell=0;cell<mesh.cells.size();++cell) {
			const auto rows=iga::NativeTetAlePetscCellRows(mesh,topology,cell);
			const auto element=builder(cell,rows);
			for(std::size_t row=0;row<rows.size();++row) {
				expected_residual[static_cast<std::size_t>(rows[row])]+=element.residual[row];
				for(std::size_t column=0;column<rows.size();++column)
					expected_matrix[static_cast<std::size_t>(rows[row])*dofs
						+static_cast<std::size_t>(rows[column])]
						+=element.jacobian[row*rows.size()+column];
			}
		}
		PetscInt begin=0,end=0;MatGetOwnershipRange(matrix,&begin,&end);
		double local_max=0.0;
		std::vector<PetscInt> columns(static_cast<std::size_t>(dofs));
		for(PetscInt column=0;column<dofs;++column) columns[column]=column;
		std::vector<PetscScalar> values(static_cast<std::size_t>(dofs));
		for(PetscInt row=begin;row<end;++row) {
			MatGetValues(matrix,1,&row,dofs,columns.data(),values.data());
			for(PetscInt column=0;column<dofs;++column) {
				const double expected=expected_matrix[static_cast<std::size_t>(row)*dofs+column];
				local_max=std::max(local_max,std::abs(PetscRealPart(values[column])-expected)
					/(1.0+std::abs(expected)));
			}
			PetscScalar value=0.0;VecGetValues(residual,1,&row,&value);
			const double expected=expected_residual[static_cast<std::size_t>(row)];
			local_max=std::max(local_max,std::abs(PetscRealPart(value)-expected)
				/(1.0+std::abs(expected)));
		}
		double global_max=0.0;MPI_Allreduce(&local_max,&global_max,1,MPI_DOUBLE,MPI_MAX,
			PETSC_COMM_WORLD);
		if(!(global_max<2e-13)) throw std::runtime_error("distributed ALE assembly mismatch");
		if(rank==0) std::cout<<"native tetrahedral ALE PETSc assembly passed ranks="
			<<ranks<<" relative_max="<<global_max<<'\n';
		VecDestroy(&residual);MatDestroy(&matrix);
	} catch(const std::exception& error) {
		int rank=0;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		if(rank==0) std::cerr<<error.what()<<'\n';
		exit_code=1;
	}
	PetscFinalize();return exit_code;
}
