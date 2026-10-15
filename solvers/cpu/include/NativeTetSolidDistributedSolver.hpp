#ifndef IGA_NATIVE_TET_SOLID_DISTRIBUTED_SOLVER_HPP
#define IGA_NATIVE_TET_SOLID_DISTRIBUTED_SOLVER_HPP

// Static hyperelastic Newton solve on a partitioned tetrahedral mesh. Each
// rank assembles the residual and the tangent of its owned cells only; the
// residual is summed into the replicated layout (one Allreduce per assembly),
// the tangent goes into one PETSc matrix over the owner-major vertex
// numbering, and the prescribed DOFs are eliminated with MatZeroRowsColumns.
// Vectors stay replicated, so the Newton loop, the line search, and the
// result are those of SolveNativeTetSolidStatic to round-off.
#include "NativeTetMeshPartition.hpp"
#include "NativeTetPetscLinearSolve.hpp"
#include "NativeTetSolidStaticSolver.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

namespace native_tet_solid_distributed_detail {

inline void Check(PetscErrorCode code,const char* context)
{
	if(code)throw std::runtime_error(std::string("native tetrahedral distributed solid: ")+context
		+" failed with PETSc error "+std::to_string(static_cast<int>(code)));
}

struct Objects
{
	Mat matrix=nullptr;Vec rhs=nullptr,update=nullptr,gathered=nullptr;VecScatter scatter=nullptr;KSP solver=nullptr;
	~Objects()
	{
		if(solver)KSPDestroy(&solver);
		if(scatter)VecScatterDestroy(&scatter);
		if(gathered)VecDestroy(&gathered);
		if(update)VecDestroy(&update);
		if(rhs)VecDestroy(&rhs);
		if(matrix)MatDestroy(&matrix);
	}
};

struct Assembly
{
	std::vector<double> residual;   // replicated
	double energy=0.0;
	double minimum_jacobian=std::numeric_limits<double>::infinity();
	bool failed=false;              // an owned element could not be evaluated on some rank
};

} // namespace native_tet_solid_distributed_detail

inline NativeTetSolidStaticResult SolveNativeTetSolidStaticDistributed(const NativeTetMeshPartition& partition,
	const NativeTetMesh& mesh,const NativeTetSolidMaterial& material,const std::vector<double>& external_force_n,
	const std::map<std::size_t,double>& prescribed_displacement_m,std::vector<double> initial_displacement_m,
	const NativeTetSolidStaticOptions& options,MPI_Comm comm,const std::string& prefix="fsi_solid_",
	const NativeTetPetscLinearSolveDefaults& defaults={})
{
	using namespace native_tet_solid_distributed_detail;
	ValidateNativeTetSolidMaterial(material);native_tet_solid_static_detail::ValidateOptions(options);
	const std::size_t nodes=mesh.points.size(),dofs=3*nodes;
	if(mesh.points.empty()||mesh.cells.empty())throw std::invalid_argument("native tetrahedral solid static mesh is empty");
	if(partition.node_owner.size()!=nodes||partition.cell_owner.size()!=mesh.cells.size())
		throw std::invalid_argument("native tetrahedral distributed solid: partition differs from the mesh");
	if(external_force_n.size()!=dofs)throw std::invalid_argument("native tetrahedral solid external-force size mismatch");
	for(double value:external_force_n)if(!std::isfinite(value))throw std::invalid_argument("native tetrahedral solid external force is nonfinite");
	if(initial_displacement_m.empty())initial_displacement_m.assign(dofs,0.0);
	if(initial_displacement_m.size()!=dofs)throw std::invalid_argument("native tetrahedral solid initial-state size mismatch");
	for(const auto& item:prescribed_displacement_m)if(item.first>=dofs||!std::isfinite(item.second))
		throw std::invalid_argument("native tetrahedral solid prescribed displacement is invalid");
	std::vector<std::size_t> free;free.reserve(dofs-prescribed_displacement_m.size());
	for(std::size_t dof=0;dof<dofs;++dof)if(!prescribed_displacement_m.count(dof))free.push_back(dof);
	if(free.empty())throw std::invalid_argument("native tetrahedral solid static solve has no free DOF");
	const auto& local=partition.local;
	const auto row_of=[&](std::uint32_t local_vertex,int component){
		return static_cast<PetscInt>(3*partition.OwnerMajorNode(partition.local_to_global_node[local_vertex])+component);};
	const PetscInt begin=3*static_cast<PetscInt>(partition.owned_node_begin[partition.rank]),
		owned=3*static_cast<PetscInt>(partition.owned_node_count),end=begin+owned,global=static_cast<PetscInt>(dofs);
	// Exact preallocation from every local cell (each owned vertex has all its cells locally).
	std::vector<std::vector<PetscInt>> columns(static_cast<std::size_t>(owned));
	for(PetscInt row=begin;row<end;++row)columns[row-begin].push_back(row);
	for(const auto& cell:local.cells){
		std::array<PetscInt,12> rows{};
		for(int a=0;a<4;++a)for(int i=0;i<3;++i)rows[3*a+i]=row_of(cell.nodes[a],i);
		for(const auto row:rows)if(row>=begin&&row<end){auto& entries=columns[row-begin];entries.insert(entries.end(),rows.begin(),rows.end());}
	}
	std::vector<PetscInt> diagonal(static_cast<std::size_t>(owned)),off_diagonal(static_cast<std::size_t>(owned));
	for(std::size_t row=0;row<columns.size();++row){
		auto& entries=columns[row];std::sort(entries.begin(),entries.end());
		entries.erase(std::unique(entries.begin(),entries.end()),entries.end());
		for(const auto column:entries)++(column>=begin&&column<end?diagonal[row]:off_diagonal[row]);
	}
	Objects petsc;
	Check(MatCreateAIJ(comm,owned,owned,global,global,0,diagonal.data(),0,off_diagonal.data(),&petsc.matrix),"MatCreateAIJ");
	Check(MatSetOption(petsc.matrix,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE),"MatSetOption");
	Check(MatCreateVecs(petsc.matrix,&petsc.update,&petsc.rhs),"MatCreateVecs");
	Check(VecScatterCreateToAll(petsc.update,&petsc.scatter,&petsc.gathered),"VecScatterCreateToAll");
	std::vector<PetscInt> constrained_rows;
	for(std::uint32_t local_vertex=0;local_vertex<partition.owned_node_count;++local_vertex)
		for(int component=0;component<3;++component)
			if(prescribed_displacement_m.count(3*partition.local_to_global_node[local_vertex]+component))
				constrained_rows.push_back(row_of(local_vertex,component));
	// Residual (replicated after the reduction) and optionally the tangent of
	// the owned cells; failures are collective so every rank takes the same branch.
	const auto assemble=[&](const std::vector<double>& displacement,bool tangent){
		Assembly result;std::vector<double> partial(dofs,0.0);double energy=0.0;
		double minimum_jacobian=std::numeric_limits<double>::infinity();int failed=0;
		if(tangent){
			Check(MatZeroEntries(petsc.matrix),"MatZeroEntries");
			// Every owned diagonal exists, including those of vertices no cell uses.
			for(PetscInt row=begin;row<end;++row)Check(MatSetValue(petsc.matrix,row,row,0.0,ADD_VALUES),"MatSetValue diagonal");
		}
		try{
			for(std::size_t cell=0;cell<local.cells.size();++cell){
				if(!partition.local_cell_owned[cell])continue;
				std::array<double,12> local_displacement{};std::array<PetscInt,12> rows{};
				for(std::size_t a=0;a<4;++a)for(int i=0;i<3;++i){
					local_displacement[3*a+i]=displacement[3*partition.local_to_global_node[local.cells[cell].nodes[a]]+i];
					rows[3*a+i]=row_of(local.cells[cell].nodes[a],i);
				}
				const auto element=BuildNativeTetHyperelasticSolidElement(local,local.cells[cell],local_displacement,material);
				energy+=element.strain_energy_j;minimum_jacobian=std::min(minimum_jacobian,element.deformation_jacobian);
				for(std::size_t a=0;a<4;++a)for(int i=0;i<3;++i)
					partial[3*partition.local_to_global_node[local.cells[cell].nodes[a]]+i]+=element.residual_n[3*a+i];
				if(tangent)Check(MatSetValues(petsc.matrix,12,rows.data(),12,rows.data(),element.tangent_n_m.data(),ADD_VALUES),"MatSetValues");
			}
		}catch(const std::runtime_error&){failed=1;}
		int any_failed=0;
		if(MPI_Allreduce(&failed,&any_failed,1,MPI_INT,MPI_MAX,comm)!=MPI_SUCCESS)
			throw std::runtime_error("native tetrahedral distributed solid: failure reduction failed");
		if(tangent){
			Check(MatAssemblyBegin(petsc.matrix,MAT_FINAL_ASSEMBLY),"MatAssemblyBegin");
			Check(MatAssemblyEnd(petsc.matrix,MAT_FINAL_ASSEMBLY),"MatAssemblyEnd");
		}
		if(any_failed){result.failed=true;return result;}
		result.residual.assign(dofs,0.0);
		if(MPI_Allreduce(partial.data(),result.residual.data(),static_cast<int>(dofs),MPI_DOUBLE,MPI_SUM,comm)!=MPI_SUCCESS
			||MPI_Bcast(result.residual.data(),static_cast<int>(dofs),MPI_DOUBLE,0,comm)!=MPI_SUCCESS
			||MPI_Allreduce(&energy,&result.energy,1,MPI_DOUBLE,MPI_SUM,comm)!=MPI_SUCCESS
			||MPI_Allreduce(&minimum_jacobian,&result.minimum_jacobian,1,MPI_DOUBLE,MPI_MIN,comm)!=MPI_SUCCESS)
			throw std::runtime_error("native tetrahedral distributed solid: residual reduction failed");
		for(std::size_t dof=0;dof<dofs;++dof)result.residual[dof]-=external_force_n[dof];
		if(tangent)Check(MatZeroRowsColumns(petsc.matrix,static_cast<PetscInt>(constrained_rows.size()),constrained_rows.data(),1.0,nullptr,nullptr),
			"MatZeroRowsColumns");
		return result;
	};
	auto displacement=std::move(initial_displacement_m);
	for(const auto& item:prescribed_displacement_m)displacement[item.first]=item.second;
	auto assembly=assemble(displacement,true);
	if(assembly.failed)throw std::runtime_error("native tetrahedral distributed solid: initial state has inverted elements");
	NativeTetSolidStaticResult result;
	result.initial_free_residual_n=native_tet_solid_static_detail::FreeNorm(assembly.residual,free);
	const double threshold=options.absolute_tolerance_n+options.relative_tolerance*result.initial_free_residual_n;
	std::vector<double> update(dofs);
	for(int iteration=0;iteration<=options.maximum_iterations;++iteration){
		const double norm=native_tet_solid_static_detail::FreeNorm(assembly.residual,free);
		if(norm<=threshold){result.iterations=iteration;result.converged=true;break;}
		if(iteration==options.maximum_iterations)break;
		// Right-hand side -residual on the owned free rows, zero on the constrained rows.
		{
			PetscScalar* values=nullptr;Check(VecGetArray(petsc.rhs,&values),"VecGetArray rhs");
			for(std::uint32_t local_vertex=0;local_vertex<partition.owned_node_count;++local_vertex)
				for(int component=0;component<3;++component){
					const std::size_t dof=3*partition.local_to_global_node[local_vertex]+component;
					values[row_of(local_vertex,component)-begin]=prescribed_displacement_m.count(dof)?0.0:-assembly.residual[dof];
				}
			Check(VecRestoreArray(petsc.rhs,&values),"VecRestoreArray rhs");
		}
		if(!petsc.solver)petsc.solver=MakeNativeTetPetscKsp(comm,prefix,petsc.matrix,defaults);
		else Check(KSPSetOperators(petsc.solver,petsc.matrix,petsc.matrix),"KSPSetOperators");
		Check(VecSet(petsc.update,0.0),"VecSet update");
		NativeTetPetscKspSolve(petsc.solver,petsc.rhs,petsc.update,prefix);
		Check(VecScatterBegin(petsc.scatter,petsc.update,petsc.gathered,INSERT_VALUES,SCATTER_FORWARD),"VecScatterBegin");
		Check(VecScatterEnd(petsc.scatter,petsc.update,petsc.gathered,INSERT_VALUES,SCATTER_FORWARD),"VecScatterEnd");
		{
			const PetscScalar* gathered=nullptr;Check(VecGetArrayRead(petsc.gathered,&gathered),"VecGetArrayRead gathered");
			for(std::size_t node=0;node<nodes;++node)for(int component=0;component<3;++component)
				update[3*node+component]=PetscRealPart(gathered[3*partition.OwnerMajorNode(static_cast<std::uint32_t>(node))+component]);
			Check(VecRestoreArrayRead(petsc.gathered,&gathered),"VecRestoreArrayRead gathered");
		}
		for(const double value:update)if(!std::isfinite(value))
			throw std::runtime_error("native tetrahedral solid static update is nonfinite");
		bool accepted=false;
		for(double scale=1.0;scale>=options.minimum_line_search;scale*=0.5){
			auto trial=displacement;for(const auto dof:free)trial[dof]+=scale*update[dof];
			auto candidate=assemble(trial,false);
			if(candidate.failed)continue;
			if(native_tet_solid_static_detail::FreeNorm(candidate.residual,free)<norm){
				displacement=std::move(trial);assembly=assemble(displacement,true);
				if(assembly.failed)throw std::runtime_error("native tetrahedral distributed solid: accepted state lost its elements");
				accepted=true;break;
			}
		}
		if(!accepted)throw std::runtime_error("native tetrahedral solid static line search failed");
	}
	if(!result.converged)throw std::runtime_error("native tetrahedral solid static Newton did not converge");
	result.displacement_m=displacement;result.reaction_n=assembly.residual;
	result.final_free_residual_n=native_tet_solid_static_detail::FreeNorm(assembly.residual,free);
	result.strain_energy_j=assembly.energy;result.minimum_deformation_jacobian=assembly.minimum_jacobian;
	for(std::size_t dof=0;dof<dofs;++dof)result.external_load_potential_j-=external_force_n[dof]*displacement[dof];
	return result;
}

} // namespace iga

#endif
