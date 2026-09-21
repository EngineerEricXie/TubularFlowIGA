#ifndef IGA_NATIVE_TET_MOVING_SPECIES_PETSC_RUNTIME_HPP
#define IGA_NATIVE_TET_MOVING_SPECIES_PETSC_RUNTIME_HPP

#include "NativeTetMovingSpeciesTransport.hpp"

#include <petscksp.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

struct NativeTetMovingSpeciesPetscResult
{
	NativeTetMovingSpeciesStep step;
	PetscInt linear_iterations=0;
	KSPConvergedReason converged_reason=KSP_CONVERGED_ITERATING;
	double residual_norm=0.;
};

namespace native_tet_moving_species_petsc_detail {

inline void Check(PetscErrorCode error,const char* context)
{
	if(error)throw std::runtime_error(std::string("native moving species PETSc failure: ")
		+context);
}

struct Objects
{
	Mat matrix=nullptr;
	Vec right=nullptr;
	Vec solution=nullptr;
	Vec gathered=nullptr;
	VecScatter scatter=nullptr;
	KSP solver=nullptr;
	~Objects()
	{
		if(scatter)VecScatterDestroy(&scatter);
		if(gathered)VecDestroy(&gathered);
		if(solver)KSPDestroy(&solver);
		if(solution)VecDestroy(&solution);
		if(right)VecDestroy(&right);
		if(matrix)MatDestroy(&matrix);
	}
};

} // namespace native_tet_moving_species_petsc_detail

// PETSc supplies only sparse algebra and MPI. Element/face integration and
// rank ownership are shared with the native dense reference implementation.
inline NativeTetMovingSpeciesPetscResult SolveNativeTetMovingSpeciesPetscStep(
	const NativeTetMesh& previous_mesh,const NativeTetMesh& current_mesh,
	const std::vector<std::array<double,3>>& fluid_velocity_nodes_m_s,
	const std::vector<std::array<double,3>>& mesh_velocity_nodes_m_s,
	const std::vector<double>& previous_concentration_mol_m3,
	const std::map<int,double>& inflow_concentration_mol_m3,
	double diffusivity_m2_s,double source_mol_m3_s,double dt_s,
	bool monotone=false,double first_order_decay_rate_s_inv=0.,
	const std::map<int,NativeTetWallExchange>& wall_exchange={},
	const std::vector<std::array<double,4>>& rt0_face_flow_m3_s={},
	const std::vector<double>& cell_source_mol_m3_s={})
{
	using namespace native_tet_moving_species_petsc_detail;
	const std::size_t nodes=current_mesh.points.size();
	if(nodes==0||nodes>static_cast<std::size_t>(std::numeric_limits<PetscInt>::max()))
		throw std::invalid_argument("native moving species PETSc node count is invalid");
	int rank=0,ranks=1;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	auto local=AssembleNativeTetMovingSpeciesStep(previous_mesh,current_mesh,
		fluid_velocity_nodes_m_s,mesh_velocity_nodes_m_s,
		previous_concentration_mol_m3,inflow_concentration_mol_m3,
		diffusivity_m2_s,source_mol_m3_s,dt_s,rank,ranks,monotone,
		first_order_decay_rate_s_inv,wall_exchange,rt0_face_flow_m3_s,
		cell_source_mol_m3_s);
	const PetscInt size=static_cast<PetscInt>(nodes);
	Objects petsc;
	Check(MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,size,size,
		32,nullptr,32,nullptr,&petsc.matrix),"MatCreateAIJ");
	// Unstructured tetrahedra can have vertex valence above the 32-entry
	// estimate. Permit PETSc to grow sparse storage without changing the
	// project-owned element assembly or silently dropping a matrix entry.
	Check(MatSetOption(petsc.matrix,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE),
		"MatSetOption dynamic sparse preallocation");
	Check(VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,size,&petsc.right),"VecCreateMPI");
	Check(VecDuplicate(petsc.right,&petsc.solution),"VecDuplicate");
	for(const auto& entry:local.matrix)
		Check(MatSetValue(petsc.matrix,static_cast<PetscInt>(entry.first.first),
			static_cast<PetscInt>(entry.first.second),entry.second,ADD_VALUES),
			"MatSetValue");
	for(std::size_t node=0;node<nodes;++node)
		if(local.rhs[node]!=0.)
			Check(VecSetValue(petsc.right,static_cast<PetscInt>(node),
				local.rhs[node],ADD_VALUES),"VecSetValue");
	Check(MatAssemblyBegin(petsc.matrix,MAT_FINAL_ASSEMBLY),"MatAssemblyBegin");
	Check(MatAssemblyEnd(petsc.matrix,MAT_FINAL_ASSEMBLY),"MatAssemblyEnd");
	Check(VecAssemblyBegin(petsc.right),"VecAssemblyBegin");
	Check(VecAssemblyEnd(petsc.right),"VecAssemblyEnd");
	Check(KSPCreate(PETSC_COMM_WORLD,&petsc.solver),"KSPCreate");
	Check(KSPSetOperators(petsc.solver,petsc.matrix,petsc.matrix),"KSPSetOperators");
	Check(KSPSetType(petsc.solver,KSPGMRES),"KSPSetType");
	PC preconditioner=nullptr;
	Check(KSPGetPC(petsc.solver,&preconditioner),"KSPGetPC");
	Check(PCSetType(preconditioner,PCBJACOBI),"PCSetType");
	Check(KSPSetTolerances(petsc.solver,1e-12,1e-14,PETSC_DEFAULT,500),
		"KSPSetTolerances");
	Check(KSPSetFromOptions(petsc.solver),"KSPSetFromOptions");
	Check(KSPSolve(petsc.solver,petsc.right,petsc.solution),"KSPSolve");
	NativeTetMovingSpeciesPetscResult result;
	Check(KSPGetConvergedReason(petsc.solver,&result.converged_reason),
		"KSPGetConvergedReason");
	Check(KSPGetIterationNumber(petsc.solver,&result.linear_iterations),
		"KSPGetIterationNumber");
	PetscReal residual=0.;
	Check(KSPGetResidualNorm(petsc.solver,&residual),"KSPGetResidualNorm");
	result.residual_norm=static_cast<double>(residual);
	if(result.converged_reason<=0||!std::isfinite(result.residual_norm))
		throw std::runtime_error("native moving species PETSc solve did not converge");
	Check(VecScatterCreateToAll(petsc.solution,&petsc.scatter,&petsc.gathered),
		"VecScatterCreateToAll");
	Check(VecScatterBegin(petsc.scatter,petsc.solution,petsc.gathered,
		INSERT_VALUES,SCATTER_FORWARD),"VecScatterBegin");
	Check(VecScatterEnd(petsc.scatter,petsc.solution,petsc.gathered,
		INSERT_VALUES,SCATTER_FORWARD),"VecScatterEnd");
	const PetscScalar* values=nullptr;
	Check(VecGetArrayRead(petsc.gathered,&values),"VecGetArrayRead");
	std::vector<double> concentration(nodes);
	for(std::size_t node=0;node<nodes;++node)
		concentration[node]=PetscRealPart(values[node]);
	Check(VecRestoreArrayRead(petsc.gathered,&values),"VecRestoreArrayRead");
	double global_source=0.,global_incoming=0.;
	if(MPI_Allreduce(&local.source_mol_s,&global_source,1,MPI_DOUBLE,MPI_SUM,
		PETSC_COMM_WORLD)!=MPI_SUCCESS
		||MPI_Allreduce(&local.incoming_flux_mol_s,&global_incoming,1,
			MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD)!=MPI_SUCCESS)
		throw std::runtime_error("native moving species MPI scalar reduction failed");
	std::vector<double> global_outward(nodes,0.);
	if(MPI_Allreduce(local.outward_coefficients.data(),global_outward.data(),size,
		MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD)!=MPI_SUCCESS)
		throw std::runtime_error("native moving species MPI flux reduction failed");
	local.source_mol_s=global_source;
	local.incoming_flux_mol_s=global_incoming;
	local.outward_coefficients=std::move(global_outward);
	for(auto& item:local.outward_coefficients_by_label){
		std::vector<double> global_labelled(nodes,0.);
		if(MPI_Allreduce(item.second.data(),global_labelled.data(),size,
			MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD)!=MPI_SUCCESS)
			throw std::runtime_error("native moving species MPI labelled flux reduction failed");
		item.second=std::move(global_labelled);
		double incoming=0.;
		if(MPI_Allreduce(&local.incoming_flux_by_label_mol_s.at(item.first),
			&incoming,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD)!=MPI_SUCCESS)
			throw std::runtime_error("native moving species MPI labelled inflow reduction failed");
		local.incoming_flux_by_label_mol_s.at(item.first)=incoming;
		double flow=0.;
		if(MPI_Allreduce(&local.outward_relative_flow_by_label_m3_s.at(item.first),
			&flow,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD)!=MPI_SUCCESS)
			throw std::runtime_error("native moving species MPI labelled flow reduction failed");
		local.outward_relative_flow_by_label_m3_s.at(item.first)=flow;
		double positive=0.,negative=0.;
		if(MPI_Allreduce(&local.positive_relative_flow_by_label_m3_s.at(item.first),
			&positive,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD)!=MPI_SUCCESS
			||MPI_Allreduce(&local.negative_relative_flow_by_label_m3_s.at(item.first),
				&negative,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD)!=MPI_SUCCESS)
			throw std::runtime_error("native moving species MPI split flow reduction failed");
		local.positive_relative_flow_by_label_m3_s.at(item.first)=positive;
		local.negative_relative_flow_by_label_m3_s.at(item.first)=negative;
	}
	for(auto& item:local.wall_exchange_coefficients_by_label){
		std::vector<double> global_labelled(nodes,0.);
		if(MPI_Allreduce(item.second.data(),global_labelled.data(),size,
			MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD)!=MPI_SUCCESS)
			throw std::runtime_error("native moving species MPI wall exchange reduction failed");
		item.second=std::move(global_labelled);
		double constant=0.;
		if(MPI_Allreduce(&local.wall_exchange_constant_by_label_mol_s.at(item.first),
			&constant,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD)!=MPI_SUCCESS)
			throw std::runtime_error("native moving species MPI wall constant reduction failed");
		local.wall_exchange_constant_by_label_mol_s.at(item.first)=constant;
	}
	result.step=FinishNativeTetMovingSpeciesStep(previous_mesh,current_mesh,
		previous_concentration_mol_m3,std::move(concentration),local,dt_s,
		monotone);
	return result;
}

} // namespace iga

#endif
