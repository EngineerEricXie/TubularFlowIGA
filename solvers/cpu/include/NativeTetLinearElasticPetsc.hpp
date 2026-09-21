#ifndef IGA_NATIVE_TET_LINEAR_ELASTIC_PETSC_HPP
#define IGA_NATIVE_TET_LINEAR_ELASTIC_PETSC_HPP

#include "NativeTetFem.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

struct NativeTetLinearElasticResult
{
	std::vector<std::array<double,3>> displacement_m;
	KSPConvergedReason converged_reason=KSP_CONVERGED_ITERATING;
	PetscInt linear_iterations=0;
	double residual_norm_n=0.;
	double maximum_displacement_m=0.;
	double minimum_deformation_jacobian=std::numeric_limits<double>::infinity();
};

namespace native_tet_linear_elastic_detail {

inline void Check(PetscErrorCode code,const char* context)
{
	if(code)throw std::runtime_error(std::string("native linear elastic PETSc failure: ")+context);
}

struct Objects
{
	Mat matrix=nullptr;
	Vec rhs=nullptr,solution=nullptr,gathered=nullptr;
	VecScatter scatter=nullptr;
	KSP solver=nullptr;
	~Objects()
	{
		if(scatter)VecScatterDestroy(&scatter);
		if(gathered)VecDestroy(&gathered);
		if(solver)KSPDestroy(&solver);
		if(solution)VecDestroy(&solution);
		if(rhs)VecDestroy(&rhs);
		if(matrix)MatDestroy(&matrix);
	}
};

} // namespace native_tet_linear_elastic_detail

// Small-strain, isotropic, total-reference P1 weak form. PETSc supplies only
// sparse algebra and MPI; the element gradients and stiffness are owned here.
inline NativeTetLinearElasticResult SolveNativeTetLinearElasticPetsc(
	const NativeTetMesh& reference_mesh,double young_pa,double poisson_ratio,
	const std::vector<double>& external_nodal_force_n,
	const std::map<std::uint32_t,std::array<double,3>>& prescribed_displacement_m)
{
	using namespace native_tet_linear_elastic_detail;
	const std::size_t nodes=reference_mesh.points.size(),dofs=3*nodes;
	if(nodes==0||reference_mesh.cells.empty()||dofs>
		static_cast<std::size_t>(std::numeric_limits<PetscInt>::max())||
		external_nodal_force_n.size()!=dofs||prescribed_displacement_m.empty()||
		!(young_pa>0.)||!std::isfinite(young_pa)||!(poisson_ratio>0.)||
		!(poisson_ratio<.45)||!std::isfinite(poisson_ratio))
		throw std::invalid_argument("native linear elastic mesh/material/load invalid");
	for(double force:external_nodal_force_n)
		if(!std::isfinite(force))throw std::invalid_argument("nonfinite solid load");
	std::map<PetscInt,double> fixed;
	for(const auto& item:prescribed_displacement_m){
		if(item.first>=nodes)throw std::invalid_argument("invalid solid support node");
		for(int axis=0;axis<3;++axis){
			if(!std::isfinite(item.second[axis]))
				throw std::invalid_argument("nonfinite solid support value");
			fixed.emplace(static_cast<PetscInt>(3*item.first+axis),item.second[axis]);
		}
	}
	if(fixed.size()==dofs)throw std::invalid_argument("all solid DOFs are fixed");
	const double mu=young_pa/(2.*(1.+poisson_ratio));
	const double lambda=young_pa*poisson_ratio/
		((1.+poisson_ratio)*(1.-2.*poisson_ratio));
	Objects petsc;
	Check(MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,
		static_cast<PetscInt>(dofs),static_cast<PetscInt>(dofs),
		60,nullptr,60,nullptr,&petsc.matrix),"MatCreateAIJ");
	Check(MatSetOption(petsc.matrix,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE),
		"dynamic preallocation");
	Check(MatCreateVecs(petsc.matrix,&petsc.solution,&petsc.rhs),"MatCreateVecs");
	int rank=0,ranks=1;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	for(std::size_t cell_index=static_cast<std::size_t>(rank);
		cell_index<reference_mesh.cells.size();cell_index+=static_cast<std::size_t>(ranks)){
		const auto& cell=reference_mesh.cells[cell_index];
		const auto geometry=EvaluateNativeTetGeometry(reference_mesh,cell);
		const double volume=geometry.determinant/6.;
		for(std::size_t a=0;a<4;++a)for(int i=0;i<3;++i){
			const PetscInt row=static_cast<PetscInt>(3*cell.nodes[a]+i);
			if(fixed.count(row))continue;
			for(std::size_t b=0;b<4;++b)for(int k=0;k<3;++k){
				const PetscInt col=static_cast<PetscInt>(3*cell.nodes[b]+k);
				double gradients_dot=0.;
				for(int axis=0;axis<3;++axis)
					gradients_dot+=geometry.barycentric_gradients[a][axis]
						*geometry.barycentric_gradients[b][axis];
				const double stiffness=volume*(
					lambda*geometry.barycentric_gradients[a][i]
						*geometry.barycentric_gradients[b][k]
					+mu*(i==k?gradients_dot:0.)
					+mu*geometry.barycentric_gradients[a][k]
						*geometry.barycentric_gradients[b][i]);
				const auto constrained=fixed.find(col);
				if(constrained!=fixed.end())
					Check(VecSetValue(petsc.rhs,row,-stiffness*constrained->second,
						ADD_VALUES),"VecSetValue constrained contribution");
				else Check(MatSetValue(petsc.matrix,row,col,stiffness,ADD_VALUES),
					"MatSetValue element stiffness");
			}
		}
	}
	PetscInt begin=0,end=0;
	Check(VecGetOwnershipRange(petsc.rhs,&begin,&end),"VecGetOwnershipRange");
	for(PetscInt dof=begin;dof<end;++dof){
		const auto constrained=fixed.find(dof);
		if(constrained!=fixed.end()){
			Check(MatSetValue(petsc.matrix,dof,dof,1.,ADD_VALUES),
				"MatSetValue prescribed diagonal");
			Check(VecSetValue(petsc.rhs,dof,constrained->second,ADD_VALUES),
				"VecSetValue prescribed displacement");
		}else Check(VecSetValue(petsc.rhs,dof,
			external_nodal_force_n[static_cast<std::size_t>(dof)],ADD_VALUES),
			"VecSetValue external load");
	}
	Check(MatAssemblyBegin(petsc.matrix,MAT_FINAL_ASSEMBLY),"MatAssemblyBegin");
	Check(MatAssemblyEnd(petsc.matrix,MAT_FINAL_ASSEMBLY),"MatAssemblyEnd");
	Check(VecAssemblyBegin(petsc.rhs),"VecAssemblyBegin");
	Check(VecAssemblyEnd(petsc.rhs),"VecAssemblyEnd");
	Check(KSPCreate(PETSC_COMM_WORLD,&petsc.solver),"KSPCreate");
	Check(KSPSetOperators(petsc.solver,petsc.matrix,petsc.matrix),"KSPSetOperators");
	Check(KSPSetType(petsc.solver,KSPCG),"KSPSetType");
	PC pc=nullptr;
	Check(KSPGetPC(petsc.solver,&pc),"KSPGetPC");
	Check(PCSetType(pc,PCJACOBI),"PCSetType");
	Check(KSPSetTolerances(petsc.solver,1e-10,1e-12,PETSC_DEFAULT,4000),
		"KSPSetTolerances");
	Check(KSPSetFromOptions(petsc.solver),"KSPSetFromOptions");
	Check(KSPSolve(petsc.solver,petsc.rhs,petsc.solution),"KSPSolve");
	NativeTetLinearElasticResult result;
	Check(KSPGetConvergedReason(petsc.solver,&result.converged_reason),
		"KSPGetConvergedReason");
	Check(KSPGetIterationNumber(petsc.solver,&result.linear_iterations),
		"KSPGetIterationNumber");
	PetscReal residual=0.;
	Check(KSPGetResidualNorm(petsc.solver,&residual),"KSPGetResidualNorm");
	result.residual_norm_n=static_cast<double>(residual);
	if(result.converged_reason<=0||!std::isfinite(result.residual_norm_n))
		throw std::runtime_error("native linear elastic solve did not converge");
	Check(VecScatterCreateToAll(petsc.solution,&petsc.scatter,&petsc.gathered),
		"VecScatterCreateToAll");
	Check(VecScatterBegin(petsc.scatter,petsc.solution,petsc.gathered,
		INSERT_VALUES,SCATTER_FORWARD),"VecScatterBegin");
	Check(VecScatterEnd(petsc.scatter,petsc.solution,petsc.gathered,
		INSERT_VALUES,SCATTER_FORWARD),"VecScatterEnd");
	const PetscScalar* values=nullptr;
	Check(VecGetArrayRead(petsc.gathered,&values),"VecGetArrayRead");
	result.displacement_m.resize(nodes);
	for(std::size_t node=0;node<nodes;++node)
		for(int axis=0;axis<3;++axis){
			const double value=PetscRealPart(values[3*node+axis]);
			if(!std::isfinite(value))throw std::runtime_error("nonfinite solid displacement");
			result.displacement_m[node][axis]=value;
		}
	for(const auto& displacement:result.displacement_m)
		result.maximum_displacement_m=std::max(result.maximum_displacement_m,
			std::sqrt(displacement[0]*displacement[0]+
				displacement[1]*displacement[1]+
				displacement[2]*displacement[2]));
	Check(VecRestoreArrayRead(petsc.gathered,&values),"VecRestoreArrayRead");
	auto current=reference_mesh;
	for(std::size_t node=0;node<nodes;++node)
		for(int axis=0;axis<3;++axis)
			current.points[node][axis]+=result.displacement_m[node][axis];
	for(const auto& cell:current.cells){
		const auto geometry=EvaluateNativeTetGeometry(current,cell);
		const auto reference=EvaluateNativeTetGeometry(reference_mesh,cell);
		result.minimum_deformation_jacobian=std::min(
			result.minimum_deformation_jacobian,
			geometry.determinant/reference.determinant);
	}
	if(!(result.minimum_deformation_jacobian>0.))
		throw std::runtime_error("solid displacement inverted a tetrahedron");
	return result;
}

} // namespace iga

#endif
