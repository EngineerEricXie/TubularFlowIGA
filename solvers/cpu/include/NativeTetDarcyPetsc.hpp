#ifndef IGA_NATIVE_TET_DARCY_PETSC_HPP
#define IGA_NATIVE_TET_DARCY_PETSC_HPP

#include "NativeTetFem.hpp"
#include "NativeTetRt0Flux.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

// q=-K grad(p), with hydraulic mobility K in m^2/(Pa s), q in m/s,
// and -div(K grad(p))=s with volumetric source s in 1/s.
struct NativeTetDarcyResult
{
	std::vector<double> pressure_pa;
	std::vector<std::array<double,3>> cell_flux_m_s;
	// Integrated outward flux on the four faces opposite each cell's local nodes.
	std::vector<std::array<double,4>> conservative_face_flow_m3_s;
	std::map<int,double> outward_boundary_flow_m3_s;
	std::map<int,double> conservative_outward_boundary_flow_m3_s;
	double maximum_cell_balance_defect_m3_s=0.;
	PetscInt flux_recovery_iterations=0;
	double volume_source_m3_s=0.;
	PetscInt linear_iterations=0;
	KSPConvergedReason converged_reason=KSP_CONVERGED_ITERATING;
	double residual_norm=0.;
};

// RT0 basis for the face opposite local vertex i is (x-x_i)/(3V).
// Each basis has unit integrated outward flow on its own face and zero
// normal flow on the other three faces. The recovered field is therefore
// H(div)-conforming across cells, though the pressure solve itself is P1.
inline std::array<double,3> EvaluateNativeTetDarcyRt0Flux(
	const NativeTetMesh& mesh,const NativeTetDarcyResult& result,
	std::size_t cell_index,const std::array<double,3>& position_m)
{
	return EvaluateNativeTetRt0Flux(mesh,result.conservative_face_flow_m3_s,
		cell_index,position_m);
}

namespace native_tet_darcy_detail {

inline void Check(PetscErrorCode error,const char* context)
{
	if(error)throw std::runtime_error(std::string("native Darcy PETSc failure: ")+context);
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

inline std::array<std::uint32_t,3> SortedFace(
	const std::array<std::uint32_t,3>& nodes)
{
	auto result=nodes;
	std::sort(result.begin(),result.end());
	return result;
}

} // namespace native_tet_darcy_detail

inline void CompleteNativeTetDarcyFlux(
	const NativeTetMesh& mesh,const std::vector<double>& mobility_by_cell_m2_pa_s,
	const std::vector<double>& source_by_cell_s_inv,
	const std::map<int,double>& pressure_by_boundary_label_pa,
	const std::map<int,double>& outward_flux_by_boundary_label_m_s,
	NativeTetDarcyResult& result)
{
	using namespace native_tet_darcy_detail;
	const bool compute_cell_flux=result.cell_flux_m_s.empty();
	for(std::size_t index=0;index<mesh.cells.size();++index){
		const auto& cell=mesh.cells[index];
		const auto geometry=EvaluateNativeTetGeometry(mesh,cell);
		if(compute_cell_flux){
			std::array<double,3> flux{};
			for(std::size_t node=0;node<4;++node)
				for(int axis=0;axis<3;++axis)
					flux[axis]-=mobility_by_cell_m2_pa_s[index]
						*result.pressure_pa[cell.nodes[node]]
						*geometry.barycentric_gradients[node][axis];
			result.cell_flux_m_s.push_back(flux);
		}
		result.volume_source_m3_s+=source_by_cell_s_inv[index]
			*geometry.determinant/6.;
	}
	// Face flux is a P1 cell-gradient diagnostic. For nonzero sources its
	// integral is not an exactly conservative recovered H(div) flux.
	std::map<std::array<std::uint32_t,3>,std::pair<std::size_t,std::uint32_t>> owners;
	for(std::size_t index=0;index<mesh.cells.size();++index)
		for(std::size_t opposite=0;opposite<4;++opposite){
			std::array<std::uint32_t,3> face{};
			std::size_t next=0;
			for(std::size_t local=0;local<4;++local)
				if(local!=opposite)face[next++]=mesh.cells[index].nodes[local];
			const auto key=SortedFace(face);
			const auto inserted=owners.emplace(key,std::make_pair(index,
				mesh.cells[index].nodes[opposite]));
			if(!inserted.second)owners.erase(inserted.first);
		}
	for(const auto& face:mesh.boundary_triangles){
		const auto found=owners.find(SortedFace(face.nodes));
		if(found==owners.end())
			throw std::invalid_argument("native Darcy labelled face is not on the boundary");
		const auto& a=mesh.points.at(face.nodes[0]);
		const auto& b=mesh.points.at(face.nodes[1]);
		const auto& c=mesh.points.at(face.nodes[2]);
		const auto& d=mesh.points.at(found->second.second);
		std::array<double,3> normal{{
			(b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1]),
			(b[2]-a[2])*(c[0]-a[0])-(b[0]-a[0])*(c[2]-a[2]),
			(b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])}};
		const double interior=normal[0]*(d[0]-a[0])+normal[1]*(d[1]-a[1])
			+normal[2]*(d[2]-a[2]);
		if(interior>0.)for(auto& value:normal)value=-value;
		const auto& flux=result.cell_flux_m_s[found->second.first];
		result.outward_boundary_flow_m3_s[face.boundary_label]
			+=0.5*(flux[0]*normal[0]+flux[1]*normal[1]+flux[2]*normal[2]);
	}
	// Recover a single integrated flux per geometric face. The least-squares
	// correction minimizes sum(delta_flux^2) while enforcing each cell's
	// integrated balance. Natural/Neumann boundary fluxes are held fixed;
	// only pressure boundaries may absorb the balancing correction.
	struct RecoveryFace
	{
		std::size_t owner=0,other=0,owner_local=0,other_local=0;
		double flow=0.;
		int label=-1;
		bool boundary=false,free=false;
	};
	std::vector<RecoveryFace> faces;
	std::map<std::array<std::uint32_t,3>,std::size_t> face_index;
	for(std::size_t cell_index=0;cell_index<mesh.cells.size();++cell_index)
		for(std::size_t opposite=0;opposite<4;++opposite){
			std::array<std::uint32_t,3> nodes_on_face{};
			std::size_t next=0;
			for(std::size_t local=0;local<4;++local)
				if(local!=opposite)
					nodes_on_face[next++]=mesh.cells[cell_index].nodes[local];
			const auto key=SortedFace(nodes_on_face);
			const auto inserted=face_index.emplace(key,faces.size());
			if(inserted.second)faces.push_back({cell_index,0,opposite,0,0.,-1,true,false});
			else{
				auto& face=faces[inserted.first->second];
				if(!face.boundary)
					throw std::invalid_argument("native Darcy nonmanifold tetrahedral face");
				face.other=cell_index;
				face.other_local=opposite;
				face.boundary=false;
				face.free=true;
			}
		}
	std::map<std::array<std::uint32_t,3>,int> boundary_labels;
	for(const auto& triangle:mesh.boundary_triangles){
		const auto key=SortedFace(triangle.nodes);
		if(!boundary_labels.emplace(key,triangle.boundary_label).second)
			throw std::invalid_argument("native Darcy duplicate boundary face");
	}
	std::vector<double> right(mesh.cells.size(),0.),diagonal(mesh.cells.size(),0.);
	for(std::size_t index=0;index<mesh.cells.size();++index)
		right[index]=source_by_cell_s_inv[index]
			*EvaluateNativeTetGeometry(mesh,mesh.cells[index]).determinant/6.;
	for(const auto& entry:face_index){
		auto& face=faces[entry.second];
		const auto& nodes_on_face=entry.first;
		const auto& a=mesh.points.at(nodes_on_face[0]);
		const auto& b=mesh.points.at(nodes_on_face[1]);
		const auto& c=mesh.points.at(nodes_on_face[2]);
		const auto& opposite=mesh.points.at(
			mesh.cells[face.owner].nodes[face.owner_local]);
		std::array<double,3> area_normal{{
			0.5*((b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1])),
			0.5*((b[2]-a[2])*(c[0]-a[0])-(b[0]-a[0])*(c[2]-a[2])),
			0.5*((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]))}};
		if(area_normal[0]*(opposite[0]-a[0])
			+area_normal[1]*(opposite[1]-a[1])
			+area_normal[2]*(opposite[2]-a[2])>0.)
			for(auto& component:area_normal)component=-component;
		const auto dot=[&](std::size_t cell){
			const auto& flux=result.cell_flux_m_s[cell];
			return flux[0]*area_normal[0]+flux[1]*area_normal[1]
				+flux[2]*area_normal[2];
		};
		if(face.boundary){
			const auto label=boundary_labels.find(nodes_on_face);
			if(label==boundary_labels.end())
				throw std::invalid_argument("native Darcy unlabelled boundary face");
			face.label=label->second;
			face.free=pressure_by_boundary_label_pa.count(face.label)!=0;
			if(face.free)face.flow=dot(face.owner);
			else{
				const auto specified=outward_flux_by_boundary_label_m_s.find(face.label);
				const double value=specified==outward_flux_by_boundary_label_m_s.end()
					?0.:specified->second;
				const double area=std::sqrt(area_normal[0]*area_normal[0]
					+area_normal[1]*area_normal[1]
					+area_normal[2]*area_normal[2]);
				face.flow=value*area;
			}
			right[face.owner]-=face.flow;
			if(face.free)diagonal[face.owner]+=1.;
		}else{
			face.flow=0.5*(dot(face.owner)+dot(face.other));
			right[face.owner]-=face.flow;
			right[face.other]+=face.flow;
			diagonal[face.owner]+=1.;
			diagonal[face.other]+=1.;
		}
	}
	if(boundary_labels.size()!=static_cast<std::size_t>(std::count_if(
		faces.begin(),faces.end(),
		[](const RecoveryFace& face){return face.boundary;})))
		throw std::invalid_argument("native Darcy boundary labels do not match tetra faces");
	const auto multiply=[&](const std::vector<double>& input,
		std::vector<double>& output){
		std::fill(output.begin(),output.end(),0.);
		for(const auto& face:faces){
			if(!face.free)continue;
			const double jump=face.boundary?input[face.owner]
				:input[face.owner]-input[face.other];
			output[face.owner]+=jump;
			if(!face.boundary)output[face.other]-=jump;
		}
	};
	std::vector<double> potential(mesh.cells.size(),0.),flux_residual=right,
		preconditioned(mesh.cells.size(),0.),direction(mesh.cells.size(),0.),
		product(mesh.cells.size(),0.);
	const auto inner=[](const std::vector<double>& left,
		const std::vector<double>& right_vector){
		double sum=0.;
		for(std::size_t index=0;index<left.size();++index)
			sum+=left[index]*right_vector[index];
		return sum;
	};
	const double initial_norm=std::sqrt(inner(right,right));
	const double target=1e-12*initial_norm;
	if(initial_norm>0.){
		for(std::size_t index=0;index<mesh.cells.size();++index){
			if(diagonal[index]<=0.)
				throw std::runtime_error("native Darcy flux recovery has isolated cell");
			preconditioned[index]=flux_residual[index]/diagonal[index];
			direction[index]=preconditioned[index];
		}
		double previous=inner(flux_residual,preconditioned);
		const std::size_t limit=std::max<std::size_t>(1000,2*mesh.cells.size());
		for(std::size_t iteration=0;iteration<limit;++iteration){
			multiply(direction,product);
			const double denominator=inner(direction,product);
			if(!(denominator>0.)||!std::isfinite(denominator))
				throw std::runtime_error("native Darcy flux recovery is singular");
			const double alpha=previous/denominator;
			for(std::size_t index=0;index<mesh.cells.size();++index){
				potential[index]+=alpha*direction[index];
				flux_residual[index]-=alpha*product[index];
			}
			result.flux_recovery_iterations=static_cast<PetscInt>(iteration+1);
			if(std::sqrt(inner(flux_residual,flux_residual))<=target)break;
			if(iteration+1==limit)
				throw std::runtime_error("native Darcy flux recovery did not converge");
			for(std::size_t index=0;index<mesh.cells.size();++index)
				preconditioned[index]=flux_residual[index]/diagonal[index];
			const double current=inner(flux_residual,preconditioned);
			const double beta=current/previous;
			for(std::size_t index=0;index<mesh.cells.size();++index)
				direction[index]=preconditioned[index]+beta*direction[index];
			previous=current;
		}
	}
	result.conservative_face_flow_m3_s.resize(mesh.cells.size());
	std::vector<double> balance(mesh.cells.size(),0.);
	for(auto& face:faces){
		if(face.free)face.flow+=face.boundary?potential[face.owner]
			:potential[face.owner]-potential[face.other];
		result.conservative_face_flow_m3_s[face.owner][face.owner_local]=face.flow;
		balance[face.owner]+=face.flow;
		if(face.boundary)
			result.conservative_outward_boundary_flow_m3_s[face.label]+=face.flow;
		else{
			result.conservative_face_flow_m3_s[face.other][face.other_local]=-face.flow;
			balance[face.other]-=face.flow;
		}
	}
	for(std::size_t index=0;index<mesh.cells.size();++index){
		const double source=source_by_cell_s_inv[index]
			*EvaluateNativeTetGeometry(mesh,mesh.cells[index]).determinant/6.;
		result.maximum_cell_balance_defect_m3_s=std::max(
			result.maximum_cell_balance_defect_m3_s,std::abs(balance[index]-source));
	}
}

inline NativeTetDarcyResult SolveNativeTetDarcyPetsc(
	const NativeTetMesh& mesh,const std::vector<double>& mobility_by_cell_m2_pa_s,
	const std::vector<double>& source_by_cell_s_inv,
	const std::map<int,double>& pressure_by_boundary_label_pa,
	const std::map<int,double>& outward_flux_by_boundary_label_m_s={})
{
	using namespace native_tet_darcy_detail;
	const std::size_t nodes=mesh.points.size();
	if(nodes==0||mesh.cells.empty()||nodes>static_cast<std::size_t>(
		std::numeric_limits<PetscInt>::max())
		||mobility_by_cell_m2_pa_s.size()!=mesh.cells.size()
		||source_by_cell_s_inv.size()!=mesh.cells.size()
		||pressure_by_boundary_label_pa.empty())
		throw std::invalid_argument("native Darcy mesh, material or pressure data is invalid");
	std::set<int> labels;
	for(const auto& face:mesh.boundary_triangles)labels.insert(face.boundary_label);
	for(const auto& item:pressure_by_boundary_label_pa)
		if(!labels.count(item.first)||!std::isfinite(item.second))
			throw std::invalid_argument("native Darcy pressure boundary is invalid");
	for(const auto& item:outward_flux_by_boundary_label_m_s)
		if(!labels.count(item.first)||pressure_by_boundary_label_pa.count(item.first)
			||!std::isfinite(item.second))
			throw std::invalid_argument("native Darcy flux boundary is invalid");
	std::map<std::uint32_t,double> fixed;
	for(const auto& face:mesh.boundary_triangles){
		const auto condition=pressure_by_boundary_label_pa.find(face.boundary_label);
		if(condition==pressure_by_boundary_label_pa.end())continue;
		for(const auto node:face.nodes){
			if(node>=nodes)throw std::invalid_argument("native Darcy boundary node is invalid");
			const auto inserted=fixed.emplace(node,condition->second);
			if(!inserted.second&&std::abs(inserted.first->second-condition->second)
				>1e-12*std::max({1.,std::abs(inserted.first->second),
					std::abs(condition->second)}))
				throw std::invalid_argument("native Darcy pressure labels conflict at a node");
		}
	}
	if(fixed.empty())throw std::invalid_argument("native Darcy has no constrained nodes");
	int rank=0,ranks=1;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	const PetscInt size=static_cast<PetscInt>(nodes);
	Objects petsc;
	Check(MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,size,size,
		16,nullptr,16,nullptr,&petsc.matrix),"MatCreateAIJ");
	// Irregular fTetWild vertex valence may exceed this initial estimate.
	Check(MatSetOption(petsc.matrix,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE),
		"MatSetOption dynamic preallocation");
	Check(VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,size,&petsc.rhs),"VecCreateMPI");
	Check(VecDuplicate(petsc.rhs,&petsc.solution),"VecDuplicate");
	for(std::size_t index=0;index<mesh.cells.size();++index){
		const auto& cell=mesh.cells[index];
		const auto geometry=EvaluateNativeTetGeometry(mesh,cell);
		const double mobility=mobility_by_cell_m2_pa_s[index];
		const double source=source_by_cell_s_inv[index];
		if(!(mobility>0.)||!std::isfinite(mobility)||!std::isfinite(source))
			throw std::invalid_argument("native Darcy cell material or source is invalid");
		if(static_cast<int>(index%static_cast<std::size_t>(ranks))!=rank)continue;
		const double volume=geometry.determinant/6.;
		for(std::size_t row=0;row<4;++row){
			const auto global_row=cell.nodes[row];
			if(global_row>=nodes)throw std::invalid_argument("native Darcy cell node is invalid");
			if(fixed.count(global_row))continue;
			double right=source*volume/4.;
			for(std::size_t column=0;column<4;++column){
				const auto global_column=cell.nodes[column];
				if(global_column>=nodes)
					throw std::invalid_argument("native Darcy cell node is invalid");
				double dot=0.;
				for(int axis=0;axis<3;++axis)
					dot+=geometry.barycentric_gradients[row][axis]
						*geometry.barycentric_gradients[column][axis];
				const double value=mobility*volume*dot;
				const auto known=fixed.find(global_column);
				if(known==fixed.end())
					Check(MatSetValue(petsc.matrix,static_cast<PetscInt>(global_row),
						static_cast<PetscInt>(global_column),value,ADD_VALUES),"MatSetValue");
				else right-=value*known->second;
			}
			Check(VecSetValue(petsc.rhs,static_cast<PetscInt>(global_row),right,
				ADD_VALUES),"VecSetValue source");
		}
	}
	for(std::size_t index=0;index<mesh.boundary_triangles.size();++index){
		if(static_cast<int>(index%static_cast<std::size_t>(ranks))!=rank)continue;
		const auto& face=mesh.boundary_triangles[index];
		const auto condition=outward_flux_by_boundary_label_m_s.find(face.boundary_label);
		if(condition==outward_flux_by_boundary_label_m_s.end())continue;
		const auto& a=mesh.points.at(face.nodes[0]);
		const auto& b=mesh.points.at(face.nodes[1]);
		const auto& c=mesh.points.at(face.nodes[2]);
		const std::array<double,3> ab{{b[0]-a[0],b[1]-a[1],b[2]-a[2]}};
		const std::array<double,3> ac{{c[0]-a[0],c[1]-a[1],c[2]-a[2]}};
		const std::array<double,3> cross{{ab[1]*ac[2]-ab[2]*ac[1],
			ab[2]*ac[0]-ab[0]*ac[2],ab[0]*ac[1]-ab[1]*ac[0]}};
		const double area=0.5*std::sqrt(cross[0]*cross[0]
			+cross[1]*cross[1]+cross[2]*cross[2]);
		if(!(area>0.)||!std::isfinite(area))
			throw std::invalid_argument("native Darcy boundary face area is invalid");
		for(const auto node:face.nodes)
			if(!fixed.count(node))
				Check(VecSetValue(petsc.rhs,static_cast<PetscInt>(node),
					-condition->second*area/3.,ADD_VALUES),"VecSetValue flux");
	}
	PetscInt owned_begin=0,owned_end=0;
	Check(VecGetOwnershipRange(petsc.rhs,&owned_begin,&owned_end),"VecGetOwnershipRange");
	for(const auto& item:fixed)
		if(static_cast<PetscInt>(item.first)>=owned_begin
			&&static_cast<PetscInt>(item.first)<owned_end){
			Check(MatSetValue(petsc.matrix,static_cast<PetscInt>(item.first),
				static_cast<PetscInt>(item.first),1.,ADD_VALUES),"MatSetValue pressure");
			Check(VecSetValue(petsc.rhs,static_cast<PetscInt>(item.first),
				item.second,ADD_VALUES),"VecSetValue pressure");
		}
	Check(MatAssemblyBegin(petsc.matrix,MAT_FINAL_ASSEMBLY),"MatAssemblyBegin");
	Check(MatAssemblyEnd(petsc.matrix,MAT_FINAL_ASSEMBLY),"MatAssemblyEnd");
	Check(VecAssemblyBegin(petsc.rhs),"VecAssemblyBegin");
	Check(VecAssemblyEnd(petsc.rhs),"VecAssemblyEnd");
	Check(KSPCreate(PETSC_COMM_WORLD,&petsc.solver),"KSPCreate");
	Check(KSPSetOperators(petsc.solver,petsc.matrix,petsc.matrix),"KSPSetOperators");
	Check(KSPSetType(petsc.solver,KSPCG),"KSPSetType");
	PC preconditioner=nullptr;
	Check(KSPGetPC(petsc.solver,&preconditioner),"KSPGetPC");
	Check(PCSetType(preconditioner,PCJACOBI),"PCSetType");
	Check(KSPSetTolerances(petsc.solver,1e-12,1e-14,PETSC_DEFAULT,1000),
		"KSPSetTolerances");
	Check(KSPSetFromOptions(petsc.solver),"KSPSetFromOptions");
	Check(KSPSolve(petsc.solver,petsc.rhs,petsc.solution),"KSPSolve");
	NativeTetDarcyResult result;
	Check(KSPGetConvergedReason(petsc.solver,&result.converged_reason),
		"KSPGetConvergedReason");
	Check(KSPGetIterationNumber(petsc.solver,&result.linear_iterations),
		"KSPGetIterationNumber");
	PetscReal residual=0.;
	Check(KSPGetResidualNorm(petsc.solver,&residual),"KSPGetResidualNorm");
	result.residual_norm=static_cast<double>(residual);
	if(result.converged_reason<=0||!std::isfinite(result.residual_norm))
		throw std::runtime_error("native Darcy PETSc solve did not converge");
	Check(VecScatterCreateToAll(petsc.solution,&petsc.scatter,&petsc.gathered),
		"VecScatterCreateToAll");
	Check(VecScatterBegin(petsc.scatter,petsc.solution,petsc.gathered,
		INSERT_VALUES,SCATTER_FORWARD),"VecScatterBegin");
	Check(VecScatterEnd(petsc.scatter,petsc.solution,petsc.gathered,
		INSERT_VALUES,SCATTER_FORWARD),"VecScatterEnd");
	const PetscScalar* values=nullptr;
	Check(VecGetArrayRead(petsc.gathered,&values),"VecGetArrayRead");
	result.pressure_pa.resize(nodes);
	for(std::size_t node=0;node<nodes;++node)
		result.pressure_pa[node]=PetscRealPart(values[node]);
	Check(VecRestoreArrayRead(petsc.gathered,&values),"VecRestoreArrayRead");
	CompleteNativeTetDarcyFlux(mesh,mobility_by_cell_m2_pa_s,
		source_by_cell_s_inv,pressure_by_boundary_label_pa,
		outward_flux_by_boundary_label_m_s,result);
	return result;
}

} // namespace iga

#endif
