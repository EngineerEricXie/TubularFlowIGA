#ifndef IGA_NATIVE_TET_MIXED_SOLID_STATIC_SOLVER_HPP
#define IGA_NATIVE_TET_MIXED_SOLID_STATIC_SOLVER_HPP

#include "NativeTetMixedHyperelasticSolid.hpp"
#include "NativeTetSolidStaticSolver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

namespace iga {

struct NativeTetMixedSolidStaticResult
{
	std::vector<double> displacement_m,pressure_pa,reaction_n;
	int iterations=0;
	double initial_free_residual=0.0,final_free_residual=0.0;
	double minimum_deformation_jacobian=0.0,maximum_absolute_volume_change=0.0;
	bool converged=false;
};

namespace native_tet_mixed_solid_static_detail {

struct Assembly
{
	std::vector<double> residual,tangent;
	double minimum_jacobian=std::numeric_limits<double>::infinity(),maximum_volume_change=0.0;
};

inline Assembly Assemble(const NativeTetMesh& mesh,const std::vector<double>& state,
	const std::vector<double>& external_force,const NativeTetMixedSolidMaterial& material)
{
	const std::size_t nodes=mesh.points.size(),dofs=4*nodes;Assembly result;
	result.residual.assign(dofs,0.);result.tangent.assign(dofs*dofs,0.);
	for(const auto& cell:mesh.cells){
		std::array<double,12> displacement{};std::array<double,4> pressure{};
		std::array<std::size_t,16> rows{};
		for(std::size_t a=0;a<4;++a){
			for(int component=0;component<3;++component){rows[3*a+component]=3*cell.nodes[a]+component;
				displacement[3*a+component]=state[rows[3*a+component]];}
			rows[12+a]=3*nodes+cell.nodes[a];pressure[a]=state[rows[12+a]];
		}
		const auto element=BuildNativeTetMixedHyperelasticSolidElement(mesh,cell,displacement,pressure,material);
		result.minimum_jacobian=std::min(result.minimum_jacobian,element.deformation_jacobian);
		result.maximum_volume_change=std::max(result.maximum_volume_change,
			std::abs(element.deformation_jacobian-1.));
		for(std::size_t i=0;i<16;++i){result.residual[rows[i]]+=element.residual[i];
			for(std::size_t j=0;j<16;++j)result.tangent[rows[i]*dofs+rows[j]]+=element.tangent[i*16+j];}
	}
	for(std::size_t dof=0;dof<3*nodes;++dof)result.residual[dof]-=external_force[dof];
	return result;
}

inline double Norm(const std::vector<double>& residual,const std::vector<std::size_t>& free)
{
	long double sum=0.;for(auto dof:free)sum+=static_cast<long double>(residual[dof])*residual[dof];
	return std::sqrt(static_cast<double>(sum));
}

} // namespace native_tet_mixed_solid_static_detail

inline NativeTetMixedSolidStaticResult SolveNativeTetMixedSolidStatic(
	const NativeTetMesh& mesh,const NativeTetMixedSolidMaterial& material,
	const std::vector<double>& external_force_n,
	const std::map<std::size_t,double>& prescribed_displacement_m,
	std::vector<double> initial_displacement_m={},std::vector<double> initial_pressure_pa={},
	const NativeTetSolidStaticOptions& options={})
{
	using namespace native_tet_mixed_solid_static_detail;
	ValidateNativeTetMixedSolidMaterial(material);
	native_tet_solid_static_detail::ValidateOptions(options);
	const std::size_t nodes=mesh.points.size(),displacement_dofs=3*nodes,dofs=4*nodes;
	if(nodes==0||mesh.cells.empty()||dofs>options.maximum_dofs
		||external_force_n.size()!=displacement_dofs)
		throw std::invalid_argument("native mixed tetrahedral solid static input shape is invalid");
	if(initial_displacement_m.empty())initial_displacement_m.assign(displacement_dofs,0.);
	if(initial_pressure_pa.empty())initial_pressure_pa.assign(nodes,0.);
	if(initial_displacement_m.size()!=displacement_dofs||initial_pressure_pa.size()!=nodes)
		throw std::invalid_argument("native mixed tetrahedral solid initial state size mismatch");
	for(double value:external_force_n)if(!std::isfinite(value))
		throw std::invalid_argument("native mixed tetrahedral solid load is nonfinite");
	for(double value:initial_displacement_m)if(!std::isfinite(value))
		throw std::invalid_argument("native mixed tetrahedral solid displacement is nonfinite");
	for(double value:initial_pressure_pa)if(!std::isfinite(value))
		throw std::invalid_argument("native mixed tetrahedral solid pressure is nonfinite");
	for(const auto& constraint:prescribed_displacement_m)
		if(constraint.first>=displacement_dofs||!std::isfinite(constraint.second))
			throw std::invalid_argument("native mixed tetrahedral solid constraint is invalid");
	std::vector<double> state(dofs);std::copy(initial_displacement_m.begin(),initial_displacement_m.end(),state.begin());
	std::copy(initial_pressure_pa.begin(),initial_pressure_pa.end(),state.begin()+displacement_dofs);
	for(const auto& constraint:prescribed_displacement_m)state[constraint.first]=constraint.second;
	std::vector<std::size_t> free;free.reserve(dofs-prescribed_displacement_m.size());
	for(std::size_t dof=0;dof<dofs;++dof)if(!prescribed_displacement_m.count(dof))free.push_back(dof);
	auto assembly=Assemble(mesh,state,external_force_n,material);
	NativeTetMixedSolidStaticResult result;result.initial_free_residual=Norm(assembly.residual,free);
	const double threshold=options.absolute_tolerance_n+options.relative_tolerance*result.initial_free_residual;
	for(int iteration=0;iteration<=options.maximum_iterations;++iteration){
		const double norm=Norm(assembly.residual,free);
		if(norm<=threshold){result.iterations=iteration;result.converged=true;break;}
		if(iteration==options.maximum_iterations)break;
		const std::size_t n=free.size();std::vector<double> reduced(n*n),rhs(n);
		for(std::size_t i=0;i<n;++i){rhs[i]=-assembly.residual[free[i]];
			for(std::size_t j=0;j<n;++j)reduced[i*n+j]=assembly.tangent[free[i]*dofs+free[j]];}
		const auto update=native_tet_solid_static_detail::SolveDense(std::move(reduced),std::move(rhs));
		bool accepted=false;
		for(double scale=1.;scale>=options.minimum_line_search;scale*=.5){
			auto trial=state;for(std::size_t i=0;i<n;++i)trial[free[i]]+=scale*update[i];
			try{auto candidate=Assemble(mesh,trial,external_force_n,material);
				if(Norm(candidate.residual,free)<norm){state=std::move(trial);assembly=std::move(candidate);
					accepted=true;break;}}
			catch(const std::runtime_error&){}
		}
		if(!accepted)throw std::runtime_error("native mixed tetrahedral solid static line search failed");
	}
	if(!result.converged)throw std::runtime_error("native mixed tetrahedral solid static Newton did not converge");
	result.displacement_m.assign(state.begin(),state.begin()+displacement_dofs);
	result.pressure_pa.assign(state.begin()+displacement_dofs,state.end());
	result.reaction_n.assign(assembly.residual.begin(),assembly.residual.begin()+displacement_dofs);
	result.final_free_residual=Norm(assembly.residual,free);
	result.minimum_deformation_jacobian=assembly.minimum_jacobian;
	result.maximum_absolute_volume_change=assembly.maximum_volume_change;return result;
}

} // namespace iga
#endif
