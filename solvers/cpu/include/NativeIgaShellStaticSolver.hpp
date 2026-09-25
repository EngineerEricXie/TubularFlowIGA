#ifndef IGA_NATIVE_IGA_SHELL_STATIC_SOLVER_HPP
#define IGA_NATIVE_IGA_SHELL_STATIC_SOLVER_HPP

#include "NativeIgaShellPatchSystem.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <vector>

namespace iga {

inline std::vector<double> NativeIgaShellPatchDeadTraction(
	const NativeNurbsShellPatch& reference,const std::array<double,3>& traction_pa)
{
	ValidateNativeNurbsShellPatch(reference);
	for(const double value:traction_pa)
		if(!std::isfinite(value)) throw std::invalid_argument("native shell traction is nonfinite");
	std::vector<double> force(3*reference.control_points.size(),0.0);
	const auto spans_u=NativeShellNonzeroSpans(reference.knots_u,reference.degree_u);
	const auto spans_v=NativeShellNonzeroSpans(reference.knots_v,reference.degree_v);
	const auto gauss=NativeShellGaussFour();
	for(const auto& span_v:spans_v) for(const auto& span_u:spans_u)
		for(const auto& gv:gauss) for(const auto& gu:gauss) {
			const double u=0.5*((1.0-gu[0])*span_u[0]+(1.0+gu[0])*span_u[1]);
			const double v=0.5*((1.0-gv[0])*span_v[0]+(1.0+gv[0])*span_v[1]);
			const auto point=EvaluateNativeNurbsShellSurface(reference,u,v);
			const double weight=0.25*(span_u[1]-span_u[0])*(span_v[1]-span_v[0])
				*gu[1]*gv[1]*point.surface_jacobian;
			for(std::size_t node=0;node<point.basis.size();++node)
				for(std::size_t component=0;component<3;++component)
					force[3*node+component]+=weight*point.basis[node]*traction_pa[component];
		}
	return force;
}

inline std::vector<double> SolveNativeShellDense(std::vector<double> matrix,
	std::vector<double> right)
{
	const std::size_t n=right.size();
	if(matrix.size()!=n*n) throw std::invalid_argument("native shell dense matrix size mismatch");
	for(std::size_t column=0;column<n;++column) {
		std::size_t pivot=column;
		for(std::size_t row=column+1;row<n;++row)
			if(std::abs(matrix[row*n+column])>std::abs(matrix[pivot*n+column])) pivot=row;
		if(!(std::abs(matrix[pivot*n+column])>1e-14)||!std::isfinite(matrix[pivot*n+column]))
			throw std::runtime_error("native shell constrained tangent is singular");
		if(pivot!=column) {
			for(std::size_t j=column;j<n;++j) std::swap(matrix[column*n+j],matrix[pivot*n+j]);
			std::swap(right[column],right[pivot]);
		}
		for(std::size_t row=column+1;row<n;++row) {
			const double factor=matrix[row*n+column]/matrix[column*n+column];
			for(std::size_t j=column;j<n;++j) matrix[row*n+j]-=factor*matrix[column*n+j];
			right[row]-=factor*right[column];
		}
	}
	std::vector<double> result(n,0.0);
	for(std::size_t reverse=0;reverse<n;++reverse) {
		const std::size_t row=n-1-reverse;
		double value=right[row];
		for(std::size_t j=row+1;j<n;++j) value-=matrix[row*n+j]*result[j];
		result[row]=value/matrix[row*n+row];
	}
	return result;
}

struct NativeIgaShellStaticOptions
{
	int maximum_newton_iterations=12;
	double residual_tolerance_n=1e-5;
	double residual_difference_step_m=2e-7;
	double tangent_difference_step_m=2e-5;
	int maximum_line_search_cuts=10;
};

struct NativeIgaShellStaticResult
{
	std::vector<double> displacement_m;
	std::vector<double> reaction_n;
	int newton_iterations=0;
	double free_residual_norm_n=0.0;
	bool converged=false;
};

inline NativeIgaShellStaticResult SolveNativeIgaShellStaticPatch(
	const NativeNurbsShellPatch& reference,const NativeKirchhoffLoveMaterial& material,
	const std::vector<double>& external_force_n,const std::map<std::size_t,double>& constraints,
	const NativeIgaShellStaticOptions& options={})
{
	const std::size_t dofs=3*reference.control_points.size();
	if(external_force_n.size()!=dofs) throw std::invalid_argument("native shell force size mismatch");
	if(options.maximum_newton_iterations<1||!(options.residual_tolerance_n>0.0)
		||options.maximum_line_search_cuts<0)
		throw std::invalid_argument("native shell static solver options are invalid");
	NativeIgaShellStaticResult result;result.displacement_m.assign(dofs,0.0);
	std::vector<std::size_t> free;
	for(const auto& entry:constraints) {
		if(entry.first>=dofs||!std::isfinite(entry.second))
			throw std::invalid_argument("native shell constraint is invalid");
		result.displacement_m[entry.first]=entry.second;
	}
	for(std::size_t dof=0;dof<dofs;++dof) if(!constraints.count(dof)) free.push_back(dof);
	if(free.empty()) throw std::invalid_argument("native shell static solve has no free DOF");
	for(int iteration=0;iteration<options.maximum_newton_iterations;++iteration) {
		const auto system=BuildNativeIgaShellPatchSystem(reference,result.displacement_m,material,
			options.residual_difference_step_m,options.tangent_difference_step_m);
		std::vector<double> residual(dofs);
		for(std::size_t i=0;i<dofs;++i) residual[i]=system.residual_n[i]-external_force_n[i];
		double norm_squared=0.0;for(const auto dof:free) norm_squared+=residual[dof]*residual[dof];
		result.free_residual_norm_n=std::sqrt(norm_squared);
		result.newton_iterations=iteration;
		if(result.free_residual_norm_n<=options.residual_tolerance_n) {
			result.converged=true;result.reaction_n=std::move(residual);return result;
		}
		std::vector<double> reduced_matrix(free.size()*free.size()),right(free.size());
		for(std::size_t i=0;i<free.size();++i) {
			right[i]=-residual[free[i]];
			for(std::size_t j=0;j<free.size();++j)
				reduced_matrix[i*free.size()+j]=system.tangent_n_m[free[i]*dofs+free[j]];
		}
		const auto update=SolveNativeShellDense(std::move(reduced_matrix),std::move(right));
		auto potential=[&](const std::vector<double>& displacement) {
			double value=NativeIgaShellPatchEnergy(reference,displacement,material);
			for(std::size_t i=0;i<dofs;++i) value-=external_force_n[i]*displacement[i];
			return value;
		};
		const double old_potential=potential(result.displacement_m);
		double scale=1.0;bool accepted=false;
		for(int cut=0;cut<=options.maximum_line_search_cuts;++cut) {
			auto trial=result.displacement_m;
			for(std::size_t i=0;i<free.size();++i) trial[free[i]]+=scale*update[i];
			if(potential(trial)<old_potential) {
				result.displacement_m=std::move(trial);accepted=true;break;
			}
			scale*=0.5;
		}
		if(!accepted) throw std::runtime_error("native shell Newton line search failed");
	}
	throw std::runtime_error("native shell Newton solve did not converge");
}

} // namespace iga

#endif
