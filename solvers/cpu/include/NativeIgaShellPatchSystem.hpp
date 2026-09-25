#ifndef IGA_NATIVE_IGA_SHELL_PATCH_SYSTEM_HPP
#define IGA_NATIVE_IGA_SHELL_PATCH_SYSTEM_HPP

#include "NativeIgaKirchhoffLove.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace iga {

struct NativeIgaShellPatchSystem
{
	double strain_energy_j=0.0;
	std::vector<double> residual_n;
	std::vector<double> tangent_n_m;
	std::vector<double> consistent_mass_kg;
};

inline std::vector<std::array<double,2>> NativeShellGaussFour()
{
	return {{{-0.8611363115940526,0.3478548451374538},
		{-0.3399810435848563,0.6521451548625461},
		{0.3399810435848563,0.6521451548625461},
		{0.8611363115940526,0.3478548451374538}}};
}

inline std::vector<std::array<double,2>> NativeShellNonzeroSpans(
	const std::vector<double>& knots,int degree)
{
	std::vector<std::array<double,2>> result;
	const std::size_t begin=static_cast<std::size_t>(degree);
	const std::size_t end=knots.size()-static_cast<std::size_t>(degree)-1;
	for(std::size_t i=begin;i<end;++i)
		if(knots[i+1]>knots[i]) result.push_back({{knots[i],knots[i+1]}});
	if(result.empty()) throw std::invalid_argument("native shell patch has no nonzero knot span");
	return result;
}

inline NativeNurbsShellPatch NativeShellCurrentPatch(const NativeNurbsShellPatch& reference,
	const std::vector<double>& displacement)
{
	if(displacement.size()!=3*reference.control_points.size())
		throw std::invalid_argument("native shell displacement size is inconsistent with controls");
	NativeNurbsShellPatch current=reference;
	for(std::size_t i=0;i<current.control_points.size();++i)
		for(std::size_t component=0;component<3;++component) {
			const double value=displacement[3*i+component];
			if(!std::isfinite(value)) throw std::invalid_argument("native shell displacement is nonfinite");
			current.control_points[i][component]+=value;
		}
	return current;
}

inline double NativeIgaShellPatchEnergy(const NativeNurbsShellPatch& reference,
	const std::vector<double>& displacement,const NativeKirchhoffLoveMaterial& material)
{
	ValidateNativeNurbsShellPatch(reference);ValidateNativeKirchhoffLoveMaterial(material);
	const auto current=NativeShellCurrentPatch(reference,displacement);
	const auto spans_u=NativeShellNonzeroSpans(reference.knots_u,reference.degree_u);
	const auto spans_v=NativeShellNonzeroSpans(reference.knots_v,reference.degree_v);
	const auto gauss=NativeShellGaussFour();
	double energy=0.0;
	for(const auto& span_v:spans_v) for(const auto& span_u:spans_u)
		for(const auto& gv:gauss) for(const auto& gu:gauss) {
			const double u=0.5*((1.0-gu[0])*span_u[0]+(1.0+gu[0])*span_u[1]);
			const double v=0.5*((1.0-gv[0])*span_v[0]+(1.0+gv[0])*span_v[1]);
			const auto reference_point=EvaluateNativeNurbsShellSurface(reference,u,v);
			const auto current_point=EvaluateNativeNurbsShellSurface(current,u,v);
			const auto response=EvaluateNativeKirchhoffLovePoint(reference_point,current_point,material);
			const double parameter_weight=0.25*(span_u[1]-span_u[0])*(span_v[1]-span_v[0])
				*gu[1]*gv[1];
			energy+=parameter_weight*reference_point.surface_jacobian*
				(response.membrane_energy_per_reference_area_j_m2
				+response.bending_energy_per_reference_area_j_m2);
		}
	if(!std::isfinite(energy)||energy<0.0)
		throw std::runtime_error("native shell patch energy is invalid");
	return energy;
}

inline std::vector<double> NativeIgaShellPatchResidual(const NativeNurbsShellPatch& reference,
	const std::vector<double>& displacement,const NativeKirchhoffLoveMaterial& material,
	double difference_step_m)
{
	if(!(difference_step_m>0.0)||!std::isfinite(difference_step_m))
		throw std::invalid_argument("native shell residual difference step must be positive");
	std::vector<double> residual(displacement.size(),0.0),plus=displacement,minus=displacement;
	for(std::size_t dof=0;dof<displacement.size();++dof) {
		const double step=difference_step_m*std::max(1.0,std::abs(displacement[dof]));
		plus[dof]+=step;minus[dof]-=step;
		residual[dof]=(NativeIgaShellPatchEnergy(reference,plus,material)
			-NativeIgaShellPatchEnergy(reference,minus,material))/(2.0*step);
		plus[dof]=minus[dof]=displacement[dof];
	}
	return residual;
}

inline std::vector<double> NativeIgaShellPatchMass(const NativeNurbsShellPatch& reference,
	const NativeKirchhoffLoveMaterial& material)
{
	ValidateNativeNurbsShellPatch(reference);ValidateNativeKirchhoffLoveMaterial(material);
	const std::size_t dofs=3*reference.control_points.size();
	std::vector<double> mass(dofs*dofs,0.0);
	const auto spans_u=NativeShellNonzeroSpans(reference.knots_u,reference.degree_u);
	const auto spans_v=NativeShellNonzeroSpans(reference.knots_v,reference.degree_v);
	const auto gauss=NativeShellGaussFour();
	for(const auto& span_v:spans_v) for(const auto& span_u:spans_u)
		for(const auto& gv:gauss) for(const auto& gu:gauss) {
			const double u=0.5*((1.0-gu[0])*span_u[0]+(1.0+gu[0])*span_u[1]);
			const double v=0.5*((1.0-gv[0])*span_v[0]+(1.0+gv[0])*span_v[1]);
			const auto point=EvaluateNativeNurbsShellSurface(reference,u,v);
			const double weight=0.25*(span_u[1]-span_u[0])*(span_v[1]-span_v[0])
				*gu[1]*gv[1]*point.surface_jacobian*material.density_kg_m3*material.thickness_m;
			for(std::size_t i=0;i<point.basis.size();++i)
				for(std::size_t j=0;j<point.basis.size();++j)
					for(std::size_t component=0;component<3;++component)
						mass[(3*i+component)*dofs+3*j+component]+=weight*point.basis[i]*point.basis[j];
		}
	return mass;
}

inline NativeIgaShellPatchSystem BuildNativeIgaShellPatchSystem(
	const NativeNurbsShellPatch& reference,const std::vector<double>& displacement,
	const NativeKirchhoffLoveMaterial& material,double residual_step_m=1.0e-6,
	double tangent_step_m=2.0e-5)
{
	if(!(tangent_step_m>0.0)||!std::isfinite(tangent_step_m))
		throw std::invalid_argument("native shell tangent difference step must be positive");
	NativeIgaShellPatchSystem result;
	result.strain_energy_j=NativeIgaShellPatchEnergy(reference,displacement,material);
	result.residual_n=NativeIgaShellPatchResidual(reference,displacement,material,residual_step_m);
	const std::size_t dofs=displacement.size();
	result.tangent_n_m.assign(dofs*dofs,0.0);
	std::vector<double> plus=displacement,minus=displacement;
	for(std::size_t column=0;column<dofs;++column) {
		const double step=tangent_step_m*std::max(1.0,std::abs(displacement[column]));
		plus[column]+=step;minus[column]-=step;
		const auto plus_residual=NativeIgaShellPatchResidual(reference,plus,material,residual_step_m);
		const auto minus_residual=NativeIgaShellPatchResidual(reference,minus,material,residual_step_m);
		for(std::size_t row=0;row<dofs;++row)
			result.tangent_n_m[row*dofs+column]=(plus_residual[row]-minus_residual[row])/(2.0*step);
		plus[column]=minus[column]=displacement[column];
	}
	result.consistent_mass_kg=NativeIgaShellPatchMass(reference,material);
	return result;
}

} // namespace iga

#endif
