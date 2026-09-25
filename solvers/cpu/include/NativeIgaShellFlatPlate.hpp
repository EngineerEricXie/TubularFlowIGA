#ifndef IGA_NATIVE_IGA_SHELL_FLAT_PLATE_HPP
#define IGA_NATIVE_IGA_SHELL_FLAT_PLATE_HPP

#include "NativeIgaShellStaticSolver.hpp"

#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <vector>

namespace iga {

struct NativeIgaFlatPlateSolution
{
	std::vector<double> transverse_control_displacement_m;
	double flexural_rigidity_n_m=0.0;
};

inline NativeIgaFlatPlateSolution SolveNativeIgaSimplySupportedFlatPlate(
	const NativeNurbsShellPatch& patch,const NativeKirchhoffLoveMaterial& material,
	const std::function<double(double,double)>& transverse_load_pa)
{
	ValidateNativeNurbsShellPatch(patch);ValidateNativeKirchhoffLoveMaterial(material);
	if(!transverse_load_pa) throw std::invalid_argument("native flat-plate load is empty");
	const std::size_t controls=patch.control_points.size();
	std::vector<double> stiffness(controls*controls,0.0),force(controls,0.0);
	const double rigidity=material.young_modulus_pa*material.thickness_m
		*material.thickness_m*material.thickness_m/
		(12.0*(1.0-material.poisson_ratio*material.poisson_ratio));
	const auto spans_u=NativeShellNonzeroSpans(patch.knots_u,patch.degree_u);
	const auto spans_v=NativeShellNonzeroSpans(patch.knots_v,patch.degree_v);
	const auto gauss=NativeShellGaussFour();
	for(const auto& span_v:spans_v) for(const auto& span_u:spans_u)
		for(const auto& gv:gauss) for(const auto& gu:gauss) {
			const double u=0.5*((1.0-gu[0])*span_u[0]+(1.0+gu[0])*span_u[1]);
			const double v=0.5*((1.0-gv[0])*span_v[0]+(1.0+gv[0])*span_v[1]);
			const auto point=EvaluateNativeNurbsShellSurface(patch,u,v);
			// This verification operator is deliberately restricted to an axis-aligned
			// unit plane, so parametric and physical second derivatives coincide.
			if(std::abs(point.tangent_u[0]-1.0)>1e-11||std::abs(point.tangent_u[1])>1e-11
				||std::abs(point.tangent_v[0])>1e-11||std::abs(point.tangent_v[1]-1.0)>1e-11
				||std::abs(point.unit_normal[2]-1.0)>1e-11)
				throw std::invalid_argument("native flat-plate verification requires the unit xy plane");
			const double weight=0.25*(span_u[1]-span_u[0])*(span_v[1]-span_v[0])
				*gu[1]*gv[1]*point.surface_jacobian;
			const double load=transverse_load_pa(point.position[0],point.position[1]);
			if(!std::isfinite(load)) throw std::runtime_error("native flat-plate load is nonfinite");
			for(std::size_t i=0;i<controls;++i) {
				force[i]+=weight*point.basis[i]*load;
				for(std::size_t j=0;j<controls;++j)
					stiffness[i*controls+j]+=weight*rigidity*(
						point.basis_uu[i]*point.basis_uu[j]
						+point.basis_vv[i]*point.basis_vv[j]
						+material.poisson_ratio*(point.basis_uu[i]*point.basis_vv[j]
							+point.basis_vv[i]*point.basis_uu[j])
						+2.0*(1.0-material.poisson_ratio)*point.basis_uv[i]*point.basis_uv[j]);
			}
		}
	std::vector<std::size_t> free;
	for(std::size_t j=1;j+1<patch.controls_v;++j)
		for(std::size_t i=1;i+1<patch.controls_u;++i) free.push_back(j*patch.controls_u+i);
	if(free.empty()) throw std::invalid_argument("native flat-plate mesh has no interior control DOF");
	std::vector<double> reduced(free.size()*free.size()),right(free.size());
	for(std::size_t i=0;i<free.size();++i) {
		right[i]=force[free[i]];
		for(std::size_t j=0;j<free.size();++j)
			reduced[i*free.size()+j]=stiffness[free[i]*controls+free[j]];
	}
	const auto interior=SolveNativeShellDense(std::move(reduced),std::move(right));
	NativeIgaFlatPlateSolution result;
	result.transverse_control_displacement_m.assign(controls,0.0);
	for(std::size_t i=0;i<free.size();++i)
		result.transverse_control_displacement_m[free[i]]=interior[i];
	result.flexural_rigidity_n_m=rigidity;
	return result;
}

} // namespace iga

#endif
