#ifndef IGA_NATIVE_TET_SOLID_PRESTRESS_HPP
#define IGA_NATIVE_TET_SOLID_PRESTRESS_HPP

// Small-case inverse-elastostatics reference initialization. A loaded image
// mesh is never treated as stress free: given material, dead loads and zero
// displacement supports, recover an unloaded mesh whose forward equilibrium
// reconstructs the observed loaded coordinates.
#include "NativeTetSolidStaticSolver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

namespace iga {

struct NativeTetSolidPrestressOptions
{
	int maximum_inverse_iterations=60;
	double reference_update_relaxation=0.5;
	double absolute_reconstruction_tolerance_m=1.e-9;
	double minimum_reference_update_scale=1.0/1024.0;
	NativeTetSolidStaticOptions forward_options{};
};

struct NativeTetSolidPrestressResult
{
	NativeTetMesh unloaded_reference_mesh;
	std::vector<double> loaded_displacement_m;
	std::vector<double> equilibrating_internal_force_n;
	std::vector<NativeTetSolidMatrix3> loaded_first_piola_pa;
	int inverse_iterations=0;
	double maximum_reconstruction_error_m=0.0;
	double minimum_loaded_deformation_jacobian=0.0;
	bool converged=false;
};

namespace native_tet_solid_prestress_detail {

inline void ValidateOptions(const NativeTetSolidPrestressOptions& options)
{
	if(options.maximum_inverse_iterations<=0
		||!std::isfinite(options.reference_update_relaxation)
		||!(options.reference_update_relaxation>0.0)||options.reference_update_relaxation>1.0
		||!std::isfinite(options.absolute_reconstruction_tolerance_m)
		||!(options.absolute_reconstruction_tolerance_m>0.0)
		||!std::isfinite(options.minimum_reference_update_scale)
		||!(options.minimum_reference_update_scale>0.0)||options.minimum_reference_update_scale>1.0)
		throw std::invalid_argument("native tetrahedral solid prestress options are invalid");
	native_tet_solid_static_detail::ValidateOptions(options.forward_options);
}

inline std::vector<double> LoadedMinusReference(const NativeTetMesh& loaded,
	const NativeTetMesh& reference)
{
	std::vector<double> result(3*loaded.points.size());
	for(std::size_t node=0;node<loaded.points.size();++node)for(int component=0;component<3;++component)
		result[3*node+component]=loaded.points[node][component]-reference.points[node][component];
	return result;
}

inline double ReconstructionError(const NativeTetMesh& loaded,const NativeTetMesh& reference,
	const std::vector<double>& displacement)
{
	double result=0.0;
	for(std::size_t node=0;node<loaded.points.size();++node)for(int component=0;component<3;++component)
		result=std::max(result,std::abs(reference.points[node][component]
			+displacement[3*node+component]-loaded.points[node][component]));
	return result;
}

inline bool HasPositiveReferenceCells(const NativeTetMesh& mesh)
{
	try{for(const auto& cell:mesh.cells)(void)EvaluateNativeTetGeometry(mesh,cell);}
	catch(const std::exception&){return false;}return true;
}

} // namespace native_tet_solid_prestress_detail

inline NativeTetSolidPrestressResult EstimateNativeTetSolidUnloadedReference(
	const NativeTetMesh& loaded_image_mesh,const NativeTetSolidMaterial& material,
	const std::vector<double>& external_force_n,
	const std::map<std::size_t,double>& zero_displacement_supports,
	const NativeTetSolidPrestressOptions& options={})
{
	using namespace native_tet_solid_prestress_detail;
	ValidateNativeTetSolidMaterial(material);ValidateOptions(options);
	const std::size_t dofs=3*loaded_image_mesh.points.size();
	if(loaded_image_mesh.points.empty()||loaded_image_mesh.cells.empty()
		||external_force_n.size()!=dofs)
		throw std::invalid_argument("native tetrahedral solid prestress input shape is invalid");
	for(double value:external_force_n)if(!std::isfinite(value))
		throw std::invalid_argument("native tetrahedral solid prestress load is nonfinite");
	for(const auto& support:zero_displacement_supports)
		if(support.first>=dofs||support.second!=0.0)
			throw std::invalid_argument("native tetrahedral solid inverse solve requires explicit zero supports");
	if(zero_displacement_supports.empty()||zero_displacement_supports.size()>=dofs)
		throw std::invalid_argument("native tetrahedral solid prestress supports are incomplete or overconstrained");
	if(!HasPositiveReferenceCells(loaded_image_mesh))
		throw std::invalid_argument("native tetrahedral solid loaded image mesh is invalid");
	NativeTetMesh reference=loaded_image_mesh;NativeTetSolidStaticResult forward;
	NativeTetSolidPrestressResult result;
	for(int iteration=0;iteration<=options.maximum_inverse_iterations;++iteration){
		const auto initial=LoadedMinusReference(loaded_image_mesh,reference);
		forward=SolveNativeTetSolidStatic(reference,material,external_force_n,
			zero_displacement_supports,initial,options.forward_options);
		const double error=ReconstructionError(loaded_image_mesh,reference,forward.displacement_m);
		result.inverse_iterations=iteration;result.maximum_reconstruction_error_m=error;
		if(error<=options.absolute_reconstruction_tolerance_m){result.converged=true;break;}
		if(iteration==options.maximum_inverse_iterations)break;
		double scale=options.reference_update_relaxation;bool accepted=false;
		while(scale>=options.minimum_reference_update_scale){
			auto candidate=reference;
			for(std::size_t dof=0;dof<dofs;++dof)if(!zero_displacement_supports.count(dof)){
				const std::size_t node=dof/3;const int component=static_cast<int>(dof%3);
				const double mismatch=reference.points[node][component]+forward.displacement_m[dof]
					-loaded_image_mesh.points[node][component];
				candidate.points[node][component]-=scale*mismatch;
			}
			if(HasPositiveReferenceCells(candidate)){reference=std::move(candidate);accepted=true;break;}
			scale*=0.5;
		}
		if(!accepted)throw std::runtime_error("native tetrahedral solid inverse reference update inverted the mesh");
	}
	if(!result.converged)
		throw std::runtime_error("native tetrahedral solid inverse elastostatics did not converge");
	result.unloaded_reference_mesh=reference;result.loaded_displacement_m=forward.displacement_m;
	result.minimum_loaded_deformation_jacobian=forward.minimum_deformation_jacobian;
	result.equilibrating_internal_force_n.resize(dofs);
	for(std::size_t dof=0;dof<dofs;++dof)
		result.equilibrating_internal_force_n[dof]=forward.reaction_n[dof]+external_force_n[dof];
	for(const auto& cell:reference.cells){
		std::array<double,12> local{};
		for(std::size_t node=0;node<4;++node)for(int component=0;component<3;++component)
			local[3*node+component]=forward.displacement_m[3*cell.nodes[node]+component];
		result.loaded_first_piola_pa.push_back(
			BuildNativeTetHyperelasticSolidElement(reference,cell,local,material).first_piola_pa);
	}
	return result;
}

} // namespace iga
#endif
