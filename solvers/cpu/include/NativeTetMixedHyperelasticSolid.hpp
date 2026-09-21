#ifndef IGA_NATIVE_TET_MIXED_HYPERELASTIC_SOLID_HPP
#define IGA_NATIVE_TET_MIXED_HYPERELASTIC_SOLID_HPP

// Stabilized equal-order P1 displacement/P1 pressure total-Lagrangian
// tetrahedron for nearly incompressible solids. This is a distinct formulation
// from the displacement-only compressible element.
#include "NativeTetHyperelasticSolid.hpp"

#include <array>
#include <cmath>
#include <stdexcept>

namespace iga {

struct NativeTetMixedSolidMaterial
{
	double shear_modulus_pa=0.0;
	double bulk_modulus_pa=0.0;
	double density_kg_m3=0.0;
	double pressure_stabilization=0.1;
};

struct NativeTetMixedSolidElementSystem
{
	static constexpr std::size_t displacement_dofs=12;
	static constexpr std::size_t pressure_dofs=4;
	static constexpr std::size_t dofs=16;
	std::array<double,dofs> residual{};
	std::array<double,dofs*dofs> tangent{};
	std::array<double,displacement_dofs*displacement_dofs> mass_kg{};
	NativeTetSolidMatrix3 deformation_gradient{},first_piola_pa{};
	double deformation_jacobian=0.0,reference_volume_m3=0.0;
	double deviatoric_energy_j=0.0,mixed_potential_j=0.0;
};

inline void ValidateNativeTetMixedSolidMaterial(const NativeTetMixedSolidMaterial& material)
{
	if(!(material.shear_modulus_pa>0.0)||!std::isfinite(material.shear_modulus_pa)
		||!(material.bulk_modulus_pa>0.0)||!std::isfinite(material.bulk_modulus_pa)
		||material.bulk_modulus_pa/material.shear_modulus_pa>1.e9
		||!(material.density_kg_m3>0.0)||!std::isfinite(material.density_kg_m3)
		||!(material.pressure_stabilization>0.0)||!std::isfinite(material.pressure_stabilization))
		throw std::invalid_argument("native mixed tetrahedral solid material is invalid");
}

inline NativeTetMixedSolidMaterial NativeTetMixedSolidMaterialFromYoungPoisson(
	double young_modulus_pa,double poisson_ratio,double density_kg_m3,
	double pressure_stabilization=0.1)
{
	if(!(young_modulus_pa>0.0)||!std::isfinite(young_modulus_pa)
		||!std::isfinite(poisson_ratio)||poisson_ratio<=-1.0||poisson_ratio>=0.5)
		throw std::invalid_argument("native mixed tetrahedral solid E/nu values are invalid");
	NativeTetMixedSolidMaterial result{young_modulus_pa/(2.*(1.+poisson_ratio)),
		young_modulus_pa/(3.*(1.-2.*poisson_ratio)),density_kg_m3,pressure_stabilization};
	ValidateNativeTetMixedSolidMaterial(result);return result;
}

inline NativeTetMixedSolidElementSystem BuildNativeTetMixedHyperelasticSolidElement(
	const NativeTetMesh& reference_mesh,const NativeTetCell& cell,
	const std::array<double,12>& displacement_m,const std::array<double,4>& pressure_pa,
	const NativeTetMixedSolidMaterial& material)
{
	ValidateNativeTetMixedSolidMaterial(material);
	for(double value:displacement_m)if(!std::isfinite(value))
		throw std::invalid_argument("native mixed tetrahedral solid displacement is nonfinite");
	for(double value:pressure_pa)if(!std::isfinite(value))
		throw std::invalid_argument("native mixed tetrahedral solid pressure is nonfinite");
	const auto geometry=EvaluateNativeTetGeometry(reference_mesh,cell);
	NativeTetMixedSolidElementSystem result;result.reference_volume_m3=geometry.determinant/6.;
	for(int i=0;i<3;++i)result.deformation_gradient[i][i]=1.;
	for(std::size_t node=0;node<4;++node)for(int i=0;i<3;++i)for(int j=0;j<3;++j)
		result.deformation_gradient[i][j]+=displacement_m[3*node+i]*geometry.barycentric_gradients[node][j];
	result.deformation_jacobian=NativeTetSolidDeterminant(result.deformation_gradient);
	if(!(result.deformation_jacobian>0.)||!std::isfinite(result.deformation_jacobian))
		throw std::runtime_error("native mixed tetrahedral solid deformation Jacobian must be positive");
	const auto inverse=NativeTetSolidInverse(result.deformation_gradient);NativeTetSolidMatrix3 inverse_transpose{};
	for(int i=0;i<3;++i)for(int j=0;j<3;++j)inverse_transpose[i][j]=inverse[j][i];
	double invariant=0.,mean_pressure=0.;
	for(int i=0;i<3;++i)for(int j=0;j<3;++j)
		invariant+=result.deformation_gradient[i][j]*result.deformation_gradient[i][j];
	for(double value:pressure_pa)mean_pressure+=value/4.;
	const double mu=material.shear_modulus_pa,kappa=material.bulk_modulus_pa;
	const double isochoric_scale=std::pow(result.deformation_jacobian,-2./3.);
	NativeTetSolidMatrix3 deviatoric_piola{};
	for(int i=0;i<3;++i)for(int j=0;j<3;++j){
		deviatoric_piola[i][j]=mu*isochoric_scale*(result.deformation_gradient[i][j]
			-invariant*inverse_transpose[i][j]/3.);
		result.first_piola_pa[i][j]=deviatoric_piola[i][j]
			+mean_pressure*result.deformation_jacobian*inverse_transpose[i][j];
	}
	result.deviatoric_energy_j=result.reference_volume_m3*.5*mu*(isochoric_scale*invariant-3.);
	const double characteristic_length=std::cbrt(6.*result.reference_volume_m3);
	const double stabilization=material.pressure_stabilization*characteristic_length*characteristic_length/mu;
	double pressure_square_integral=0.,pressure_gradient_square=0.;
	for(std::size_t a=0;a<4;++a)for(std::size_t b=0;b<4;++b){
		const double mass=result.reference_volume_m3*(a==b?2.:1.)/20.;
		pressure_square_integral+=mass*pressure_pa[a]*pressure_pa[b];
		for(int j=0;j<3;++j)pressure_gradient_square+=result.reference_volume_m3
			*geometry.barycentric_gradients[a][j]*geometry.barycentric_gradients[b][j]
			*pressure_pa[a]*pressure_pa[b];
	}
	result.mixed_potential_j=result.deviatoric_energy_j+result.reference_volume_m3
		*mean_pressure*(result.deformation_jacobian-1.)-.5*pressure_square_integral/kappa
		-.5*stabilization*pressure_gradient_square;
	for(std::size_t a=0;a<4;++a)for(int i=0;i<3;++i){
		const std::size_t row=3*a+i;
		for(int j=0;j<3;++j)result.residual[row]+=result.reference_volume_m3
			*result.first_piola_pa[i][j]*geometry.barycentric_gradients[a][j];
		for(std::size_t b=0;b<4;++b)for(int k=0;k<3;++k){
			const std::size_t column=3*b+k;
			for(int j=0;j<3;++j)for(int l=0;l<3;++l){
				const double g=result.deformation_gradient[i][j]-invariant*inverse_transpose[i][j]/3.;
				double elasticity=mu*isochoric_scale*((i==k&&j==l?1.:0.)
					-2.*inverse_transpose[k][l]*g/3.-2.*result.deformation_gradient[k][l]
						*inverse_transpose[i][j]/3.+invariant*inverse_transpose[i][l]
						*inverse_transpose[k][j]/3.);
				elasticity+=mean_pressure*result.deformation_jacobian*(inverse_transpose[k][l]
					*inverse_transpose[i][j]-inverse_transpose[i][l]*inverse_transpose[k][j]);
				result.tangent[row*16+column]+=result.reference_volume_m3
					*geometry.barycentric_gradients[a][j]*elasticity*geometry.barycentric_gradients[b][l];
			}
		}
		for(std::size_t b=0;b<4;++b){
			double coupling=0.;for(int j=0;j<3;++j)coupling+=result.reference_volume_m3/4.
				*result.deformation_jacobian*inverse_transpose[i][j]*geometry.barycentric_gradients[a][j];
			result.tangent[row*16+12+b]+=coupling;result.tangent[(12+b)*16+row]+=coupling;
		}
	}
	for(std::size_t a=0;a<4;++a){
		const std::size_t row=12+a;result.residual[row]+=result.reference_volume_m3/4.
			*(result.deformation_jacobian-1.);
		for(std::size_t b=0;b<4;++b){
			const double mass=result.reference_volume_m3*(a==b?2.:1.)/20.;double gradient=0.;
			for(int j=0;j<3;++j)gradient+=result.reference_volume_m3
				*geometry.barycentric_gradients[a][j]*geometry.barycentric_gradients[b][j];
			result.residual[row]-=(mass/kappa+stabilization*gradient)*pressure_pa[b];
			result.tangent[row*16+12+b]-=mass/kappa+stabilization*gradient;
		}
	}
	for(std::size_t a=0;a<4;++a)for(std::size_t b=0;b<4;++b){
		const double scalar=material.density_kg_m3*result.reference_volume_m3*(a==b?2.:1.)/20.;
		for(int component=0;component<3;++component)
			result.mass_kg[(3*a+component)*12+3*b+component]=scalar;
	}
	return result;
}

} // namespace iga
#endif
