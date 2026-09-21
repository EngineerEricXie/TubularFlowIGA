#ifndef IGA_NATIVE_TET_HYPERELASTIC_SOLID_HPP
#define IGA_NATIVE_TET_HYPERELASTIC_SOLID_HPP

#include "NativeTetFem.hpp"

#include <array>
#include <cmath>
#include <stdexcept>

namespace iga {

struct NativeTetSolidMaterial
{
	double young_modulus_pa=0.0;
	double poisson_ratio=0.0;
	double density_kg_m3=0.0;
};

struct NativeTetSolidElementSystem
{
	static constexpr std::size_t nodes=4;
	static constexpr std::size_t dofs=12;
	std::array<double,dofs> residual_n{};
	std::array<double,dofs*dofs> tangent_n_m{};
	std::array<double,dofs*dofs> mass_kg{};
	std::array<std::array<double,3>,3> deformation_gradient{};
	std::array<std::array<double,3>,3> first_piola_pa{};
	double deformation_jacobian=0.0;
	double strain_energy_j=0.0;
	double reference_volume_m3=0.0;
};

using NativeTetSolidMatrix3=std::array<std::array<double,3>,3>;

inline double NativeTetSolidDeterminant(const NativeTetSolidMatrix3& a)
{
	return a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
		-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
		+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);
}

inline NativeTetSolidMatrix3 NativeTetSolidInverse(const NativeTetSolidMatrix3& a)
{
	const double determinant=NativeTetSolidDeterminant(a);
	if(!(determinant>0.0)||!std::isfinite(determinant))
		throw std::runtime_error("native tetrahedral solid deformation Jacobian must be finite and positive");
	NativeTetSolidMatrix3 result{};
	result[0][0]=(a[1][1]*a[2][2]-a[1][2]*a[2][1])/determinant;
	result[0][1]=(a[0][2]*a[2][1]-a[0][1]*a[2][2])/determinant;
	result[0][2]=(a[0][1]*a[1][2]-a[0][2]*a[1][1])/determinant;
	result[1][0]=(a[1][2]*a[2][0]-a[1][0]*a[2][2])/determinant;
	result[1][1]=(a[0][0]*a[2][2]-a[0][2]*a[2][0])/determinant;
	result[1][2]=(a[0][2]*a[1][0]-a[0][0]*a[1][2])/determinant;
	result[2][0]=(a[1][0]*a[2][1]-a[1][1]*a[2][0])/determinant;
	result[2][1]=(a[0][1]*a[2][0]-a[0][0]*a[2][1])/determinant;
	result[2][2]=(a[0][0]*a[1][1]-a[0][1]*a[1][0])/determinant;
	return result;
}

inline void ValidateNativeTetSolidMaterial(const NativeTetSolidMaterial& material)
{
	if(!(material.young_modulus_pa>0.0)||!std::isfinite(material.young_modulus_pa))
		throw std::invalid_argument("native tetrahedral solid Young modulus must be finite and positive");
	if(!std::isfinite(material.poisson_ratio)||material.poisson_ratio<=-1.0
		||material.poisson_ratio>0.45)
		throw std::invalid_argument("native displacement-only tetrahedral solid Poisson ratio must lie in (-1,0.45]");
	if(!(material.density_kg_m3>0.0)||!std::isfinite(material.density_kg_m3))
		throw std::invalid_argument("native tetrahedral solid density must be finite and positive");
}

inline NativeTetSolidElementSystem BuildNativeTetHyperelasticSolidElement(
	const NativeTetMesh& reference_mesh,const NativeTetCell& cell,
	const std::array<double,NativeTetSolidElementSystem::dofs>& displacement_m,
	const NativeTetSolidMaterial& material)
{
	ValidateNativeTetSolidMaterial(material);
	for(double value:displacement_m) if(!std::isfinite(value))
		throw std::invalid_argument("native tetrahedral solid displacement is nonfinite");
	const auto geometry=EvaluateNativeTetGeometry(reference_mesh,cell);
	NativeTetSolidElementSystem result;
	result.reference_volume_m3=geometry.determinant/6.0;
	for(int i=0;i<3;++i) result.deformation_gradient[i][i]=1.0;
	for(std::size_t node=0;node<4;++node)
		for(int i=0;i<3;++i) for(int j=0;j<3;++j)
			result.deformation_gradient[i][j]+=displacement_m[3*node+i]
				*geometry.barycentric_gradients[node][j];
	result.deformation_jacobian=NativeTetSolidDeterminant(result.deformation_gradient);
	if(!(result.deformation_jacobian>0.0)||!std::isfinite(result.deformation_jacobian))
		throw std::runtime_error("native tetrahedral solid deformation Jacobian must be finite and positive");
	const auto inverse=NativeTetSolidInverse(result.deformation_gradient);
	NativeTetSolidMatrix3 inverse_transpose{};
	for(int i=0;i<3;++i) for(int j=0;j<3;++j) inverse_transpose[i][j]=inverse[j][i];
	const double young=material.young_modulus_pa,poisson=material.poisson_ratio;
	const double mu=young/(2.0*(1.0+poisson));
	const double lambda=young*poisson/((1.0+poisson)*(1.0-2.0*poisson));
	const double log_j=std::log(result.deformation_jacobian);
	double trace_c=0.0;
	for(int i=0;i<3;++i) for(int j=0;j<3;++j)
		trace_c+=result.deformation_gradient[i][j]*result.deformation_gradient[i][j];
	const double energy_density=0.5*mu*(trace_c-3.0)-mu*log_j+0.5*lambda*log_j*log_j;
	result.strain_energy_j=result.reference_volume_m3*energy_density;
	if(!std::isfinite(result.strain_energy_j))
		throw std::runtime_error("native tetrahedral solid strain energy is nonfinite");
	for(int i=0;i<3;++i) for(int j=0;j<3;++j)
		result.first_piola_pa[i][j]=mu*result.deformation_gradient[i][j]
			+(lambda*log_j-mu)*inverse_transpose[i][j];
	for(std::size_t a=0;a<4;++a) for(int i=0;i<3;++i) {
		const std::size_t row=3*a+i;
		for(int j=0;j<3;++j) result.residual_n[row]+=result.reference_volume_m3
			*result.first_piola_pa[i][j]*geometry.barycentric_gradients[a][j];
		for(std::size_t b=0;b<4;++b) for(int k=0;k<3;++k) {
			const std::size_t column=3*b+k;
			for(int j=0;j<3;++j) for(int l=0;l<3;++l) {
				double elasticity=mu*(i==k&&j==l?1.0:0.0)
					+lambda*inverse_transpose[i][j]*inverse_transpose[k][l]
					+(mu-lambda*log_j)*inverse_transpose[i][l]*inverse_transpose[k][j];
				result.tangent_n_m[row*12+column]+=result.reference_volume_m3
					*geometry.barycentric_gradients[a][j]*elasticity
					*geometry.barycentric_gradients[b][l];
			}
		}
	}
	for(std::size_t a=0;a<4;++a) for(std::size_t b=0;b<4;++b) {
		const double scalar=material.density_kg_m3*result.reference_volume_m3
			*(a==b?2.0:1.0)/20.0;
		for(int component=0;component<3;++component)
			result.mass_kg[(3*a+component)*12+3*b+component]=scalar;
	}
	return result;
}

} // namespace iga

#endif
