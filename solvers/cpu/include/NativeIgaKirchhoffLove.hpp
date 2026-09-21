#ifndef IGA_NATIVE_IGA_KIRCHHOFF_LOVE_HPP
#define IGA_NATIVE_IGA_KIRCHHOFF_LOVE_HPP

#include "NativeIgaShellSurface.hpp"

#include <array>
#include <cmath>
#include <stdexcept>

namespace iga {

struct NativeKirchhoffLoveMaterial
{
	double young_modulus_pa=0.0;
	double poisson_ratio=0.0;
	double thickness_m=0.0;
	double density_kg_m3=0.0;
};

using NativeShellTensor2=std::array<std::array<double,2>,2>;

struct NativeKirchhoffLovePointResponse
{
	NativeShellTensor2 membrane_strain{};
	NativeShellTensor2 curvature_change{};
	NativeShellTensor2 membrane_resultant_contravariant{};
	NativeShellTensor2 bending_moment_contravariant{};
	double membrane_energy_per_reference_area_j_m2=0.0;
	double bending_energy_per_reference_area_j_m2=0.0;
	double mass_per_reference_area_kg_m2=0.0;
};

inline double NativeShellDot(const std::array<double,3>& a,const std::array<double,3>& b)
{
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

inline NativeShellTensor2 NativeShellMetric(const NativeNurbsShellSurfacePoint& point)
{
	return {{{{NativeShellDot(point.tangent_u,point.tangent_u),
		NativeShellDot(point.tangent_u,point.tangent_v)}},
		{{NativeShellDot(point.tangent_v,point.tangent_u),
		NativeShellDot(point.tangent_v,point.tangent_v)}}}};
}

inline NativeShellTensor2 NativeShellCurvature(const NativeNurbsShellSurfacePoint& point)
{
	return {{{{NativeShellDot(point.unit_normal,point.second_uu),
		NativeShellDot(point.unit_normal,point.second_uv)}},
		{{NativeShellDot(point.unit_normal,point.second_uv),
		NativeShellDot(point.unit_normal,point.second_vv)}}}};
}

inline NativeShellTensor2 NativeShellInverseMetric(const NativeShellTensor2& metric)
{
	const double determinant=metric[0][0]*metric[1][1]-metric[0][1]*metric[1][0];
	if(!(determinant>0.0)||!std::isfinite(determinant))
		throw std::runtime_error("native shell reference metric is singular");
	return {{{{metric[1][1]/determinant,-metric[0][1]/determinant}},
		{{-metric[1][0]/determinant,metric[0][0]/determinant}}}};
}

inline void ValidateNativeKirchhoffLoveMaterial(const NativeKirchhoffLoveMaterial& material)
{
	if(!(material.young_modulus_pa>0.0)||!std::isfinite(material.young_modulus_pa))
		throw std::invalid_argument("native shell Young modulus must be finite and positive");
	if(!std::isfinite(material.poisson_ratio)||material.poisson_ratio<=-1.0
		||material.poisson_ratio>=0.5)
		throw std::invalid_argument("native shell Poisson ratio must lie in (-1,0.5)");
	if(!(material.thickness_m>0.0)||!std::isfinite(material.thickness_m))
		throw std::invalid_argument("native shell thickness must be finite and positive");
	if(!(material.density_kg_m3>0.0)||!std::isfinite(material.density_kg_m3))
		throw std::invalid_argument("native shell density must be finite and positive");
}

inline NativeShellTensor2 NativeShellPlaneStressResultant(const NativeShellTensor2& strain,
	const NativeShellTensor2& inverse_metric,double lambda,double mu,double scale)
{
	double trace=0.0;
	for(int a=0;a<2;++a) for(int b=0;b<2;++b)
		trace+=inverse_metric[a][b]*strain[a][b];
	NativeShellTensor2 result{};
	for(int a=0;a<2;++a) for(int b=0;b<2;++b) {
		result[a][b]=lambda*trace*inverse_metric[a][b];
		for(int c=0;c<2;++c) for(int d=0;d<2;++d)
			result[a][b]+=2.0*mu*inverse_metric[a][c]*inverse_metric[b][d]*strain[c][d];
		result[a][b]*=scale;
	}
	return result;
}

inline double NativeShellEnergyDensity(const NativeShellTensor2& strain,
	const NativeShellTensor2& resultant)
{
	double result=0.0;
	for(int a=0;a<2;++a) for(int b=0;b<2;++b) result+=strain[a][b]*resultant[a][b];
	return 0.5*result;
}

inline NativeKirchhoffLovePointResponse EvaluateNativeKirchhoffLovePoint(
	const NativeNurbsShellSurfacePoint& reference,
	const NativeNurbsShellSurfacePoint& current,
	const NativeKirchhoffLoveMaterial& material)
{
	ValidateNativeKirchhoffLoveMaterial(material);
	const auto reference_metric=NativeShellMetric(reference);
	const auto current_metric=NativeShellMetric(current);
	const auto inverse_metric=NativeShellInverseMetric(reference_metric);
	const auto reference_curvature=NativeShellCurvature(reference);
	const auto current_curvature=NativeShellCurvature(current);
	NativeKirchhoffLovePointResponse result;
	for(int a=0;a<2;++a) for(int b=0;b<2;++b) {
		result.membrane_strain[a][b]=0.5*(current_metric[a][b]-reference_metric[a][b]);
		result.curvature_change[a][b]=current_curvature[a][b]-reference_curvature[a][b];
	}
	const double lambda=material.young_modulus_pa*material.poisson_ratio/
		(1.0-material.poisson_ratio*material.poisson_ratio);
	const double mu=material.young_modulus_pa/(2.0*(1.0+material.poisson_ratio));
	result.membrane_resultant_contravariant=NativeShellPlaneStressResultant(
		result.membrane_strain,inverse_metric,lambda,mu,material.thickness_m);
	result.bending_moment_contravariant=NativeShellPlaneStressResultant(
		result.curvature_change,inverse_metric,lambda,mu,
		material.thickness_m*material.thickness_m*material.thickness_m/12.0);
	result.membrane_energy_per_reference_area_j_m2=NativeShellEnergyDensity(
		result.membrane_strain,result.membrane_resultant_contravariant);
	result.bending_energy_per_reference_area_j_m2=NativeShellEnergyDensity(
		result.curvature_change,result.bending_moment_contravariant);
	result.mass_per_reference_area_kg_m2=material.density_kg_m3*material.thickness_m;
	return result;
}

} // namespace iga

#endif
