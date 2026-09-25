#include "NativeIgaKirchhoffLove.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

namespace {

bool Near(double a,double b,double tolerance=3.0e-11)
{
	return std::abs(a-b)<=tolerance*std::max({1.0,std::abs(a),std::abs(b)});
}

iga::NativeNurbsShellPatch Plane()
{
	iga::NativeNurbsShellPatch patch;
	patch.degree_u=patch.degree_v=2;patch.knots_u=patch.knots_v={0,0,0,1,1,1};
	patch.controls_u=patch.controls_v=3;
	for(std::size_t j=0;j<3;++j) for(std::size_t i=0;i<3;++i) {
		patch.control_points.push_back({{0.5*i,0.5*j,0}});patch.weights.push_back(1);
	}
	return patch;
}

}

int main()
{
	const iga::NativeKirchhoffLoveMaterial material{2.0e6,0.3,0.01,1000.0};
	const auto reference_patch=Plane();
	const auto reference=iga::EvaluateNativeNurbsShellSurface(reference_patch,0.0,0.37);
	for(const double angle:{0.43,2.40}) {
		auto rigid_patch=reference_patch;
		const double c=std::cos(angle),s=std::sin(angle);
		for(auto& point:rigid_patch.control_points) {
			const auto original=point;
			point={{c*original[0]-s*original[1]+2.0,
				s*original[0]+c*original[1]-3.0,original[2]+0.7}};
		}
		const auto rigid=iga::EvaluateNativeKirchhoffLovePoint(reference,
			iga::EvaluateNativeNurbsShellSurface(rigid_patch,0.0,0.37),material);
		assert(Near(rigid.membrane_energy_per_reference_area_j_m2,0.0,1e-9));
		assert(Near(rigid.bending_energy_per_reference_area_j_m2,0.0,1e-9));
	}
	auto stretched_patch=reference_patch;
	const double stretch=1.02;
	for(auto& point:stretched_patch.control_points) point[0]*=stretch;
	const auto stretched=iga::EvaluateNativeKirchhoffLovePoint(reference,
		iga::EvaluateNativeNurbsShellSurface(stretched_patch,0.0,0.37),material);
	const double strain=0.5*(stretch*stretch-1.0);
	const double lambda=material.young_modulus_pa*material.poisson_ratio/
		(1.0-material.poisson_ratio*material.poisson_ratio);
	const double mu=material.young_modulus_pa/(2.0*(1.0+material.poisson_ratio));
	const double expected_membrane=0.5*material.thickness_m*(lambda+2.0*mu)*strain*strain;
	assert(Near(stretched.membrane_strain[0][0],strain));
	assert(Near(stretched.membrane_energy_per_reference_area_j_m2,expected_membrane));
	assert(Near(stretched.bending_energy_per_reference_area_j_m2,0.0));
	auto bent_patch=reference_patch;
	const double curvature=0.08;
	for(std::size_t j=0;j<3;++j) {
		bent_patch.control_points[j*3+0][2]=0.0;
		bent_patch.control_points[j*3+1][2]=0.0;
		bent_patch.control_points[j*3+2][2]=0.5*curvature;
	}
	const auto bent=iga::EvaluateNativeKirchhoffLovePoint(reference,
		iga::EvaluateNativeNurbsShellSurface(bent_patch,0.0,0.37),material);
	const double rigidity_scale=material.thickness_m*material.thickness_m
		*material.thickness_m/12.0;
	const double expected_bending=0.5*rigidity_scale*(lambda+2.0*mu)
		*curvature*curvature;
	assert(Near(bent.membrane_energy_per_reference_area_j_m2,0.0));
	assert(Near(bent.curvature_change[0][0],curvature));
	assert(Near(bent.bending_energy_per_reference_area_j_m2,expected_bending));
	assert(Near(bent.mass_per_reference_area_kg_m2,10.0));
	auto rotated_bent=bent_patch;
	const double large_angle=2.10,large_c=std::cos(large_angle),large_s=std::sin(large_angle);
	for(auto& point:rotated_bent.control_points) {
		const auto original=point;
		point={{large_c*original[0]+large_s*original[2]+1.7,original[1]-0.8,
			-large_s*original[0]+large_c*original[2]+2.4}};
	}
	const auto objective_bent=iga::EvaluateNativeKirchhoffLovePoint(reference,
		iga::EvaluateNativeNurbsShellSurface(rotated_bent,0.0,0.37),material);
	assert(Near(objective_bent.membrane_energy_per_reference_area_j_m2,
		bent.membrane_energy_per_reference_area_j_m2,2e-9));
	assert(Near(objective_bent.bending_energy_per_reference_area_j_m2,
		bent.bending_energy_per_reference_area_j_m2,2e-9));
	auto invalid=material;invalid.poisson_ratio=0.5;
	try { iga::EvaluateNativeKirchhoffLovePoint(reference,reference,invalid);assert(false); }
	catch(const std::invalid_argument&) {}
	std::cout<<"native IGA Kirchhoff-Love point tests passed\n";
	return 0;
}
