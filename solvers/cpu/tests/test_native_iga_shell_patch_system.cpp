#include "NativeIgaShellPatchSystem.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>

namespace {

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

double Norm(const std::vector<double>& values)
{
	double sum=0.0;for(const double value:values) sum+=value*value;return std::sqrt(sum);
}

}

int main()
{
	const auto patch=Plane();
	const iga::NativeKirchhoffLoveMaterial material{2.0e5,0.3,0.01,1000.0};
	const std::size_t dofs=3*patch.control_points.size();
	std::vector<double> zero(dofs,0.0),translation(dofs,0.0),stretch(dofs,0.0);
	for(std::size_t i=0;i<patch.control_points.size();++i) {
		translation[3*i]=0.07;translation[3*i+1]=-0.04;translation[3*i+2]=0.02;
		stretch[3*i]=0.01*patch.control_points[i][0];
	}
	assert(iga::NativeIgaShellPatchEnergy(patch,zero,material)<1e-20);
	assert(iga::NativeIgaShellPatchEnergy(patch,translation,material)<1e-20);
	const auto residual=iga::NativeIgaShellPatchResidual(patch,stretch,material,2e-7);
	assert(Norm(residual)>0.0);
	const auto system=iga::BuildNativeIgaShellPatchSystem(patch,stretch,material,2e-7,2e-5);
	assert(system.strain_energy_j>0.0&&system.residual_n.size()==dofs);
	assert(system.tangent_n_m.size()==dofs*dofs&&system.consistent_mass_kg.size()==dofs*dofs);
	double maximum_asymmetry=0.0,maximum_tangent=0.0;
	for(std::size_t i=0;i<dofs;++i) for(std::size_t j=0;j<dofs;++j) {
		maximum_asymmetry=std::max(maximum_asymmetry,std::abs(
			system.tangent_n_m[i*dofs+j]-system.tangent_n_m[j*dofs+i]));
		maximum_tangent=std::max(maximum_tangent,std::abs(system.tangent_n_m[i*dofs+j]));
		assert(std::abs(system.consistent_mass_kg[i*dofs+j]
			-system.consistent_mass_kg[j*dofs+i])<1e-13);
	}
	assert(maximum_asymmetry<2e-4*std::max(1.0,maximum_tangent));
	std::vector<double> direction(dofs,0.0);
	for(std::size_t i=0;i<dofs;++i) direction[i]=std::sin(static_cast<double>(i+1));
	const double h=1e-5;
	std::vector<double> plus=stretch,minus=stretch;
	for(std::size_t i=0;i<dofs;++i) { plus[i]+=h*direction[i];minus[i]-=h*direction[i]; }
	const auto rp=iga::NativeIgaShellPatchResidual(patch,plus,material,2e-7);
	const auto rm=iga::NativeIgaShellPatchResidual(patch,minus,material,2e-7);
	double error=0.0,reference=0.0;
	for(std::size_t row=0;row<dofs;++row) {
		double tangent_action=0.0;
		for(std::size_t column=0;column<dofs;++column)
			tangent_action+=system.tangent_n_m[row*dofs+column]*direction[column];
		const double finite_difference=(rp[row]-rm[row])/(2*h);
		error+=(tangent_action-finite_difference)*(tangent_action-finite_difference);
		reference+=finite_difference*finite_difference;
	}
	assert(std::sqrt(error/reference)<1e-5);
	double total_scalar_mass=0.0;
	for(std::size_t i=0;i<patch.control_points.size();++i)
		for(std::size_t j=0;j<patch.control_points.size();++j)
			total_scalar_mass+=system.consistent_mass_kg[(3*i)*dofs+3*j];
	assert(std::abs(total_scalar_mass-10.0)<1e-11);
	std::cout<<"native IGA shell patch system tests passed\n";
	return 0;
}
