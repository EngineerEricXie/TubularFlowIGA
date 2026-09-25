#include "NativeIgaShellPatchInterface.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

namespace {

iga::NativeNurbsShellPatch Patch(double begin,double end,double stretch,double curvature)
{
	iga::NativeNurbsShellPatch patch;patch.degree_u=patch.degree_v=2;
	patch.knots_u=patch.knots_v={0,0,0,1,1,1};patch.controls_u=patch.controls_v=3;
	for(int j=0;j<3;++j) for(int i=0;i<3;++i) {
		const double fraction=0.5*i;
		const double x=begin+(end-begin)*fraction;
		const double z=0.5*curvature*x*x;
		patch.control_points.push_back({{stretch*x,0.5*j,z}});patch.weights.push_back(1.0);
	}
	// Exact quadratic Bezier middle coefficient for z(x)=k*x^2/2.
	for(int j=0;j<3;++j) {
		const double width=end-begin;
		patch.control_points[static_cast<std::size_t>(3*j+1)][2]=
			0.5*curvature*begin*begin+0.5*curvature*begin*width;
	}
	return patch;
}

}

int main()
{
	const iga::NativeKirchhoffLoveMaterial material{2e6,0.3,0.01,1000};
	const auto left_reference=Patch(0.0,0.5,1.0,0.0);
	const auto right_reference=Patch(0.5,1.0,1.0,0.0);
	const auto left_current=Patch(0.0,0.5,1.01,0.06);
	const auto right_current=Patch(0.5,1.0,1.01,0.06);
	const auto diagnostics=iga::CheckNativeIgaShellCompatibleUInterface(left_reference,
		right_reference,left_current,right_current,material,1e-11,1e-10,1e-6);
	assert(diagnostics.maximum_reference_gap_m<1e-13);
	assert(diagnostics.maximum_current_gap_m<1e-13);
	assert(diagnostics.maximum_membrane_resultant_jump_n_m<1e-8);
	assert(diagnostics.maximum_bending_moment_jump_n<1e-8);
	const auto constraints=iga::BuildNativeIgaShellCompatibleUConstraints(
		left_reference,right_reference);
	assert(constraints.size()==6*left_reference.controls_v);
	std::vector<double> affine(3*(left_reference.control_points.size()
		+right_reference.control_points.size()),0.0);
	std::size_t offset=0;
	for(const auto* patch:{&left_reference,&right_reference}) {
		for(const auto& point:patch->control_points) {
			affine[3*offset]=0.02*point[0]-0.01*point[1];
			affine[3*offset+1]=0.03*point[1];
			affine[3*offset+2]=0.04*point[0]+0.01*point[1];
			++offset;
		}
	}
	for(const auto& constraint:constraints)
		assert(std::abs(iga::EvaluateNativeIgaShellConstraint(constraint,affine))<1e-14);
	auto c1_broken=affine;
	const std::size_t right_offset=3*left_reference.control_points.size();
	c1_broken[right_offset+3]+=1e-3;
	double maximum_constraint_violation=0.0;
	for(const auto& constraint:constraints)
		maximum_constraint_violation=std::max(maximum_constraint_violation,
			std::abs(iga::EvaluateNativeIgaShellConstraint(constraint,c1_broken)));
	assert(maximum_constraint_violation>1e-3);
	const auto transformation=iga::BuildNativeIgaShellCompatibleUTransformation(
		left_reference,right_reference);
	assert(transformation.local_dofs==affine.size());
	assert(transformation.independent_dofs+6*right_reference.controls_v==affine.size());
	std::vector<double> independent(transformation.independent_dofs);
	for(std::size_t i=0;i<independent.size();++i) independent[i]=std::sin(0.37*(i+1));
	const auto constrained=iga::ApplyNativeIgaShellInterfaceTransformation(
		transformation,independent);
	for(const auto& constraint:constraints)
		assert(std::abs(iga::EvaluateNativeIgaShellConstraint(constraint,constrained))<1e-12);
	std::vector<double> local_residual(transformation.local_dofs);
	for(std::size_t i=0;i<local_residual.size();++i) local_residual[i]=std::cos(0.23*(i+1));
	const auto reduced_residual=iga::ReduceNativeIgaShellInterfaceResidual(
		transformation,local_residual);
	double local_work=0.0,reduced_work=0.0;
	for(std::size_t i=0;i<constrained.size();++i) local_work+=local_residual[i]*constrained[i];
	for(std::size_t i=0;i<independent.size();++i) reduced_work+=reduced_residual[i]*independent[i];
	assert(std::abs(local_work-reduced_work)<1e-11);
	std::vector<double> local_tangent(transformation.local_dofs*transformation.local_dofs,0.0);
	for(std::size_t i=0;i<transformation.local_dofs;++i) {
		local_tangent[i*transformation.local_dofs+i]=2.0;
		if(i+1<transformation.local_dofs) {
			local_tangent[i*transformation.local_dofs+i+1]=0.1;
			local_tangent[(i+1)*transformation.local_dofs+i]=0.1;
		}
	}
	const auto reduced_tangent=iga::ReduceNativeIgaShellInterfaceTangent(
		transformation,local_tangent);
	for(std::size_t i=0;i<transformation.independent_dofs;++i)
		for(std::size_t j=0;j<transformation.independent_dofs;++j)
			assert(std::abs(reduced_tangent[i*transformation.independent_dofs+j]
				-reduced_tangent[j*transformation.independent_dofs+i])<1e-12);
	auto gapped=right_current;for(auto& point:gapped.control_points) point[0]+=1e-4;
	try {
		iga::CheckNativeIgaShellCompatibleUInterface(left_reference,right_reference,
			left_current,gapped,material);assert(false);
	} catch(const std::runtime_error&) {}
	auto kinked=right_current;
	for(int j=0;j<3;++j) kinked.control_points[static_cast<std::size_t>(3*j+1)][2]+=0.02;
	try {
		iga::CheckNativeIgaShellCompatibleUInterface(left_reference,right_reference,
			left_current,kinked,material);assert(false);
	} catch(const std::runtime_error&) {}
	std::cout<<"native IGA shell compatible-interface diagnostics passed\n";
	return 0;
}
