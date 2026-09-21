#include "NativeIgaShellStaticSolver.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <map>

int main()
{
	iga::NativeNurbsShellPatch patch;
	patch.degree_u=patch.degree_v=2;patch.knots_u=patch.knots_v={0,0,0,1,1,1};
	patch.controls_u=patch.controls_v=3;
	for(std::size_t j=0;j<3;++j) for(std::size_t i=0;i<3;++i) {
		patch.control_points.push_back({{0.5*i,0.5*j,0}});patch.weights.push_back(1);
	}
	const iga::NativeKirchhoffLoveMaterial material{2e5,0.3,0.01,1000};
	const auto force=iga::NativeIgaShellPatchDeadTraction(patch,{{20.0,0,0}});
	double total_force=0.0;for(std::size_t i=0;i<patch.control_points.size();++i) total_force+=force[3*i];
	assert(std::abs(total_force-20.0)<1e-11);
	std::map<std::size_t,double> constraints;
	for(std::size_t node=0;node<patch.control_points.size();++node) {
		constraints[3*node+1]=0.0;constraints[3*node+2]=0.0;
		if(node%3==0) constraints[3*node]=0.0;
	}
	iga::NativeIgaShellStaticOptions options;options.residual_tolerance_n=2e-4;
	const auto result=iga::SolveNativeIgaShellStaticPatch(patch,material,force,constraints,options);
	assert(result.converged&&result.newton_iterations>=1);
	assert(result.free_residual_norm_n<=options.residual_tolerance_n);
	for(std::size_t node=0;node<patch.control_points.size();++node) {
		assert(result.displacement_m[3*node+1]==0.0&&result.displacement_m[3*node+2]==0.0);
		if(node%3==0) assert(result.displacement_m[3*node]==0.0);
	}
	assert(result.displacement_m[3*2]>0.0&&result.displacement_m[3*5]>0.0
		&&result.displacement_m[3*8]>0.0);
	double reaction=0.0;for(std::size_t node=0;node<patch.control_points.size();node+=3)
		reaction+=result.reaction_n[3*node];
	assert(std::abs(reaction+20.0)<2e-3);
	std::cout<<"native IGA shell static solver tests passed iterations="
		<<result.newton_iterations<<" residual="<<result.free_residual_norm_n<<'\n';
	return 0;
}
