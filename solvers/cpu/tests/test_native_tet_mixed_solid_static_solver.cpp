#include "NativeTetMixedSolidStaticSolver.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

namespace {
iga::NativeTetMesh Mesh(){iga::NativeTetMesh mesh;mesh.points={{{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}}}};
	mesh.cells={iga::NativeTetCell{1,{{0,1,2,3}}}};return mesh;}
}

int main()
{
	const auto mesh=Mesh();const auto material=iga::NativeTetMixedSolidMaterialFromYoungPoisson(1500.,.4999,980.);
	const double axial=1.04,lateral=1./std::sqrt(axial);
	const iga::NativeTetSolidMatrix3 target_f{{{{axial,0,0}},{{0,lateral,0}},{{0,0,lateral}}}};
	std::array<double,12> target{};
	for(std::size_t node=0;node<4;++node)for(int i=0;i<3;++i){target[3*node+i]=-mesh.points[node][i];
		for(int j=0;j<3;++j)target[3*node+i]+=target_f[i][j]*mesh.points[node][j];}
	std::array<double,4> pressure{};const auto target_element=iga::BuildNativeTetMixedHyperelasticSolidElement(
		mesh,mesh.cells[0],target,pressure,material);
	const std::map<std::size_t,double> constraints{{0,0},{1,0},{2,0},{4,0},{5,0},{6,0},{8,0},{9,0},{10,0}};
	std::vector<double> external(12,0.);for(std::size_t dof=0;dof<12;++dof)
		if(!constraints.count(dof))external[dof]=target_element.residual[dof];
	iga::NativeTetSolidStaticOptions options;options.relative_tolerance=1e-10;
	options.absolute_tolerance_n=1e-11;
	const auto result=iga::SolveNativeTetMixedSolidStatic(mesh,material,external,constraints,{}, {},options);
	assert(result.converged&&result.iterations>0&&result.minimum_deformation_jacobian>0.);
	for(std::size_t dof=0;dof<12;++dof)
		assert(std::abs(result.displacement_m[dof]-target[dof])<2e-8);
	for(double value:result.pressure_pa)assert(std::abs(value)<2e-7);
	assert(result.final_free_residual<1e-8);
	std::cout<<"native mixed tetrahedral solid static solver passed iterations="<<result.iterations
		<<" residual="<<result.final_free_residual<<'\n';return 0;
}
