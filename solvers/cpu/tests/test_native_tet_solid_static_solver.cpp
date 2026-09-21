#include "NativeTetSolidStaticSolver.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

namespace {
iga::NativeTetMesh Mesh(){iga::NativeTetMesh m;m.points={{{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}}}};m.cells={iga::NativeTetCell{1,{{0,1,2,3}}}};return m;}
bool Close(double a,double b,double tolerance=2e-8){return std::abs(a-b)<=tolerance*std::max({1.0,std::abs(a),std::abs(b)});}
}

int main()
{
	const auto mesh=Mesh();const iga::NativeTetSolidMaterial material{1500.0,0.3,980.0};
	iga::NativeTetSolidMatrix3 target_f{{{{1.08,0,0}},{{0,0.96,0}},{{0,0,1.03}}}};
	std::array<double,12> target{};
	for(std::size_t node=0;node<4;++node)for(int i=0;i<3;++i){target[3*node+i]=-mesh.points[node][i];for(int j=0;j<3;++j)target[3*node+i]+=target_f[i][j]*mesh.points[node][j];}
	const auto target_element=iga::BuildNativeTetHyperelasticSolidElement(mesh,mesh.cells[0],target,material);
	const std::map<std::size_t,double> constraints{{0,0},{1,0},{2,0},{4,0},{5,0},{6,0},{8,0},{9,0},{10,0}};
	std::vector<double> external(12,0.0);
	for(std::size_t dof=0;dof<12;++dof)if(!constraints.count(dof))external[dof]=target_element.residual_n[dof];
	iga::NativeTetSolidStaticOptions options;options.relative_tolerance=1e-10;options.absolute_tolerance_n=1e-11;
	const auto result=iga::SolveNativeTetSolidStatic(mesh,material,external,constraints,{},options);
	assert(result.converged&&result.iterations>0&&result.iterations<10);
	for(std::size_t dof=0;dof<12;++dof)assert(Close(result.displacement_m[dof],target[dof]));
	assert(result.final_free_residual_n<1e-8&&result.minimum_deformation_jacobian>0.0);
	// reaction_n is internal-external on constrained DOFs, and approximately zero on free DOFs.
	for(int component=0;component<3;++component){
		double external_total=0.0,reaction_total=0.0;
		for(std::size_t node=0;node<4;++node){external_total+=external[3*node+component];if(constraints.count(3*node+component))reaction_total+=result.reaction_n[3*node+component];}
		assert(Close(reaction_total+external_total,0.0));
	}
	std::cout<<"native tetrahedral solid static solver test passed iterations="<<result.iterations
		<<" residual="<<result.final_free_residual_n<<'\n';
	return 0;
}
