#include "NativeTetAleDenseRuntime.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <vector>

namespace {
iga::NativeTetMesh Mesh()
{
	iga::NativeTetMesh mesh;
	mesh.points={{{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}},{{.25,.25,.25}}}};
	mesh.cells={iga::NativeTetCell{1,{{4,1,2,3}}},iga::NativeTetCell{2,{{0,4,2,3}}},
		iga::NativeTetCell{3,{{0,1,4,3}}},iga::NativeTetCell{4,{{0,1,2,4}}}};
	mesh.boundary_triangles={iga::NativeTetTriangle{1,{{1,2,3}},0},
		iga::NativeTetTriangle{2,{{0,3,2}},0},iga::NativeTetTriangle{3,{{0,1,3}},0},
		iga::NativeTetTriangle{4,{{0,2,1}},0}};return mesh;
}
double Velocity(double time){return .01*std::sin(2.*time);}
double Run(double dt)
{
	const auto mesh=Mesh();const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
	const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
	std::vector<double> committed(3*velocity_nodes+mesh.points.size(),0.0);
	const std::vector<std::array<double,3>> stationary_grid(mesh.points.size(),{{0,0,0}});
	const double final_time=.4;const std::size_t steps=static_cast<std::size_t>(std::llround(final_time/dt));
	for(std::size_t step=1;step<=steps;++step){
		const double value=Velocity(step*dt);std::map<std::uint32_t,std::array<double,3>> boundary;
		for(auto node:topology.boundary_velocity_nodes.at(0))boundary[node]={{value,0,0}};
		const auto result=iga::SolveNativeTetAleDenseTransient(mesh,stationary_grid,committed,
			committed,boundary,0,{1000.,.004},dt,1e-10,8);committed=result.state;
	}
	const double numerical_gradient=committed[3*velocity_nodes+1]-committed[3*velocity_nodes];
	const double exact_gradient=-1000.*.01*2.*std::cos(2.*final_time);
	return std::abs(numerical_gradient-exact_gradient);
}
}

int main()
{
	const std::array<double,4> dt{{.1,.05,.025,.0125}};std::array<double,4> error{};
	for(std::size_t level=0;level<dt.size();++level)error[level]=Run(dt[level]);
	for(std::size_t level=1;level<error.size();++level){
		const double order=std::log(error[level-1]/error[level])/std::log(2.);
		assert(error[level]<error[level-1]);assert(order>.85&&order<1.15);
	}
	std::cout<<"native fixed-domain backward-Euler temporal convergence passed errors="
		<<error[0]<<','<<error[1]<<','<<error[2]<<','<<error[3]<<'\n';return 0;
}
