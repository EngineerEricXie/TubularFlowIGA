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
	mesh.points={{{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}},{{0.25,0.25,0.25}}}};
	mesh.cells={iga::NativeTetCell{1,{{4,1,2,3}}},iga::NativeTetCell{2,{{0,4,2,3}}},
		iga::NativeTetCell{3,{{0,1,4,3}}},iga::NativeTetCell{4,{{0,1,2,4}}}};
	mesh.boundary_triangles={iga::NativeTetTriangle{1,{{1,2,3}},0},
		iga::NativeTetTriangle{2,{{0,3,2}},0},iga::NativeTetTriangle{3,{{0,1,3}},0},
		iga::NativeTetTriangle{4,{{0,2,1}},0}};
	return mesh;
}

double Position(double time)
{
	return 0.01*std::sin(2.0*time);
}

double Run(double dt)
{
	auto reference=Mesh();
	const auto topology=iga::BuildNativeTaylorHoodTopology(reference);
	const std::size_t velocity_nodes=reference.points.size()+topology.edges.size();
	const std::size_t dofs=3*velocity_nodes+reference.points.size();
	std::vector<double> committed(dofs,0.0);
	double previous_velocity=(Position(0.0)-Position(-dt))/dt;
	for(std::size_t node=0;node<velocity_nodes;++node) committed[3*node]=previous_velocity;
	const double final_time=0.4;
	const std::size_t steps=static_cast<std::size_t>(std::llround(final_time/dt));
	for(std::size_t step=1;step<=steps;++step) {
		const double time=step*dt;
		auto current=reference;
		for(auto& point:current.points) point[0]+=Position(time);
		const double velocity=(Position(time)-Position(time-dt))/dt;
		const std::array<double,3> value{{velocity,0,0}};
		const std::vector<std::array<double,3>> grid(current.points.size(),value);
		std::map<std::uint32_t,std::array<double,3>> boundary;
		for(const auto node:topology.boundary_velocity_nodes.at(0)) boundary[node]=value;
		const auto result=iga::SolveNativeTetAleDenseTransient(current,grid,committed,
			committed,boundary,0,{1000.0,0.004},dt,1e-10,6);
		committed=result.state;previous_velocity=velocity;
	}
	(void)previous_velocity;
	const double numerical_gradient=committed[3*velocity_nodes+1]
		-committed[3*velocity_nodes+0];
	const double exact_gradient=1000.0*0.01*4.0*std::sin(2.0*final_time);
	return std::abs(numerical_gradient-exact_gradient);
}

} // namespace

int main()
{
	const std::array<double,4> dt{{0.1,0.05,0.025,0.0125}};
	std::array<double,4> error{};
	for(std::size_t level=0;level<dt.size();++level) error[level]=Run(dt[level]);
	for(std::size_t level=1;level<error.size();++level) {
		const double order=std::log(error[level-1]/error[level])/std::log(2.0);
		assert(error[level]<error[level-1]);
		assert(order>0.85&&order<1.15);
	}
	std::cout<<"native tetrahedral ALE backward-Euler temporal convergence passed errors="
		<<error[0]<<','<<error[1]<<','<<error[2]<<','<<error[3]<<'\n';
	return 0;
}
