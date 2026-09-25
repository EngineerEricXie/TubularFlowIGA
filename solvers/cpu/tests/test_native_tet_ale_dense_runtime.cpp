#include "NativeTetAleDenseRuntime.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

int main()
{
	iga::NativeTetMesh mesh;
	mesh.points={{{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}},{{0.25,0.25,0.25}}}};
	mesh.cells={iga::NativeTetCell{1,{{4,1,2,3}}},iga::NativeTetCell{2,{{0,4,2,3}}},
		iga::NativeTetCell{3,{{0,1,4,3}}},iga::NativeTetCell{4,{{0,1,2,4}}}};
	mesh.boundary_triangles={iga::NativeTetTriangle{1,{{1,2,3}},0},
		iga::NativeTetTriangle{2,{{0,3,2}},0},iga::NativeTetTriangle{3,{{0,1,3}},0},
		iga::NativeTetTriangle{4,{{0,2,1}},0}};
	const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
	const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
	const std::size_t dofs=3*velocity_nodes+mesh.points.size();
	const std::array<double,3> translation{{0.2,-0.1,0.04}};
	const double dt=0.05;
	for(auto& point:mesh.points)
		for(int component=0;component<3;++component)
			point[component]+=dt*translation[component];
	std::vector<std::array<double,3>> mesh_velocity(mesh.points.size(),translation);
	std::vector<double> committed(dofs,0.0),trial(dofs,0.0);
	for(std::size_t node=0;node<velocity_nodes;++node)
		for(int component=0;component<3;++component) {
			committed[3*node+component]=translation[component];
			trial[3*node+component]=translation[component];
		}
	for(std::size_t node=0;node<velocity_nodes;++node)
		if(topology.boundary_velocity_nodes.at(0).count(static_cast<std::uint32_t>(node))==0)
			trial[3*node]+=0.03;
	std::map<std::uint32_t,std::array<double,3>> boundary_velocity;
	for(const auto node:topology.boundary_velocity_nodes.at(0))
		boundary_velocity[node]=translation;
	const auto result=iga::SolveNativeTetAleDenseTransient(mesh,mesh_velocity,committed,
		trial,boundary_velocity,0,{1000.0,0.004},dt,1e-10,8);
	assert(result.newton_iterations>0&&result.newton_iterations<=4);
	assert(result.final_free_residual_l2<1e-10);
	for(std::size_t node=0;node<velocity_nodes;++node)
		for(int component=0;component<3;++component)
			assert(std::abs(result.state[3*node+component]-translation[component])<2e-11);
	for(std::size_t node=0;node<mesh.points.size();++node)
		assert(std::abs(result.state[3*velocity_nodes+node])<2e-9);
	std::cout<<"native tetrahedral ALE dense global runtime test passed\n";
	return 0;
}
