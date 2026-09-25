#include "NativeTetAleConservation.hpp"
#include "NativeTetAleDenseRuntime.hpp"
#include "NativeTetAleKinematics.hpp"
#include "NativeTetAleMeshMotion.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <vector>

int main()
{
	iga::NativeTetMesh reference;
	reference.points={{{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}},{{0.25,0.25,0.25}}}};
	reference.cells={iga::NativeTetCell{1,{{4,1,2,3}}},
		iga::NativeTetCell{2,{{0,4,2,3}}},iga::NativeTetCell{3,{{0,1,4,3}}},
		iga::NativeTetCell{4,{{0,1,2,4}}}};
	reference.boundary_triangles={iga::NativeTetTriangle{1,{{1,2,3}},0},
		iga::NativeTetTriangle{2,{{0,3,2}},0},iga::NativeTetTriangle{3,{{0,1,3}},0},
		iga::NativeTetTriangle{4,{{0,2,1}},0}};
	const auto topology=iga::BuildNativeTaylorHoodTopology(reference);
	const std::size_t velocity_nodes=reference.points.size()+topology.edges.size();
	const std::size_t dofs=3*velocity_nodes+reference.points.size();
	const std::array<double,3> wall_velocity{{0.12,-0.03,0.05}};
	std::vector<double> committed_fluid(dofs,0.0);
	for(std::size_t node=0;node<velocity_nodes;++node)
		for(int component=0;component<3;++component)
			committed_fluid[3*node+component]=wall_velocity[component];
	std::map<std::uint32_t,std::array<double,3>> fluid_boundary;
	for(const auto node:topology.boundary_velocity_nodes.at(0))
		fluid_boundary[node]=wall_velocity;
	iga::NativeTetAleKinematics kinematics(reference.points,reference.cells);
	auto previous_mesh=reference;
	const double dt=0.04;
	for(std::size_t step=1;step<=5;++step) {
		const double next_time=step*dt;
		std::map<std::uint32_t,std::array<double,3>> wall_displacement;
		for(std::uint32_t node=0;node<4;++node)
			for(int component=0;component<3;++component)
				wall_displacement[node][component]=next_time*wall_velocity[component];
		const auto harmonic=iga::SolveNativeTetAleHarmonicMotion(reference,wall_displacement);
		assert(harmonic.maximum_free_residual<1e-14);
		const auto& trial_geometry=kinematics.BeginTrial(harmonic.displacement_m,next_time);
		const auto current_mesh=iga::BuildNativeTetAleCurrentMesh(reference,trial_geometry);
		auto fluid_trial=committed_fluid;
		for(std::size_t node=0;node<velocity_nodes;++node)
			if(topology.boundary_velocity_nodes.at(0).count(static_cast<std::uint32_t>(node))==0)
				fluid_trial[3*node]+=0.002;
		const auto fluid=iga::SolveNativeTetAleDenseTransient(current_mesh,
			trial_geometry.mesh_velocity_m_s,committed_fluid,fluid_trial,fluid_boundary,
			0,{1000.0,0.004},dt,1e-10,8);
		assert(fluid.newton_iterations>0&&fluid.final_free_residual_l2<1e-10);
		for(std::size_t node=0;node<velocity_nodes;++node)
			for(int component=0;component<3;++component)
				assert(std::abs(fluid.state[3*node+component]-wall_velocity[component])<2e-11);
		const std::vector<std::array<double,3>> nodal_fluid(reference.points.size(),wall_velocity);
		const auto conservation=iga::EvaluateNativeTetAleConservation(previous_mesh,current_mesh,
			trial_geometry.mesh_velocity_m_s,nodal_fluid,nodal_fluid,dt);
		assert(std::abs(conservation.gcl_residual_m3_s)<1e-13);
		assert(std::abs(conservation.moving_domain_balance_residual_m3_s)<1e-13);
		committed_fluid=fluid.state;
		kinematics.CommitTrial();
		previous_mesh=current_mesh;
		assert(std::abs(kinematics.CommittedTime()-next_time)<1e-14);
	}
	std::cout<<"native tetrahedral ALE prescribed-motion multi-step runtime test passed\n";
	return 0;
}
