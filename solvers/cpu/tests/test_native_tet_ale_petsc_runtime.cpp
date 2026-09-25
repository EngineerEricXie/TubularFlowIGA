#include "NativeTetAlePetscRuntime.hpp"

#include <petscksp.h>

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);int exit_code=0;
	try {
		int rank=0,ranks=1;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
		iga::NativeTetMesh mesh;
		mesh.points={{{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}},{{0.25,0.25,0.25}}}};
		mesh.cells={iga::NativeTetCell{1,{{4,1,2,3}}},
			iga::NativeTetCell{2,{{0,4,2,3}}},iga::NativeTetCell{3,{{0,1,4,3}}},
			iga::NativeTetCell{4,{{0,1,2,4}}}};
		mesh.boundary_triangles={iga::NativeTetTriangle{1,{{1,2,3}},0},
			iga::NativeTetTriangle{2,{{0,3,2}},0},iga::NativeTetTriangle{3,{{0,1,3}},0},
			iga::NativeTetTriangle{4,{{0,2,1}},0}};
		const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
		const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
		const std::size_t dofs=3*velocity_nodes+mesh.points.size();
		const std::array<double,3> velocity{{0.2,-0.1,0.04}};const double dt=0.05;
		for(auto& point:mesh.points)
			for(int component=0;component<3;++component) point[component]+=dt*velocity[component];
		std::vector<std::array<double,3>> grid(mesh.points.size(),velocity);
		std::vector<double> committed(dofs,0.0),trial(dofs,0.0);
		for(std::size_t node=0;node<velocity_nodes;++node)
			for(int component=0;component<3;++component)
				committed[3*node+component]=trial[3*node+component]=velocity[component];
		for(std::size_t node=0;node<velocity_nodes;++node)
			if(topology.boundary_velocity_nodes.at(0).count(static_cast<std::uint32_t>(node))==0)
				trial[3*node]+=0.03;
		std::map<std::uint32_t,std::array<double,3>> boundary;
		for(const auto node:topology.boundary_velocity_nodes.at(0)) boundary[node]=velocity;
		const auto result=iga::SolveNativeTetAlePetscTransient(mesh,grid,committed,trial,
			boundary,0,{1000.0,0.004},dt,1e-10,8);
		if(!(result.newton_iterations>0&&result.newton_iterations<=4
			&&result.final_residual_l2<1e-10))
			throw std::runtime_error("native ALE PETSc nonlinear result failed");
		for(std::size_t node=0;node<velocity_nodes;++node)
			for(int component=0;component<3;++component)
				if(std::abs(result.replicated_state[3*node+component]-velocity[component])>2e-11)
					throw std::runtime_error("native ALE PETSc velocity differs from exact translation");
		auto steady_trial=committed;
		for(std::size_t node=0;node<velocity_nodes;++node)
			if(topology.boundary_velocity_nodes.at(0).count(static_cast<std::uint32_t>(node))==0)
				steady_trial[3*node]+=1e-4;
		const auto steady=iga::SolveNativeTetAlePetscSteady(mesh,grid,steady_trial,boundary,0,
			{1000.0,0.004},1e-10,8);
		if(!(steady.newton_iterations>0&&steady.newton_iterations<=4
			&&steady.final_residual_l2<1e-10))
			throw std::runtime_error("native steady PETSc nonlinear result failed");
		for(std::size_t node=0;node<velocity_nodes;++node)
			for(int component=0;component<3;++component)
				if(std::abs(steady.replicated_state[3*node+component]-velocity[component])>2e-11)
					throw std::runtime_error("native steady PETSc velocity differs from exact translation");
		auto port_mesh=mesh;
		port_mesh.boundary_triangles[0].boundary_label=1;
		for(std::size_t triangle=1;triangle<port_mesh.boundary_triangles.size();++triangle)
			port_mesh.boundary_triangles[triangle].boundary_label=2;
		const auto port_topology=iga::BuildNativeTaylorHoodTopology(port_mesh);
		const auto port_operators=iga::BuildNativeTetAleBoundaryFluxOperators(
			port_mesh,port_topology);
		iga::NativeTetAleBoundaryConditions port_conditions;
		port_conditions.prescribed_pressure_pa[2]=0.0;
		port_conditions.flow_rate_controls.push_back({1,
			port_operators.at(1).OutwardFlowM3S(committed)});
		auto port_trial=committed;
		for(std::size_t node=0;node<velocity_nodes;++node)
			port_trial[3*node]+=(node%2?1.0:-1.0)*0.003;
		const auto port_result=iga::SolveNativeTetAlePetscTransient(port_mesh,grid,
			committed,port_trial,{},std::numeric_limits<std::uint32_t>::max(),
			{1000.0,0.004},dt,1e-10,8,port_conditions);
		const double controlled_flow=port_operators.at(1).OutwardFlowM3S(
			port_result.replicated_state);
		if(!(port_result.newton_iterations>0&&port_result.final_residual_l2<1e-10
			&&std::abs(controlled_flow-port_conditions.flow_rate_controls[0]
				.target_outward_flow_m3_s)<1e-13
			&&port_result.flow_controller_multipliers_pa.size()==1
			&&std::isfinite(port_result.flow_controller_multipliers_pa[0])))
			throw std::runtime_error("native ALE PETSc pressure/flow boundary solve failed");
		auto pressurized=port_conditions;
		pressurized.prescribed_pressure_pa[2]=5.0;
		auto pressurized_trial=port_result.replicated_state;
		for(std::size_t node=0;node<port_mesh.points.size();++node)
			pressurized_trial[3*velocity_nodes+node]+=5.0;
		const auto shifted=iga::SolveNativeTetAlePetscTransient(port_mesh,grid,
			committed,pressurized_trial,{},std::numeric_limits<std::uint32_t>::max(),
			{1000.0,0.004},dt,1e-10,8,pressurized);
		if(!(shifted.final_residual_l2<1e-10
			&&std::abs(port_operators.at(1).OutwardFlowM3S(shifted.replicated_state)
				-pressurized.flow_rate_controls[0].target_outward_flow_m3_s)<1e-13
			&&std::abs(shifted.flow_controller_multipliers_pa.at(0)
				-port_result.flow_controller_multipliers_pa.at(0)-5.0)<1e-8))
			throw std::runtime_error("native ALE natural-pressure shift invariance failed");
		for(std::size_t dof=0;dof<3*velocity_nodes;++dof)
			if(std::abs(shifted.replicated_state[dof]-port_result.replicated_state[dof])>1e-8)
				throw std::runtime_error("native ALE outlet pressure shift changed velocity");
		if(rank==0) std::cout<<"native tetrahedral ALE PETSc runtime passed ranks="<<ranks
			<<" transient_newton="<<result.newton_iterations
			<<" transient_residual="<<result.final_residual_l2
			<<" steady_newton="<<steady.newton_iterations
			<<" steady_residual="<<steady.final_residual_l2
			<<" controlled_flow="<<controlled_flow
			<<" controller_pa="<<port_result.flow_controller_multipliers_pa[0]
			<<" pressure_shift_controller_pa="
			<<shifted.flow_controller_multipliers_pa[0]<<'\n';
	} catch(const std::exception& error) {
		int rank=0;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		if(rank==0) std::cerr<<error.what()<<'\n';
		exit_code=1;
	}
	PetscFinalize();return exit_code;
}
