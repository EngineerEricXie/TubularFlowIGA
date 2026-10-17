#include "NativeTetAlePetscRuntime.hpp"

#include <petscksp.h>

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
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
		// Flow-rate controllers do not support -native_ale_scale yet, and their exact
		// flow and pressure-shift checks assume a direct solve; the iterative
		// fieldsplit run turns them off with -native_ale_test_controllers false.
		PetscBool controllers=PETSC_TRUE;
		PetscOptionsGetBool(nullptr,nullptr,"-native_ale_test_controllers",&controllers,nullptr);
		std::ostringstream controller_summary;
		if(controllers){
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
			controller_summary<<" controlled_flow="<<controlled_flow
				<<" controller_pa="<<port_result.flow_controller_multipliers_pa[0]
				<<" pressure_shift_controller_pa="<<shifted.flow_controller_multipliers_pa[0];
			// Implicit resistance ports: label 1 keeps a pressure-dependent flow
			// control, labels 2 and 3 obey pi=offset+R*(relative outward flow).
			// Prescribing the returned pressures must reproduce the same state.
			auto split_mesh=port_mesh;
			split_mesh.boundary_triangles[1].boundary_label=2;
			for(std::size_t triangle=2;triangle<split_mesh.boundary_triangles.size();++triangle)
				split_mesh.boundary_triangles[triangle].boundary_label=3;
			const auto split_operators=iga::BuildNativeTetAleBoundaryFluxOperators(split_mesh,
				iga::BuildNativeTaylorHoodTopology(split_mesh));
			double weight_sum=0.0;
			for(const double weight:split_operators.at(3).pressure_weights)weight_sum+=weight;
			if(std::abs(weight_sum-1.0)>1e-14)
				throw std::runtime_error("native ALE mean-pressure weights do not sum to one");
			const auto mesh_flow=[&](int label){return split_operators.at(label).OutwardFlowM3S(committed);};
			const auto mean_pressure=[&](const std::vector<double>& state,int label){
				const auto& operator_value=split_operators.at(label);
				double value=0.0;
				for(std::size_t entry=0;entry<operator_value.pressure_nodes.size();++entry)
					value+=operator_value.pressure_weights[entry]
						*state[3*velocity_nodes+operator_value.pressure_nodes[entry]];
				return value;
			};
			iga::NativeTetAleBoundaryConditions implicit;
			iga::NativeTetAleFlowRateControl inflow;
			inflow.boundary_label=1;inflow.target_outward_flow_m3_s=mesh_flow(1)-2e-3;
			inflow.target_pressure_slope_m3_s_pa=1e-4;
			implicit.flow_rate_controls.push_back(inflow);
			implicit.resistance_controls.push_back({2,20.0,1.0,mesh_flow(2)});
			implicit.resistance_controls.push_back({3,60.0,-2.0,mesh_flow(3)});
			iga::NativeTetAlePetscSolveControls tight;
			tight.control_absolute_tolerance_m3_s=1e-15;tight.scaled_line_search=true;
			const auto implicit_result=iga::SolveNativeTetAlePetscTransient(split_mesh,grid,
				committed,committed,{},std::numeric_limits<std::uint32_t>::max(),
				{1000.0,0.004},dt,1e-10,8,implicit,tight);
			if(implicit_result.resistance_pressures_pa.size()!=2
				||implicit_result.flow_controller_multipliers_pa.size()!=1)
				throw std::runtime_error("native ALE resistance controls returned the wrong unknowns");
			iga::NativeTetAleBoundaryConditions explicit_conditions;
			for(std::size_t control=0;control<2;++control){
				const auto& resistance=implicit.resistance_controls[control];
				const double relative=split_operators.at(resistance.boundary_label)
					.OutwardFlowM3S(implicit_result.replicated_state)-resistance.mesh_flow_m3_s;
				const double pressure=implicit_result.resistance_pressures_pa[control];
				if(std::abs(pressure-resistance.offset_pa-resistance.resistance_pa_s_m3*relative)
					>1e-10*std::max(1.0,std::abs(pressure)))
					throw std::runtime_error("native ALE resistance relation is not satisfied");
				explicit_conditions.prescribed_pressure_pa[resistance.boundary_label]=pressure;
			}
			const double inlet_pressure=mean_pressure(implicit_result.replicated_state,1);
			auto fixed_inflow=inflow;
			fixed_inflow.target_outward_flow_m3_s+=inflow.target_pressure_slope_m3_s_pa*inlet_pressure;
			fixed_inflow.target_pressure_slope_m3_s_pa=0.0;
			explicit_conditions.flow_rate_controls.push_back(fixed_inflow);
			iga::NativeTetAlePetscSolveControls warm;
			warm.initial_control_values=implicit_result.flow_controller_multipliers_pa;
			const auto explicit_result=iga::SolveNativeTetAlePetscTransient(split_mesh,grid,
				committed,implicit_result.replicated_state,{},std::numeric_limits<std::uint32_t>::max(),
				{1000.0,0.004},dt,1e-10,8,explicit_conditions,warm);
			if(explicit_result.newton_iterations!=0)
				throw std::runtime_error("native ALE implicit state is not a fixed-pressure solution");
			const auto cold=iga::SolveNativeTetAlePetscTransient(split_mesh,grid,
				committed,committed,{},std::numeric_limits<std::uint32_t>::max(),
				{1000.0,0.004},dt,1e-10,8,explicit_conditions);
			for(std::size_t dof=0;dof<dofs;++dof)
				if(std::abs(cold.replicated_state[dof]-implicit_result.replicated_state[dof])
					>1e-9*std::max(1.0,std::abs(implicit_result.replicated_state[dof])))
					throw std::runtime_error("native ALE implicit and prescribed-pressure states differ");
			controller_summary<<" resistance_pa="<<implicit_result.resistance_pressures_pa[0]
				<<","<<implicit_result.resistance_pressures_pa[1];
		}
		if(rank==0) std::cout<<"native tetrahedral ALE PETSc runtime passed ranks="<<ranks
			<<" transient_newton="<<result.newton_iterations
			<<" transient_residual="<<result.final_residual_l2
			<<" steady_newton="<<steady.newton_iterations
			<<" steady_residual="<<steady.final_residual_l2
			<<controller_summary.str()<<'\n';
	} catch(const std::exception& error) {
		int rank=0;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		if(rank==0) std::cerr<<error.what()<<'\n';
		exit_code=1;
	}
	PetscFinalize();return exit_code;
}
