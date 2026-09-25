#include "NativeTetAleGraphPorts.hpp"
#include "NativeTetAlePetscRuntime.hpp"
#include "NativeTetMovingSpeciesPetscRuntime.hpp"
#include "NativeTetMovingSpeciesPorts.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

iga::CouplingPort Port(const char* id,int label,iga::PortQuantity input)
{
	iga::CouplingPort port;
	port.id=id;port.subsystem_id="tet";
	port.locator_kind="boundary_label";port.locator=std::to_string(label);
	port.provides={iga::PortQuantity::Area,iga::PortQuantity::FlowRate,
		iga::PortQuantity::MeanPressure,
		iga::PortQuantity::SpeciesConcentration,
		iga::PortQuantity::SpeciesFlux};
	port.requires={input,iga::PortQuantity::SpeciesConcentration,
		iga::PortQuantity::SpeciesFlux};
	port.species={"tracer"};return port;
}

void Close(double actual,double expected,double tolerance=1e-9)
{
	if(!std::isfinite(actual)
		||std::abs(actual-expected)>tolerance*std::max(1.,std::abs(expected)))
		throw std::runtime_error("native ALE graph PETSc port value differs");
}

} // namespace

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int exit_code=0;
	try{
		int rank=0,ranks=1;
		MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
		iga::NativeTetMesh mesh;
		mesh.points={{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}},{{0.25,0.25,0.25}}};
		mesh.cells={{1,{{4,1,2,3}}},{2,{{0,4,2,3}}},
			{3,{{0,1,4,3}}},{4,{{0,1,2,4}}}};
		mesh.boundary_triangles={{1,{{1,2,3}},1},
			{2,{{0,3,2}},2},{3,{{0,1,3}},2},{4,{{0,2,1}},2}};
		const std::vector<iga::CouplingPort> ports{
			Port("controlled",1,iga::PortQuantity::FlowRate),
			Port("pressure",2,iga::PortQuantity::MeanPressure)};
		const iga::NativeTetAleGraphPorts bridge("tet",mesh,ports);
		const iga::NativeTetMovingSpeciesPorts species_ports(
			"tet","tracer",mesh,ports);
		const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
		const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
		const std::size_t dofs=3*velocity_nodes+mesh.points.size();
		const std::array<double,3> velocity{{0.2,0.,0.}};
		const std::vector<std::array<double,3>> grid(mesh.points.size(),{{0,0,0}});
		std::vector<double> committed(dofs,0.);
		for(std::size_t node=0;node<velocity_nodes;++node)
			for(int component=0;component<3;++component)
				committed[3*node+component]=velocity[component];
		const double dt=0.05;
		const auto initial_ports=bridge.Observe(mesh,grid,committed,dt);
		const double target=*initial_ports.at("controlled").outward_flow_m3_s;
		iga::PortBoundaryData flow;flow.time_s=dt;flow.outward_flow_m3_s=target;
		iga::PortBoundaryData pressure;pressure.time_s=dt;pressure.mean_pressure_pa=0.;
		const auto conditions=bridge.BoundaryConditions(mesh,grid,
			{{"controlled",flow},{"pressure",pressure}},dt);
		auto trial=committed;
		for(std::size_t node=0;node<velocity_nodes;++node)
			trial[3*node]+=(node%2?1.:-1.)*0.003;
		const auto solved=iga::SolveNativeTetAlePetscTransient(mesh,grid,
			committed,trial,{},std::numeric_limits<std::uint32_t>::max(),
			{1000.,0.004},dt,1e-10,8,conditions);
		if(solved.newton_iterations==0||solved.final_residual_l2>=1e-10)
			throw std::runtime_error("native ALE graph PETSc flow solve did not converge");
		const auto measured=bridge.Observe(mesh,grid,solved.replicated_state,dt);
		Close(*measured.at("controlled").outward_flow_m3_s,target,1e-11);
		Close(*measured.at("pressure").outward_flow_m3_s,-target,1e-10);
		if(!measured.at("controlled").mean_pressure_pa
			||!measured.at("pressure").mean_pressure_pa)
			throw std::runtime_error("native ALE graph PETSc pressure observation is absent");
		const std::vector<double> initial_concentration(mesh.points.size(),2.);
		auto solve_species=[&](const iga::NativeTetMesh& previous_mesh,
			const iga::NativeTetMesh& current_mesh,
			const std::vector<std::array<double,3>>& mesh_velocity,
			const std::vector<double>& flow_state,
			const std::map<std::string,iga::PortState>& hydraulic){
			std::vector<std::array<double,3>> nodal_velocity(velocity_nodes);
			for(std::size_t node=0;node<velocity_nodes;++node)
				for(int axis=0;axis<3;++axis)
					nodal_velocity[node][axis]=flow_state[3*node+axis];
			iga::PortBoundaryData species_input;
			species_input.time_s=dt;
			species_input.concentration={{"tracer",2.}};
			const auto inflow=species_ports.InflowByLabel(
				{{"controlled",*hydraulic.at("controlled").outward_flow_m3_s},
					{"pressure",*hydraulic.at("pressure").outward_flow_m3_s}},
				{{"pressure",species_input}},dt,1e-10);
			const auto scalar=iga::SolveNativeTetMovingSpeciesPetscStep(
				previous_mesh,current_mesh,nodal_velocity,mesh_velocity,
				initial_concentration,inflow,0.,0.,dt);
			const auto scalar_ports=species_ports.Observe(current_mesh,scalar.step,dt,1e-10);
			for(const auto& port:ports){
				const double q=*hydraulic.at(port.id).outward_flow_m3_s;
				Close(*scalar_ports.at(port.id).outward_flow_m3_s,q,1e-10);
				Close(scalar_ports.at(port.id).outward_species_flux.at("tracer"),
					2.*q,1e-8);
			}
			const auto balance=species_ports.Accounting(scalar.step,dt);
			Close(balance.residual,0.,1e-10);
			return scalar;
		};
		const auto fixed_species=solve_species(mesh,mesh,grid,
			solved.replicated_state,measured);
		auto moved=mesh;
		const std::array<double,3> grid_velocity{{0.05,0.,0.}};
		for(auto& point:moved.points)
			for(int axis=0;axis<3;++axis)point[axis]+=dt*grid_velocity[axis];
		const std::vector<std::array<double,3>> moving_grid(
			moved.points.size(),grid_velocity);
		const auto moving_initial=bridge.Observe(moved,moving_grid,committed,dt);
		const double relative_target=
			*moving_initial.at("controlled").outward_flow_m3_s;
		flow.outward_flow_m3_s=relative_target;
		const auto moving_conditions=bridge.BoundaryConditions(moved,moving_grid,
			{{"controlled",flow},{"pressure",pressure}},dt);
		Close(moving_conditions.flow_rate_controls[0].target_outward_flow_m3_s,
			target,1e-12);
		const auto moving_solved=iga::SolveNativeTetAlePetscTransient(moved,moving_grid,
			committed,trial,{},std::numeric_limits<std::uint32_t>::max(),
			{1000.,0.004},dt,1e-10,8,moving_conditions);
		if(moving_solved.newton_iterations==0||moving_solved.final_residual_l2>=1e-10)
			throw std::runtime_error("native moving ALE graph PETSc solve did not converge");
		const auto moving_measured=bridge.Observe(moved,moving_grid,
			moving_solved.replicated_state,dt);
		Close(*moving_measured.at("controlled").outward_flow_m3_s,
			relative_target,1e-11);
		Close(*moving_measured.at("pressure").outward_flow_m3_s,
			-relative_target,1e-10);
		const auto moving_species=solve_species(mesh,moved,moving_grid,
			moving_solved.replicated_state,moving_measured);
		Close(moving_species.step.current_inventory_mol,
			fixed_species.step.current_inventory_mol,1e-8);
		if(rank==0)std::cout<<"native tetra ALE graph PETSc ports passed ranks="
			<<ranks<<" target_m3_s="<<target
			<<" measured_m3_s="
			<<*measured.at("controlled").outward_flow_m3_s
			<<" moving_relative_m3_s="<<relative_target
			<<" moving_inventory_mol="
			<<moving_species.step.current_inventory_mol
			<<" newton="<<solved.newton_iterations<<'\n';
	}catch(const std::exception& error){
		int rank=0;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		if(rank==0)std::cerr<<error.what()<<'\n';
		exit_code=1;
	}
	PetscFinalize();return exit_code;
}
