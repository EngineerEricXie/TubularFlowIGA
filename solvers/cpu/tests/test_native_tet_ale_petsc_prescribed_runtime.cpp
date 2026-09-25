#include "NativeTetAleConservation.hpp"
#include "NativeTetAleKinematics.hpp"
#include "NativeTetAleMeshMotion.hpp"
#include "NativeTetAlePetscRuntime.hpp"

#include <petscksp.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);int exit_code=0;
	try {
		int rank=0,ranks=1;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
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
		const std::array<double,3> velocity{{0.12,-0.03,0.05}};
		std::vector<double> committed_fluid(dofs,0.0);
		for(std::size_t node=0;node<velocity_nodes;++node)
			for(int component=0;component<3;++component)
				committed_fluid[3*node+component]=velocity[component];
		std::map<std::uint32_t,std::array<double,3>> fluid_boundary;
		for(const auto node:topology.boundary_velocity_nodes.at(0)) fluid_boundary[node]=velocity;
		iga::NativeTetAleKinematics kinematics(reference.points,reference.cells);
		auto previous_mesh=reference;const double dt=0.04;
		PetscInt total_linear_iterations=0;
		for(std::size_t step=1;step<=5;++step) {
			const double next_time=step*dt;
			std::map<std::uint32_t,std::array<double,3>> boundary_displacement;
			for(std::uint32_t node=0;node<4;++node)
				for(int component=0;component<3;++component)
					boundary_displacement[node][component]=next_time*velocity[component];
			const auto motion=iga::SolveNativeTetAleHarmonicMotion(reference,boundary_displacement);
			const auto& geometry_trial=kinematics.BeginTrial(motion.displacement_m,next_time);
			const auto current_mesh=iga::BuildNativeTetAleCurrentMesh(reference,geometry_trial);
			auto fluid_trial=committed_fluid;
			for(std::size_t node=0;node<velocity_nodes;++node)
				if(topology.boundary_velocity_nodes.at(0).count(static_cast<std::uint32_t>(node))==0)
					fluid_trial[3*node]+=0.002;
			const auto fluid=iga::SolveNativeTetAlePetscTransient(current_mesh,
				geometry_trial.mesh_velocity_m_s,committed_fluid,fluid_trial,fluid_boundary,
				0,{1000.0,0.004},dt,1e-10,8);
			if(!(fluid.newton_iterations>0&&fluid.final_residual_l2<1e-10))
				throw std::runtime_error("distributed prescribed ALE step did not converge");
			const std::vector<std::array<double,3>> nodal_fluid(reference.points.size(),velocity);
			const auto conservation=iga::EvaluateNativeTetAleConservation(previous_mesh,current_mesh,
				geometry_trial.mesh_velocity_m_s,nodal_fluid,nodal_fluid,dt);
			if(std::abs(conservation.gcl_residual_m3_s)>1e-13
				||std::abs(conservation.moving_domain_balance_residual_m3_s)>1e-13)
				throw std::runtime_error("distributed prescribed ALE conservation gate failed");
			committed_fluid=fluid.replicated_state;total_linear_iterations+=fluid.total_linear_iterations;
			kinematics.CommitTrial();previous_mesh=current_mesh;
		}
		if(std::abs(kinematics.CommittedTime()-0.2)>1e-14)
			throw std::runtime_error("distributed prescribed ALE history did not advance");
		if(rank==0) std::cout<<"native tetrahedral ALE PETSc five-step runtime passed ranks="
			<<ranks<<" linear_iterations="<<total_linear_iterations<<'\n';
	} catch(const std::exception& error) {
		int rank=0;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		if(rank==0) std::cerr<<error.what()<<'\n';
		exit_code=1;
	}
	PetscFinalize();return exit_code;
}
