#include "NativeTetAleConservation.hpp"
#include "NativeTetAleKinematics.hpp"
#include "NativeTetAlePetscRuntime.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);int exit_code=0;
	try {
		int rank=0,ranks=1;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
		if(argc<2||argc>4)
			throw std::invalid_argument("usage: native_tet_ale_petsc_mesh_runtime_test mesh.msh [result.json] [translation|shear]");
		std::ifstream input(argv[1]);
		if(!input) throw std::runtime_error("cannot open native ALE mesh");
		const auto reference=iga::ReadNativeTetMeshGmsh41(input);
		const auto topology=iga::BuildNativeTaylorHoodTopology(reference);
		const std::size_t velocity_nodes=reference.points.size()+topology.edges.size();
		const std::size_t dofs=3*velocity_nodes+reference.points.size();
		const std::string mode=argc==4?argv[3]:"translation";
		if(mode!="translation"&&mode!="shear")
			throw std::invalid_argument("native ALE mesh runtime mode must be translation or shear");
		const std::array<double,3> translation{{0.01,-0.004,0.002}};
		const double shear_rate=0.2,dt=0.01;
		std::vector<std::array<double,3>> displacement(reference.points.size());
		for(std::size_t node=0;node<displacement.size();++node)
			if(mode=="translation")
				for(int component=0;component<3;++component)
					displacement[node][component]=dt*translation[component];
			else displacement[node][0]=dt*shear_rate*reference.points[node][1];
		iga::NativeTetAleKinematics kinematics(reference.points,reference.cells);
		const auto& geometry=kinematics.BeginTrial(displacement,dt);
		const auto current=iga::BuildNativeTetAleCurrentMesh(reference,geometry);
		std::set<std::uint32_t> boundary_nodes;
		for(const auto& labelled:topology.boundary_velocity_nodes)
			boundary_nodes.insert(labelled.second.begin(),labelled.second.end());
		std::map<std::uint32_t,std::array<double,3>> boundary;
		auto exact_velocity=[&](std::size_t node) {
			if(mode=="translation") return translation;
			double y=0.0;
			if(node<reference.points.size()) y=reference.points[node][1];
			else {
				const auto& edge=topology.edges[node-reference.points.size()];
				y=0.5*(reference.points[edge[0]][1]+reference.points[edge[1]][1]);
			}
			return std::array<double,3>{{shear_rate*y,0,0}};
		};
		for(const auto node:boundary_nodes) boundary[node]=exact_velocity(node);
		std::vector<double> committed(dofs,0.0),trial(dofs,0.0);
		for(std::size_t node=0;node<velocity_nodes;++node) {
			const auto velocity=exact_velocity(node);
			for(int component=0;component<3;++component)
				committed[3*node+component]=trial[3*node+component]=velocity[component];
		}
		for(std::size_t node=0;node<velocity_nodes;++node)
			if(boundary_nodes.count(static_cast<std::uint32_t>(node))==0)
				trial[3*node]+=1e-4;
		const auto result=iga::SolveNativeTetAlePetscTransient(current,
			geometry.mesh_velocity_m_s,committed,trial,boundary,0,{1000.0,0.004},dt,1e-8,8);
		double maximum_velocity_error=0.0;
		for(std::size_t node=0;node<velocity_nodes;++node)
			for(int component=0;component<3;++component) {
				const auto velocity=exact_velocity(node);
				maximum_velocity_error=std::max(maximum_velocity_error,
					std::abs(result.replicated_state[3*node+component]-velocity[component]));
			}
		std::vector<std::array<double,3>> nodal_fluid(reference.points.size());
		for(std::size_t node=0;node<nodal_fluid.size();++node)
			nodal_fluid[node]=exact_velocity(node);
		const auto conservation=iga::EvaluateNativeTetAleConservation(reference,current,
			geometry.mesh_velocity_m_s,nodal_fluid,nodal_fluid,dt);
		if(rank==0) std::cout<<"native ALE mesh diagnostics residual="
			<<result.final_residual_l2<<" velocity_error="<<maximum_velocity_error
			<<" gcl="<<conservation.gcl_residual_m3_s
			<<" moving_balance="<<conservation.moving_domain_balance_residual_m3_s<<'\n';
		if(result.final_residual_l2>1e-8||maximum_velocity_error>1e-6
			||std::abs(conservation.gcl_residual_m3_s)>1e-12
			||std::abs(conservation.moving_domain_balance_residual_m3_s)>1e-12)
			throw std::runtime_error("native ALE mesh-scale runtime failed its numerical gates");
		const auto accepted_quality=geometry.quality;
		kinematics.CommitTrial();
		if(rank==0) {
			std::cout<<"native tetrahedral ALE mesh runtime passed ranks="<<ranks
				<<" cells="<<current.cells.size()<<" dofs="<<dofs
				<<" newton="<<result.newton_iterations
				<<" residual="<<result.final_residual_l2
				<<" velocity_error="<<maximum_velocity_error<<'\n';
			if(argc>=3) {
				std::ofstream output(argv[2]);
				if(!output) throw std::runtime_error("cannot open native ALE mesh result");
				output<<std::setprecision(17)<<"{\n"
					<<"  \"schema_version\": 1,\n"
					<<"  \"backend\": \"native_cpp_petsc_tetrahedral_ALE\",\n"
					<<"  \"motion_mode\": \""<<mode<<"\",\n"
					<<"  \"mpi_ranks\": "<<ranks<<",\n"
					<<"  \"mesh_vertices\": "<<current.points.size()<<",\n"
					<<"  \"velocity_edges\": "<<topology.edges.size()<<",\n"
					<<"  \"tetrahedra\": "<<current.cells.size()<<",\n"
					<<"  \"mixed_dofs\": "<<dofs<<",\n"
					<<"  \"newton_iterations\": "<<result.newton_iterations<<",\n"
					<<"  \"linear_iterations\": "<<result.total_linear_iterations<<",\n"
					<<"  \"final_residual_l2\": "<<result.final_residual_l2<<",\n"
					<<"  \"maximum_velocity_error\": "<<maximum_velocity_error<<",\n"
					<<"  \"minimum_determinant_ratio\": "
					<<accepted_quality.minimum_determinant_ratio<<",\n"
					<<"  \"minimum_scaled_jacobian\": "
					<<accepted_quality.minimum_scaled_jacobian<<",\n"
					<<"  \"gcl_residual_m3_s\": "<<conservation.gcl_residual_m3_s<<",\n"
					<<"  \"moving_domain_balance_residual_m3_s\": "
					<<conservation.moving_domain_balance_residual_m3_s<<"\n}\n";
				if(!output) throw std::runtime_error("cannot write native ALE mesh result");
			}
		}
	} catch(const std::exception& error) {
		int rank=0;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		if(rank==0) std::cerr<<error.what()<<'\n';
		exit_code=1;
	}
	PetscFinalize();return exit_code;
}
