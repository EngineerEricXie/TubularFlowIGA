#include "IdealizedLeftVentricleFixture.hpp"
#include "NativeTetAleKinematics.hpp"
#include "NativeTetAleConservation.hpp"
#include "NativeTetAleMeshMotion.hpp"
#include "NativeTetAlePetscRuntime.hpp"
#include "NativeTetBoundaryFlow.hpp"
#include "NativeTetStarRadialRefinement.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

iga::NativeTetMesh StarMesh(const iga::RawSurfaceSoup& surface)
{
	iga::NativeTetMesh mesh;mesh.points=surface.vertices;
	const auto center=static_cast<std::uint32_t>(mesh.points.size());
	mesh.points.push_back({{0.,0.,.045}});std::uint64_t id=1;
	for(const auto& triangle:surface.triangles){
		auto nodes=std::array<std::uint32_t,4>{{center,
			static_cast<std::uint32_t>(triangle.indices[0]),
			static_cast<std::uint32_t>(triangle.indices[1]),
			static_cast<std::uint32_t>(triangle.indices[2])}};
		if(iga::NativeAleDeterminant(mesh.points[nodes[0]],mesh.points[nodes[1]],
			mesh.points[nodes[2]],mesh.points[nodes[3]])<0.)std::swap(nodes[2],nodes[3]);
		mesh.cells.push_back({id,nodes});
		mesh.boundary_triangles.push_back({id,{{
			static_cast<std::uint32_t>(triangle.indices[0]),
			static_cast<std::uint32_t>(triangle.indices[1]),
			static_cast<std::uint32_t>(triangle.indices[2])}},
			static_cast<int>(triangle.boundary_id)});++id;
	}
	return mesh;
}

std::array<double,3> EdgeAverage(const std::vector<std::array<double,3>>& values,
	const std::array<std::uint32_t,2>& edge)
{
	std::array<double,3> result{};
	for(int component=0;component<3;++component)
		result[component]=0.5*(values[edge[0]][component]+values[edge[1]][component]);
	return result;
}

std::pair<double,double> VolumeMeanPressure(const iga::NativeTetMesh& mesh,
	const iga::NativeTaylorHoodTopology& topology,const std::vector<double>& state)
{
	const std::size_t pressure_offset=3*(mesh.points.size()+topology.edges.size());
	if(state.size()!=pressure_offset+mesh.points.size())
		throw std::invalid_argument("native LV pressure state size is invalid");
	double pressure_volume=0.,volume=0.;
	for(const auto& cell:mesh.cells){
		const double cell_volume=iga::EvaluateNativeTetGeometry(mesh,cell).determinant/6.;
		double cell_pressure=0.;
		for(const auto node:cell.nodes)cell_pressure+=state[pressure_offset+node]/4.;
		volume+=cell_volume;
		pressure_volume+=cell_volume*cell_pressure;
	}
	if(!(volume>0.)||!std::isfinite(volume)||!std::isfinite(pressure_volume))
		throw std::runtime_error("native LV volume-mean pressure is invalid");
	return {pressure_volume/volume,volume};
}

}

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);int status=0;
	try {
		int rank=0,ranks=1;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
		std::cout.precision(17);
		std::cerr.precision(17);
		const auto motion=iga::IdealizedLeftVentricleFixture::Motion();
		auto reference=StarMesh(motion.Frames().front().surface);
#ifdef IGA_NATIVE_LV_RADIAL_REFINEMENT
		reference=iga::RefineNativeTetStarRadially(reference,
			static_cast<std::uint32_t>(reference.points.size()-1),0.5);
#endif
#ifdef IGA_NATIVE_LV_TWO_RADIAL_SHELLS
		reference=iga::RefineNativeTetStarRadialLayers(reference,
			static_cast<std::uint32_t>(reference.points.size()-1),{1./3.,2./3.});
#endif
		const auto topology=iga::BuildNativeTaylorHoodTopology(reference);
		const std::size_t velocity_nodes=reference.points.size()+topology.edges.size();
		const std::size_t dofs=3*velocity_nodes+reference.points.size();
		std::vector<double> committed(dofs,0.0);constexpr double omega=0.35;
		for(std::size_t node=0;node<velocity_nodes;++node){
			const auto coordinate=node<reference.points.size()?reference.points[node]:
				EdgeAverage(reference.points,topology.edges[node-reference.points.size()]);
			committed[3*node]=-omega*coordinate[1];
			committed[3*node+1]=omega*coordinate[0];
		}
		iga::NativeTetAleKinematics kinematics(reference.points,reference.cells);
		const auto initial_pressure_volume=VolumeMeanPressure(reference,topology,committed);
		const auto initial_velocity=iga::EvaluateNativeTetVolumeVelocity(reference,topology,
			committed,1050.0);
		double prior_pressure_pa=initial_pressure_volume.first;
		double prior_volume_m3=initial_pressure_volume.second;
		double pressure_volume_integral_j=0.;
		double final_mean_pressure_pa=0.,final_kinetic_energy_j=0.;
		double maximum_absolute_mean_pressure_pa=0.;
		iga::NativeTetAleBoundaryConditions ports;
		ports.prescribed_pressure_pa[2]=0.0;
		ports.flow_rate_controls.push_back({1,-1e-6});
#ifdef IGA_NATIVE_LV_BACKFLOW_STABILIZATION
		ports.backflow_stabilization_beta[2]=0.5;
#endif
		constexpr int steps_per_cycle=16;
#ifdef IGA_NATIVE_LV_TWO_CYCLES
		constexpr int cycles=2;
#else
		constexpr int cycles=1;
#endif
#ifdef IGA_NATIVE_LV_SINGLE_STEP
		constexpr int steps=1;
#else
		constexpr int steps=steps_per_cycle*cycles;
#endif
		const double dt=iga::IdealizedLeftVentricleFixture::PeriodS/steps_per_cycle;
		double maximum_flow_relative_error=0.0,maximum_mass_defect=0.0;
		double es_rms=0.0,final_rms=0.0,es_outlet=0.0,final_outlet=0.0;
		double cycle_return_velocity_relative_l2=0.;
		std::vector<double> cycle_reference_state=committed;
		std::size_t total_newton=0;PetscInt total_linear=0;
		auto previous_mesh=reference;
		double maximum_gcl_relative_error=0.0;
		for(int step=1;step<=steps;++step){
			if(rank==0)std::cerr<<"native ALE LV attempting step="<<step<<'\n';
			const int local_step=(step-1)%steps_per_cycle;
			const double source=local_step*dt,target=(local_step+1)*dt;
			const double absolute_target=step*dt;
			const auto evaluated=motion.Evaluate(target,source,target);
			std::map<std::uint32_t,std::array<double,3>> displacement;
			for(std::uint32_t node=0;node<motion.Frames().front().surface.vertices.size();++node)
				for(int component=0;component<3;++component)
					displacement[node][component]=evaluated.SourceVerticesM()[node][component]
						-reference.points[node][component];
			const auto harmonic=iga::SolveNativeTetAleHarmonicMotion(reference,displacement);
			const auto& trial_geometry=kinematics.BeginTrial(harmonic.displacement_m,absolute_target);
			if(rank==0)std::cerr<<"native ALE LV quality step="<<step
				<<" minimum_det_ratio="<<trial_geometry.quality.minimum_determinant_ratio
				<<" minimum_scaled_jacobian="<<trial_geometry.quality.minimum_scaled_jacobian<<'\n';
			const auto current=iga::BuildNativeTetAleCurrentMesh(reference,trial_geometry);
			std::map<std::uint32_t,std::array<double,3>> wall_velocity;
			for(const auto node:topology.boundary_velocity_nodes.at(0))
				wall_velocity[node]=node<reference.points.size()
					?evaluated.SourceVertexVelocitiesMPerS()[node]
					:EdgeAverage(evaluated.SourceVertexVelocitiesMPerS(),
						topology.edges[node-reference.points.size()]);
			const auto solved=iga::SolveNativeTetAlePetscTransient(current,
				trial_geometry.mesh_velocity_m_s,committed,committed,wall_velocity,
				std::numeric_limits<std::uint32_t>::max(),{1050.0,0.012},dt,
				1e-12,40,ports);
			const auto flows=iga::EvaluateNativeTetBoundaryFlows(current,topology,
				solved.replicated_state);
			const auto volume=iga::EvaluateNativeTetVolumeVelocity(current,topology,
				solved.replicated_state,1050.0);
			const auto mean_pressure_volume=VolumeMeanPressure(current,topology,
				solved.replicated_state);
			const double mean_pressure_pa=mean_pressure_volume.first;
			if(std::abs(mean_pressure_volume.second-volume.volume_m3)
				>1e-12*volume.volume_m3)
				throw std::runtime_error("native LV pressure and kinetic quadrature volumes disagree");
			pressure_volume_integral_j+=0.5*(prior_pressure_pa+mean_pressure_pa)
				*(volume.volume_m3-prior_volume_m3);
			prior_pressure_pa=mean_pressure_pa;
			prior_volume_m3=volume.volume_m3;
			maximum_absolute_mean_pressure_pa=std::max(maximum_absolute_mean_pressure_pa,
				std::abs(mean_pressure_pa));
			final_mean_pressure_pa=mean_pressure_pa;
			final_kinetic_energy_j=volume.kinetic_energy_j;
			std::vector<std::array<double,3>> old_vertex_velocity(reference.points.size());
			std::vector<std::array<double,3>> new_vertex_velocity(reference.points.size());
			for(std::size_t node=0;node<reference.points.size();++node)
				for(int component=0;component<3;++component){
					old_vertex_velocity[node][component]=committed[3*node+component];
					new_vertex_velocity[node][component]=solved.replicated_state[3*node+component];
				}
			const auto conservation=iga::EvaluateNativeTetAleConservation(previous_mesh,
				current,trial_geometry.mesh_velocity_m_s,old_vertex_velocity,
				new_vertex_velocity,dt);
			const double gcl_relative_error=std::abs(conservation.gcl_residual_m3_s)
				/std::max(1e-30,std::abs(conservation.volume_rate_m3_s));
			maximum_gcl_relative_error=std::max(maximum_gcl_relative_error,gcl_relative_error);
			const double inlet_error=std::abs(flows.at(1).outward_flow_m3_s+1e-6)/1e-6;
			const double mass=std::abs(flows.at(0).outward_flow_m3_s
				+flows.at(1).outward_flow_m3_s+flows.at(2).outward_flow_m3_s)/1e-6;
			maximum_flow_relative_error=std::max(maximum_flow_relative_error,inlet_error);
			maximum_mass_defect=std::max(maximum_mass_defect,mass);
			if(rank==0)std::cerr<<"native ALE LV step="<<step<<" newton="
				<<solved.newton_iterations<<" residual="<<solved.final_residual_l2
				<<" inlet_error="<<inlet_error<<" mass="<<mass
				<<" rms="<<volume.rms_speed_m_s
				<<" outlet_flow_m3_s="<<flows.at(2).outward_flow_m3_s
				<<" wall_flow_m3_s="<<flows.at(0).outward_flow_m3_s
				<<" volume_m3="<<volume.volume_m3
				<<" mean_pressure_pa="<<mean_pressure_pa
				<<" kinetic_energy_j="<<volume.kinetic_energy_j<<'\n';
			total_newton+=solved.newton_iterations;total_linear+=solved.total_linear_iterations;
			if(step==8){es_rms=volume.rms_speed_m_s;es_outlet=flows.at(2).outward_flow_m3_s;}
			if(step==steps){final_rms=volume.rms_speed_m_s;final_outlet=flows.at(2).outward_flow_m3_s;}
			if(step==steps){
				auto difference=solved.replicated_state;
				for(std::size_t dof=0;dof<3*velocity_nodes;++dof)
					difference[dof]-=cycle_reference_state[dof];
				for(std::size_t dof=3*velocity_nodes;dof<difference.size();++dof)
					difference[dof]=0.;
				const auto difference_norm=iga::EvaluateNativeTetVolumeVelocity(current,
					topology,difference,1050.0);
				cycle_return_velocity_relative_l2=difference_norm.rms_speed_m_s
					/volume.rms_speed_m_s;
			}
			if(step==steps_per_cycle&&cycles>1)cycle_reference_state=solved.replicated_state;
			if(!(solved.final_residual_l2<1e-9&&inlet_error<1e-3&&mass<2e-8
				&&gcl_relative_error<2e-12
				&&std::isfinite(volume.rms_speed_m_s)))
				throw std::runtime_error("native ALE idealized LV step failed its functional gate");
			committed=solved.replicated_state;kinematics.CommitTrial();previous_mesh=current;
		}
		const double periodic_volume_relative_error=
			std::abs(prior_volume_m3-initial_pressure_volume.second)
				/initial_pressure_volume.second;
		if(!std::isfinite(final_mean_pressure_pa)
			||!std::isfinite(maximum_absolute_mean_pressure_pa)
			||!std::isfinite(pressure_volume_integral_j)
			||!std::isfinite(final_kinetic_energy_j)
			||(steps%steps_per_cycle==0&&periodic_volume_relative_error>1e-12)
			||!std::isfinite(cycle_return_velocity_relative_l2))
			throw std::runtime_error("native ALE LV cycle QoI is nonfinite or geometric volume is not periodic");
		if(rank==0)std::cout<<"native ALE idealized LV flow passed ranks="<<ranks
			<<" tetrahedra="<<reference.cells.size()
			<<" cycles="<<cycles
			<<" steps="<<steps<<" max_flow_relative_error="<<maximum_flow_relative_error
			<<" max_mass_defect="<<maximum_mass_defect
			<<" max_gcl_relative_error="<<maximum_gcl_relative_error
			<<" es_rms_m_s="<<es_rms
			<<" final_rms_m_s="<<final_rms<<" es_outlet_m3_s="<<es_outlet
			<<" final_outlet_m3_s="<<final_outlet<<" total_newton="<<total_newton
			<<" total_linear="<<total_linear
			<<" final_mean_pressure_pa="<<final_mean_pressure_pa
			<<" maximum_absolute_mean_pressure_pa="<<maximum_absolute_mean_pressure_pa
			<<" pressure_volume_integral_j="<<pressure_volume_integral_j
			<<" initial_kinetic_energy_j="<<initial_velocity.kinetic_energy_j
			<<" final_kinetic_energy_j="<<final_kinetic_energy_j
			<<" periodic_volume_relative_error="<<periodic_volume_relative_error
			<<" cycle_return_velocity_relative_l2="<<cycle_return_velocity_relative_l2<<'\n';
	} catch(const std::exception& error){
		int rank=0;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		if(rank==0)std::cerr<<error.what()<<'\n';
		status=1;
	}
	PetscFinalize();return status;
}
