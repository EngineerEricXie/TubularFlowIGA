#include "NativeTetAleFlowTransportDomainAdapter.hpp"

#include <petscksp.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {

iga::CouplingPort Port(const char* id,int label,iga::PortQuantity input)
{
	iga::CouplingPort port;
	port.id=id;port.subsystem_id="tet";
	port.locator_kind="boundary_label";port.locator=std::to_string(label);
	port.provides={iga::PortQuantity::Area,iga::PortQuantity::FlowRate,
		iga::PortQuantity::MeanPressure,iga::PortQuantity::SpeciesConcentration,
		iga::PortQuantity::SpeciesFlux};
	port.requires={input,iga::PortQuantity::SpeciesConcentration,
		iga::PortQuantity::SpeciesFlux};
	port.species={"tracer"};return port;
}

void Close(double actual,double expected,double tolerance=1e-8)
{
	if(!std::isfinite(actual)
		||std::abs(actual-expected)>tolerance*std::max(1.,std::abs(expected)))
		throw std::runtime_error("native staged ALE transport value differs");
}

} // namespace

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int exit_code=0;
	try{
		int rank=0,ranks=1;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
		iga::NativeTetMesh mesh;
		mesh.points={{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}},{{0.25,0.25,0.25}}};
		mesh.cells={{1,{{4,1,2,3}}},{2,{{0,4,2,3}}},
			{3,{{0,1,4,3}}},{4,{{0,1,2,4}}}};
		mesh.boundary_triangles={{1,{{1,2,3}},1},
			{2,{{0,3,2}},2},{3,{{0,1,3}},2},{4,{{0,2,1}},2}};
		const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
		const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
		std::vector<double> initial_flow(3*velocity_nodes+mesh.points.size(),0.);
		for(std::size_t node=0;node<velocity_nodes;++node)initial_flow[3*node]=0.2;
		const std::vector<iga::CouplingPort> ports{
			Port("controlled",1,iga::PortQuantity::FlowRate),
			Port("pressure",2,iga::PortQuantity::MeanPressure)};
		iga::NativeTetAleFlowTransportDomainAdapter runtime("tet",mesh,ports,
			"tracer",{1000.,0.004},initial_flow,
			std::vector<double>(mesh.points.size(),2.),0.,0.,
			[](const iga::DomainStepContext& step,const iga::NativeTetMesh& previous){
				auto current=previous;
				if(step.step_index>0)
					for(auto& point:current.points)point[0]+=0.05*step.dt_s;
				return current;
			},{},1e-8,12,"rigid-plus-x-0.05");
		const double dt=0.05;
		auto hydraulic_inputs=[&](double target,double time){
			iga::PortBoundaryData flow;flow.time_s=time;
			flow.outward_flow_m3_s=target;
			iga::PortBoundaryData pressure;pressure.time_s=time;
			pressure.mean_pressure_pa=0.;
			runtime.SetPortInput("controlled",flow);
			runtime.SetPortInput("pressure",pressure);
		};
		auto solve_step=[&](int index,double start,double target){
			runtime.BeginStep({index,start,dt});
			hydraulic_inputs(target,start+dt);
			runtime.SolveHydraulicTrial();
			if(index==0){
				runtime.RollbackHydraulicTrial();
				hydraulic_inputs(target,start+dt);
				runtime.SolveHydraulicTrial();
			}
			const auto hydraulic=runtime.GetHydraulicPortState("controlled");
			Close(*hydraulic.outward_flow_m3_s,target,1e-10);
			runtime.SetTransportConcentration("pressure",start+dt,
				{{"tracer",2.}});
			runtime.SolveTransportTrial();
			if(index==0){
				runtime.RollbackTransportTrial();
				runtime.SetTransportConcentration("pressure",start+dt,
					{{"tracer",2.}});
				runtime.SolveTransportTrial();
			}
			const auto scalar=runtime.GetTransportPortState("controlled");
			Close(*scalar.outward_flow_m3_s,target,1e-10);
			Close(scalar.outward_species_flux.at("tracer"),2.*target,1e-8);
			const auto accounting=runtime.GetSpeciesStepAccounting().at("tracer");
			Close(accounting.residual,0.,1e-10);
			Close(accounting.final_mass,1./3.,1e-8);
			runtime.PrepareCommitStep();runtime.FinalizeCommitStep();
			Close(runtime.CommittedSpeciesState().time_s,start+dt);
			if(runtime.CommittedSpeciesState().accepted_steps
				!=static_cast<std::uint64_t>(index+1))
				throw std::runtime_error("native staged ALE accepted step count differs");
		};
		solve_step(0,0.,0.1);
		const auto pair=runtime.CaptureCheckpoints();
		const auto flow_bytes=iga::SerializeNativeTetAleFlowCheckpoint(pair.flow);
		const auto species_bytes=iga::SerializeNativeTetMovingSpeciesCheckpoint(pair.species);
		iga::NativeTetAleFlowTransportDomainAdapter restored("tet",mesh,ports,
			"tracer",{1000.,0.004},initial_flow,
			std::vector<double>(mesh.points.size(),2.),0.,0.,
			[](const iga::DomainStepContext& step,const iga::NativeTetMesh& previous){
				auto current=previous;
				if(step.step_index>0)
					for(auto& point:current.points)point[0]+=0.05*step.dt_s;
				return current;
			},{},1e-8,12,"rigid-plus-x-0.05");
		const iga::NativeTetAleFlowSpeciesCheckpoint decoded{
			iga::ParseNativeTetAleFlowCheckpoint(flow_bytes),
			iga::ParseNativeTetMovingSpeciesCheckpoint(species_bytes)};
		iga::NativeTetAleFlowTransportDomainAdapter monotone_runtime(
			"tet",mesh,ports,"tracer",{1000.,0.004},initial_flow,
			std::vector<double>(mesh.points.size(),2.),0.,0.,{}, {},1e-8,12,
			"rigid-plus-x-0.05",true);
		if(monotone_runtime.CheckpointModelIdentitySha256()
			==restored.CheckpointModelIdentitySha256())
			throw std::runtime_error("native staged ALE stabilization mode is absent from identity");
		bool mode_rejected=false;
		try{monotone_runtime.RestoreCheckpoints(decoded);}
		catch(const std::invalid_argument&){mode_rejected=true;}
		if(!mode_rejected||monotone_runtime.CommittedSpeciesState().accepted_steps!=0)
			throw std::runtime_error("native staged ALE mismatched stabilization checkpoint was accepted");
		auto mismatched=decoded;mismatched.flow.time_s+=dt;
		mismatched.flow.state_identity_sha256=
			iga::native_tet_ale_flow_checkpoint_detail::StateIdentity(mismatched.flow);
		bool mismatch_rejected=false;
		try{restored.RestoreCheckpoints(mismatched);}
		catch(const std::exception&){mismatch_rejected=true;}
		if(!mismatch_rejected||restored.CommittedSpeciesState().accepted_steps!=0)
			throw std::runtime_error("native staged ALE mismatched pair was accepted");
		restored.RestoreCheckpoints(decoded);
		if(iga::SerializeNativeTetAleFlowCheckpoint(restored.CaptureCheckpoints().flow)
			!=flow_bytes
			||iga::SerializeNativeTetMovingSpeciesCheckpoint(restored.CaptureCheckpoints().species)
			!=species_bytes)
			throw std::runtime_error("native staged ALE pair restore changed accepted bytes");
		const auto committed_before=iga::SerializeNativeTetMovingSpeciesCheckpoint(
			iga::CaptureNativeTetMovingSpeciesCheckpoint(runtime.CommittedSpeciesState()));
		runtime.BeginStep({1,dt,dt});
		hydraulic_inputs(0.075,2*dt);
		runtime.SolveHydraulicTrial();
		bool rejected=false;
		try{runtime.SolveTransportTrial();}
		catch(const std::invalid_argument&){rejected=true;}
		if(!rejected)throw std::runtime_error("native staged ALE missing inlet was accepted");
		runtime.AbortStep();
		if(committed_before!=iga::SerializeNativeTetMovingSpeciesCheckpoint(
			iga::CaptureNativeTetMovingSpeciesCheckpoint(runtime.CommittedSpeciesState())))
			throw std::runtime_error("native staged ALE rejected step changed committed species");
		solve_step(1,dt,0.075);
		restored.BeginStep({1,dt,dt});
		iga::PortBoundaryData flow_input;flow_input.time_s=2*dt;
		flow_input.outward_flow_m3_s=0.075;
		iga::PortBoundaryData pressure_input;pressure_input.time_s=2*dt;
		pressure_input.mean_pressure_pa=0.;
		restored.SetPortInput("controlled",flow_input);
		restored.SetPortInput("pressure",pressure_input);
		restored.SolveHydraulicTrial();
		restored.SetTransportConcentration("pressure",2*dt,{{"tracer",2.}});
		restored.SolveTransportTrial();
		restored.PrepareCommitStep();restored.FinalizeCommitStep();
		const auto expected_pair=runtime.CaptureCheckpoints();
		const auto actual_pair=restored.CaptureCheckpoints();
		{
			double flow_difference=0.,species_difference=0.;
			for(std::size_t i=0;i<expected_pair.flow.flow_state.size();++i)
				flow_difference=std::max(flow_difference,std::abs(
					expected_pair.flow.flow_state[i]-actual_pair.flow.flow_state[i]));
			for(std::size_t i=0;i<expected_pair.species.concentration_mol_m3.size();++i)
				species_difference=std::max(species_difference,std::abs(
					expected_pair.species.concentration_mol_m3[i]
					-actual_pair.species.concentration_mol_m3[i]));
			if(expected_pair.flow.model_identity_sha256!=actual_pair.flow.model_identity_sha256
				||expected_pair.species.model_identity_sha256
					!=actual_pair.species.model_identity_sha256
				||expected_pair.flow.accepted_steps!=actual_pair.flow.accepted_steps
				||expected_pair.species.accepted_steps!=actual_pair.species.accepted_steps
				||expected_pair.flow.time_s!=actual_pair.flow.time_s
				||expected_pair.species.time_s!=actual_pair.species.time_s
				||expected_pair.flow.current_points_m!=actual_pair.flow.current_points_m
				||expected_pair.species.current_points_m!=actual_pair.species.current_points_m
				||flow_difference>1e-12||species_difference>1e-12){
				std::ostringstream details;details<<std::setprecision(17)
					<<"native staged ALE pair restart continuation differs: flow="
					<<flow_difference<<" species="<<species_difference;
				throw std::runtime_error(details.str());
			}
		}
		if(rank==0)std::cout<<"native staged ALE flow/transport passed ranks="
			<<ranks<<" accepted_steps="
			<<runtime.CommittedSpeciesState().accepted_steps<<'\n';
	}catch(const std::exception& error){
		int rank=0;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		if(rank==0)std::cerr<<error.what()<<'\n';
		exit_code=1;
	}
	PetscFinalize();return exit_code;
}
