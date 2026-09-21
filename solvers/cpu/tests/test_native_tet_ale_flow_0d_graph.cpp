#include "CollectivePressureFlowExecution.hpp"
#include "NativeTetAleFlowDomainAdapter.hpp"
#include "NativeTetHydraulicGraphCheckpoint.hpp"
#include "PressureFlowComponentExecutor.hpp"
#include "ZeroDFlowDomain.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

iga::NativeTetMesh Mesh()
{
	iga::NativeTetMesh mesh;
	mesh.points={{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}},{{0.25,0.25,0.25}}};
	mesh.cells={{1,{{4,1,2,3}}},{2,{{0,4,2,3}}},
		{3,{{0,1,4,3}}},{4,{{0,1,2,4}}}};
	mesh.boundary_triangles={{1,{{1,2,3}},1},
		{2,{{0,3,2}},2},{3,{{0,1,3}},2},{4,{{0,2,1}},2}};
	return mesh;
}

iga::CouplingPort TetPort(const std::string& id,int label,
	iga::PortQuantity required)
{
	iga::CouplingPort port;
	port.id=id;port.subsystem_id="tet";
	port.locator_kind="boundary_label";port.locator=std::to_string(label);
	port.provides={iga::PortQuantity::Area,iga::PortQuantity::FlowRate,
		iga::PortQuantity::MeanPressure};
	port.requires={required};
	return port;
}

void Close(double actual,double expected,double tolerance=1e-7)
{
	if(!std::isfinite(actual)||!std::isfinite(expected)
		||std::abs(actual-expected)>tolerance*std::max({1.,std::abs(expected)}))
		throw std::runtime_error("native tetra/0D graph quantity differs");
}

} // namespace

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int exit_code=0;
	try{
		int rank=0,ranks=1;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
		const std::string argument=argc>1?argv[1]:"";
		const bool process_save=argument=="--cross-process-save";
		const bool process_resume=argument=="--cross-process-resume";
		if((process_save||process_resume)?argc!=3:argc>2)
			throw std::invalid_argument("native tetra/0D graph test arguments differ");
		const std::filesystem::path process_root=
			process_save||process_resume?argv[2]:"";
		const std::string mode=process_save||process_resume
			?"--balanced-strong":argument;
		const bool balanced=mode=="--balanced-strong"||mode=="--balanced-explicit";
		const bool strong=mode=="--strong"||mode=="--balanced-strong";
		const auto mesh=Mesh();
		const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
		const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
		std::vector<double> initial_flow(3*velocity_nodes+mesh.points.size(),0.);
		for(std::size_t node=0;node<velocity_nodes;++node)
			initial_flow[3*node]=-0.2;
		const std::vector<iga::CouplingPort> tet_ports{
			TetPort("inlet",1,iga::PortQuantity::FlowRate),
			TetPort("outlet",2,iga::PortQuantity::MeanPressure)};
		iga::ZeroDFlowModel source_model;
		source_model.role=iga::ZeroDFlowRole::SourceReservoir;
		source_model.source=balanced
			?iga::ZeroDSourceReservoirModel{1e-9,500.,0.1}
			:iga::ZeroDSourceReservoirModel{0.001,1.,0.1};
		iga::ZeroDFlowModel terminal_model;
		terminal_model.role=iga::ZeroDFlowRole::TerminalRcr;
		terminal_model.terminal=balanced
			?iga::ZeroDTerminalRcrModel{10.,50.,1e-9,0.}
			:iga::ZeroDTerminalRcrModel{0.1,1.,0.001,0.};
		const auto source_port=iga::MakeZeroDFlowPort("source",
			iga::ZeroDFlowRole::SourceReservoir);
		const auto terminal_port=iga::MakeZeroDFlowPort("terminal",
			iga::ZeroDFlowRole::TerminalRcr);
		iga::SimulationGraph graph(
			{{"source",iga::DomainKind::ZeroDFlow,{source_port}},
			 {"tet",iga::DomainKind::ThreeDBodyFittedFlow,tet_ports},
			 {"terminal",iga::DomainKind::ZeroDFlow,{terminal_port}}},
			{{"source_tet",{"source","port"},{"tet","inlet"},
				iga::CouplingLaw::PressureFlow},
			 {"tet_rcr",{"tet","outlet"},{"terminal","port"},
				iga::CouplingLaw::PressureFlow}});
		auto source=std::make_unique<iga::ZeroDFlowDomainRuntime>(
			"source",source_model,iga::ZeroDFlowState{0.},
			std::vector<iga::CouplingPort>{source_port});
		auto terminal=std::make_unique<iga::ZeroDFlowDomainRuntime>(
			"terminal",terminal_model,iga::ZeroDFlowState{0.},
			std::vector<iga::CouplingPort>{terminal_port});
		const auto motion=[](const iga::DomainStepContext& step,
			const iga::NativeTetMesh& previous){
			auto current=previous;
			if(step.step_index>0)
				for(auto& point:current.points)point[0]+=0.001*step.dt_s;
			return current;
		};
		const std::string motion_id="rigid-translation-x-0.001-v1";
		auto native=std::make_unique<iga::NativeTetAleFlowDomainAdapter>(
			"tet",mesh,tet_ports,
			iga::NativeNavierStokesParameters{strong&&!balanced?1.:1000.,0.004},
			initial_flow,motion,motion_id);
		auto* source_ptr=source.get();
		auto* terminal_ptr=terminal.get();
		auto* native_ptr=native.get();
		std::vector<std::unique_ptr<iga::CoupledDomainRuntime>> runtimes;
		runtimes.push_back(std::move(source));
		runtimes.push_back(std::move(native));
		runtimes.push_back(std::move(terminal));
		iga::DomainRuntimeRegistry registry(graph,std::move(runtimes));
		iga::PressureFlowExecutionControls controls;
		controls.method=strong?iga::PressureFlowIterationMethod::Fixed
			:iga::PressureFlowIterationMethod::Explicit;
		controls.maximum_iterations=strong?20:1;
		controls.pressure_relative_tolerance=strong?1e-4:1e-6;
		controls.flow_relative_tolerance=1e-6;
		iga::PressureFlowComponentExecutor executor(registry,"source",controls,
			iga::CollectivePressureFlowExecution(PETSC_COMM_WORLD));
		const auto case_hash=iga::coupled_checkpoint_detail::Hash(
			"native-tet-0d-graph-test-case/v1");
		const auto execution_hash=iga::coupled_checkpoint_detail::Hash(
			"native-tet-0d-graph-test-execution/v1:"+std::to_string(ranks));
		if(process_resume){
			iga::NativeTetHydraulicGraphCheckpoint reader(PETSC_COMM_WORLD,
				graph,controls,*native_ptr,
				{{"source",source_ptr},{"terminal",terminal_ptr}});
			const auto image=reader.LoadLatest(process_root,
				reader.Compatibility(case_hash,execution_hash));
			if(image.epoch.accepted_steps!=1||image.epoch.time_s!=0.05
				||image.epoch.dt_s!=0.05)
				throw std::runtime_error("native cross-process checkpoint clock differs");
			reader.RestoreFresh(image);
		}
		std::size_t final_iterations=0;
		double final_outlet_m3_s=0.;
		std::vector<double> restarted_flow;
		std::vector<std::array<double,3>> restarted_points;
		double restarted_source_pa=0.,restarted_terminal_pa=0.;
		for(int step=process_resume?1:0;step<2;++step){
			if(step==1&&!process_save&&!process_resume){
				const auto checkpoint_bytes=iga::SerializeNativeTetAleFlowCheckpoint(
					native_ptr->CaptureCheckpoint());
				const auto checkpoint=iga::ParseNativeTetAleFlowCheckpoint(
					checkpoint_bytes);
				iga::NativeTetAleFlowDomainAdapter restored("tet",mesh,tet_ports,
					iga::NativeNavierStokesParameters{strong&&!balanced?1.:1000.,0.004},
					initial_flow,motion,motion_id);
				restored.RestoreCheckpoint(checkpoint);
				if(restored.CommittedFlowState()!=native_ptr->CommittedFlowState()
					||restored.CommittedMesh().points!=native_ptr->CommittedMesh().points
					||restored.CommittedTime()!=native_ptr->CommittedTime()
					||restored.AcceptedSteps()!=native_ptr->AcceptedSteps())
					throw std::runtime_error("native tetra flow checkpoint restart differs");
				auto damaged=checkpoint_bytes;
				damaged[damaged.size()/2]^=1;
				bool rejected_damage=false;
				try{iga::ParseNativeTetAleFlowCheckpoint(damaged);}
				catch(const std::exception&){rejected_damage=true;}
				iga::NativeTetAleFlowDomainAdapter wrong_motion("tet",mesh,tet_ports,
					iga::NativeNavierStokesParameters{strong&&!balanced?1.:1000.,0.004},
					initial_flow,motion,"different-motion-v1");
				bool rejected_model=false;
				try{wrong_motion.RestoreCheckpoint(checkpoint);}
				catch(const std::exception&){rejected_model=true;}
				if(!rejected_damage||!rejected_model)
					throw std::runtime_error("native tetra flow checkpoint validation missed damage");
				auto restarted_source=std::make_unique<iga::ZeroDFlowDomainRuntime>(
					"source",source_model,iga::ZeroDFlowState{0.},
					std::vector<iga::CouplingPort>{source_port});
				auto restarted_terminal=std::make_unique<iga::ZeroDFlowDomainRuntime>(
					"terminal",terminal_model,iga::ZeroDFlowState{0.},
					std::vector<iga::CouplingPort>{terminal_port});
				auto restarted_native=std::make_unique<iga::NativeTetAleFlowDomainAdapter>(
					"tet",mesh,tet_ports,
					iga::NativeNavierStokesParameters{strong&&!balanced?1.:1000.,0.004},
					initial_flow,motion,motion_id);
				auto* restarted_source_ptr=restarted_source.get();
				auto* restarted_terminal_ptr=restarted_terminal.get();
				auto* restarted_native_ptr=restarted_native.get();
				iga::NativeTetHydraulicGraphCheckpoint writer(PETSC_COMM_WORLD,
					graph,controls,
					*native_ptr,{{"source",source_ptr},{"terminal",terminal_ptr}});
				iga::NativeTetHydraulicGraphCheckpoint reader(PETSC_COMM_WORLD,
					graph,controls,
					*restarted_native_ptr,
					{{"source",restarted_source_ptr},{"terminal",restarted_terminal_ptr}});
				auto changed_controls=controls;
				changed_controls.flow_relative_tolerance*=2.;
				iga::NativeTetHydraulicGraphCheckpoint changed_reader(PETSC_COMM_WORLD,
					graph,changed_controls,*restarted_native_ptr,
					{{"source",restarted_source_ptr},{"terminal",restarted_terminal_ptr}});
				const iga::CoupledCheckpointEpoch epoch{
					iga::coupled_checkpoint_detail::Hash("native-tet-0d-graph-test-epoch/v1:"+mode),
					{},writer.Compatibility(case_hash,execution_hash),1,0.05,0.05};
				if(changed_reader.Compatibility(case_hash,execution_hash)
					.configuration_sha256==epoch.compatibility.configuration_sha256)
					throw std::runtime_error("native graph changed solver controls share checkpoint identity");
				std::array<char,64> temp_path{};
				iga::CollectiveLocalStage(PETSC_COMM_WORLD,"native graph temp directory",[&]{
					if(rank!=0)return;
					const std::string pattern="/tmp/iga-native-hydraulic-XXXXXX";
					std::copy(pattern.begin(),pattern.end(),temp_path.begin());
					if(!mkdtemp(temp_path.data()))
						throw std::runtime_error("native graph temp directory creation failed");
				});
				MPI_Bcast(temp_path.data(),static_cast<int>(temp_path.size()),MPI_CHAR,
					0,PETSC_COMM_WORLD);
				const std::filesystem::path checkpoint_root(temp_path.data());
				writer.Save(checkpoint_root,epoch);
				std::string original_shard;
				const auto shard_path=checkpoint_root/epoch.id/
					(writer.Catalog().front().id+".shard");
				iga::CollectiveLocalStage(PETSC_COMM_WORLD,
					"native hydraulic shard corruption",[&]{
						if(rank!=0)return;
						std::ifstream input(shard_path,std::ios::binary);
						original_shard.assign(std::istreambuf_iterator<char>(input),{});
						if(!input.is_open()||input.bad()||original_shard.empty())
							throw std::runtime_error("native hydraulic shard read failed");
						auto damaged=original_shard;
						damaged.back()^=1;
						std::ofstream output(shard_path,std::ios::binary|std::ios::trunc);
						output.write(damaged.data(),damaged.size());
						if(!output)throw std::runtime_error("native hydraulic shard damage failed");
					});
				bool rejected_shard=false;
				try{reader.LoadLatest(checkpoint_root,
					reader.Compatibility(case_hash,execution_hash));}
				catch(const std::exception&){rejected_shard=true;}
				if(!rejected_shard||restarted_native_ptr->AcceptedSteps()!=0
					||restarted_source_ptr->CommittedStepCount()!=0
					||restarted_terminal_ptr->CommittedStepCount()!=0)
					throw std::runtime_error("damaged native graph checkpoint was accepted");
				iga::CollectiveLocalStage(PETSC_COMM_WORLD,
					"native hydraulic shard repair",[&]{
						if(rank!=0)return;
						std::ofstream output(shard_path,std::ios::binary|std::ios::trunc);
						output.write(original_shard.data(),original_shard.size());
						if(!output)throw std::runtime_error("native hydraulic shard repair failed");
					});
				const auto image=reader.LoadLatest(checkpoint_root,
					reader.Compatibility(case_hash,execution_hash));
				if(image.epoch.id!=epoch.id)
					throw std::runtime_error("native graph file checkpoint epoch differs");
				reader.RestoreFresh(image);
				iga::CollectiveLocalStage(PETSC_COMM_WORLD,"native graph temp cleanup",[&]{
					if(rank==0)std::filesystem::remove_all(checkpoint_root);
				});
				std::vector<std::unique_ptr<iga::CoupledDomainRuntime>> restarted_runtimes;
				restarted_runtimes.push_back(std::move(restarted_source));
				restarted_runtimes.push_back(std::move(restarted_native));
				restarted_runtimes.push_back(std::move(restarted_terminal));
				iga::DomainRuntimeRegistry restarted_registry(graph,
					std::move(restarted_runtimes));
				iga::PressureFlowComponentExecutor restarted_executor(restarted_registry,
					"source",controls,
					iga::CollectivePressureFlowExecution(PETSC_COMM_WORLD));
				restarted_executor.Advance({step,0.05*step,0.05},
					{{"source_tet",0.},{"tet_rcr",0.}});
				restarted_flow=restarted_native_ptr->CommittedFlowState();
				restarted_points=restarted_native_ptr->CommittedMesh().points;
				restarted_source_pa=
					restarted_source_ptr->CommittedState().stored_pressure_pa;
				restarted_terminal_pa=
					restarted_terminal_ptr->CommittedState().stored_pressure_pa;
				const auto native_flow=native_ptr->CommittedFlowState();
				const auto source_state=source_ptr->CommittedState();
				const auto terminal_state=terminal_ptr->CommittedState();
				bool rejected=false;
				try{
					executor.Advance({step,0.05*step,0.05},
						{{"source_tet",0.},{"tet_rcr",0.}},
						[rank](const iga::PressureFlowStepResult&){
							if(rank==0)
								throw std::runtime_error("injected source/RCR precommit failure");
						});
				}catch(const std::exception& error){
					rejected=std::string(error.what()).find(
						"injected source/RCR precommit failure")!=std::string::npos;
				}
				if(!rejected||native_ptr->CommittedFlowState()!=native_flow
					||native_ptr->AcceptedSteps()!=1
					||source_ptr->CommittedState().stored_pressure_pa
						!=source_state.stored_pressure_pa
					||terminal_ptr->CommittedState().stored_pressure_pa
						!=terminal_state.stored_pressure_pa
					||source_ptr->CommittedStepCount()!=1
					||terminal_ptr->CommittedStepCount()!=1)
					throw std::runtime_error("source/RCR graph one-rank failure did not roll back");
			}
			const auto result=executor.Advance({step,0.05*step,0.05},
				{{"source_tet",0.},{"tet_rcr",0.}});
			if(result.iterations.empty()||(!strong&&result.iterations.size()!=1)
				||(strong&&(result.iterations.size()<2
					||!result.iterations.back().converged))
				||result.iterations.back().edges.size()!=2)
				throw std::runtime_error("native tetra/0D graph hydraulic result is incomplete");
			final_iterations=result.iterations.size();
			for(const auto& edge:result.iterations.back().edges)
				Close(edge.flow_residual_m3_s,0.,1e-6);
			const auto inlet=*result.accepted_ports.at({"tet","inlet"}).outward_flow_m3_s;
			const auto outlet=*result.accepted_ports.at({"tet","outlet"}).outward_flow_m3_s;
			Close(inlet+outlet,0.,1e-6);
			if(!(inlet<-1e-4&&outlet>1e-4))
				throw std::runtime_error("native tetra/0D graph has no throughflow");
			final_outlet_m3_s=outlet;
			if(source_ptr->CommittedStepCount()!=static_cast<std::size_t>(step+1)
				||terminal_ptr->CommittedStepCount()!=static_cast<std::size_t>(step+1)
				||native_ptr->AcceptedSteps()!=static_cast<std::uint64_t>(step+1))
				throw std::runtime_error("native tetra/0D graph commit count differs");
			Close(source_ptr->CommittedStepAccounting()->residual_m3,0.,1e-8);
			Close(terminal_ptr->CommittedStepAccounting()->residual_m3,0.,1e-8);
			Close(native_ptr->CommittedTime(),0.05*(step+1));
			if(step==0&&process_save){
				iga::NativeTetHydraulicGraphCheckpoint writer(PETSC_COMM_WORLD,
					graph,controls,*native_ptr,
					{{"source",source_ptr},{"terminal",terminal_ptr}});
				const iga::CoupledCheckpointEpoch epoch{
					iga::coupled_checkpoint_detail::Hash(
						"native-tet-0d-graph-test-epoch/v1:--balanced-strong"),
					{},writer.Compatibility(case_hash,execution_hash),1,0.05,0.05};
				writer.Save(process_root,epoch);
			}
			if(step==1&&!process_save&&!process_resume){
				if(restarted_flow.size()!=native_ptr->CommittedFlowState().size()
					||restarted_points!=native_ptr->CommittedMesh().points)
					throw std::runtime_error("restarted native tetra/0D graph shape differs");
				for(std::size_t i=0;i<restarted_flow.size();++i)
					Close(restarted_flow[i],native_ptr->CommittedFlowState()[i],1e-10);
				Close(restarted_source_pa,
					source_ptr->CommittedState().stored_pressure_pa,1e-10);
				Close(restarted_terminal_pa,
					terminal_ptr->CommittedState().stored_pressure_pa,1e-10);
			}
		}
		Close(native_ptr->CommittedMesh().points[0][0],0.00005);
		if(process_save)
			iga::CollectiveLocalStage(PETSC_COMM_WORLD,
				"native cross-process reference write",[&]{
					if(rank!=0)return;
					std::ofstream output(process_root/"reference.txt");
					if(!output)throw std::runtime_error("native reference open failed");
					output<<std::setprecision(17);
					output<<native_ptr->CommittedFlowState().size()<<'\n';
					for(const double value:native_ptr->CommittedFlowState())
						output<<value<<'\n';
					output<<native_ptr->CommittedMesh().points.size()<<'\n';
					for(const auto& point:native_ptr->CommittedMesh().points)
						for(const double value:point)output<<value<<'\n';
					output<<source_ptr->CommittedState().stored_pressure_pa<<'\n'
						<<terminal_ptr->CommittedState().stored_pressure_pa<<'\n'
						<<final_outlet_m3_s<<'\n';
					output.close();
					if(!output)throw std::runtime_error("native reference write failed");
				});
		if(process_resume)
			iga::CollectiveLocalStage(PETSC_COMM_WORLD,
				"native cross-process reference compare",[&]{
					std::ifstream input(process_root/"reference.txt");
					if(!input)throw std::runtime_error("native reference open failed");
					std::size_t count=0;
					input>>count;
					if(!input||count!=native_ptr->CommittedFlowState().size())
						throw std::runtime_error("native reference flow shape differs");
					for(const double value:native_ptr->CommittedFlowState()){
						double expected=0.;input>>expected;
						if(!input)throw std::runtime_error("native reference flow is truncated");
						Close(value,expected,1e-10);
					}
					input>>count;
					if(!input||count!=native_ptr->CommittedMesh().points.size())
						throw std::runtime_error("native reference geometry shape differs");
					for(const auto& point:native_ptr->CommittedMesh().points)
						for(const double value:point){
							double expected=0.;input>>expected;
							if(!input||value!=expected)
								throw std::runtime_error("native reference geometry differs");
						}
					for(const double value:{source_ptr->CommittedState().stored_pressure_pa,
						terminal_ptr->CommittedState().stored_pressure_pa,final_outlet_m3_s}){
						double expected=0.;input>>expected;
						if(!input)throw std::runtime_error("native reference pressure is truncated");
						Close(value,expected,1e-10);
					}
					std::string extra;
					if(input>>extra)throw std::runtime_error("native reference has extra data");
				});
		if(rank==0)std::cout<<"native tetra ALE source/RCR graph passed ranks="
			<<ranks<<" steps=2 strong="<<strong<<" balanced_parameters="<<balanced
			<<" iterations="<<final_iterations
			<<" outlet_m3_s="<<final_outlet_m3_s<<" source_pa="
			<<source_ptr->CommittedState().stored_pressure_pa
			<<" terminal_pa="<<terminal_ptr->CommittedState().stored_pressure_pa<<'\n';
	}catch(const std::exception& error){
		int rank=0;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		if(rank==0)std::cerr<<error.what()<<'\n';
		exit_code=1;
	}
	PetscFinalize();return exit_code;
}
