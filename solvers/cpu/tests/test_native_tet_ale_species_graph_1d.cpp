#include "NativeTetAleFlowTransportDomainAdapter.hpp"
#include "NativeTetAleKinematics.hpp"
#include "CollectiveSpeciesPressureFlowExecution.hpp"
#include "OneDFlowDomainAdapter.hpp"
#include "SpeciesPressureFlowComponentExecutor.hpp"
#include "PressureFlowCheckpointControls.hpp"
#include "SpeciesGraphCheckpointIdentity.hpp"
#include "CoupledCheckpointBundle.hpp"
#include "ZeroDSourceReservoirSpeciesDomainRuntime.hpp"
#include "ZeroDTerminalRcrSpeciesDomainRuntime.hpp"
#include "NativeTetSpeciesVisualization.hpp"
#include "ParallelVtkOutput.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

iga::CouplingPort Port(const std::string& domain,const std::string& id,
	const std::string& locator_kind,const std::string& locator,
	std::optional<iga::PortQuantity> hydraulic_input)
{
	iga::CouplingPort port;
	port.id=id;port.subsystem_id=domain;
	port.locator_kind=locator_kind;port.locator=locator;
	port.provides={iga::PortQuantity::Area,iga::PortQuantity::FlowRate,
		iga::PortQuantity::MeanPressure,iga::PortQuantity::SpeciesConcentration,
		iga::PortQuantity::SpeciesFlux};
	if(hydraulic_input)port.requires.insert(*hydraulic_input);
	if(hydraulic_input){
		port.requires.insert(iga::PortQuantity::SpeciesConcentration);
		port.requires.insert(iga::PortQuantity::SpeciesFlux);
	}
	port.species={"tracer"};return port;
}

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

iga::NativeTetMesh ChannelMesh()
{
	iga::NativeTetMesh mesh;
	auto node=[](int i,int j,int k){
		return static_cast<std::uint32_t>((i*3+j)*3+k);
	};
	for(int i=0;i<2;++i)for(int j=0;j<3;++j)for(int k=0;k<3;++k)
		mesh.points.push_back({{double(i),0.5*j,0.5*k}});
	std::uint64_t id=1;
	for(int j=0;j<2;++j)for(int k=0;k<2;++k){
		const auto a=node(0,j,k),b=node(1,j,k),c=node(0,j+1,k),d=node(1,j+1,k);
		const auto e=node(0,j,k+1),f=node(1,j,k+1),g=node(0,j+1,k+1),h=node(1,j+1,k+1);
		const std::array<std::array<std::uint32_t,4>,6> cells{{
			{{a,b,d,h}},{{a,d,c,h}},{{a,c,g,h}},
			{{a,g,e,h}},{{a,e,f,h}},{{a,f,b,h}}}};
		for(auto nodes:cells){
			if(iga::NativeAleDeterminant(mesh.points[nodes[0]],
				mesh.points[nodes[1]],mesh.points[nodes[2]],
				mesh.points[nodes[3]])<0.)std::swap(nodes[1],nodes[2]);
			mesh.cells.push_back({id++,nodes});
		}
	}
	using Face=std::array<std::uint32_t,3>;
	std::map<Face,std::pair<int,Face>> faces;
	for(const auto& cell:mesh.cells)for(std::size_t opposite=0;opposite<4;++opposite){
		Face oriented{};std::size_t n=0;
		for(std::size_t local=0;local<4;++local)
			if(local!=opposite)oriented[n++]=cell.nodes[local];
		auto key=oriented;std::sort(key.begin(),key.end());
		auto& face=faces[key];++face.first;face.second=oriented;
	}
	id=1;
	for(const auto& item:faces)if(item.second.first==1){
		const auto face=item.second.second;
		bool x0=true,x1=true;
		for(const auto vertex:face){
			x0=x0&&mesh.points[vertex][0]==0.;
			x1=x1&&mesh.points[vertex][0]==1.;
		}
		mesh.boundary_triangles.push_back({id++,face,x1?1:(x0?2:3)});
	}
	return mesh;
}

std::string OneDConfiguration()
{
	return R"json({
"schema_version":3,"dimension":"1d",
"geometry":{"kind":"swc_network","file":"tree.swc","length_scale_to_m":1.0},
"fields":[{"name":"area","kind":"scalar"},
 {"name":"flow_rate","kind":"scalar"},{"name":"pressure","kind":"pressure"},
 {"name":"signal","kind":"scalar","initial_value":2.0}],
"time":{"dt":0.05,"steps":2,"output_every":1},
"equation_systems":[
 {"name":"flow","kind":"network_flow_1d","unknowns":["area","flow_rate","pressure"],
  "model":"rigid","scheme":"steady_poiseuille","dynamic_viscosity":0.004,
  "density":1060.0,"discretization":{"cells_per_segment":1}},
 {"name":"transport","kind":"network_transport_1d","unknowns":["signal"],
  "flow_system":"flow","species":[{"field":"signal","diffusivity":0.0}]}],
"boundaries":[
 {"name":"inlet","role":"inlet","conditions":[
  {"field":"flow_rate","type":"dirichlet","value":0.1},
  {"field":"signal","type":"dirichlet","value":2.0}]},
 {"name":"outlet","role":"outlet","conditions":[
  {"field":"pressure","type":"pressure","value":0.0}]}],
"physiology":{"enabled":false}
})json";
}

void Close(double actual,double expected,double tolerance=1e-7)
{
	if(!std::isfinite(actual)
		||std::abs(actual-expected)>tolerance*std::max(1.,std::abs(expected)))
		throw std::runtime_error("native tetra/1D graph quantity differs");
}

std::string Digest(const std::string& value)
{
	iga::Sha256 hash;hash.Append(value.data(),value.size());return hash.Hex();
}

iga::CoupledCheckpointEpoch GraphEpoch(int ranks,
	const iga::SimulationGraph& graph,
	const std::map<std::string,std::string>& domain_models,
	const iga::SpeciesPressureFlowExecutionControls& controls)
{
	return {Digest("native-zero-d-source-species-step-0"),{},
		iga::BuildSpeciesGraphCheckpointCompatibility(graph,domain_models,
			controls,static_cast<std::uint32_t>(ranks)),1,0.05,0.05};
}

std::vector<iga::CoupledCheckpointShardSpec> GraphCatalog()
{
	return {{"controls","pressure-flow-controls-v1",0},
		{"source","zero-d-flow-species-v1",0},
		{"terminal","zero-d-flow-species-v1",0},
		{"tet-flow","native-tet-ale-flow-v1",0},
		{"tet-species","native-moving-species-v1",0}};
}

std::vector<double> GraphFingerprint(
	const iga::ZeroDSourceReservoirSpeciesDomainRuntime& source,
	const iga::NativeTetAleFlowTransportDomainAdapter& native,
	const iga::ZeroDTerminalRcrSpeciesDomainRuntime& terminal,
	const iga::SpeciesPressureFlowStepResult& result)
{
	std::vector<double> values{
		source.CommittedHydraulicState().stored_pressure_pa,
		source.CommittedSpeciesState().volume_m3,
		source.CommittedSpeciesState().amount_mol.at("tracer"),
		terminal.CommittedHydraulicState().stored_pressure_pa,
		terminal.CommittedSpeciesState().volume_m3,
		terminal.CommittedSpeciesState().amount_mol.at("tracer"),
		result.global_balances.at("tracer").residual};
	for(const auto& edge:result.edge_amounts){
		values.push_back(edge.first_outward_amount);
		values.push_back(edge.second_outward_amount);
	}
	const auto& flow=native.CommittedFlowState();
	values.insert(values.end(),flow.begin(),flow.end());
	const auto& species=native.CommittedSpeciesState();
	values.insert(values.end(),species.concentration_mol_m3.begin(),
		species.concentration_mol_m3.end());
	for(const auto& point:species.current_mesh.points)
		for(const double component:point)values.push_back(component);
	return values;
}

} // namespace

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int exit_code=0;
	try{
		int rank=0,ranks=1;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
		const bool rcr_strong=argc>1&&std::string(argv[1])=="--channel-rcr-strong";
		const bool zero_d_source_strong=argc>1
			&&std::string(argv[1])=="--channel-zero-d-source-rcr-strong";
		const bool checkpoint_save=argc>1
			&&(std::string(argv[1])=="--channel-zero-d-source-rcr-checkpoint-save"
				||std::string(argv[1])=="--channel-zero-d-source-rcr-monotone-front-checkpoint-save");
		const bool checkpoint_resume=argc>1
			&&(std::string(argv[1])=="--channel-zero-d-source-rcr-checkpoint-resume"
				||std::string(argv[1])=="--channel-zero-d-source-rcr-monotone-front-checkpoint-resume");
		const bool checkpoint=checkpoint_save||checkpoint_resume;
		const bool monotone_checkpoint=checkpoint&&std::string(argv[1]).find(
			"--channel-zero-d-source-rcr-monotone-front-")==0;
		const bool ftetwild_pipe=argc>1&&std::string(argv[1])=="--ftetwild-pipe";
		if(checkpoint&&argc!=3)
			throw std::invalid_argument("native graph checkpoint mode requires a root directory");
		std::string transport_variant;
		bool emit_concentration=false;
		bool moving_two_step=false;
		std::filesystem::path species_output_dir;
		bool monotone_species=monotone_checkpoint,sharp_front=monotone_checkpoint;
		if(ftetwild_pipe){
			if(argc<3||argc>9)
				throw std::invalid_argument("native species fTetWild mode requires a Gmsh mesh and optional test flags");
			for(int index=3;index<argc;++index){
				const std::string flag(argv[index]);
				if(flag=="--emit-concentration"){
					if(emit_concentration)throw std::invalid_argument("duplicate concentration flag");
					emit_concentration=true;
				}else if(flag=="--monotone-species"){
					if(monotone_species)throw std::invalid_argument("duplicate monotone species flag");
					monotone_species=true;
				}else if(flag=="--sharp-front"){
					if(sharp_front)throw std::invalid_argument("duplicate sharp front flag");
					sharp_front=true;
				}else if(flag=="--moving-two-step"){
					if(moving_two_step)throw std::invalid_argument("duplicate moving two-step flag");
					moving_two_step=true;
				}else if(flag.rfind("--species-output-dir=",0)==0){
					if(!species_output_dir.empty()||flag.size()==21)
						throw std::invalid_argument("native species output directory is invalid or duplicated");
					species_output_dir=flag.substr(21);
				}else if(flag=="--material-transport"
					||flag=="--transport-variant=diffusion"
					||flag=="--transport-variant=source"
					||flag=="--transport-variant=both"){
					if(!transport_variant.empty())throw std::invalid_argument("duplicate transport variant");
					transport_variant=flag;
				}else throw std::invalid_argument("unknown native species fTetWild test flag");
			}
		}
		const bool diffusion_variant=transport_variant=="--material-transport"
			||transport_variant=="--transport-variant=diffusion"
			||transport_variant=="--transport-variant=both";
		const bool source_variant=transport_variant=="--material-transport"
			||transport_variant=="--transport-variant=source"
			||transport_variant=="--transport-variant=both";
		const bool rcr_distal_strong=argc>1
			&&std::string(argv[1])=="--channel-rcr-distal-backflow-strong";
		const bool rcr_distal_reverse=rcr_distal_strong||(argc>1
			&&std::string(argv[1])=="--channel-rcr-distal-backflow");
		const bool zero_d_source=ftetwild_pipe||checkpoint||zero_d_source_strong||(argc>1
			&&std::string(argv[1])=="--channel-zero-d-source-rcr");
		const bool strong=zero_d_source_strong||rcr_strong||rcr_distal_strong
			||(argc>1&&std::string(argv[1])=="--channel-strong");
		const bool rcr=zero_d_source||rcr_strong||rcr_distal_reverse
			||(argc>1&&std::string(argv[1])=="--channel-rcr");
		const bool channel=strong||rcr||(argc>1&&std::string(argv[1])=="--channel");
		const auto directory=std::filesystem::temp_directory_path()
			/"tubularflowiga-native-t7-graph-1d";
		if(rank==0){
			std::filesystem::create_directories(directory);
			std::ofstream file(directory/"tree.swc");
			file<<"1 2 0 0 0 0.5 -1\n2 2 1 0 0 0.5 1\n";
			if(!file)throw std::runtime_error("native T7 1D fixture write failed");
		}
		MPI_Barrier(PETSC_COMM_WORLD);
		const auto configuration=iga::ParseOneDConfiguration(OneDConfiguration());
		const auto flow=configuration.flow_systems.front();
		const auto network=iga::ReadOneDNetwork(directory/"tree.swc",1.,1,
			flow.dynamic_viscosity);
		iga::OneDFlowRuntime upstream_native(configuration,flow,network,
			iga::ResolveOneDInlet(configuration),directory);
		iga::OneDFlowRuntime downstream_native(configuration,flow,network,
			iga::ResolveOneDInlet(configuration),directory);
		upstream_native.InitializeOpenLoop(0.1);
		downstream_native.InitializeOpenLoop(0.1);
		const std::vector<iga::CouplingPort> up_ports=zero_d_source
			?std::vector<iga::CouplingPort>{
				iga::MakeZeroDSourceReservoirSpeciesPort("up",{"tracer"})}
			:std::vector<iga::CouplingPort>{
			Port("up","root_observation","runtime_port","root",std::nullopt),
			Port("up","outlet","runtime_port","outlet:2",
				iga::PortQuantity::MeanPressure)};
		const std::vector<iga::CouplingPort> tet_ports{
			Port("tet","inlet","boundary_label","1",iga::PortQuantity::FlowRate),
			Port("tet","outlet","boundary_label","2",iga::PortQuantity::MeanPressure)};
		const std::vector<iga::CouplingPort> down_ports=rcr
			?std::vector<iga::CouplingPort>{
				iga::MakeZeroDTerminalRcrSpeciesPort("down",{"tracer"})}
			:std::vector<iga::CouplingPort>{
			Port("down","root","runtime_port","root",iga::PortQuantity::FlowRate),
			Port("down","outlet_observation","runtime_port","outlet:2",std::nullopt)};
		iga::SimulationGraph graph(
			{{"up",zero_d_source?iga::DomainKind::ZeroDFlow:iga::DomainKind::OneDFlow,
				up_ports,{{"tracer",zero_d_source?"tracer":"signal"}}},
			 {"tet",iga::DomainKind::ThreeDBodyFittedFlow,tet_ports,{{"tracer","tracer"}}},
			 {"down",rcr?iga::DomainKind::ZeroDFlow:iga::DomainKind::OneDFlow,
				down_ports,{{"tracer",rcr?"tracer":"signal"}}}},
			{{"up_tet",{"up",zero_d_source?"port":"outlet"},{"tet","inlet"},
				iga::CouplingLaw::PressureFlow,{"tracer"}},
			 {"tet_down",{"tet","outlet"},{"down",rcr?"port":"root"},
				iga::CouplingLaw::PressureFlow,{"tracer"}}},
			{{"tracer","mol/m^3"}});
		iga::NativeTetMesh mesh;
		if(ftetwild_pipe){
			std::ifstream input(argv[2]);
			if(!input)throw std::runtime_error("native species fTetWild mesh is unavailable");
			mesh=iga::ReadNativeTetMeshGmsh41(input);
		}else mesh=channel?ChannelMesh():Mesh();
		const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
		const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
		std::vector<double> initial_flow(3*velocity_nodes+mesh.points.size(),0.);
		for(std::size_t node=0;node<velocity_nodes;++node)
			initial_flow[3*node]=ftetwild_pipe?-0.01:(channel?-0.1:-0.2);
		bool enable_motion=!ftetwild_pipe||moving_two_step;
		std::unique_ptr<iga::CoupledDomainRuntime> up;
		iga::ZeroDSourceReservoirSpeciesDomainRuntime* source_ptr=nullptr;
		if(zero_d_source){
			iga::ZeroDFlowModel hydraulic;
			hydraulic.role=iga::ZeroDFlowRole::SourceReservoir;
			hydraulic.source=ftetwild_pipe
				?iga::ZeroDSourceReservoirModel{0.01,10.,1e-6}
				:(zero_d_source_strong||checkpoint)
				?iga::ZeroDSourceReservoirModel{1e-4,1e5,0.1}
				:iga::ZeroDSourceReservoirModel{0.01,0.1,0.1};
			iga::ZeroDSpeciesReservoirModel transport;
			transport.port_ids={"graph","pump"};
			transport.species_ids={"tracer"};
			auto source=std::make_unique<iga::ZeroDSourceReservoirSpeciesDomainRuntime>(
				"up",hydraulic,iga::ZeroDFlowState{
					(zero_d_source_strong||checkpoint)?1e4:0.},transport,
				iga::ZeroDSpeciesReservoirState{ftetwild_pipe?1e-5:1.,
					{{"tracer",ftetwild_pipe?2e-5:2.}}},
				std::map<std::string,double>{{"tracer",ftetwild_pipe?5.:2.}});
			source_ptr=source.get();up=std::move(source);
		}else up=std::make_unique<iga::OneDFlowTransportDomainAdapter>(
			"up",upstream_native,up_ports,iga::OneDInletPolicy::ConfiguredOpenLoop,
			std::map<std::string,std::string>{{"tracer","signal"}});
		std::unique_ptr<iga::CoupledDomainRuntime> down;
		iga::ZeroDTerminalRcrSpeciesDomainRuntime* rcr_ptr=nullptr;
		if(rcr){
			iga::ZeroDFlowModel hydraulic;
			hydraulic.role=iga::ZeroDFlowRole::TerminalRcr;
			hydraulic.terminal=ftetwild_pipe
				?iga::ZeroDTerminalRcrModel{1.,10.,0.01,0.}
				:iga::ZeroDTerminalRcrModel{1.,10.,1.,rcr_distal_reverse?20.:0.};
			iga::ZeroDSpeciesReservoirModel transport;
			transport.port_ids={"graph","distal"};
			transport.species_ids={"tracer"};
			auto terminal=std::make_unique<iga::ZeroDTerminalRcrSpeciesDomainRuntime>(
				"down",hydraulic,iga::ZeroDFlowState{0.},transport,
				iga::ZeroDSpeciesReservoirState{ftetwild_pipe?1e-5:1.,
					{{"tracer",ftetwild_pipe?2e-5:2.}}},
				rcr_distal_reverse?std::map<std::string,double>{{"tracer",3.}}
					:std::map<std::string,double>{});
			rcr_ptr=terminal.get();down=std::move(terminal);
		}else down=std::make_unique<iga::OneDFlowTransportDomainAdapter>(
			"down",downstream_native,down_ports,iga::OneDInletPolicy::CoupledRoot,
			std::map<std::string,std::string>{{"tracer","signal"}});
		auto native=std::make_unique<iga::NativeTetAleFlowTransportDomainAdapter>(
			"tet",mesh,tet_ports,"tracer",ftetwild_pipe
				?iga::NativeNavierStokesParameters{1.,0.001}
				:iga::NativeNavierStokesParameters{1000.,0.004},
			initial_flow,std::vector<double>(mesh.points.size(),sharp_front?0.:2.),
			diffusion_variant?1e-5:0.,source_variant?0.1:0.,
			[&enable_motion](const iga::DomainStepContext& step,
				const iga::NativeTetMesh& previous){
				auto current=previous;
				if(enable_motion&&step.step_index>0)
					for(auto& point:current.points)point[0]-=0.001*step.dt_s;
				return current;
			},ftetwild_pipe?std::set<int>{0}:(channel?std::set<int>{3}:std::set<int>{}),
			ftetwild_pipe?1e-8:1e-6,20,
			zero_d_source?"translation-x-minus-0.001-v1":"",monotone_species);
		auto* native_ptr=native.get();
		std::vector<std::unique_ptr<iga::CoupledDomainRuntime>> runtimes;
		runtimes.push_back(std::move(up));
		runtimes.push_back(std::move(native));
		runtimes.push_back(std::move(down));
		iga::DomainRuntimeRegistry registry(graph,std::move(runtimes));
		iga::SpeciesPressureFlowExecutionControls controls;
		controls.amount_tolerances.emplace("tracer",
			ftetwild_pipe?iga::SpeciesAmountTolerance{1e-12,1e-6,1e-6}
				:iga::SpeciesAmountTolerance{1e-7,1e-2,1e-6});
		controls.routing.flow_switch_m3_s=1e-10;
		controls.routing.flow_absolute_tolerance_m3_s=ftetwild_pipe?1e-10:1e-7;
		controls.hydraulic.method=strong
			?iga::PressureFlowIterationMethod::Fixed
			:iga::PressureFlowIterationMethod::Explicit;
		controls.hydraulic.maximum_iterations=strong?16:1;
		controls.hydraulic.pressure_relative_tolerance=strong?1e-4:1e-6;
		controls.hydraulic.flow_relative_tolerance=1e-6;
		iga::SpeciesPressureFlowComponentExecutor executor(registry,"up",controls,
			iga::CollectiveSpeciesPressureFlowExecution(PETSC_COMM_WORLD));
		if(ftetwild_pipe){
			const auto publish_species=[&](int step){
				if(species_output_dir.empty())return;
				const auto& accepted=native_ptr->CommittedSpeciesState();
				iga::VtkPartition piece;
				iga::CollectiveLocalStage(PETSC_COMM_WORLD,
					"native species visualization build",[&]{
						piece=iga::BuildNativeTetSpeciesVtkPartition(mesh,
							accepted.current_mesh,accepted.concentration_mol_m3,rank,ranks);
					});
				iga::WriteParallelVtkSnapshot(PETSC_COMM_WORLD,
					species_output_dir/("step_"+std::to_string(step)),piece,0.01*step);
			};
			auto result=executor.Advance({0,0.,0.01},
				{{"up_tet",0.},{"tet_down",0.}});
			publish_species(1);
			if(result.edge_amounts.size()!=2||result.domain_balances.size()!=3
				||native_ptr->CommittedSpeciesState().accepted_steps!=1
				||source_ptr->CommittedStepIndex()!=0||rcr_ptr->CommittedStepIndex()!=0)
				throw std::runtime_error("fTetWild species graph result is incomplete");
			for(const auto& edge:result.edge_amounts){
				Close(edge.residual,0.,1e-12);
				const std::string donor=edge.edge_id=="up_tet"?"up":"tet";
				const std::string receiver=edge.edge_id=="up_tet"?"tet":"down";
				const double donor_amount=edge.first.domain_id==donor
					?edge.first_outward_amount:edge.second_outward_amount;
				const double receiver_amount=edge.first.domain_id==receiver
					?edge.first_outward_amount:edge.second_outward_amount;
				if((!sharp_front||edge.edge_id=="up_tet")
					&&(!(donor_amount>1e-11)||!(receiver_amount<-1e-11))){
					std::ostringstream message;message.precision(17);
					message<<"fTetWild species graph edge has no directed throughflow: "
						<<edge.edge_id<<" first="<<edge.first_outward_amount
						<<" second="<<edge.second_outward_amount;
					throw std::runtime_error(message.str());
				}
			}
			Close(result.global_balances.at("tracer").residual,0.,1e-12);
			if(moving_two_step){
				result=executor.Advance({1,0.01,0.01},
					{{"up_tet",0.},{"tet_down",0.}});
				publish_species(2);
				if(result.edge_amounts.size()!=2||result.domain_balances.size()!=3
					||native_ptr->CommittedSpeciesState().accepted_steps!=2
					||source_ptr->CommittedStepIndex()!=1||rcr_ptr->CommittedStepIndex()!=1)
					throw std::runtime_error("fTetWild moving species graph result is incomplete");
				for(const auto& edge:result.edge_amounts)Close(edge.residual,0.,1e-12);
				Close(result.global_balances.at("tracer").residual,0.,1e-12);
				Close(native_ptr->CommittedSpeciesState().current_mesh.points[0][0]
					-mesh.points[0][0],-1e-5,1e-12);
			}
			const auto native_balance=std::find_if(result.domain_balances.begin(),
				result.domain_balances.end(),[](const auto& item){
					return item.domain_id=="tet";});
			if(native_balance==result.domain_balances.end())
				throw std::runtime_error("fTetWild species graph native balance is missing");
			const double source_amount=result.global_balances.at("tracer").source_amount;
			if(source_variant?!(source_amount>1e-10):std::abs(source_amount)>1e-12)
				throw std::runtime_error("fTetWild species graph material source amount differs");
			const auto& concentration=native_ptr->CommittedSpeciesState().concentration_mol_m3;
			const double minimum_concentration=*std::min_element(
				concentration.begin(),concentration.end());
			if(monotone_species&&minimum_concentration<-1e-12)
				throw std::runtime_error("fTetWild monotone species front is negative");
			double squared_delta=0.;
			for(const double value:concentration)
				squared_delta+=(value-(sharp_front?0.:2.))
					*(value-(sharp_front?0.:2.));
			const double rms_delta=std::sqrt(squared_delta/concentration.size());
			if(rank==0)std::cout<<"fTetWild native source/FEM/RCR species graph passed"
				<<" ranks="<<ranks<<" tetrahedra="<<mesh.cells.size()
				<<" steps="<<(moving_two_step?2:1)
				<<" diffusion="<<diffusion_variant<<" source="<<source_variant
				<<" monotone="<<monotone_species<<" sharp_front="<<sharp_front
				<<std::setprecision(17)
				<<" native_final_mol="<<native_balance->accounting.final_mass
				<<" native_minimum_mol_m3="<<minimum_concentration
				<<" native_rms_delta_mol_m3="<<rms_delta
				<<" source_mol="<<source_amount
				<<" global_residual_mol="
				<<result.global_balances.at("tracer").residual<<'\n';
			if(rank==0&&emit_concentration){
				std::cout<<"native_concentration_mol_m3=";
				for(std::size_t index=0;index<concentration.size();++index){
					if(index>0)std::cout<<',';
					std::cout<<concentration[index];
				}
				std::cout<<'\n';
			}
			PetscFinalize();return 0;
		}
		if(checkpoint){
			const auto root=std::filesystem::path(argv[2]);
			const std::map<std::string,std::string> domain_models{
				{"up",source_ptr->CheckpointModelIdentitySha256()},
				{"tet",native_ptr->CheckpointModelIdentitySha256()},
				{"down",rcr_ptr->CheckpointModelIdentitySha256()}};
			const auto epoch=GraphEpoch(ranks,graph,domain_models,controls);
			const auto catalog=GraphCatalog();
			iga::PressureFlowCheckpointControls graph_controls;
			if(checkpoint_save){
				const iga::DomainStepContext accepted{0,0.,0.05};
				if(monotone_checkpoint){
					bool rejected=false;
					try{executor.Advance(accepted,
						{{"up_tet",0.},{"tet_down",0.}},[rank](const auto&){
							if(rank==0)throw std::runtime_error("injected monotone front precommit failure");
						});}
					catch(const std::exception&){rejected=true;}
					const auto& species=native_ptr->CommittedSpeciesState();
					if(!rejected||species.accepted_steps!=0
						||species.current_mesh.points!=mesh.points
						||std::any_of(species.concentration_mol_m3.begin(),
							species.concentration_mol_m3.end(),[](double value){return value!=0.;})
						||native_ptr->CommittedFlowState()!=initial_flow
						||source_ptr->CommittedStepIndex()!=-1
						||source_ptr->CommittedSpeciesState().volume_m3!=1.
						||source_ptr->CommittedSpeciesState().amount_mol.at("tracer")!=2.
						||source_ptr->CommittedHydraulicState().stored_pressure_pa!=1e4
						||rcr_ptr->CommittedStepIndex()!=-1
						||rcr_ptr->CommittedSpeciesState().volume_m3!=1.
						||rcr_ptr->CommittedSpeciesState().amount_mol.at("tracer")!=2.
						||rcr_ptr->CommittedHydraulicState().stored_pressure_pa!=0.
						||!executor.CommittedDonorOwnership().empty())
						throw std::runtime_error("native monotone front rollback changed committed graph");
				}
				const auto first=executor.Advance(accepted,
					{{"up_tet",0.},{"tet_down",0.}});
				graph_controls=iga::MakePressureFlowCheckpointControls(accepted,
					first.hydraulic_iterations.back(),executor.CaptureCheckpointDonors());
				iga::ValidatePressureFlowCheckpointControls(graph_controls,graph,epoch);
				const auto tet=native_ptr->CaptureCheckpoints();
				const std::vector<std::string> payloads{
					iga::SerializePressureFlowCheckpointControls(graph_controls),
					iga::SerializeZeroDFlowSpeciesCheckpoint(source_ptr->CaptureCheckpointState()),
					iga::SerializeZeroDFlowSpeciesCheckpoint(rcr_ptr->CaptureCheckpointState()),
					iga::SerializeNativeTetAleFlowCheckpoint(tet.flow),
					iga::SerializeNativeTetMovingSpeciesCheckpoint(tet.species)};
				if(rank==0){
					iga::CreateCoupledCheckpointEpoch(root,epoch);
					std::vector<iga::CoupledCheckpointShard> receipts;
					for(std::size_t index=0;index<catalog.size();++index){
						const auto& bytes=payloads[index];
						receipts.push_back(iga::WriteCoupledCheckpointShard(root,epoch,
							catalog[index],bytes.size(),[&](auto& output){
								output.Write(bytes.data(),bytes.size());
							}));
					}
					iga::PublishCoupledCheckpoint(root,epoch,catalog,receipts);
				}
				MPI_Barrier(PETSC_COMM_WORLD);
			}else{
				std::vector<iga::DomainNode> altered_domains;
				for(const auto& domain:graph.Domains())
					altered_domains.push_back(domain.second);
				iga::SimulationGraph altered_graph(altered_domains,graph.Edges(),
					{{"tracer","mol/L"}});
				bool rejected_changed_graph=false;
				try{
					iga::LoadCoupledCheckpoint(root,epoch.id,
						iga::BuildSpeciesGraphCheckpointCompatibility(altered_graph,
							domain_models,controls,static_cast<std::uint32_t>(ranks)),
						catalog);
				}catch(const std::exception&){rejected_changed_graph=true;}
				if(!rejected_changed_graph)
					throw std::runtime_error("native graph changed species unit checkpoint was accepted");
				auto altered_models=domain_models;
				altered_models["up"]=Digest("changed-source-model");
				bool rejected_changed_model=false;
				try{
					iga::LoadCoupledCheckpoint(root,epoch.id,
						iga::BuildSpeciesGraphCheckpointCompatibility(graph,
							altered_models,controls,static_cast<std::uint32_t>(ranks)),
						catalog);
				}catch(const std::exception&){rejected_changed_model=true;}
				if(!rejected_changed_model)
					throw std::runtime_error("native graph changed model checkpoint was accepted");
				auto altered_controls=controls;
				altered_controls.routing.flow_switch_m3_s+=1e-9;
				bool rejected_changed_execution=false;
				try{
					iga::LoadCoupledCheckpoint(root,epoch.id,
						iga::BuildSpeciesGraphCheckpointCompatibility(graph,
							domain_models,altered_controls,static_cast<std::uint32_t>(ranks)),
						catalog);
				}catch(const std::exception&){rejected_changed_execution=true;}
				if(!rejected_changed_execution)
					throw std::runtime_error("native graph changed execution checkpoint was accepted");
				const auto manifest=iga::LoadCoupledCheckpoint(root,epoch.id,
					epoch.compatibility,catalog);
				std::vector<std::string> payloads(catalog.size());
				for(std::size_t index=0;index<catalog.size();++index)
					iga::ReadCoupledCheckpointShard(root,manifest,index,
						[&](const void* data,std::size_t count){
							payloads[index].append(static_cast<const char*>(data),count);
						});
				graph_controls=iga::ParsePressureFlowCheckpointControls(payloads[0],
					graph,manifest.epoch);
				const auto source=iga::ParseZeroDFlowSpeciesCheckpoint(payloads[1],
					manifest.epoch);
				const auto terminal=iga::ParseZeroDFlowSpeciesCheckpoint(payloads[2],
					manifest.epoch);
				iga::NativeTetAleFlowSpeciesCheckpoint tet{
					iga::ParseNativeTetAleFlowCheckpoint(payloads[3]),
					iga::ParseNativeTetMovingSpeciesCheckpoint(payloads[4])};
				if(tet.flow.accepted_steps!=epoch.accepted_steps
					||tet.flow.time_s!=epoch.time_s
					||tet.species.accepted_steps!=epoch.accepted_steps
					||tet.species.time_s!=epoch.time_s)
					throw std::runtime_error("native graph checkpoint FEM clock differs");
				source_ptr->RestoreCheckpointState(source);
				rcr_ptr->RestoreCheckpointState(terminal);
				native_ptr->RestoreCheckpoints(tet);
				executor.RestoreCheckpointDonors(graph_controls.donors);
				if(source_ptr->CommittedStepIndex()!=0
					||rcr_ptr->CommittedStepIndex()!=0
					||native_ptr->CommittedSpeciesState().accepted_steps!=1
					||executor.CommittedDonorOwnership()!=graph_controls.donors)
					throw std::runtime_error("native graph restored accepted state differs");
			}
			const auto second=executor.Advance(graph_controls.NextStep(),
				graph_controls.next_pressure_pa);
			if(monotone_checkpoint){
				const auto& concentration=native_ptr->CommittedSpeciesState().concentration_mol_m3;
				if(native_ptr->CommittedSpeciesState().accepted_steps!=2
					||*std::min_element(concentration.begin(),concentration.end())<-1e-12
					||std::abs(second.global_balances.at("tracer").residual)>1e-8)
					throw std::runtime_error("native monotone moving checkpoint front is invalid");
			}
			const auto values=GraphFingerprint(*source_ptr,*native_ptr,*rcr_ptr,second);
			if(rank==0){
				const auto reference=root/"expected-step-1.bin";
				if(checkpoint_save){
					iga::checkpoint_metadata::Writer output;
					output.Unsigned(values.size());
					for(const double value:values)output.Real(value);
					std::ofstream file(reference,std::ios::binary);
					file.write(output.Bytes().data(),output.Bytes().size());
					if(!file)throw std::runtime_error("native graph reference write failed");
				}else{
					std::ifstream file(reference,std::ios::binary);
					const std::string bytes((std::istreambuf_iterator<char>(file)),
						std::istreambuf_iterator<char>());
					iga::checkpoint_metadata::Reader input(bytes);
					if(input.Unsigned()!=values.size())
						throw std::runtime_error("native graph restart field shape differs");
					for(const double value:values)Close(value,input.Real(),1e-10);
					input.Finish();
				}
				std::cout<<"native source/FEM/RCR graph checkpoint "
					<<(checkpoint_save?"saved":"resumed")
					<<" ranks="<<ranks<<" fields="<<values.size();
				if(monotone_checkpoint)std::cout<<std::setprecision(17)
					<<" minimum_concentration_mol_m3="<<*std::min_element(
						native_ptr->CommittedSpeciesState().concentration_mol_m3.begin(),
						native_ptr->CommittedSpeciesState().concentration_mol_m3.end())
					<<" global_residual_mol="
					<<second.global_balances.at("tracer").residual;
				std::cout<<'\n';
			}
			PetscFinalize();return 0;
		}
		if(rcr){
			bool reached_precommit=false,rejected=false;
			std::string rejected_reason;
			try{
				executor.Advance({0,0.,0.05},{{"up_tet",0.},{"tet_down",0.}},
					[&reached_precommit,rank](const auto&){
						reached_precommit=true;
						if(rank==0)
							throw std::runtime_error("injected native FEM/RCR precommit failure");
					});
			}catch(const std::exception& error){
				rejected=true;rejected_reason=error.what();
			}
			if(!reached_precommit||!rejected
				||native_ptr->CommittedSpeciesState().accepted_steps!=0
				||(zero_d_source&&(source_ptr->CommittedStepIndex()!=-1
					||source_ptr->CommittedSpeciesState().volume_m3!=1.
					||source_ptr->CommittedSpeciesState().amount_mol.at("tracer")!=2.
					||source_ptr->CommittedHydraulicState().stored_pressure_pa
						!=((zero_d_source_strong||checkpoint)?1e4:0.)))
				||rcr_ptr->CommittedStepIndex()!=-1
				||rcr_ptr->CommittedSpeciesState().volume_m3!=1.
				||rcr_ptr->CommittedSpeciesState().amount_mol.at("tracer")!=2.
				||!executor.CommittedDonorOwnership().empty())
				throw std::runtime_error("native FEM/RCR rejected trial changed committed state: "
					+rejected_reason);
		}
		const auto result=executor.Advance({0,0.,0.05},
			{{"up_tet",0.},{"tet_down",0.}});
		if(result.hydraulic_iterations.empty()
			||(!strong&&result.hydraulic_iterations.size()!=1)
			||(strong&&!result.hydraulic_iterations.back().converged)
			||result.edge_amounts.size()!=2
			||result.domain_balances.size()!=3)
			throw std::runtime_error("native tetra/1D graph result is incomplete");
		for(const auto& edge:result.edge_amounts)Close(edge.residual,0.,1e-6);
		Close(result.global_balances.at("tracer").residual,0.,1e-6);
		Close(native_ptr->CommittedSpeciesState().time_s,0.05);
		if(native_ptr->CommittedSpeciesState().accepted_steps!=1)
			throw std::runtime_error("native tetra/1D graph did not commit");
		if(channel){
			const auto moving=executor.Advance({1,0.05,0.05},
				{{"up_tet",0.},{"tet_down",0.}});
			if(moving.hydraulic_iterations.empty()
				||(!strong&&moving.hydraulic_iterations.size()!=1)
				||(strong&&(moving.hydraulic_iterations.size()<2
					||!moving.hydraulic_iterations.back().converged))
				||moving.edge_amounts.size()!=2
				||moving.domain_balances.size()!=3)
				throw std::runtime_error("native channel moving graph result is incomplete");
			for(const auto& edge:moving.edge_amounts)Close(edge.residual,0.,1e-6);
			Close(moving.global_balances.at("tracer").residual,0.,1e-6);
			Close(native_ptr->CommittedSpeciesState().time_s,0.1);
			Close(native_ptr->CommittedSpeciesState().current_mesh.points[0][0],-0.00005);
			if(native_ptr->CommittedSpeciesState().accepted_steps!=2)
				throw std::runtime_error("native channel moving graph did not commit");
			if(rcr&&(!rcr_ptr||rcr_ptr->CommittedStepIndex()!=1
				||!(rcr_ptr->CommittedSpeciesState().volume_m3>0.)
				||std::abs(rcr_ptr->CommittedSpeciesState().volume_m3-1.)<1e-8
				||!(rcr_ptr->CommittedSpeciesState().amount_mol.at("tracer")>2.)
				||!rcr_ptr->CommittedHydraulicState().stored_pressure_pa))
				throw std::runtime_error("native channel RCR species did not commit");
			if(zero_d_source&&(!source_ptr||source_ptr->CommittedStepIndex()!=1
				||!(source_ptr->CommittedSpeciesState().volume_m3>0.)
				||!(source_ptr->CommittedSpeciesState().amount_mol.at("tracer")>0.)))
				throw std::runtime_error("native channel source species did not commit");
			if(zero_d_source){
				const auto source=std::find_if(moving.domain_balances.begin(),
					moving.domain_balances.end(),[](const auto& value){
						return value.domain_id=="up";});
				if(source==moving.domain_balances.end()
					||!(source->accounting.outward_port_amount.at("port")>0.)
					||std::abs(source->accounting.outward_port_amount.at("pump")+0.01)>1e-8)
					throw std::runtime_error("native channel source pump/graph species flux is absent");
			}
			if(rcr){
				const auto terminal=std::find_if(moving.domain_balances.begin(),
					moving.domain_balances.end(),[](const auto& value){
						return value.domain_id=="down";});
				if(terminal==moving.domain_balances.end()
					||!(terminal->accounting.outward_port_amount.at("port")<0.)
					||!(rcr_distal_reverse
						?terminal->accounting.outward_port_amount.at("distal")<0.
						:terminal->accounting.outward_port_amount.at("distal")>0.))
					throw std::runtime_error("native channel RCR graph/distal species flux is absent");
			}
			if(rank==0)std::cout<<"native channel moving ALE/species graph passed ranks="
				<<ranks<<" steps=2 strong="<<strong<<" rcr="<<rcr
				<<" zero_d_source="<<zero_d_source
				<<" distal_reverse="<<rcr_distal_reverse
				<<" iterations="<<moving.hydraulic_iterations.size()
				<<" global_residual_mol="
				<<moving.global_balances.at("tracer").residual<<'\n';
		}else{
		const auto native_before=iga::SerializeNativeTetMovingSpeciesCheckpoint(
			iga::CaptureNativeTetMovingSpeciesCheckpoint(
				native_ptr->CommittedSpeciesState()));
		const auto flow_before=native_ptr->CommittedFlowState();
		const auto up_before=upstream_native.FlowState();
		const auto down_before=downstream_native.FlowState();
		const auto donors_before=executor.CommittedDonorOwnership();
		bool rejected_mixed_flow=false;
		try{
			executor.Advance({1,0.05,0.05},
				{{"up_tet",0.},{"tet_down",0.}});
		}catch(const std::exception& error){
			rejected_mixed_flow=std::string(error.what()).find(
				"boundary label 2 has no concentration")!=std::string::npos;
		}
		if(!rejected_mixed_flow)
			throw std::runtime_error("moving graph did not reject unresolved local backflow");
		if(iga::SerializeNativeTetMovingSpeciesCheckpoint(
			iga::CaptureNativeTetMovingSpeciesCheckpoint(
				native_ptr->CommittedSpeciesState()))!=native_before
			||native_ptr->CommittedFlowState()!=flow_before
			||upstream_native.FlowState().completed_step!=up_before.completed_step
			||upstream_native.FlowState().physical_time!=up_before.physical_time
			||upstream_native.FlowState().flow!=up_before.flow
			||downstream_native.FlowState().completed_step!=down_before.completed_step
			||downstream_native.FlowState().physical_time!=down_before.physical_time
			||downstream_native.FlowState().flow!=down_before.flow
			||executor.CommittedDonorOwnership()!=donors_before)
			throw std::runtime_error("failed moving graph altered committed state");
		enable_motion=false;
		bool rejected_single_rank=false;
		try{
			executor.Advance({1,0.05,0.05},
				{{"up_tet",0.},{"tet_down",0.}},
				[rank](const iga::SpeciesPressureFlowStepResult&){
					if(rank==0)throw std::runtime_error("injected one-rank precommit failure");
				});
		}catch(const std::exception& error){
			rejected_single_rank=std::string(error.what()).find(
				"injected one-rank precommit failure")!=std::string::npos;
		}
		if(!rejected_single_rank)
			throw std::runtime_error("one-rank graph failure was not propagated");
		if(iga::SerializeNativeTetMovingSpeciesCheckpoint(
			iga::CaptureNativeTetMovingSpeciesCheckpoint(
				native_ptr->CommittedSpeciesState()))!=native_before
			||native_ptr->CommittedFlowState()!=flow_before
			||upstream_native.FlowState().completed_step!=up_before.completed_step
			||downstream_native.FlowState().completed_step!=down_before.completed_step
			||executor.CommittedDonorOwnership()!=donors_before)
			throw std::runtime_error("one-rank graph failure altered committed state");
		const auto second=executor.Advance({1,0.05,0.05},
			{{"up_tet",0.},{"tet_down",0.}});
		if(second.hydraulic_iterations.size()!=1
			||second.edge_amounts.size()!=2
			||second.domain_balances.size()!=3)
			throw std::runtime_error("native tetra/1D retry graph result is incomplete");
		for(const auto& edge:second.edge_amounts)Close(edge.residual,0.,1e-6);
		Close(second.global_balances.at("tracer").residual,0.,1e-6);
		Close(native_ptr->CommittedSpeciesState().time_s,0.1);
		Close(native_ptr->CommittedSpeciesState().current_mesh.points[0][0],0.);
		if(native_ptr->CommittedSpeciesState().accepted_steps!=2)
			throw std::runtime_error("native tetra/1D retry graph did not commit");
		if(rank==0)std::cout<<"native tetra ALE/1D species graph passed ranks="
			<<ranks<<" steps=2 rejected_moving_backflow=1"
			<<" rejected_one_rank=1 edges="
			<<second.edge_amounts.size()
			<<" global_residual_mol="
			<<second.global_balances.at("tracer").residual<<'\n';
		}
	}catch(const std::exception& error){
		int rank=0;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		if(rank==0)std::cerr<<error.what()<<'\n';
		exit_code=1;
	}
	PetscFinalize();return exit_code;
}
