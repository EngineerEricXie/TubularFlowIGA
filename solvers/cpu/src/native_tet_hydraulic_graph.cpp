#include "CaseConfig.hpp"
#include "CollectivePressureFlowExecution.hpp"
#include "CollectiveSpeciesPressureFlowExecution.hpp"
#include "CollectivePetscOptions.hpp"
#include "NativeTetHydraulicGraphCheckpoint.hpp"
#include "NativeTetHydraulicVisualization.hpp"
#include "NativeTetAleFlowTransportDomainAdapter.hpp"
#include "NativeTetMovingSpeciesPetscRuntime.hpp"
#include "NativeTetSpeciesVisualization.hpp"
#include "ParallelVtkOutput.hpp"
#include "OneDFlowDomainAdapter.hpp"
#include "SpeciesPressureFlowComponentExecutor.hpp"
#include "ZeroDTerminalRcrSpeciesDomainRuntime.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using Object=std::map<std::string,iga::config_detail::JsonValue>;

const iga::config_detail::JsonValue& Member(const Object& object,
	const std::string& key,const std::string& context)
{
	const auto* value=iga::config_detail::Find(object,key);
	if(!value)throw std::invalid_argument(context+" requires "+key);
	return *value;
}

const Object& Section(const Object& object,const std::string& key)
{
	return iga::config_detail::RequireObject(Member(object,key,"native case"),key);
}

double Number(const Object& object,const std::string& key,const std::string& context)
{
	return iga::config_detail::RequireNumber(Member(object,key,context),context+"."+key);
}

int Integer(const Object& object,const std::string& key,const std::string& context)
{
	return iga::config_detail::RequireInteger(Member(object,key,context),context+"."+key);
}

struct Case
{
	struct Species
	{
		std::string id;
		double initial_concentration=0.,diffusivity=0.,source=0.,decay=0.;
		bool monotone=false;
		std::map<int,double> inflow;
		std::map<int,iga::NativeTetWallExchange> wall_exchange;
	};
	struct Outlet
	{
		int label=-1;
		std::string domain_id,port_id,edge_id;
		iga::ZeroDFlowModel model;
		iga::ZeroDFlowState initial;
		double species_volume_m3=0.,species_initial_concentration=0.;
	};
	std::filesystem::path mesh_file;
	int inlet_label=-1,steps=0;
	int schema_version=1;
	std::set<int> wall_labels;
	std::vector<Outlet> outlets;
	std::vector<Species> species;
	double dt_s=0.,initial_velocity_x_m_s=0.,motion_speed_x_m_s=0.;
	std::string motion_kind="fixed";
	iga::NativeNavierStokesParameters fluid;
	iga::ZeroDFlowModel source;
	iga::ZeroDFlowState initial_source;
	std::filesystem::path source_1d_dir;
	int source_1d_outlet_node=-1;
	std::string source_1d_system;
	std::map<std::string,std::string> source_1d_species_bindings;
	iga::PressureFlowExecutionControls coupling;
	double nonlinear_tolerance=1e-8;
	std::size_t nonlinear_iterations=12;
};

Case ParseCase(const std::string& text,const std::filesystem::path& file)
{
	using namespace iga::config_detail;
	const auto parsed=JsonParser(text).Parse();
	const auto& root=RequireObject(parsed,"native case");
	RequireKnownKeys(root,{"schema_version","mesh_file","boundaries","fluid",
		"source","source_1d","terminal_rcr","terminal_rcrs","motion","coupling","time",
		"species"},"native case");
	const int schema=Integer(root,"schema_version","native case");
	if(schema!=1&&schema!=2)
		throw std::invalid_argument("unsupported native hydraulic case schema");
	Case value;
	value.schema_version=schema;
	value.mesh_file=file.parent_path()/RequireString(Member(root,"mesh_file","native case"),
		"native case.mesh_file");
	const auto& boundaries=Section(root,"boundaries");
	RequireKnownKeys(boundaries,{"inlet_label","outlet_label","wall_labels"},"boundaries");
	value.inlet_label=Integer(boundaries,"inlet_label","boundaries");
	for(const auto& item:RequireArray(Member(boundaries,"wall_labels","boundaries"),
		"boundaries.wall_labels")){
		const int label=RequireInteger(item,"boundaries.wall_labels[]");
		if(!value.wall_labels.insert(label).second)
			throw std::invalid_argument("duplicate native wall label");
	}
	if(value.inlet_label<0||value.wall_labels.count(value.inlet_label))
		throw std::invalid_argument("native inlet and wall labels overlap or are negative");
	const auto& fluid=Section(root,"fluid");
	RequireKnownKeys(fluid,{"density_kg_m3","dynamic_viscosity_pa_s",
		"initial_velocity_x_m_s","nonlinear_tolerance","maximum_iterations"},"fluid");
	value.fluid.density=Number(fluid,"density_kg_m3","fluid");
	value.fluid.dynamic_viscosity=Number(fluid,"dynamic_viscosity_pa_s","fluid");
	value.initial_velocity_x_m_s=Number(fluid,"initial_velocity_x_m_s","fluid");
	value.nonlinear_tolerance=Number(fluid,"nonlinear_tolerance","fluid");
	value.nonlinear_iterations=Integer(fluid,"maximum_iterations","fluid");
	if(const auto* source_1d=Find(root,"source_1d")){
		if(Find(root,"source"))throw std::invalid_argument("native hydraulic case has two sources");
		const auto& source=RequireObject(*source_1d,"source_1d");
		RequireKnownKeys(source,{"case_dir","outlet_node_id","system","species_bindings"},"source_1d");
		value.source_1d_dir=file.parent_path()/RequireString(
			Member(source,"case_dir","source_1d"),"source_1d.case_dir");
		value.source_1d_outlet_node=Integer(source,"outlet_node_id","source_1d");
		if(const auto* system=Find(source,"system"))
			value.source_1d_system=RequireString(*system,"source_1d.system");
		if(const auto* bindings=Find(source,"species_bindings"))
			for(const auto& binding:RequireObject(*bindings,"source_1d.species_bindings"))
				value.source_1d_species_bindings.emplace(binding.first,
					RequireString(binding.second,"source_1d species field"));
	}else{
		const auto& source=Section(root,"source");
		RequireKnownKeys(source,{"capacitance_m3_pa","resistance_pa_s_m3",
			"prescribed_flow_m3_s","initial_pressure_pa"},"source");
		value.source.role=iga::ZeroDFlowRole::SourceReservoir;
		value.source.source={Number(source,"capacitance_m3_pa","source"),
			Number(source,"resistance_pa_s_m3","source"),
			Number(source,"prescribed_flow_m3_s","source")};
		value.initial_source.stored_pressure_pa=Number(source,"initial_pressure_pa","source");
	}
	const auto parse_terminal=[&](const Object& terminal,const std::string& context,
		int label,std::string domain,std::string port,std::string edge){
		Case::Outlet outlet;
		outlet.label=label;outlet.domain_id=std::move(domain);
		outlet.port_id=std::move(port);outlet.edge_id=std::move(edge);
		outlet.model.role=iga::ZeroDFlowRole::TerminalRcr;
		outlet.model.terminal={
			Number(terminal,"proximal_resistance_pa_s_m3",context),
			Number(terminal,"distal_resistance_pa_s_m3",context),
			Number(terminal,"capacitance_m3_pa",context),
			Number(terminal,"distal_pressure_pa",context)};
			outlet.initial.stored_pressure_pa=Number(terminal,"initial_pressure_pa",context);
			if(const auto* volume=Find(terminal,"species_volume_m3"))
				outlet.species_volume_m3=RequireNumber(*volume,context+".species_volume_m3");
			if(const auto* concentration=Find(terminal,
				"species_initial_concentration_mol_m3"))
				outlet.species_initial_concentration=RequireNumber(*concentration,
					context+".species_initial_concentration_mol_m3");
		if(label<0||label==value.inlet_label||value.wall_labels.count(label))
			throw std::invalid_argument("native outlet label overlaps inlet/wall or is negative");
		value.outlets.push_back(std::move(outlet));
	};
	if(schema==1){
		if(Find(root,"terminal_rcrs")||!Find(boundaries,"outlet_label"))
			throw std::invalid_argument("schema 1 requires one outlet_label and terminal_rcr");
		const auto& terminal=Section(root,"terminal_rcr");
			RequireKnownKeys(terminal,{"proximal_resistance_pa_s_m3",
				"distal_resistance_pa_s_m3","capacitance_m3_pa","distal_pressure_pa",
				"initial_pressure_pa","species_volume_m3",
				"species_initial_concentration_mol_m3"},"terminal_rcr");
		parse_terminal(terminal,"terminal_rcr",
			Integer(boundaries,"outlet_label","boundaries"),
			"terminal","outlet","tet_rcr");
	}else{
		if(Find(root,"terminal_rcr")||Find(boundaries,"outlet_label"))
			throw std::invalid_argument("schema 2 uses terminal_rcrs, not singular outlet fields");
		const auto& terminals=RequireArray(Member(root,"terminal_rcrs","native case"),
			"terminal_rcrs");
		if(terminals.empty()||terminals.size()>16)
			throw std::invalid_argument("schema 2 requires 1 to 16 terminal RCRs");
		std::set<int> labels;
		for(const auto& item:terminals){
			const auto& terminal=RequireObject(item,"terminal_rcrs[]");
				RequireKnownKeys(terminal,{"boundary_label","proximal_resistance_pa_s_m3",
					"distal_resistance_pa_s_m3","capacitance_m3_pa",
					"distal_pressure_pa","initial_pressure_pa","species_volume_m3",
					"species_initial_concentration_mol_m3"},"terminal_rcrs[]");
			const int label=Integer(terminal,"boundary_label","terminal_rcrs[]");
			if(!labels.insert(label).second)
				throw std::invalid_argument("duplicate native terminal boundary label");
			const auto suffix=std::to_string(label);
			parse_terminal(terminal,"terminal_rcrs[]",label,
				"terminal_"+suffix,"outlet_"+suffix,"tet_rcr_"+suffix);
		}
		std::sort(value.outlets.begin(),value.outlets.end(),
			[](const auto& a,const auto& b){return a.label<b.label;});
	}
	const auto& motion=Section(root,"motion");
	RequireKnownKeys(motion,{"kind","speed_x_m_s"},"motion");
	value.motion_kind=RequireString(Member(motion,"kind","motion"),"motion.kind");
	value.motion_speed_x_m_s=Number(motion,"speed_x_m_s","motion");
	if((value.motion_kind!="fixed"&&value.motion_kind!="rigid_translation_x")
		||(value.motion_kind=="fixed"&&value.motion_speed_x_m_s!=0.))
		throw std::invalid_argument("native motion kind or speed is unsupported");
	const auto& coupling=Section(root,"coupling");
	RequireKnownKeys(coupling,{"method","maximum_iterations",
		"pressure_relative_tolerance","flow_relative_tolerance",
		"relaxation_factor"},"coupling");
	const auto method=RequireString(Member(coupling,"method","coupling"),
		"coupling.method");
	if(method=="explicit")value.coupling.method=iga::PressureFlowIterationMethod::Explicit;
	else if(method=="fixed")value.coupling.method=iga::PressureFlowIterationMethod::Fixed;
	else throw std::invalid_argument("native hydraulic coupling method is unsupported");
	value.coupling.maximum_iterations=Integer(coupling,"maximum_iterations","coupling");
	value.coupling.pressure_relative_tolerance=
		Number(coupling,"pressure_relative_tolerance","coupling");
	value.coupling.flow_relative_tolerance=
		Number(coupling,"flow_relative_tolerance","coupling");
	value.coupling.relaxation_factor=Number(coupling,"relaxation_factor","coupling");
	iga::ValidatePressureFlowExecutionControls(value.coupling);
	const auto& time=Section(root,"time");
	RequireKnownKeys(time,{"dt_s","steps"},"time");
	value.dt_s=Number(time,"dt_s","time");
	value.steps=Integer(time,"steps","time");
	if(!(value.dt_s>0.)||value.steps<1||value.nonlinear_iterations==0
		||!(value.nonlinear_tolerance>0.)
		||!(value.fluid.density>0.)||!(value.fluid.dynamic_viscosity>0.))
		throw std::invalid_argument("native hydraulic time or fluid parameters are invalid");
	if(const auto* species=Find(root,"species")){
		std::set<std::string> ids;
		for(const auto& item:RequireArray(*species,"species")){
			const auto& object=RequireObject(item,"species[]");
			RequireKnownKeys(object,{"id","initial_concentration_mol_m3",
				"diffusivity_m2_s","source_mol_m3_s",
				"first_order_decay_rate_s_inv","monotone",
				"inflow_concentration_by_label_mol_m3","wall_exchange_by_label"},
				"species[]");
			Case::Species model;
			model.id=RequireString(Member(object,"id","species[]"),"species.id");
			if(model.id.empty()||!ids.insert(model.id).second
				||model.id.find_first_not_of("abcdefghijklmnopqrstuvwxyz"
					"ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_")!=std::string::npos)
				throw std::invalid_argument("species id must be a unique ASCII identifier");
			model.initial_concentration=Number(object,"initial_concentration_mol_m3","species[]");
			model.diffusivity=Number(object,"diffusivity_m2_s","species[]");
			model.source=Number(object,"source_mol_m3_s","species[]");
			if(const auto* decay=Find(object,"first_order_decay_rate_s_inv"))
				model.decay=RequireNumber(*decay,"species decay");
			if(const auto* monotone=Find(object,"monotone"))
				model.monotone=RequireBoolean(*monotone,"species monotone");
			if(model.initial_concentration<0.||model.diffusivity<0.||model.decay<0.)
				throw std::invalid_argument("species material or initial value is negative");
			const auto& inflow=RequireObject(Member(object,
				"inflow_concentration_by_label_mol_m3","species[]"),"species inflow");
			for(const auto& entry:inflow){
				const int label=std::stoi(entry.first);
				model.inflow.emplace(label,RequireNumber(entry.second,"species inflow"));
			}
			if(const auto* exchange=Find(object,"wall_exchange_by_label"))
				for(const auto& entry:RequireObject(*exchange,"species wall exchange")){
					const auto& condition=RequireObject(entry.second,"species wall exchange label");
					model.wall_exchange.emplace(std::stoi(entry.first),iga::NativeTetWallExchange{
						Number(condition,"transfer_coefficient_m_s","species wall exchange"),
						Number(condition,"external_concentration_mol_m3","species wall exchange")});
				}
			value.species.push_back(std::move(model));
		}
	}
	return value;
}

struct Options
{
	std::filesystem::path case_file,checkpoint_root,restart_root,output_root;
	int stop_after_step=0;
	bool check_input=false;
};

Options ParseOptions(int argc,char** argv)
{
	if(argc<2)throw std::invalid_argument(
		"usage: native_tet_hydraulic_graph case.json [--checkpoint-dir DIR] "
		"[--restart-dir DIR] [--output-dir DIR] [--stop-after-step N] "
		"[--check-input]");
	Options result;result.case_file=argv[1];
	for(int i=2;i<argc;){
		const std::string name=argv[i];
		if(name=="--check-input"){
			if(result.check_input)
				throw std::invalid_argument("duplicate native hydraulic check-input option");
			result.check_input=true;++i;continue;
		}
		if(i+1>=argc)throw std::invalid_argument("native hydraulic option needs a value");
		if(name=="--checkpoint-dir"&&result.checkpoint_root.empty())
			result.checkpoint_root=argv[i+1];
		else if(name=="--restart-dir"&&result.restart_root.empty())
			result.restart_root=argv[i+1];
		else if(name=="--output-dir"&&result.output_root.empty())
			result.output_root=argv[i+1];
		else if(name=="--stop-after-step"&&result.stop_after_step==0){
			const std::string value=argv[i+1];
			std::size_t consumed=0;
			result.stop_after_step=std::stoi(value,&consumed);
			if(consumed!=value.size()||result.stop_after_step<1)
				throw std::invalid_argument("native hydraulic stop step is invalid");
		}
		else throw std::invalid_argument("unknown or duplicate native hydraulic option: "+name);
		i+=2;
	}
	if(result.check_input&&(!result.checkpoint_root.empty()
		||!result.restart_root.empty()||!result.output_root.empty()
		||result.stop_after_step!=0))
		throw std::invalid_argument("check-input cannot write or resume native hydraulic output");
	return result;
}

std::string FileText(const std::filesystem::path& path)
{
	std::ifstream input(path,std::ios::binary);
	if(!input)throw std::runtime_error("cannot open native hydraulic input: "+path.string());
	return iga::ReadCheckedText(input);
}

std::string ExecutionSha256(int ranks)
{
	iga::Sha256 hash;
	const auto binary=FileText("/proc/self/exe");
	hash.AppendLittleEndian64(binary.size());hash.Append(binary.data(),binary.size());
	char petsc[256]{};
	if(PetscGetVersion(petsc,sizeof(petsc)))
		throw std::runtime_error("cannot inspect PETSc version");
	const std::string version(petsc);
	hash.AppendLittleEndian64(version.size());hash.Append(version.data(),version.size());
	hash.AppendLittleEndian64(static_cast<std::uint64_t>(ranks));
	const auto options=iga::CapturePetscOptions(nullptr,
		{"--checkpoint-dir","--restart-dir","--output-dir",
			"--stop-after-step","--check-input"});
	hash.AppendLittleEndian64(options.size());hash.Append(options.data(),options.size());
	return hash.Hex();
}

iga::CouplingPort TetPort(const std::string& id,int label,
	iga::PortQuantity required)
{
	iga::CouplingPort port;
	port.id=id;port.subsystem_id="tet";
	port.locator_kind="boundary_label";port.locator=std::to_string(label);
	port.provides={iga::PortQuantity::Area,iga::PortQuantity::FlowRate,
		iga::PortQuantity::MeanPressure};
	port.requires={required};return port;
}

iga::CouplingPort SpeciesPort(iga::CouplingPort port,const std::string& species)
{
	port.provides.insert(iga::PortQuantity::SpeciesConcentration);
	port.provides.insert(iga::PortQuantity::SpeciesFlux);
	port.requires.insert(iga::PortQuantity::SpeciesConcentration);
	port.requires.insert(iga::PortQuantity::SpeciesFlux);
	port.species.insert(species);
	return port;
}

iga::CouplingPort OneDObservationSpeciesPort(const std::string& species)
{
	iga::CouplingPort port;
	port.id="root_observation";port.subsystem_id="source";
	port.locator_kind="runtime_port";port.locator="root";
	port.provides={iga::PortQuantity::Area,iga::PortQuantity::FlowRate,
		iga::PortQuantity::MeanPressure,iga::PortQuantity::SpeciesConcentration,
		iga::PortQuantity::SpeciesFlux};
	port.species={species};return port;
}

void RunT7SpeciesGraph(const Case& scenario,const iga::NativeTetMesh& mesh,
	const std::vector<double>& initial_flow,const Options& options,int rank,int ranks)
{
	const auto& species=scenario.species.front();
	const auto source_observation=OneDObservationSpeciesPort(species.id);
	auto source_port=SpeciesPort(TetPort("outlet",scenario.source_1d_outlet_node,
		iga::PortQuantity::MeanPressure),species.id);
	source_port.subsystem_id="source";
	source_port.locator_kind="runtime_port";
	source_port.locator="outlet:"+std::to_string(scenario.source_1d_outlet_node);
	auto inlet_port=SpeciesPort(TetPort("inlet",scenario.inlet_label,
		iga::PortQuantity::FlowRate),species.id);
	std::vector<iga::CouplingPort> tet_ports{inlet_port};
	std::vector<iga::DomainNode> nodes{{"source",iga::DomainKind::OneDFlow,
		{source_observation,source_port},
		{{species.id,scenario.source_1d_species_bindings.at(species.id)}}}};
	std::vector<iga::CouplingEdge> edges{{"source_tet",{"source","outlet"},
		{"tet","inlet"},iga::CouplingLaw::PressureFlow,{species.id}}};
	for(const auto& outlet:scenario.outlets){
		auto port=SpeciesPort(TetPort(outlet.port_id,outlet.label,
			iga::PortQuantity::MeanPressure),species.id);
		tet_ports.push_back(port);
		const auto terminal=iga::MakeZeroDTerminalRcrSpeciesPort(
			outlet.domain_id,{species.id});
		nodes.push_back({outlet.domain_id,iga::DomainKind::ZeroDFlow,{terminal},
			{{species.id,species.id}}});
		edges.push_back({outlet.edge_id,{"tet",outlet.port_id},
			{outlet.domain_id,"port"},iga::CouplingLaw::PressureFlow,{species.id}});
	}
	nodes.push_back({"tet",iga::DomainKind::ThreeDBodyFittedFlow,tet_ports,
		{{species.id,species.id}}});
	iga::SimulationGraph graph(std::move(nodes),std::move(edges),
		{{species.id,"mol/m^3"}});
	const auto config_text=FileText(scenario.source_1d_dir/"simulation_config.json");
	const auto configuration=iga::ParseOneDConfiguration(config_text);
	const iga::OneDFlowSystemDefinition* flow=nullptr;
	for(const auto& system:configuration.flow_systems)
		if(scenario.source_1d_system.empty()||system.name==scenario.source_1d_system){
			flow=&system;break;
		}
	auto network=iga::ReadOneDNetwork(scenario.source_1d_dir/
		configuration.geometry.file,configuration.geometry.length_scale_to_m,
		flow->discretization.cells_per_segment,flow->dynamic_viscosity,
		configuration.geometry.root_node_id);
	const auto inlet=iga::ResolveOneDInlet(configuration);
	const double initial_inlet=iga::EvaluateOneDInlet(configuration,inlet,
		scenario.source_1d_dir,0.,network.segments.front().area0);
	auto one_d=std::make_unique<iga::OneDFlowRuntime>(configuration,*flow,
		std::move(network),inlet,scenario.source_1d_dir);
	one_d->InitializeOpenLoop(one_d->OpenLoopInlet(0.,initial_inlet));
	iga::NativeTetAleFlowTransportDomainAdapter::MeshProvider provider;
	std::string motion_id;
	if(scenario.motion_kind=="rigid_translation_x"){
		const double speed=scenario.motion_speed_x_m_s;
		provider=[speed](const iga::DomainStepContext& step,
			const iga::NativeTetMesh& previous){
			auto current=previous;
			for(auto& point:current.points)point[0]+=speed*step.dt_s;
			return current;
		};
		std::ostringstream identity;
		identity<<"rigid_translation_x:"<<std::hexfloat<<speed;
		motion_id=identity.str();
	}
	auto source=std::make_unique<iga::OneDFlowTransportDomainAdapter>("source",
		*one_d,std::vector<iga::CouplingPort>{source_observation,source_port},
		iga::OneDInletPolicy::ConfiguredOpenLoop,scenario.source_1d_species_bindings);
	auto native=std::make_unique<iga::NativeTetAleFlowTransportDomainAdapter>(
		"tet",mesh,tet_ports,species.id,scenario.fluid,initial_flow,
		std::vector<double>(mesh.points.size(),species.initial_concentration),
		species.diffusivity,species.source,provider,scenario.wall_labels,
		scenario.nonlinear_tolerance,scenario.nonlinear_iterations,motion_id,
		species.monotone,true);
	auto* native_ptr=native.get();
	std::vector<std::unique_ptr<iga::CoupledDomainRuntime>> domains;
	domains.push_back(std::move(source));domains.push_back(std::move(native));
	for(const auto& outlet:scenario.outlets){
		iga::ZeroDSpeciesReservoirModel transport;
		transport.port_ids={"graph","distal"};transport.species_ids={species.id};
		auto terminal=std::make_unique<iga::ZeroDTerminalRcrSpeciesDomainRuntime>(
			outlet.domain_id,outlet.model,outlet.initial,transport,
			iga::ZeroDSpeciesReservoirState{outlet.species_volume_m3,
				{{species.id,outlet.species_volume_m3
					*outlet.species_initial_concentration}}});
		domains.push_back(std::move(terminal));
	}
	iga::DomainRuntimeRegistry registry(graph,std::move(domains));
	iga::SpeciesPressureFlowExecutionControls controls;
	controls.hydraulic=scenario.coupling;
	controls.routing.flow_switch_m3_s=0.;
	controls.routing.flow_absolute_tolerance_m3_s=1e-12;
	controls.routing.flow_relative_tolerance=
		std::max(1e-6,scenario.coupling.flow_relative_tolerance);
	controls.amount_tolerances.emplace(species.id,
		iga::SpeciesAmountTolerance{1e-7,1e-2,1e-6});
	iga::SpeciesPressureFlowComponentExecutor executor(registry,"source",controls,
		iga::CollectiveSpeciesPressureFlowExecution(PETSC_COMM_WORLD));
	std::vector<std::pair<double,std::filesystem::path>> flow_series,species_series;
	for(int step=0;step<scenario.steps;++step){
		std::map<std::string,double> guesses{{"source_tet",0.}};
		for(const auto& outlet:scenario.outlets)guesses.emplace(outlet.edge_id,0.);
		const auto result=executor.Advance({step,step*scenario.dt_s,scenario.dt_s},guesses);
		const auto& accepted=native_ptr->CommittedSpeciesState();
		if(!options.output_root.empty()){
			iga::VtkPartition flow_piece,species_piece;
			iga::CollectiveLocalStage(PETSC_COMM_WORLD,"T7 FEM visualization",[&]{
				flow_piece=iga::BuildNativeTetHydraulicVtkPartition(mesh,
					accepted.current_mesh,native_ptr->CommittedFlowState(),
					scenario.fluid.dynamic_viscosity,rank,ranks);
				species_piece=iga::BuildNativeTetSpeciesVtkPartition(mesh,
					accepted.current_mesh,accepted.concentration_mol_m3,rank,ranks);
			});
			const auto flow_snapshot=options.output_root/("step_"+std::to_string(step+1));
			const auto species_snapshot=options.output_root/("species_"+species.id
				+"_step_"+std::to_string(step+1));
			iga::WriteParallelVtkSnapshot(PETSC_COMM_WORLD,flow_snapshot,flow_piece,
				accepted.time_s);
			iga::WriteParallelVtkSnapshot(PETSC_COMM_WORLD,species_snapshot,
				species_piece,accepted.time_s);
			flow_series.push_back({accepted.time_s,flow_snapshot/"snapshot.pvtu"});
			species_series.push_back({accepted.time_s,species_snapshot/"snapshot.pvtu"});
			iga::WriteParallelVtkSeries(PETSC_COMM_WORLD,
				options.output_root/"flow.pvd",flow_series);
			iga::WriteParallelVtkSeries(PETSC_COMM_WORLD,
				options.output_root/("species_"+species.id+".pvd"),species_series);
		}
		if(rank==0)std::cout<<std::setprecision(17)
			<<"native_t7_species_step id="<<species.id<<" step="<<step+1
			<<" inventory_mol="
			<<iga::native_tet_moving_species_detail::Inventory(accepted.current_mesh,
				accepted.concentration_mol_m3)
			<<" global_balance_residual_mol="
			<<result.global_balances.at(species.id).residual
			<<" domains="<<result.domain_balances.size()
			<<" edges="<<result.edge_amounts.size();
		if(rank==0)for(const auto& edge:result.edge_amounts)
			std::cout<<" edge_"<<edge.edge_id<<"_residual_mol="<<edge.residual;
		if(rank==0)std::cout<<'\n';
	}
	if(rank==0)std::cout<<"native_t7_species_complete accepted_steps="
		<<native_ptr->CommittedSpeciesState().accepted_steps<<'\n';
}

} // namespace

int main(int argc,char** argv)
{
	if(PetscInitialize(&argc,&argv,nullptr,nullptr))return 1;
	int status=0,rank=0,ranks=1;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	try{
		Options options;Case scenario;iga::NativeTetMesh mesh;
		std::string case_text,mesh_text,case_sha256,execution_sha256;
		iga::CollectiveLocalStage(PETSC_COMM_WORLD,"native hydraulic inputs",[&]{
			options=ParseOptions(argc,argv);
			case_text=FileText(options.case_file);
			scenario=ParseCase(case_text,options.case_file);
			mesh_text=FileText(scenario.mesh_file);
			std::istringstream input(mesh_text);
			mesh=iga::ReadNativeTetMeshGmsh41(input);
			for(const auto& cell:mesh.cells)
				iga::EvaluateNativeTetGeometry(mesh,cell);
			std::set<int> observed_labels;
			for(const auto& face:mesh.boundary_triangles)
				observed_labels.insert(face.boundary_label);
			auto expected_labels=scenario.wall_labels;
			expected_labels.insert(scenario.inlet_label);
			for(const auto& outlet:scenario.outlets)
				expected_labels.insert(outlet.label);
			if(observed_labels!=expected_labels)
				throw std::invalid_argument("native hydraulic mesh boundary labels differ from case");
			if(options.stop_after_step<0||options.stop_after_step>scenario.steps)
				throw std::invalid_argument("native hydraulic stop step is invalid");
			iga::Sha256 hash;hash.AppendLittleEndian64(case_text.size());
			hash.Append(case_text.data(),case_text.size());
			hash.AppendLittleEndian64(mesh_text.size());
			hash.Append(mesh_text.data(),mesh_text.size());
			if(!scenario.source_1d_dir.empty()){
				const auto text=FileText(scenario.source_1d_dir/"simulation_config.json");
				const auto config=iga::ParseOneDConfiguration(text);
				if(config.time.dt!=scenario.dt_s)
					throw std::invalid_argument("1D source and FEM graph time steps differ");
				const iga::OneDFlowSystemDefinition* flow=nullptr;
				for(const auto& system:config.flow_systems)
					if(scenario.source_1d_system.empty()
						||system.name==scenario.source_1d_system){flow=&system;break;}
				if(!flow)throw std::invalid_argument("1D source flow system is unavailable");
				const auto swc=FileText(scenario.source_1d_dir/config.geometry.file);
				const auto network=iga::ReadOneDNetwork(scenario.source_1d_dir/
					config.geometry.file,config.geometry.length_scale_to_m,
					flow->discretization.cells_per_segment,flow->dynamic_viscosity,
					config.geometry.root_node_id);
				iga::ValidateOneDTopologyReferences(config,network);
				hash.AppendLittleEndian64(text.size());hash.Append(text.data(),text.size());
				hash.AppendLittleEndian64(swc.size());hash.Append(swc.data(),swc.size());
			}
			case_sha256=hash.Hex();
			execution_sha256=ExecutionSha256(ranks);
		});
		iga::RequireCollectiveSameText(PETSC_COMM_WORLD,"native hydraulic case",case_sha256);
		iga::RequireCollectiveSameText(PETSC_COMM_WORLD,"native hydraulic execution",
			execution_sha256);
		const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
		if(!scenario.species.empty()&&(!options.checkpoint_root.empty()
			||!options.restart_root.empty()))
			throw std::invalid_argument("coupled flow/species checkpoint is not available");
		if(options.check_input){
			if(rank==0)std::cout<<"native_hydraulic_input: PASS tetrahedra="
				<<mesh.cells.size()<<" boundary_labels="
				<<scenario.wall_labels.size()+1+scenario.outlets.size()<<'\n';
			PetscFinalize();return 0;
		}
		const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
			std::vector<double> initial_flow(3*velocity_nodes+mesh.points.size(),0.);
			for(std::size_t node=0;node<velocity_nodes;++node)
				initial_flow[3*node]=scenario.initial_velocity_x_m_s;
			const bool t7_species=scenario.source_1d_dir.empty()==false
				&&scenario.species.size()==1
				&&scenario.source_1d_species_bindings.count(scenario.species.front().id)
				&&std::all_of(scenario.outlets.begin(),scenario.outlets.end(),
					[](const auto& outlet){return outlet.species_volume_m3>0.;});
			if(t7_species){
				RunT7SpeciesGraph(scenario,mesh,initial_flow,options,rank,ranks);
				PetscFinalize();return 0;
			}
			std::vector<iga::CouplingPort> tet_ports{
			TetPort("inlet",scenario.inlet_label,iga::PortQuantity::FlowRate)};
		for(const auto& outlet:scenario.outlets)
			tet_ports.push_back(TetPort(outlet.port_id,outlet.label,
				iga::PortQuantity::MeanPressure));
		iga::CouplingPort source_port;
		const bool one_d_source=!scenario.source_1d_dir.empty();
		if(one_d_source){
			source_port.id="outlet";source_port.subsystem_id="source";
			source_port.locator_kind="runtime_port";
			source_port.locator="outlet:"+std::to_string(scenario.source_1d_outlet_node);
			source_port.provides={iga::PortQuantity::Area,iga::PortQuantity::FlowRate,
				iga::PortQuantity::MeanPressure};
			source_port.requires={iga::PortQuantity::MeanPressure};
		}else source_port=iga::MakeZeroDFlowPort("source",
			iga::ZeroDFlowRole::SourceReservoir);
		std::vector<iga::DomainNode> nodes{
			{"source",one_d_source?iga::DomainKind::OneDFlow:
				iga::DomainKind::ZeroDFlow,{source_port}},
			{"tet",iga::DomainKind::ThreeDBodyFittedFlow,tet_ports}};
		std::vector<iga::CouplingEdge> edges{
			{"source_tet",{"source",source_port.id},{"tet","inlet"},
				iga::CouplingLaw::PressureFlow}};
		for(const auto& outlet:scenario.outlets){
			const auto terminal_port=iga::MakeZeroDFlowPort(outlet.domain_id,
				iga::ZeroDFlowRole::TerminalRcr);
			nodes.push_back({outlet.domain_id,iga::DomainKind::ZeroDFlow,
				{terminal_port}});
			edges.push_back({outlet.edge_id,{"tet",outlet.port_id},
				{outlet.domain_id,"port"},iga::CouplingLaw::PressureFlow});
		}
		iga::SimulationGraph graph(std::move(nodes),std::move(edges));
		std::unique_ptr<iga::OneDFlowRuntime> one_d_runtime;
		std::unique_ptr<iga::CoupledDomainRuntime> source;
		iga::ZeroDFlowDomainRuntime* source_ptr=nullptr;
		if(one_d_source){
			const auto config_text=FileText(scenario.source_1d_dir/"simulation_config.json");
			const auto configuration=iga::ParseOneDConfiguration(config_text);
			if(configuration.time.dt!=scenario.dt_s)
				throw std::invalid_argument("1D source and FEM graph time steps differ");
			const iga::OneDFlowSystemDefinition* flow=nullptr;
			for(const auto& system:configuration.flow_systems)
				if(scenario.source_1d_system.empty()
					||system.name==scenario.source_1d_system){flow=&system;break;}
			if(!flow)throw std::invalid_argument("1D source flow system is unavailable");
			auto network=iga::ReadOneDNetwork(scenario.source_1d_dir/
				configuration.geometry.file,configuration.geometry.length_scale_to_m,
				flow->discretization.cells_per_segment,flow->dynamic_viscosity,
				configuration.geometry.root_node_id);
			const auto inlet=iga::ResolveOneDInlet(configuration);
			const double initial_inlet=iga::EvaluateOneDInlet(configuration,inlet,
				scenario.source_1d_dir,0.,network.segments.front().area0);
			one_d_runtime=std::make_unique<iga::OneDFlowRuntime>(configuration,*flow,
				std::move(network),inlet,scenario.source_1d_dir);
			one_d_runtime->InitializeOpenLoop(
				one_d_runtime->OpenLoopInlet(0.,initial_inlet));
			source=std::make_unique<iga::OneDFlowDomainAdapter>("source",
				*one_d_runtime,std::vector<iga::CouplingPort>{source_port},
				iga::OneDInletPolicy::ConfiguredOpenLoop);
		}else{
			auto zero_d=std::make_unique<iga::ZeroDFlowDomainRuntime>("source",
				scenario.source,scenario.initial_source,
				std::vector<iga::CouplingPort>{source_port});
			source_ptr=zero_d.get();source=std::move(zero_d);
		}
		iga::NativeTetAleFlowDomainAdapter::MeshProvider provider;
		std::string motion_id;
		if(scenario.motion_kind=="rigid_translation_x"){
			const double speed=scenario.motion_speed_x_m_s;
			provider=[speed](const iga::DomainStepContext& step,
				const iga::NativeTetMesh& previous){
				auto current=previous;
				for(auto& point:current.points)point[0]+=speed*step.dt_s;
				return current;
			};
			std::ostringstream identity;
			identity<<"rigid_translation_x:"<<std::hexfloat<<speed;
			motion_id=identity.str();
		}
		auto native=std::make_unique<iga::NativeTetAleFlowDomainAdapter>(
			"tet",mesh,tet_ports,scenario.fluid,initial_flow,provider,motion_id,
			scenario.wall_labels,scenario.nonlinear_tolerance,
			scenario.nonlinear_iterations);
		auto* native_ptr=native.get();
		std::map<std::string,iga::ZeroDFlowDomainRuntime*> zero_d_ptrs;
		if(source_ptr)zero_d_ptrs.emplace("source",source_ptr);
		std::map<int,iga::ZeroDFlowDomainRuntime*> outlet_ptrs;
		std::vector<std::unique_ptr<iga::CoupledDomainRuntime>> domains;
		domains.push_back(std::move(source));
		domains.push_back(std::move(native));
		for(const auto& outlet:scenario.outlets){
			const auto port=iga::MakeZeroDFlowPort(outlet.domain_id,
				iga::ZeroDFlowRole::TerminalRcr);
			auto terminal=std::make_unique<iga::ZeroDFlowDomainRuntime>(
				outlet.domain_id,outlet.model,outlet.initial,
				std::vector<iga::CouplingPort>{port});
			zero_d_ptrs.emplace(outlet.domain_id,terminal.get());
			outlet_ptrs.emplace(outlet.label,terminal.get());
			domains.push_back(std::move(terminal));
		}
		iga::DomainRuntimeRegistry registry(graph,std::move(domains));
		iga::PressureFlowComponentExecutor executor(registry,"source",scenario.coupling,
			iga::CollectivePressureFlowExecution(PETSC_COMM_WORLD));
		std::unique_ptr<iga::NativeTetHydraulicGraphCheckpoint> checkpoint;
		iga::CoupledCheckpointCompatibility identity;
		if(!one_d_source){
			checkpoint=std::make_unique<iga::NativeTetHydraulicGraphCheckpoint>(
				PETSC_COMM_WORLD,graph,scenario.coupling,*native_ptr,zero_d_ptrs);
			identity=checkpoint->Compatibility(case_sha256,execution_sha256);
		}else if(!options.checkpoint_root.empty()||!options.restart_root.empty())
			throw std::invalid_argument("1D/FEM graph checkpoint is unavailable");
		std::string previous_epoch;
		if(!options.restart_root.empty()){
			const auto image=checkpoint->LoadLatest(options.restart_root,identity);
			if(image.epoch.dt_s!=scenario.dt_s
				||image.epoch.accepted_steps>static_cast<std::uint64_t>(scenario.steps))
				throw std::runtime_error("native hydraulic restart horizon differs");
			checkpoint->RestoreFresh(image);previous_epoch=image.epoch.id;
		}
		const int final_step=options.stop_after_step>0
			?options.stop_after_step:scenario.steps;
		if(native_ptr->AcceptedSteps()>static_cast<std::uint64_t>(final_step))
			throw std::runtime_error("native hydraulic checkpoint exceeds stop step");
		if(!options.checkpoint_root.empty())
			iga::CollectiveLocalStage(PETSC_COMM_WORLD,"native hydraulic output root",[&]{
				if(rank==0)std::filesystem::create_directories(options.checkpoint_root);
			});
		std::vector<std::vector<double>> concentrations;
		std::vector<std::vector<std::pair<double,std::filesystem::path>>> species_series;
		for(const auto& species:scenario.species)
			concentrations.emplace_back(mesh.points.size(),species.initial_concentration);
		species_series.resize(scenario.species.size());
		for(int step=static_cast<int>(native_ptr->AcceptedSteps());step<final_step;++step){
			const auto previous_mesh=native_ptr->CommittedMesh();
			const iga::DomainStepContext context{step,native_ptr->CommittedTime(),
				scenario.dt_s};
			std::map<std::string,double> edge_guesses{{"source_tet",0.}};
			for(const auto& outlet:scenario.outlets)
				edge_guesses.emplace(outlet.edge_id,0.);
			iga::PressureFlowStepResult result;
			try{result=executor.Advance(context,edge_guesses);}
			catch(const iga::PressureFlowConvergenceError& error){
				std::ostringstream diagnostic;
				diagnostic<<error.what();
				if(!error.Result().iterations.empty()){
					const auto& last=error.Result().iterations.back();
					diagnostic<<" iteration="<<last.iteration;
					for(const auto& edge:last.edges)
						diagnostic<<" edge="<<edge.edge_id
						<<" pressure_residual_pa="<<edge.pressure_residual_pa
						<<" normalized_pressure_residual="
							<<edge.normalized_pressure_residual
						<<" normalized_flow_residual="
							<<edge.normalized_flow_residual;
				}
				throw std::runtime_error(diagnostic.str());
			}
			if(result.iterations.empty()
				||result.iterations.back().edges.size()!=1+scenario.outlets.size())
				throw std::runtime_error("native hydraulic graph result is incomplete");
			for(const auto& edge:result.iterations.back().edges)
				if(!std::isfinite(edge.flow_residual_m3_s)
					||std::abs(edge.flow_residual_m3_s)>
						std::max(1e-12,1e-4*(std::abs(edge.first_outward_flow_m3_s)
							+std::abs(edge.second_outward_flow_m3_s))))
					throw std::runtime_error("native hydraulic edge flow does not balance");
			const double inlet=*result.accepted_ports.at({"tet","inlet"})
				.outward_flow_m3_s;
			std::map<int,double> outlet_flows;
			double boundary_sum=inlet,flow_scale=std::abs(inlet);
			for(const auto& outlet:scenario.outlets){
				const double flow=*result.accepted_ports.at({"tet",outlet.port_id})
					.outward_flow_m3_s;
				outlet_flows.emplace(outlet.label,flow);
				boundary_sum+=flow;
				flow_scale+=std::abs(flow);
			}
			if(std::abs(boundary_sum)>std::max(1e-12,1e-4*flow_scale))
				throw std::runtime_error("native hydraulic 3D flow does not balance");
			for(const auto& owner:zero_d_ptrs){
				const auto& accounting=owner.second->CommittedStepAccounting();
				if(!accounting||std::abs(accounting->residual_m3)>1e-8)
					throw std::runtime_error("native hydraulic 0D storage does not balance");
			}
			if(!scenario.species.empty()){
				const auto& current_mesh=native_ptr->CommittedMesh();
				const auto& state=native_ptr->CommittedFlowState();
				std::vector<std::array<double,3>> fluid_velocity(velocity_nodes);
				for(std::size_t node=0;node<velocity_nodes;++node)
					for(int axis=0;axis<3;++axis)
						fluid_velocity[node][axis]=state[3*node+axis];
				std::vector<std::array<double,3>> mesh_velocity(mesh.points.size());
				for(std::size_t node=0;node<mesh.points.size();++node)
					for(int axis=0;axis<3;++axis)
						mesh_velocity[node][axis]=(current_mesh.points[node][axis]
							-previous_mesh.points[node][axis])/scenario.dt_s;
				std::map<std::string,double> one_d_donor;
				if(one_d_runtime){
					const auto port=one_d_runtime->GetPortState(
						"outlet:"+std::to_string(scenario.source_1d_outlet_node));
					for(const auto& binding:scenario.source_1d_species_bindings)
						one_d_donor.emplace(binding.first,
							port.concentration.at(binding.second));
				}
				for(std::size_t index=0;index<scenario.species.size();++index){
					const auto& species=scenario.species[index];
					auto inflow=species.inflow;
					const auto donor=one_d_donor.find(species.id);
					if(donor!=one_d_donor.end())inflow[scenario.inlet_label]=donor->second;
					const auto result=iga::SolveNativeTetMovingSpeciesPetscStep(
						previous_mesh,current_mesh,fluid_velocity,mesh_velocity,
						concentrations[index],inflow,species.diffusivity,
						species.source,scenario.dt_s,species.monotone,species.decay,
						species.wall_exchange);
					const double scale=std::max({1e-12,
						std::abs(result.step.current_inventory_mol/scenario.dt_s),
						std::abs(result.step.outward_advective_flux_mol_s),
						std::abs(result.step.source_mol_s)});
					if(result.converged_reason<=0
						||std::abs(result.step.balance_defect_mol_s)>1e-8*scale)
						throw std::runtime_error("native coupled species step did not conserve or converge");
					concentrations[index]=result.step.concentration_mol_m3;
					if(!options.output_root.empty()){
						iga::VtkPartition piece;
						iga::CollectiveLocalStage(PETSC_COMM_WORLD,
							"native species visualization build",[&]{
								piece=iga::BuildNativeTetSpeciesVtkPartition(mesh,
									current_mesh,concentrations[index],rank,ranks);
							});
						const auto snapshot=options.output_root/("species_"+species.id
							+"_step_"+std::to_string(step+1));
						iga::WriteParallelVtkSnapshot(PETSC_COMM_WORLD,snapshot,piece,
							native_ptr->CommittedTime());
						species_series[index].push_back({native_ptr->CommittedTime(),
							snapshot/"snapshot.pvtu"});
						iga::WriteParallelVtkSeries(PETSC_COMM_WORLD,
							options.output_root/("species_"+species.id+".pvd"),
							species_series[index]);
					}
					if(rank==0)std::cout<<std::setprecision(17)
						<<"native_flow_species_step id="<<species.id
						<<" step="<<step+1
						<<" inventory_mol="<<result.step.current_inventory_mol
						<<" balance_defect_mol_s="<<result.step.balance_defect_mol_s
						<<" inlet_donor_mol_m3="<<inflow.at(scenario.inlet_label)
						<<" linear_iterations="<<result.linear_iterations<<'\n';
				}
			}
			if(!options.output_root.empty()){
				iga::VtkPartition piece;
				iga::CollectiveLocalStage(PETSC_COMM_WORLD,
					"native hydraulic visualization build",[&]{
						piece=iga::BuildNativeTetHydraulicVtkPartition(
							native_ptr->ReferenceMesh(),
							native_ptr->CommittedMesh(),
							native_ptr->CommittedFlowState(),
							scenario.fluid.dynamic_viscosity,rank,ranks);
					});
				iga::WriteParallelVtkSnapshot(PETSC_COMM_WORLD,
					options.output_root/("step_"+std::to_string(step+1)),
					piece,native_ptr->CommittedTime());
			}
			if(rank==0){
				std::cout<<std::setprecision(17)
					<<"native_hydraulic_step step="<<step+1
					<<" time_s="<<native_ptr->CommittedTime()
					<<" iterations="<<result.iterations.size()
					<<" inlet_m3_s="<<inlet;
				if(scenario.schema_version==1){
					const auto label=scenario.outlets.front().label;
					std::cout<<" outlet_m3_s="<<outlet_flows.at(label)
						<<" source_pa="<<(one_d_source
							?*result.accepted_ports.at({"source",source_port.id}).mean_pressure_pa
							:source_ptr->CommittedState().stored_pressure_pa)
						<<" terminal_pa="<<outlet_ptrs.at(label)->CommittedState()
							.stored_pressure_pa;
				}else{
					std::cout<<" source_pa="<<(one_d_source
						?*result.accepted_ports.at({"source",source_port.id}).mean_pressure_pa
						:source_ptr->CommittedState().stored_pressure_pa);
					for(const auto& outlet:scenario.outlets)
						std::cout<<" outlet_"<<outlet.label<<"_m3_s="
							<<outlet_flows.at(outlet.label)
							<<" terminal_"<<outlet.label<<"_pa="
							<<outlet_ptrs.at(outlet.label)->CommittedState()
								.stored_pressure_pa;
				}
				std::cout<<'\n';
			}
			if(!options.checkpoint_root.empty()){
				iga::CoupledCheckpointEpoch epoch{
					iga::coupled_checkpoint_detail::Hash(case_sha256+":"
						+std::to_string(step+1)),previous_epoch,identity,
					static_cast<std::uint64_t>(step+1),native_ptr->CommittedTime(),
					scenario.dt_s};
				checkpoint->Save(options.checkpoint_root,epoch);
				previous_epoch=epoch.id;
			}
		}
		if(rank==0)std::cout<<"native_hydraulic_complete accepted_steps="
			<<native_ptr->AcceptedSteps()<<" model_sha256="
			<<native_ptr->ModelIdentitySha256()<<'\n';
		iga::CollectiveLocalStage(PETSC_COMM_WORLD,"native hydraulic text output",[&]{
			if(rank==0)iga::FlushCheckedText(std::cout);
		});
	}catch(const std::exception& error){
		if(rank==0)std::cerr<<"native_tet_hydraulic_graph: "<<error.what()<<'\n';
		status=1;
	}
	PetscFinalize();return status;
}
