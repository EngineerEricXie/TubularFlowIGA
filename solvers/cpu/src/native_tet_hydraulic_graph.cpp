#include "CaseConfig.hpp"
#include "CollectivePressureFlowExecution.hpp"
#include "CollectiveSpeciesPressureFlowExecution.hpp"
#include "CollectivePetscOptions.hpp"
#include "NativeTetHydraulicGraphCheckpoint.hpp"
#include "NativeTetHydraulicVisualization.hpp"
#include "NativeTetAleFlowTransportDomainAdapter.hpp"
#include "NativeTetMovingSpeciesPetscRuntime.hpp"
#include "NativeTetSpeciesVisualization.hpp"
#include "NativeTetTemporalOutput.hpp"
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
	// Opt-in monolithic 0D predictor for the coupling iteration's initial
	// edge pressures ("none" keeps the stateless zero guesses).
	std::string predictor="none";
	bool steady=false;
	double nonlinear_tolerance=1e-8;
	std::size_t nonlinear_iterations=12;
	// Opt-in controlled restart: case paths whose values may change between
	// a checkpoint and its restart (see MaskedCaseText).
	std::set<std::string> mutable_boundary_values;
};

const std::set<std::string>& AllowedMutableBoundaryValues()
{
	static const std::set<std::string> allowed{"source.prescribed_flow_m3_s",
		"terminal_rcrs[].proximal_resistance_pa_s_m3","terminal_rcrs[].distal_resistance_pa_s_m3",
		"terminal_rcrs[].distal_pressure_pa","terminal_rcr.proximal_resistance_pa_s_m3",
		"terminal_rcr.distal_resistance_pa_s_m3","terminal_rcr.distal_pressure_pa"};
	return allowed;
}

// Canonical, formatting-independent serialization of a parsed case. Leaves
// on a masked path (object keys joined by '.', array elements as '[]') are
// written to `values` instead of `structure`. Returns the masked leaf count.
std::size_t MaskedCaseText(const iga::config_detail::JsonValue& value,const std::string& path,
	const std::set<std::string>& masked,std::string& structure,std::string& values)
{
	using Type=iga::config_detail::JsonValue::Type;
	const auto number=[](double item){
		std::ostringstream text;text<<std::hexfloat<<item;return text.str();
	};
	if(masked.count(path)){
		if(value.type!=Type::Number)
			throw std::invalid_argument("native hydraulic mutable boundary value is not a number: "+path);
		structure+="m;";values+=path+"="+number(value.number)+";";
		return 1;
	}
	std::size_t count=0;
	switch(value.type){
	case Type::Null:structure+="n;";break;
	case Type::Boolean:structure+=value.boolean?"t;":"f;";break;
	case Type::Number:structure+="d"+number(value.number)+";";break;
	case Type::String:structure+="s"+std::to_string(value.string.size())+":"+value.string+";";break;
	case Type::Array:
		structure+="a"+std::to_string(value.array.size())+"[";
		for(const auto& item:value.array)count+=MaskedCaseText(item,path+"[]",masked,structure,values);
		structure+="]";break;
	case Type::Object:
		structure+="o"+std::to_string(value.object.size())+"{";
		for(const auto& item:value.object){
			structure+="k"+std::to_string(item.first.size())+":"+item.first;
			count+=MaskedCaseText(item.second,path.empty()?item.first:path+"."+item.first,
				masked,structure,values);
		}
		structure+="}";break;
	}
	return count;
}

// 0D model identity without the declared mutable values; it binds the
// role, capacitance and every value that must stay fixed across a restart.
std::string StructuralZeroDIdentity(const iga::ZeroDFlowModel& model,
	const std::set<std::string>& mutable_values,const std::string& prefix)
{
	iga::Sha256 hash;
	const auto text=[&](const std::string& value){
		hash.AppendLittleEndian64(value.size());hash.Append(value.data(),value.size());
	};
	const auto field=[&](const std::string& name,double value){
		text(name);
		if(mutable_values.count(prefix+name))text("mutable");
		else hash.AppendNormalizedDouble(value);
	};
	text("ZeroDFlowStructure/v1");
	hash.AppendLittleEndian32(static_cast<std::uint32_t>(model.role));
	if(model.role==iga::ZeroDFlowRole::SourceReservoir){
		field("capacitance_m3_pa",model.source.capacitance_m3_pa);
		field("resistance_pa_s_m3",model.source.resistance_pa_s_m3);
		field("prescribed_flow_m3_s",model.source.prescribed_flow_m3_s);
	}else{
		field("proximal_resistance_pa_s_m3",model.terminal.proximal_resistance_pa_s_m3);
		field("distal_resistance_pa_s_m3",model.terminal.distal_resistance_pa_s_m3);
		field("capacitance_m3_pa",model.terminal.capacitance_m3_pa);
		field("distal_pressure_pa",model.terminal.distal_pressure_pa);
	}
	return hash.Hex();
}

Case ParseCase(const std::string& text,const std::filesystem::path& file)
{
	using namespace iga::config_detail;
	const auto parsed=JsonParser(text).Parse();
	const auto& root=RequireObject(parsed,"native case");
	RequireKnownKeys(root,{"schema_version","mesh_file","boundaries","fluid",
		"source","source_1d","terminal_rcr","terminal_rcrs","motion","coupling","time",
		"species","boundary_updates"},"native case");
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
		"initial_velocity_x_m_s","nonlinear_tolerance","maximum_iterations","steady"},"fluid");
	if(const auto* steady=Find(fluid,"steady"))
		value.steady=RequireBoolean(*steady,"fluid.steady");
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
		// A sanity bound only: flux operators are sparse, the checkpoint
		// catalog allows 4096 shards, and the implicit 0D predictor makes the
		// coupling cost independent of the outlet count.
		if(terminals.empty()||terminals.size()>1024)
			throw std::invalid_argument("schema 2 requires 1 to 1024 terminal RCRs");
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
		"relaxation_factor","predictor"},"coupling");
	if(const auto* predictor=Find(coupling,"predictor")){
		value.predictor=RequireString(*predictor,"coupling.predictor");
		if(value.predictor!="none"&&value.predictor!="implicit_0d")
			throw std::invalid_argument("native hydraulic coupling predictor is unsupported");
	}
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
	if(const auto* updates=Find(root,"boundary_updates")){
		const auto& object=RequireObject(*updates,"boundary_updates");
		RequireKnownKeys(object,{"mutable"},"boundary_updates");
		for(const auto& item:RequireArray(Member(object,"mutable","boundary_updates"),
			"boundary_updates.mutable")){
			const auto path=RequireString(item,"boundary_updates.mutable[]");
			if(!AllowedMutableBoundaryValues().count(path)
				||!value.mutable_boundary_values.insert(path).second)
				throw std::invalid_argument("native hydraulic mutable boundary value is unknown "
					"or duplicate: "+path);
		}
		if(value.mutable_boundary_values.empty()||!value.source_1d_dir.empty()
			||Find(root,"species"))
			throw std::invalid_argument("native hydraulic boundary updates need a nonempty "
				"mutable list, a 0D source and no species");
	}
	if(value.steady&&value.motion_kind!="fixed")
		throw std::invalid_argument("native hydraulic steady momentum needs fixed motion");
	if((value.steady||value.predictor!="none")
		&&(!value.source_1d_dir.empty()||Find(root,"species")))
		throw std::invalid_argument("native hydraulic steady momentum and the implicit 0D "
			"predictor need a 0D source and no species");
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
	std::filesystem::path case_file,checkpoint_root,restart_root,output_root,summary_json;
	int stop_after_step=0;
	bool check_input=false;
	bool accept_boundary_update=false;
	std::string visualization_format="auto";
	bool format_set=false;
};

Options ParseOptions(int argc,char** argv)
{
	if(argc<2)throw std::invalid_argument(
		"usage: native_tet_hydraulic_graph case.json [--checkpoint-dir DIR] "
		"[--restart-dir DIR] [--output-dir DIR] [--stop-after-step N] "
		"[--visualization-format auto|vtkhdf|pvtu] [--summary-json FILE] [--check-input]");
	Options result;result.case_file=argv[1];
	for(int i=2;i<argc;){
		const std::string name=argv[i];
		if(name=="--check-input"){
			if(result.check_input)
				throw std::invalid_argument("duplicate native hydraulic check-input option");
			result.check_input=true;++i;continue;
		}
		if(name=="--accept-boundary-update"){
			if(result.accept_boundary_update)
				throw std::invalid_argument("duplicate native hydraulic boundary-update option");
			result.accept_boundary_update=true;++i;continue;
		}
		if(i+1>=argc)throw std::invalid_argument("native hydraulic option needs a value");
		if(name=="--checkpoint-dir"&&result.checkpoint_root.empty())
			result.checkpoint_root=argv[i+1];
		else if(name=="--restart-dir"&&result.restart_root.empty())
			result.restart_root=argv[i+1];
		else if(name=="--output-dir"&&result.output_root.empty())
			result.output_root=argv[i+1];
		else if(name=="--summary-json"&&result.summary_json.empty())
			result.summary_json=argv[i+1];
		else if(name=="--visualization-format"&&!result.format_set){
			result.visualization_format=argv[i+1];result.format_set=true;
			iga::NativeTetUseVtkHdf(result.visualization_format);
		}
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
		||!result.summary_json.empty()||result.stop_after_step!=0))
		throw std::invalid_argument("check-input cannot write or resume native hydraulic output");
	if(result.accept_boundary_update&&result.restart_root.empty())
		throw std::invalid_argument("--accept-boundary-update needs --restart-dir");
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
			"--stop-after-step","--check-input","--visualization-format","--summary-json",
			"--accept-boundary-update","-native_ale_monitor"});

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

// Affine relations of the 0D neighbors at the coming step, probed through the
// same trial evaluation the coupling executor uses, then one monolithic native
// solve. Returns initial edge pressures that the executor's first iteration
// reproduces when the prediction is consistent.
std::map<std::string,double> PredictImplicit0D(const Case& scenario,
	iga::NativeTetAleFlowDomainAdapter& native,const iga::ZeroDFlowDomainRuntime& source,
	const std::map<int,iga::ZeroDFlowDomainRuntime*>& terminals,
	const iga::DomainStepContext& step)
{
	const double flow_scale=std::max(1e-12,
		std::abs(source.Model().source.prescribed_flow_m3_s));
	const auto source_flow=[&](double pressure){
		return *iga::EvaluateZeroDFlowTrial(source.Model(),source.CommittedState(),pressure,
			step.dt_s,step.start_time_s).port.outward_flow_m3_s;
	};
	const auto terminal_pressure=[&](const iga::ZeroDFlowDomainRuntime& terminal,
		double native_outward_flow){
		return *iga::EvaluateZeroDFlowTrial(terminal.Model(),terminal.CommittedState(),
			-native_outward_flow,step.dt_s,step.start_time_s).port.mean_pressure_pa;
	};
	// The native inlet receives the negated source outflow.
	const double reference_pressure=1000.0;
	const double source_zero=source_flow(0.0);
	std::map<std::string,iga::NativeTetAleImplicitPortRelation> relations{{"inlet",
		{-source_zero,-(source_flow(reference_pressure)-source_zero)/reference_pressure}}};
	for(const auto& outlet:scenario.outlets){
		const auto& terminal=*terminals.at(outlet.label);
		const double zero=terminal_pressure(terminal,0.0);
		relations.emplace(outlet.port_id,iga::NativeTetAleImplicitPortRelation{zero,
			(terminal_pressure(terminal,flow_scale)-zero)/flow_scale});
	}
	const auto states=native.PredictImplicit(step,relations,1e-9*flow_scale);
	std::map<std::string,double> guesses{{"source_tet",*states.at("inlet").mean_pressure_pa}};
	for(const auto& outlet:scenario.outlets)
		guesses.emplace(outlet.edge_id,terminal_pressure(*terminals.at(outlet.label),
			*states.at(outlet.port_id).outward_flow_m3_s));
	return guesses;
}

// Peak resident set of this process (VmHWM), in bytes; 0 if unavailable.
std::uint64_t PeakResidentBytes()
{
	std::ifstream status("/proc/self/status");
	std::string line;
	while(std::getline(status,line))
		if(line.rfind("VmHWM:",0)==0){
			std::istringstream fields(line.substr(6));
			std::uint64_t kilobytes=0;fields>>kilobytes;return 1024*kilobytes;
		}
	return 0;
}

std::string JsonNumber(double value)
{
	if(!std::isfinite(value))throw std::runtime_error("native hydraulic summary value is nonfinite");
	std::ostringstream text;text<<std::setprecision(17)<<value;return text.str();
}

void RunT7SpeciesGraph(const Case& scenario,const iga::NativeTetMesh& mesh,
	const std::vector<double>& initial_flow,const Options& options,int rank)
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
	iga::NativeTetTemporalOutput flow_output(PETSC_COMM_WORLD,options.output_root,
		"flow",options.visualization_format,scenario.motion_kind!="fixed");
	iga::NativeTetTemporalOutput species_output(PETSC_COMM_WORLD,options.output_root,
		"species_"+species.id,options.visualization_format,scenario.motion_kind!="fixed");
	for(int step=0;step<scenario.steps;++step){
		std::map<std::string,double> guesses{{"source_tet",0.}};
		for(const auto& outlet:scenario.outlets)guesses.emplace(outlet.edge_id,0.);
		const auto result=executor.Advance({step,step*scenario.dt_s,scenario.dt_s},guesses);
		const auto& accepted=native_ptr->CommittedSpeciesState();
		flow_output.Append(accepted.time_s,"step_"+std::to_string(step+1),[&](int owner,int peers){
			return iga::BuildNativeTetHydraulicVtkPartition(mesh,accepted.current_mesh,
				native_ptr->CommittedFlowState(),scenario.fluid.dynamic_viscosity,owner,peers);
		});
		species_output.Append(accepted.time_s,"species_"+species.id+"_step_"+std::to_string(step+1),
			[&](int owner,int peers){return iga::BuildNativeTetSpeciesVtkPartition(mesh,
				accepted.current_mesh,accepted.concentration_mol_m3,owner,peers);});
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
	flow_output.Close();species_output.Close();
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
		std::string full_case_sha256,boundary_values_text;
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
			if(options.accept_boundary_update&&scenario.mutable_boundary_values.empty())
				throw std::invalid_argument("--accept-boundary-update needs boundary_updates in the case");
			{
				iga::Sha256 full;full.AppendLittleEndian64(case_text.size());
				full.Append(case_text.data(),case_text.size());full_case_sha256=full.Hex();
			}
			iga::Sha256 hash;
			if(scenario.mutable_boundary_values.empty()){
				hash.AppendLittleEndian64(case_text.size());
				hash.Append(case_text.data(),case_text.size());
			}else{
				// Controlled restart: the checkpoint binds the canonical case with
				// the declared boundary values masked out.
				std::string structure;
				MaskedCaseText(iga::config_detail::JsonParser(case_text).Parse(),"",
					scenario.mutable_boundary_values,structure,boundary_values_text);
				for(const auto& path:scenario.mutable_boundary_values)
					if(boundary_values_text.find(path+"=")==std::string::npos)
						throw std::invalid_argument("native hydraulic mutable boundary value is absent: "+path);
				const std::string tag="native-hydraulic-structure/v1";
				hash.AppendLittleEndian64(tag.size());hash.Append(tag.data(),tag.size());
				hash.AppendLittleEndian64(structure.size());hash.Append(structure.data(),structure.size());
			}
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
		// P2/P1 with no-slip walls: a vertex whose every tetrahedron has all ten
		// velocity nodes on a wall leaves its pressure decoupled and the system
		// singular. Fail before solving instead of inside the factorization.
		iga::CollectiveLocalStage(PETSC_COMM_WORLD,"native hydraulic wall-locked pressure",[&]{
			std::set<std::uint32_t> wall_nodes;
			for(const int label:scenario.wall_labels){
				const auto found=topology.boundary_velocity_nodes.find(label);
				if(found!=topology.boundary_velocity_nodes.end())
					wall_nodes.insert(found->second.begin(),found->second.end());
			}
			std::vector<char> free_vertex(mesh.points.size(),0);
			for(std::size_t cell=0;cell<mesh.cells.size();++cell){
				bool free=false;
				for(const auto node:topology.cell_velocity_nodes[cell])free=free||!wall_nodes.count(node);
				if(free)for(const auto vertex:mesh.cells[cell].nodes)free_vertex[vertex]=1;
			}
			std::size_t locked=0,first=0;
			for(std::size_t vertex=0;vertex<free_vertex.size();++vertex)
				if(!free_vertex[vertex]&&!locked++)first=vertex;
			if(locked){
				std::ostringstream message;
				message<<std::setprecision(9)<<"native hydraulic mesh has "<<locked
					<<" pressure vertices whose tetrahedra have every P2 velocity node on a no-slip "
					"wall (first at "<<mesh.points[first][0]<<","<<mesh.points[first][1]<<","
					<<mesh.points[first][2]<<" m); split those tetrahedra, for example with "
					"surface_to_fem_volume.py --split-boundary-locked-tetrahedra";
				throw std::invalid_argument(message.str());
			}
		});
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
				RunT7SpeciesGraph(scenario,mesh,initial_flow,options,rank);
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
			scenario.nonlinear_iterations,scenario.steady);
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
		std::optional<iga::NativeTetHydraulicBoundaryPolicy> boundary_policy;
		if(!scenario.mutable_boundary_values.empty()){
			iga::NativeTetHydraulicBoundaryPolicy policy;
			const std::string terminal_prefix=scenario.schema_version==1
				?"terminal_rcr.":"terminal_rcrs[].";
			for(const auto& item:zero_d_ptrs)
				policy.zero_d_structural_identity.emplace(item.first,StructuralZeroDIdentity(
					item.second->Model(),scenario.mutable_boundary_values,
					item.first=="source"?"source.":terminal_prefix));
			boundary_policy=std::move(policy);
		}
		if(!one_d_source){
			checkpoint=std::make_unique<iga::NativeTetHydraulicGraphCheckpoint>(
				PETSC_COMM_WORLD,graph,scenario.coupling,*native_ptr,zero_d_ptrs,boundary_policy);
			identity=checkpoint->Compatibility(case_sha256,execution_sha256);
		}else if(!options.checkpoint_root.empty()||!options.restart_root.empty())
			throw std::invalid_argument("1D/FEM graph checkpoint is unavailable");
		std::string previous_epoch;
		std::vector<std::string> updated_domains;
		std::uint64_t restored_steps=0;
		if(!options.restart_root.empty()){
			const auto image=checkpoint->LoadLatest(options.restart_root,identity);
			if(image.epoch.dt_s!=scenario.dt_s
				||image.epoch.accepted_steps>static_cast<std::uint64_t>(scenario.steps))
				throw std::runtime_error("native hydraulic restart horizon differs");
			checkpoint->RestoreFresh(image);previous_epoch=image.epoch.id;
			restored_steps=image.epoch.accepted_steps;
			for(const auto& item:zero_d_ptrs)
				if(item.second->ModelIdentitySha256()
					!=iga::BuildZeroDFlowModelIdentitySha256(item.second->AcceptedModel()))
					updated_domains.push_back(item.first);
			if(!updated_domains.empty()&&!options.accept_boundary_update){
				std::string names;
				for(const auto& name:updated_domains)names+=(names.empty()?"":",")+name;
				throw std::runtime_error("native hydraulic restart changes declared boundary values of "
					+names+"; pass --accept-boundary-update to continue from them");
			}
		}
		const auto current_boundary_values=[&]{
			std::map<std::string,iga::ZeroDFlowModel> models;
			for(const auto& item:zero_d_ptrs)models.emplace(item.first,item.second->Model());
			return models;
		};
		const std::string boundary_values_sha256=
			iga::native_tet_hydraulic_checkpoint_detail::BoundaryValuesSha256(current_boundary_values());
		const int final_step=options.stop_after_step>0
			?options.stop_after_step:scenario.steps;
		if(native_ptr->AcceptedSteps()>static_cast<std::uint64_t>(final_step))
			throw std::runtime_error("native hydraulic checkpoint exceeds stop step");
		if(!options.checkpoint_root.empty())
			iga::CollectiveLocalStage(PETSC_COMM_WORLD,"native hydraulic output root",[&]{
				if(rank==0)std::filesystem::create_directories(options.checkpoint_root);
			});
		std::vector<std::vector<double>> concentrations;
		iga::NativeTetTemporalOutput flow_output(PETSC_COMM_WORLD,options.output_root,
			"flow",options.visualization_format,scenario.motion_kind!="fixed",
			native_ptr->AcceptedSteps(),native_ptr->CommittedTime());
		std::vector<std::unique_ptr<iga::NativeTetTemporalOutput>> species_outputs;
		for(const auto& species:scenario.species)
			concentrations.emplace_back(mesh.points.size(),species.initial_concentration);
		for(const auto& species:scenario.species)
			species_outputs.push_back(std::make_unique<iga::NativeTetTemporalOutput>(
				PETSC_COMM_WORLD,options.output_root,"species_"+species.id,
				options.visualization_format,scenario.motion_kind!="fixed"));
		const std::string previous_epoch_at_start=previous_epoch;
		std::vector<std::string> summary_steps;
		const auto write_summary=[&]{
			iga::CollectiveLocalStage(PETSC_COMM_WORLD,"native hydraulic summary",[&]{
				if(rank!=0)return;
				std::ostringstream text;
				text<<"{\n \"schema_version\": 1,\n \"kind\": \"native_tet_hydraulic_summary\",\n"
					<<" \"case_sha256\": \""<<case_sha256<<"\",\n"
					<<" \"full_case_sha256\": \""<<full_case_sha256<<"\",\n"
					<<" \"boundary_values_sha256\": \""<<boundary_values_sha256<<"\",\n"
					<<" \"mutable_boundary_values\": [";
				std::size_t mutable_index=0;
				for(const auto& path:scenario.mutable_boundary_values)
					text<<(mutable_index++?", ":"")<<"\""<<path<<"\"";
				text<<"],\n \"restart\": {\"restored_accepted_steps\": "<<restored_steps
					<<", \"previous_epoch\": \""<<(options.restart_root.empty()?"":previous_epoch_at_start)
					<<"\", \"boundary_updated_domains\": [";
				for(std::size_t index=0;index<updated_domains.size();++index)
					text<<(index?", ":"")<<"\""<<updated_domains[index]<<"\"";
				text<<"]},\n \"boundary_history\": [";
				if(checkpoint&&boundary_policy){
					const auto& history=checkpoint->BoundaryHistory();
					for(std::size_t index=0;index<history.size();++index)
						text<<(index?", ":"")<<"{\"first_step\": "<<history[index].first_step
							<<", \"last_step\": "<<history[index].last_step
							<<", \"values_sha256\": \""<<history[index].values_sha256<<"\"}";
				}
				text<<"],\n"
					<<" \"execution_sha256\": \""<<execution_sha256<<"\",\n"
					<<" \"model_sha256\": \""<<native_ptr->ModelIdentitySha256()<<"\",\n"
					<<" \"ranks\": "<<ranks<<",\n \"tetrahedra\": "<<mesh.cells.size()
					<<",\n \"physical_dofs\": "<<initial_flow.size()
					<<",\n \"steady_momentum\": "<<(scenario.steady?"true":"false")
					<<",\n \"coupling_predictor\": \""<<scenario.predictor<<"\",\n \"steps\": [";
				for(std::size_t index=0;index<summary_steps.size();++index)
					text<<(index?",":"")<<"\n  "<<summary_steps[index];
				text<<"\n ]\n}\n";
				const auto temporary=options.summary_json.string()+".tmp";
				{
					std::ofstream output(temporary,std::ios::binary|std::ios::trunc);
					output<<text.str();
					if(!output.flush())throw std::runtime_error("cannot write native hydraulic summary");
				}
				std::filesystem::rename(temporary,options.summary_json);
			});
		};
		for(int step=static_cast<int>(native_ptr->AcceptedSteps());step<final_step;++step){
			const double step_start=MPI_Wtime();
			const auto previous_mesh=native_ptr->CommittedMesh();
			const iga::DomainStepContext context{step,native_ptr->CommittedTime(),
				scenario.dt_s};
			std::map<std::string,double> edge_guesses{{"source_tet",0.}};
			for(const auto& outlet:scenario.outlets)
				edge_guesses.emplace(outlet.edge_id,0.);
			if(scenario.predictor=="implicit_0d")
				edge_guesses=PredictImplicit0D(scenario,*native_ptr,*source_ptr,
					outlet_ptrs,context);
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
					species_outputs[index]->Append(native_ptr->CommittedTime(),
						"species_"+species.id+"_step_"+std::to_string(step+1),
						[&](int owner,int peers){return iga::BuildNativeTetSpeciesVtkPartition(
							mesh,current_mesh,concentrations[index],owner,peers);});
					if(rank==0)std::cout<<std::setprecision(17)
						<<"native_flow_species_step id="<<species.id
						<<" step="<<step+1
						<<" inventory_mol="<<result.step.current_inventory_mol
						<<" balance_defect_mol_s="<<result.step.balance_defect_mol_s
						<<" inlet_donor_mol_m3="<<inflow.at(scenario.inlet_label)
						<<" linear_iterations="<<result.linear_iterations<<'\n';
				}
			}
			flow_output.Append(native_ptr->CommittedTime(),"step_"+std::to_string(step+1),
				[&](int owner,int peers){return iga::BuildNativeTetHydraulicVtkPartition(
					native_ptr->ReferenceMesh(),native_ptr->CommittedMesh(),
					native_ptr->CommittedFlowState(),scenario.fluid.dynamic_viscosity,owner,peers);});
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
			const auto& statistics=native_ptr->CommittedStepStatistics();
			const double step_wall=MPI_Wtime()-step_start;
			std::uint64_t local_peak=PeakResidentBytes(),peak_max=0,peak_sum=0;
			MPI_Allreduce(&local_peak,&peak_max,1,MPI_UINT64_T,MPI_MAX,PETSC_COMM_WORLD);
			MPI_Allreduce(&local_peak,&peak_sum,1,MPI_UINT64_T,MPI_SUM,PETSC_COMM_WORLD);
			const auto& inlet_state=result.accepted_ports.at({"tet","inlet"});
			if(rank==0)std::cout<<std::setprecision(17)
				<<"native_hydraulic_solver step="<<step+1
				<<" coupling_iterations="<<result.iterations.size()
				<<" trial_solves="<<statistics.trial_solves
				<<" newton_iterations="<<statistics.newton_iterations
				<<" linear_iterations="<<statistics.linear_iterations
				<<" predictor="<<(statistics.predictor?1:0)
				<<" predictor_newton_iterations="<<statistics.predictor_newton_iterations
				<<" predictor_linear_iterations="<<statistics.predictor_linear_iterations
				<<" assembly_s="<<statistics.assembly_seconds+statistics.predictor_assembly_seconds
				<<" linear_solve_s="<<statistics.linear_solve_seconds
					+statistics.predictor_linear_solve_seconds
				<<" wall_s="<<step_wall
				<<" peak_rss_max_bytes="<<peak_max
				<<" inlet_mean_pressure_pa="<<*inlet_state.mean_pressure_pa<<'\n';
			if(!options.summary_json.empty()){
				std::ostringstream record;
				record<<std::setprecision(17)<<"{\"step\": "<<step+1
					<<", \"time_s\": "<<JsonNumber(native_ptr->CommittedTime())
					<<", \"coupling_iterations\": "<<result.iterations.size()
					<<", \"inlet\": {\"label\": "<<scenario.inlet_label
					<<", \"outward_flow_m3_s\": "<<JsonNumber(inlet)
					<<", \"total_inflow_m3_s\": "<<JsonNumber(-inlet)
					<<", \"mean_pressure_pa\": "<<JsonNumber(*inlet_state.mean_pressure_pa)
					<<", \"source_stored_pressure_pa\": "<<JsonNumber(one_d_source
						?*result.accepted_ports.at({"source",source_port.id}).mean_pressure_pa
						:source_ptr->CommittedState().stored_pressure_pa)
					<<"}, \"outlets\": [";
				double outlet_total=0.;
				for(std::size_t index=0;index<scenario.outlets.size();++index){
					const auto& outlet=scenario.outlets[index];
					const auto& state=result.accepted_ports.at({"tet",outlet.port_id});
					const auto& terminal=result.accepted_ports.at({outlet.domain_id,"port"});
					outlet_total+=outlet_flows.at(outlet.label);
					record<<(index?", ":"")<<"{\"label\": "<<outlet.label
						<<", \"outward_flow_m3_s\": "<<JsonNumber(outlet_flows.at(outlet.label))
						<<", \"mean_pressure_pa\": "<<JsonNumber(*state.mean_pressure_pa)
						<<", \"terminal_port_pressure_pa\": "<<JsonNumber(*terminal.mean_pressure_pa)
						<<", \"terminal_stored_pressure_pa\": "<<JsonNumber(outlet_ptrs.at(outlet.label)
							->CommittedState().stored_pressure_pa)<<"}";
				}
				record<<"], \"total_outlet_flow_m3_s\": "<<JsonNumber(outlet_total)
					<<", \"solver\": {\"trial_solves\": "<<statistics.trial_solves
					<<", \"newton_iterations\": "<<statistics.newton_iterations
					<<", \"linear_iterations\": "<<statistics.linear_iterations
					<<", \"assembly_s\": "<<JsonNumber(statistics.assembly_seconds)
					<<", \"linear_solve_s\": "<<JsonNumber(statistics.linear_solve_seconds)
					<<", \"ksp_type\": \""<<statistics.linear_solver_type
					<<"\", \"pc_type\": \""<<statistics.preconditioner_type
					<<"\", \"predictor\": "<<(statistics.predictor?"true":"false")
					<<", \"predictor_newton_iterations\": "<<statistics.predictor_newton_iterations
					<<", \"predictor_linear_iterations\": "<<statistics.predictor_linear_iterations
					<<", \"predictor_assembly_s\": "<<JsonNumber(statistics.predictor_assembly_seconds)
					<<", \"predictor_linear_solve_s\": "<<JsonNumber(statistics.predictor_linear_solve_seconds)
					<<"}, \"wall_time_s\": "<<JsonNumber(step_wall)
					<<", \"peak_rss_bytes_max_rank\": "<<peak_max
					<<", \"peak_rss_bytes_sum\": "<<peak_sum
					<<", \"boundary_values\": {\"values_sha256\": \""<<boundary_values_sha256
					<<"\", \"updated_at_this_step\": "<<(static_cast<std::uint64_t>(step)==restored_steps
						&&!updated_domains.empty()?"true":"false")
					<<", \"source_prescribed_flow_m3_s\": "
					<<JsonNumber(one_d_source?0.0:source_ptr->Model().source.prescribed_flow_m3_s)
					<<", \"terminals\": [";
				for(std::size_t index=0;index<scenario.outlets.size();++index){
					const auto& model=outlet_ptrs.at(scenario.outlets[index].label)->Model().terminal;
					record<<(index?", ":"")<<"{\"label\": "<<scenario.outlets[index].label
						<<", \"proximal_resistance_pa_s_m3\": "<<JsonNumber(model.proximal_resistance_pa_s_m3)
						<<", \"distal_resistance_pa_s_m3\": "<<JsonNumber(model.distal_resistance_pa_s_m3)
						<<", \"capacitance_m3_pa\": "<<JsonNumber(model.capacitance_m3_pa)
						<<", \"distal_pressure_pa\": "<<JsonNumber(model.distal_pressure_pa)<<"}";
				}
				record<<"]}}";
				summary_steps.push_back(record.str());
			}
			if(!options.checkpoint_root.empty()){
				// Controlled restarts may branch from one epoch with different
				// boundary values, so their epoch ids also bind those values.
				iga::CoupledCheckpointEpoch epoch{
					iga::coupled_checkpoint_detail::Hash(case_sha256+":"
						+(boundary_policy?boundary_values_sha256+":":std::string())
						+std::to_string(step+1)),previous_epoch,identity,
					static_cast<std::uint64_t>(step+1),native_ptr->CommittedTime(),
					scenario.dt_s};
				checkpoint->Save(options.checkpoint_root,epoch);
				previous_epoch=epoch.id;
			}
			if(!options.summary_json.empty())write_summary();

		}
		flow_output.Close();
		for(auto& output:species_outputs)output->Close();
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
