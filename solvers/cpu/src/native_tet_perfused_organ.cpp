// Perfused-organ step solver: steady tree/Darcy hemodynamics plus quasi-steady
// multispecies tissue transport with hepatic kinetics. See
// docs/PERFUSED_ORGAN.md for the case, request and field-file contracts.

#include "CaseConfig.hpp"
#include "CheckedText.hpp"
#include "PerfusedOrganHemodynamics.hpp"
#include "PerfusedOrganKinetics.hpp"
#include "PerfusedOrganTransport.hpp"
#include "Sha256.hpp"

#include <petscsys.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using iga::config_detail::JsonValue;
using Object=std::map<std::string,JsonValue>;

const JsonValue& Field(const Object& object,const std::string& key,const std::string& context)
{
	const auto* value=iga::config_detail::Find(object,key);
	if(!value)throw std::invalid_argument(context+" requires "+key);
	return *value;
}

double Number(const Object& object,const std::string& key,const std::string& context)
{
	return iga::config_detail::RequireNumber(Field(object,key,context),context+"."+key);
}

double NumberOr(const Object& object,const std::string& key,double fallback,const std::string& context)
{
	const auto* value=iga::config_detail::Find(object,key);
	return value?iga::config_detail::RequireNumber(*value,context+"."+key):fallback;
}

std::string Hash(const std::string& value)
{
	iga::Sha256 digest;
	digest.Append(value.data(),value.size());
	return digest.Hex();
}

std::string ReadText(const std::filesystem::path& path,const char* what)
{
	std::ifstream stream(path);
	if(!stream)throw std::runtime_error(std::string("cannot open ")+what+": "+path.string());
	return iga::ReadCheckedText(stream);
}

std::string Escape(const std::string& text)
{
	std::string result;
	for(const char c:text){
		if(c=='"'||c=='\\'){result+='\\';result+=c;}
		else if(c=='\n')result+="\\n";
		else if(static_cast<unsigned char>(c)<0x20)result+=' ';
		else result+=c;
	}
	return result;
}

// Minimal JSON writer with full double precision.
class Json
{
public:
	Json(){out_<<std::setprecision(17);}
	Json& Open(){Separate();out_<<'{';first_=true;return *this;}
	Json& Close(){out_<<'}';first_=false;return *this;}
	Json& OpenArray(){Separate();out_<<'[';first_=true;return *this;}
	Json& CloseArray(){out_<<']';first_=false;return *this;}
	Json& Key(const std::string& key){Separate();out_<<'"'<<Escape(key)<<"\":";first_=true;return *this;}
	Json& Value(double value)
	{
		Separate();
		if(std::isfinite(value))out_<<value;else out_<<"null";
		return *this;
	}
	Json& Value(std::size_t value){Separate();out_<<value;return *this;}
	Json& Value(long value){Separate();out_<<value;return *this;}
	Json& Value(bool value){Separate();out_<<(value?"true":"false");return *this;}
	Json& Value(const std::string& value){Separate();out_<<'"'<<Escape(value)<<'"';return *this;}
	Json& Value(const char* value){return Value(std::string(value));}
	template<typename T>Json& Entry(const std::string& key,const T& value){Key(key);return Value(value);}
	std::string Text() const{return out_.str();}
private:
	void Separate(){if(!first_)out_<<',';first_=false;}
	std::ostringstream out_;
	bool first_=true;
};

struct TreeSpec
{
	std::string id;
	bool inflow=true;
	std::filesystem::path skeleton;
	double length_scale=1.;
};

struct Organ
{
	std::filesystem::path case_path;
	iga::NativeTetMesh mesh;
	std::vector<iga::PerfusedOrganTree> trees;
	std::vector<std::map<int,int>> terminal_swc_id; // per tree: node index -> SWC id
	std::vector<double> volume;
	std::vector<std::array<double,3>> centroid;
	double total_volume=0.;
	std::string case_hash,mesh_hash,territories_hash;
	std::vector<std::string> skeleton_hashes;
	std::string reaction_model;
	std::size_t zones=1; // lobular zones per cell
};

Organ LoadOrgan(const std::filesystem::path& case_path)
{
	using namespace iga::config_detail;
	Organ organ;
	organ.case_path=std::filesystem::absolute(case_path);
	const auto case_text=ReadText(organ.case_path,"perfused organ case");
	organ.case_hash=Hash(case_text);
	const auto parsed=JsonParser(case_text).Parse();
	const auto& root=RequireObject(parsed,"perfused organ case");
	RequireKnownKeys(root,{"schema_version","kind","tissue_mesh_file","mesh_length_scale_to_m",
		"trees","territories_file","reaction_model","lobular_zones","description"},"perfused organ case");
	if(RequireInteger(Field(root,"schema_version","perfused organ case"),"schema_version")!=1)
		throw std::invalid_argument("unsupported perfused organ case schema");
	if(RequireString(Field(root,"kind","perfused organ case"),"kind")!="perfused_organ")
		throw std::invalid_argument("perfused organ case kind must be perfused_organ");
	organ.reaction_model=RequireString(Field(root,"reaction_model","perfused organ case"),"reaction_model");
	if(organ.reaction_model!="hepatic_v1")
		throw std::invalid_argument("unsupported reaction_model "+organ.reaction_model);
	if(const auto* zones=Find(root,"lobular_zones"))
		organ.zones=static_cast<std::size_t>(RequireInteger(*zones,"lobular_zones"));
	if(organ.zones<1||organ.zones>16)throw std::invalid_argument("lobular_zones must be 1..16");
	const auto base=organ.case_path.parent_path();
	const auto mesh_path=base/RequireString(Field(root,"tissue_mesh_file","perfused organ case"),"tissue_mesh_file");
	const double mesh_scale=NumberOr(root,"mesh_length_scale_to_m",1.,"perfused organ case");
	if(!(mesh_scale>0.))throw std::invalid_argument("mesh_length_scale_to_m must be positive");
	const auto mesh_text=ReadText(mesh_path,"tissue mesh");
	organ.mesh_hash=Hash(mesh_text);
	{
		std::istringstream input(mesh_text);
		organ.mesh=iga::ReadNativeTetMeshGmsh41(input);
	}
	for(auto& point:organ.mesh.points)for(auto& x:point)x*=mesh_scale;
	if(organ.mesh.cells.empty())throw std::invalid_argument("tissue mesh has no tetrahedra");
	std::unordered_map<std::uint64_t,std::size_t> cell_index;
	organ.volume.resize(organ.mesh.cells.size());
	organ.centroid.resize(organ.mesh.cells.size());
	for(std::size_t i=0;i<organ.mesh.cells.size();++i){
		if(!cell_index.emplace(organ.mesh.cells[i].id,i).second)
			throw std::invalid_argument("tissue mesh repeats a cell id");
		organ.volume[i]=iga::PerfusedOrganCellVolume(organ.mesh,i);
		organ.total_volume+=organ.volume[i];
		for(const auto node:organ.mesh.cells[i].nodes)
			for(int axis=0;axis<3;++axis)organ.centroid[i][axis]+=0.25*organ.mesh.points[node][axis];
	}

	std::vector<TreeSpec> specs;
	for(const auto& item:RequireArray(Field(root,"trees","perfused organ case"),"trees")){
		const auto& tree=RequireObject(item,"tree");
		RequireKnownKeys(tree,{"id","role","skeleton_file","length_scale_to_m"},"tree");
		TreeSpec spec;
		spec.id=RequireString(Field(tree,"id","tree"),"tree.id");
		const auto role=RequireString(Field(tree,"role","tree"),"tree.role");
		if(role!="inflow"&&role!="outflow")throw std::invalid_argument("tree role must be inflow or outflow");
		spec.inflow=role=="inflow";
		spec.skeleton=base/RequireString(Field(tree,"skeleton_file","tree"),"tree.skeleton_file");
		spec.length_scale=NumberOr(tree,"length_scale_to_m",1.,"tree");
		for(const auto& other:specs)
			if(other.id==spec.id)throw std::invalid_argument("tree ids must be unique");
		specs.push_back(spec);
	}
	if(specs.empty())throw std::invalid_argument("perfused organ case needs at least one tree");
	if(std::none_of(specs.begin(),specs.end(),[](const TreeSpec& s){return s.inflow;})
		||std::all_of(specs.begin(),specs.end(),[](const TreeSpec& s){return s.inflow;}))
		throw std::invalid_argument("perfused organ case needs inflow and outflow trees");

	const auto territories_path=base/RequireString(Field(root,"territories_file","perfused organ case"),
		"territories_file");
	const auto territories_text=ReadText(territories_path,"territories");
	organ.territories_hash=Hash(territories_text);
	const auto territories_parsed=JsonParser(territories_text).Parse();
	const auto& territories_root=RequireObject(territories_parsed,"territories");
	if(RequireInteger(Field(territories_root,"schema_version","territories"),"schema_version")!=1)
		throw std::invalid_argument("unsupported territories schema");
	const auto& by_tree=RequireObject(Field(territories_root,"trees","territories"),"territories.trees");
	for(const auto& spec:specs){
		organ.skeleton_hashes.push_back(Hash(ReadText(spec.skeleton,"tree skeleton")));
		iga::PerfusedOrganTree tree;
		tree.id=spec.id;
		tree.inflow=spec.inflow;
		// Viscosity is applied per request; the reader needs a placeholder.
		tree.network=iga::ReadOneDNetwork(spec.skeleton,spec.length_scale,1,1e-3);
		const auto found=by_tree.find(spec.id);
		if(found==by_tree.end())throw std::invalid_argument("territories missing tree "+spec.id);
		std::map<int,int> swc_ids;
		for(const auto& item:RequireObject(found->second,"territories tree")){
			std::size_t used=0;
			int swc_id=0;
			try{swc_id=std::stoi(item.first,&used);}catch(const std::exception&){used=0;}
			if(used!=item.first.size())throw std::invalid_argument("territory key must be an SWC id");
			const auto node=tree.network.node_index.find(swc_id);
			if(node==tree.network.node_index.end())
				throw std::invalid_argument("territory names unknown SWC node "+item.first);
			const auto& entry=RequireObject(item.second,"territory");
			RequireKnownKeys(entry,{"cell_ids","weights"},"territory");
			const auto& ids=RequireArray(Field(entry,"cell_ids","territory"),"cell_ids");
			const auto& weights=RequireArray(Field(entry,"weights","territory"),"weights");
			if(ids.size()!=weights.size())throw std::invalid_argument("territory cell_ids and weights differ");
			auto& list=tree.territories[node->second];
			for(std::size_t k=0;k<ids.size();++k){
				const double id=RequireNumber(ids[k],"cell_ids");
				const auto cell=cell_index.find(static_cast<std::uint64_t>(id));
				if(id<1.||std::floor(id)!=id||cell==cell_index.end())
					throw std::invalid_argument("territory cell id is not in the tissue mesh");
				list.push_back({cell->second,RequireNumber(weights[k],"weights")});
			}
			swc_ids[node->second]=swc_id;
		}
		iga::ValidatePerfusedOrganTree(organ.mesh,tree);
		organ.trees.push_back(std::move(tree));
		organ.terminal_swc_id.push_back(swc_ids);
	}
	return organ;
}

struct Fields
{
	std::size_t cells=0;
	std::map<std::string,std::vector<double>> values;
};

// Per-cell fields have one value per tetrahedron; per-node fields one value per
// lobular zone node (cell-major: node = cell*zones + zone).
bool PerCellField(const std::string& name)
{
	static const std::set<std::string> names={"mobility_scale","exchange_scale","tissue_pressure_pa",
		"exchange_in_m3_s","exchange_out_m3_s","perfusion_s_inv","darcy_flux_m_s"};
	return names.count(name)!=0;
}

Fields ReadFields(const Object& spec,std::size_t cells,std::size_t zones,const std::filesystem::path& base)
{
	using namespace iga::config_detail;
	RequireKnownKeys(spec,{"file","names"},"field_input");
	Fields fields;
	fields.cells=cells;
	auto path=std::filesystem::path(RequireString(Field(spec,"file","field_input"),"field_input.file"));
	if(path.is_relative())path=base/path;
	std::vector<std::string> names;
	for(const auto& item:RequireArray(Field(spec,"names","field_input"),"field_input.names"))
		names.push_back(RequireString(item,"field_input.names"));
	std::ifstream input(path,std::ios::binary);
	if(!input)throw std::runtime_error("cannot open field input "+path.string());
	for(const auto& name:names){
		std::vector<double> values(PerCellField(name)?cells:cells*zones);
		input.read(reinterpret_cast<char*>(values.data()),
			static_cast<std::streamsize>(values.size()*sizeof(double)));
		if(!input)throw std::runtime_error("field input is shorter than its named fields");
		for(const double value:values)
			if(!std::isfinite(value))throw std::invalid_argument("field input "+name+" is not finite");
		if(!fields.values.emplace(name,std::move(values)).second)
			throw std::invalid_argument("field input repeats "+name);
	}
	char extra=0;
	if(input.read(&extra,1))throw std::invalid_argument("field input is longer than its named fields");
	return fields;
}

iga::OrganConcentration SpeciesObject(const Object& object,const std::string& context)
{
	iga::OrganConcentration result{};
	const auto& names=iga::PerfusedOrganSpeciesNames();
	std::set<std::string> allowed(names.begin(),names.end());
	iga::config_detail::RequireKnownKeys(object,allowed,context);
	for(std::size_t s=0;s<iga::OrganSpeciesCount;++s)result[s]=Number(object,names[s],context);
	return result;
}

iga::HepaticKinetics ParseKinetics(const Object& object)
{
	iga::config_detail::RequireKnownKeys(object,{"vmax_o2_mol_m3_s","km_o2_mmhg","atp_per_o2",
		"glycolysis_capacity_mol_m3_s","glycogen_half_saturation_mol_m3","km_glucose_mol_m3",
		"lactate_vmax_mol_m3_s","lactate_km_mol_m3","gluconeogenic_fraction",
		"glycogenolysis_vmax_mol_m3_s","glucose_uptake_vmax_mol_m3_s","insulin_half_effect_mU_L",
		"insulin_clearance_s_inv","respiratory_quotient","q10","temperature_c"},"kinetics");
	iga::HepaticKinetics k;
	const std::string c="kinetics";
	k.vmax_o2_mol_m3_s=Number(object,"vmax_o2_mol_m3_s",c);
	k.glycolysis_capacity_mol_m3_s=Number(object,"glycolysis_capacity_mol_m3_s",c);
	k.lactate_vmax_mol_m3_s=Number(object,"lactate_vmax_mol_m3_s",c);
	k.glycogenolysis_vmax_mol_m3_s=Number(object,"glycogenolysis_vmax_mol_m3_s",c);
	k.glucose_uptake_vmax_mol_m3_s=Number(object,"glucose_uptake_vmax_mol_m3_s",c);
	k.insulin_clearance_s_inv=Number(object,"insulin_clearance_s_inv",c);
	k.temperature_c=Number(object,"temperature_c",c);
	k.km_o2_mmhg=NumberOr(object,"km_o2_mmhg",k.km_o2_mmhg,c);
	k.atp_per_o2=NumberOr(object,"atp_per_o2",k.atp_per_o2,c);
	k.glycogen_half_saturation_mol_m3=NumberOr(object,"glycogen_half_saturation_mol_m3",
		k.glycogen_half_saturation_mol_m3,c);
	k.km_glucose_mol_m3=NumberOr(object,"km_glucose_mol_m3",k.km_glucose_mol_m3,c);
	k.lactate_km_mol_m3=NumberOr(object,"lactate_km_mol_m3",k.lactate_km_mol_m3,c);
	k.gluconeogenic_fraction=NumberOr(object,"gluconeogenic_fraction",k.gluconeogenic_fraction,c);
	k.insulin_half_effect_mU_L=NumberOr(object,"insulin_half_effect_mU_L",k.insulin_half_effect_mU_L,c);
	k.respiratory_quotient=NumberOr(object,"respiratory_quotient",k.respiratory_quotient,c);
	k.q10=NumberOr(object,"q10",k.q10,c);
	return k;
}

iga::PerfusateChemistry ParseChemistry(const Object& object)
{
	iga::config_detail::RequireKnownKeys(object,{"haemoglobin_g_dl","o2_solubility_mol_m3_mmhg",
		"p50_mmhg","hill_coefficient","bohr_coefficient","co2_solubility_mol_m3_mmhg",
		"bicarbonate_pk","non_bicarbonate_buffer_mol_m3_ph","buffer_anion_at_7_4_mol_m3"},"chemistry");
	iga::PerfusateChemistry chemistry;
	const std::string c="chemistry";
	chemistry.haemoglobin_g_dl=Number(object,"haemoglobin_g_dl",c);
	chemistry.o2_solubility_mol_m3_mmhg=NumberOr(object,"o2_solubility_mol_m3_mmhg",
		chemistry.o2_solubility_mol_m3_mmhg,c);
	chemistry.p50_mmhg=NumberOr(object,"p50_mmhg",chemistry.p50_mmhg,c);
	chemistry.hill_coefficient=NumberOr(object,"hill_coefficient",chemistry.hill_coefficient,c);
	chemistry.bohr_coefficient=NumberOr(object,"bohr_coefficient",chemistry.bohr_coefficient,c);
	chemistry.co2_solubility_mol_m3_mmhg=NumberOr(object,"co2_solubility_mol_m3_mmhg",
		chemistry.co2_solubility_mol_m3_mmhg,c);
	chemistry.bicarbonate_pk=NumberOr(object,"bicarbonate_pk",chemistry.bicarbonate_pk,c);
	chemistry.non_bicarbonate_buffer_mol_m3_ph=NumberOr(object,"non_bicarbonate_buffer_mol_m3_ph",
		chemistry.non_bicarbonate_buffer_mol_m3_ph,c);
	chemistry.buffer_anion_at_7_4_mol_m3=NumberOr(object,"buffer_anion_at_7_4_mol_m3",
		chemistry.buffer_anion_at_7_4_mol_m3,c);
	return chemistry;
}

struct Quantiles{double mean=0.,min=0.,p05=0.,p50=0.,p95=0.,max=0.;};

Quantiles VolumeQuantiles(const std::vector<double>& values,const std::vector<double>& volume)
{
	Quantiles q;
	std::vector<std::size_t> order(values.size());
	std::iota(order.begin(),order.end(),0);
	std::sort(order.begin(),order.end(),[&](std::size_t a,std::size_t b){return values[a]<values[b];});
	double total=0.,weighted=0.;
	for(std::size_t i=0;i<values.size();++i){total+=volume[i];weighted+=volume[i]*values[i];}
	q.mean=weighted/total;
	q.min=values[order.front()];
	q.max=values[order.back()];
	const auto at=[&](double fraction){
		double running=0.;
		for(const auto i:order){running+=volume[i];if(running>=fraction*total)return values[i];}
		return values[order.back()];
	};
	q.p05=at(0.05);q.p50=at(0.5);q.p95=at(0.95);
	return q;
}

void WriteQuantiles(Json& json,const std::string& key,const Quantiles& q)
{
	json.Key(key).Open().Entry("mean",q.mean).Entry("min",q.min).Entry("p05",q.p05)
		.Entry("p50",q.p50).Entry("p95",q.p95).Entry("max",q.max).Close();
}

void WriteLegacyVtk(const std::filesystem::path& path,const Organ& organ,
	const std::map<std::string,std::vector<double>>& fields)
{
	std::ofstream out(path);
	if(!out)throw std::runtime_error("cannot write VTK "+path.string());
	out<<std::setprecision(9);
	out<<"# vtk DataFile Version 3.0\nperfused organ\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	out<<"POINTS "<<organ.mesh.points.size()<<" double\n";
	for(const auto& p:organ.mesh.points)out<<p[0]<<' '<<p[1]<<' '<<p[2]<<'\n';
	const auto cells=organ.mesh.cells.size();
	out<<"CELLS "<<cells<<' '<<cells*5<<'\n';
	for(const auto& cell:organ.mesh.cells)
		out<<"4 "<<cell.nodes[0]<<' '<<cell.nodes[1]<<' '<<cell.nodes[2]<<' '<<cell.nodes[3]<<'\n';
	out<<"CELL_TYPES "<<cells<<'\n';
	for(std::size_t i=0;i<cells;++i)out<<"10\n";
	out<<"CELL_DATA "<<cells<<'\n';
	for(const auto& field:fields){
		out<<"SCALARS "<<field.first<<" double 1\nLOOKUP_TABLE default\n";
		for(const double value:field.second)out<<value<<'\n';
	}
	if(!out)throw std::runtime_error("VTK write failed "+path.string());
}

std::string Solve(const Organ& organ,const std::string& request_text,const std::filesystem::path& base)
{
	using namespace iga::config_detail;
	const auto started=std::chrono::steady_clock::now();
	const auto parsed=JsonParser(request_text).Parse();
	const auto& request=RequireObject(parsed,"request");
	RequireKnownKeys(request,{"dt_s","viscosity_pa_s","root_pressure_pa",
		"specific_conductance_inv_pa_s","mobility_m2_pa_s","inlet","initial",
		"distribution_volume_fraction","chemistry","kinetics","viability","glycogen_mol_m3",
		"field_input","field_output","vtk_file","hypoxic_po2_mmhg","tag","resolution"},"request");
	const std::size_t cells=organ.mesh.cells.size(),zones=organ.zones,nodes=cells*zones;
	Fields fields;
	if(const auto* spec=Find(request,"field_input"))
		fields=ReadFields(RequireObject(*spec,"field_input"),cells,zones,base);
	const auto field_or=[&](const std::string& name,double fallback){
		const auto found=fields.values.find(name);
		return found==fields.values.end()
			?std::vector<double>(PerCellField(name)?cells:nodes,fallback):found->second;
	};
	const auto& species_names=iga::PerfusedOrganSpeciesNames();
	{
		std::set<std::string> known={"mobility_scale","exchange_scale","viability","glycogen_mol_m3"};
		known.insert(species_names.begin(),species_names.end());
		for(const auto& item:fields.values)
			if(!known.count(item.first))throw std::invalid_argument("unknown input field "+item.first);
	}

	iga::PerfusedOrganFlowInput flow;
	flow.viscosity_pa_s=Number(request,"viscosity_pa_s","request");
	const auto& pressures=RequireObject(Field(request,"root_pressure_pa","request"),"root_pressure_pa");
	const auto& conductances=RequireObject(Field(request,"specific_conductance_inv_pa_s","request"),
		"specific_conductance_inv_pa_s");
	for(const auto& tree:organ.trees){
		flow.root_pressure_pa.push_back(Number(pressures,tree.id,"root_pressure_pa"));
		flow.specific_conductance_inv_pa_s.push_back(Number(conductances,tree.id,
			"specific_conductance_inv_pa_s"));
	}
	if(pressures.size()!=organ.trees.size()||conductances.size()!=organ.trees.size())
		throw std::invalid_argument("root pressures and conductances must name every tree once");
	const double mobility=Number(request,"mobility_m2_pa_s","request");
	const auto scale=field_or("mobility_scale",1.);
	flow.mobility_by_cell_m2_pa_s.resize(cells);
	for(std::size_t i=0;i<cells;++i)flow.mobility_by_cell_m2_pa_s[i]=mobility*scale[i];
	if(fields.values.count("exchange_scale"))flow.exchange_scale_by_cell=fields.values.at("exchange_scale");
	const auto flow_started=std::chrono::steady_clock::now();
	const auto hemodynamics=iga::SolvePerfusedOrganFlow(organ.mesh,organ.trees,flow);
	const auto flow_finished=std::chrono::steady_clock::now();

	iga::OrganTransportInput transport;
	transport.dt_s=Number(request,"dt_s","request");
	transport.inlet=SpeciesObject(RequireObject(Field(request,"inlet","request"),"inlet"),"inlet");
	transport.distribution_volume_fraction=SpeciesObject(RequireObject(
		Field(request,"distribution_volume_fraction","request"),"distribution_volume_fraction"),
		"distribution_volume_fraction");
	transport.kinetics=ParseKinetics(RequireObject(Field(request,"kinetics","request"),"kinetics"));
	transport.chemistry=ParseChemistry(RequireObject(Field(request,"chemistry","request"),"chemistry"));
	iga::OrganConcentration initial=transport.inlet;
	if(const auto* value=Find(request,"initial"))initial=SpeciesObject(RequireObject(*value,"initial"),"initial");
	transport.previous.assign(nodes,initial);
	for(std::size_t s=0;s<iga::OrganSpeciesCount;++s){
		const auto found=fields.values.find(species_names[s]);
		if(found==fields.values.end())continue;
		for(std::size_t i=0;i<nodes;++i){
			if(found->second[i]<0.)throw std::invalid_argument("previous concentration is negative");
			transport.previous[i][s]=found->second[i];
		}
	}
	const auto viability=field_or("viability",NumberOr(request,"viability",1.,"request"));
	const auto glycogen=field_or("glycogen_mol_m3",NumberOr(request,"glycogen_mol_m3",0.,"request"));
	transport.state.resize(nodes);
	for(std::size_t i=0;i<nodes;++i){
		if(viability[i]<0.||viability[i]>1.||glycogen[i]<0.)
			throw std::invalid_argument("viability must be in [0,1] and glycogen nonnegative");
		transport.state[i]={viability[i],glycogen[i]};
	}
	// Transport resolution: every cell, one well-mixed pool per territory of the
	// first inflow tree, or one pool for the organ. Per-node inputs are
	// volume-averaged into pools and pool results are broadcast back.
	const std::string resolution=Find(request,"resolution")
		?RequireString(*Find(request,"resolution"),"resolution"):"cell";
	std::vector<std::size_t> group(cells,0);
	std::size_t groups=cells;
	if(resolution=="cell")std::iota(group.begin(),group.end(),0);
	else if(resolution=="organ")groups=1;
	else if(resolution=="portal_territory"){
		const auto inflow=std::find_if(organ.trees.begin(),organ.trees.end(),
			[](const iga::PerfusedOrganTree& tree){return tree.inflow;});
		std::vector<char> assigned(cells,0);
		groups=0;
		for(const auto& territory:inflow->territories){
			for(const auto& entry:territory.second){
				if(assigned[entry.first])throw std::invalid_argument("inflow territories overlap");
				assigned[entry.first]=1;
				group[entry.first]=groups;
			}
			++groups;
		}
		if(std::find(assigned.begin(),assigned.end(),0)!=assigned.end())
			throw std::invalid_argument("inflow territories do not cover every cell");
	}else throw std::invalid_argument("resolution must be cell, portal_territory or organ");
	const auto cell_graph=iga::BuildTetTransportGraph(organ.mesh,
		hemodynamics.darcy.conservative_face_flow_m3_s,hemodynamics.cell_exchange_in_m3_s,
		hemodynamics.cell_exchange_out_m3_s);
	const bool lumped=resolution!="cell";
	const auto graph=iga::SplitIntoLobularZones(
		lumped?iga::AggregateTransportGraph(cell_graph,group,groups):cell_graph,zones);
	if(lumped){
		iga::OrganTransportInput pooled=transport;
		pooled.previous.assign(groups*zones,iga::OrganConcentration{});
		pooled.state.assign(groups*zones,{0.,0.});
		std::vector<double> pool_volume(groups*zones,0.);
		for(std::size_t i=0;i<cells;++i)
			for(std::size_t z=0;z<zones;++z){
				const std::size_t node=i*zones+z,pool=group[i]*zones+z;
				const double v=organ.volume[i];
				pool_volume[pool]+=v;
				for(std::size_t s=0;s<iga::OrganSpeciesCount;++s)
					pooled.previous[pool][s]+=v*transport.previous[node][s];
				pooled.state[pool].viability+=v*transport.state[node].viability;
				pooled.state[pool].glycogen_mol_m3+=v*transport.state[node].glycogen_mol_m3;
			}
		for(std::size_t pool=0;pool<groups*zones;++pool){
			for(auto& value:pooled.previous[pool])value/=pool_volume[pool];
			pooled.state[pool].viability/=pool_volume[pool];
			pooled.state[pool].glycogen_mol_m3/=pool_volume[pool];
		}
		transport=pooled;
	}
	const auto pooled_species=iga::SolveOrganTransport(graph,transport);
	// Per-node view (cell-major, zone-minor) for outputs and summaries.
	iga::OrganTransportResult species=pooled_species;
	std::vector<double> node_volume(nodes);
	for(std::size_t i=0;i<cells;++i)
		for(std::size_t z=0;z<zones;++z)node_volume[i*zones+z]=organ.volume[i]/double(zones);
	if(lumped){
		species.concentration.resize(nodes);
		species.rates.resize(nodes);
		for(std::size_t i=0;i<cells;++i)
			for(std::size_t z=0;z<zones;++z){
				species.concentration[i*zones+z]=pooled_species.concentration[group[i]*zones+z];
				species.rates[i*zones+z]=pooled_species.rates[group[i]*zones+z];
			}
	}
	const auto transport_finished=std::chrono::steady_clock::now();

	// Cell output fields.
	std::map<std::string,std::vector<double>> out;
	for(std::size_t s=0;s<iga::OrganSpeciesCount;++s){
		auto& values=out[species_names[s]];
		values.resize(nodes);
		for(std::size_t i=0;i<nodes;++i)values[i]=species.concentration[i][s];
	}
	const auto rate_field=[&](const char* name,auto getter){
		auto& values=out[name];
		values.resize(nodes);
		for(std::size_t i=0;i<nodes;++i)values[i]=getter(species.rates[i]);
	};
	rate_field("po2_mmHg",[](const iga::HepaticRates& r){return r.po2_mmhg;});
	rate_field("ph",[](const iga::HepaticRates& r){return r.ph;});
	rate_field("pco2_mmHg",[](const iga::HepaticRates& r){return r.pco2_mmhg;});
	rate_field("bicarbonate_mol_m3",[](const iga::HepaticRates& r){return r.bicarbonate;});
	rate_field("o2_uptake_mol_m3_s",[](const iga::HepaticRates& r){return r.o2_uptake;});
	rate_field("lactate_production_mol_m3_s",[](const iga::HepaticRates& r){return r.lactate_production;});
	rate_field("lactate_clearance_mol_m3_s",[](const iga::HepaticRates& r){return r.lactate_clearance;});
	rate_field("glycogen_rate_mol_m3_s",[](const iga::HepaticRates& r){return r.glycogen_rate;});
	rate_field("energy_deficit",[](const iga::HepaticRates& r){return r.energy_deficit;});
	out["tissue_pressure_pa"]=hemodynamics.cell_mean_pressure_pa;
	out["exchange_in_m3_s"]=hemodynamics.cell_exchange_in_m3_s;
	out["exchange_out_m3_s"]=hemodynamics.cell_exchange_out_m3_s;
	{
		auto& perfusion=out["perfusion_s_inv"];
		auto& flux=out["darcy_flux_m_s"];
		perfusion.resize(cells);flux.resize(cells);
		for(std::size_t i=0;i<cells;++i){
			perfusion[i]=hemodynamics.cell_exchange_in_m3_s[i]/organ.volume[i];
			const auto& q=hemodynamics.darcy.cell_flux_m_s[i];
			flux[i]=std::sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]);
		}
	}

	const double hypoxic_threshold=NumberOr(request,"hypoxic_po2_mmhg",5.,"request");
	double hypoxic=0.,o2_uptake=0.;
	for(std::size_t i=0;i<nodes;++i){
		if(species.rates[i].po2_mmhg<hypoxic_threshold)hypoxic+=node_volume[i];
		o2_uptake+=node_volume[i]*species.rates[i].o2_uptake;
	}
	// Last-zone (pericentral) values, one per cell.
	const auto last_zone=[&](const std::string& name){
		std::vector<double> values(cells);
		for(std::size_t i=0;i<cells;++i)values[i]=out.at(name)[i*zones+zones-1];
		return values;
	};

	Json json;
	json.Open().Entry("ok",true);
	if(const auto* tag=Find(request,"tag"))json.Entry("tag",RequireString(*tag,"tag"));
	json.Entry("cells",cells).Entry("lobular_zones",zones).Entry("resolution",resolution)
		.Entry("transport_pools",groups).Entry("volume_m3",organ.total_volume);
	json.Key("trees").Open();
	for(std::size_t t=0;t<organ.trees.size();++t){
		const auto& tree=organ.trees[t];
		json.Key(tree.id).Open().Entry("role",tree.inflow?"inflow":"outflow")
			.Entry("root_pressure_pa",flow.root_pressure_pa[t])
			.Entry("root_inflow_m3_s",hemodynamics.root_inflow_m3_s[t]);
		json.Key("terminal_flow_m3_s").Open();
		for(const auto& item:hemodynamics.terminal_flow_m3_s[t])
			json.Entry(std::to_string(organ.terminal_swc_id[t].at(item.first)),item.second);
		json.Close();
		json.Key("terminal_pressure_pa").Open();
		for(const auto& item:hemodynamics.terminal_flow_m3_s[t])
			json.Entry(std::to_string(organ.terminal_swc_id[t].at(item.first)),
				hemodynamics.tree_pressure_pa[t][static_cast<std::size_t>(item.first)]);
		json.Close();
		// Concentration of the flow leaving tissue into this tree, mixed.
		iga::OrganConcentration mixed{};
		double leaving=0.;
		for(const auto& exchange:hemodynamics.exchanges){
			if(exchange.tree!=t||exchange.flow_m3_s>=0.)continue;
			leaving-=exchange.flow_m3_s;
			for(std::size_t s=0;s<iga::OrganSpeciesCount;++s)
				mixed[s]-=exchange.flow_m3_s*species.concentration[exchange.cell*zones+zones-1][s];
		}
		json.Entry("tissue_to_tree_flow_m3_s",leaving);
		json.Key("tissue_to_tree_concentration").Open();
		for(std::size_t s=0;s<iga::OrganSpeciesCount;++s)
			json.Entry(species_names[s],leaving>0.?mixed[s]/leaving:std::numeric_limits<double>::quiet_NaN());
		json.Close().Close();
	}
	json.Close();
	json.Key("species_totals").Open();
	for(std::size_t s=0;s<iga::OrganSpeciesCount;++s)
		json.Key(species_names[s]).Open().Entry("inflow",species.inflow[s])
			.Entry("outflow",species.outflow[s]).Entry("storage_rate",species.storage_rate[s])
			.Entry("reaction",species.reaction[s]).Entry("balance_defect",species.balance_defect[s]).Close();
	json.Close();
	json.Key("tissue").Open();
	for(const char* name:{"po2_mmHg","ph","pco2_mmHg","energy_deficit","lactate_mol_m3","glucose_mol_m3"})
		WriteQuantiles(json,name,VolumeQuantiles(out.at(name),node_volume));
	for(const char* name:{"perfusion_s_inv","tissue_pressure_pa"})
		WriteQuantiles(json,name,VolumeQuantiles(out.at(name),organ.volume));
	WriteQuantiles(json,"last_zone_po2_mmHg",VolumeQuantiles(last_zone("po2_mmHg"),organ.volume));
	WriteQuantiles(json,"last_zone_energy_deficit",VolumeQuantiles(last_zone("energy_deficit"),organ.volume));
	json.Entry("hypoxic_po2_mmhg",hypoxic_threshold)
		.Entry("hypoxic_volume_fraction",hypoxic/organ.total_volume)
		.Entry("o2_uptake_mol_s",o2_uptake).Close();
	const auto seconds=[](auto a,auto b){return std::chrono::duration<double>(b-a).count();};
	json.Key("diagnostics").Open()
		.Entry("flow_residual_norm",hemodynamics.residual_norm)
		.Entry("flow_linear_iterations",static_cast<long>(hemodynamics.linear_iterations))
		.Entry("exchange_imbalance_m3_s",hemodynamics.exchange_imbalance_m3_s)
		.Entry("darcy_cell_balance_defect_m3_s",hemodynamics.darcy.maximum_cell_balance_defect_m3_s)
		.Entry("transport_components",species.components)
		.Entry("transport_largest_component",species.largest_component)
		.Entry("gauss_seidel_sweeps",species.gauss_seidel_sweeps)
		.Entry("newton_failures",species.newton_failures)
		.Entry("flow_seconds",seconds(flow_started,flow_finished))
		.Entry("transport_seconds",seconds(flow_finished,transport_finished))
		.Entry("total_seconds",seconds(started,std::chrono::steady_clock::now())).Close();

	if(const auto* spec=Find(request,"field_output")){
		const auto& object=RequireObject(*spec,"field_output");
		RequireKnownKeys(object,{"file","names"},"field_output");
		auto path=std::filesystem::path(RequireString(Field(object,"file","field_output"),"field_output.file"));
		if(path.is_relative())path=base/path;
		std::vector<std::string> names;
		for(const auto& item:RequireArray(Field(object,"names","field_output"),"field_output.names")){
			names.push_back(RequireString(item,"field_output.names"));
			if(!out.count(names.back()))throw std::invalid_argument("unknown output field "+names.back());
		}
		const auto temporary=path.string()+".partial";
		{
			std::ofstream file(temporary,std::ios::binary|std::ios::trunc);
			if(!file)throw std::runtime_error("cannot write field output "+path.string());
			for(const auto& name:names)
				file.write(reinterpret_cast<const char*>(out.at(name).data()),
					static_cast<std::streamsize>(out.at(name).size()*sizeof(double)));
			if(!file)throw std::runtime_error("field output write failed");
		}
		std::filesystem::rename(temporary,path);
		json.Key("field_output").Open().Entry("file",path.string()).Key("names").OpenArray();
		for(const auto& name:names)json.Value(name);
		json.CloseArray().Key("lengths").OpenArray();
		for(const auto& name:names)json.Value(out.at(name).size());
		json.CloseArray().Close();
	}
	if(const auto* vtk=Find(request,"vtk_file")){
		auto path=std::filesystem::path(RequireString(*vtk,"vtk_file"));
		if(path.is_relative())path=base/path;
		std::map<std::string,std::vector<double>> cell_fields;
		for(const auto& item:out){
			if(item.second.size()==cells){cell_fields[item.first]=item.second;continue;}
			auto& mean=cell_fields[item.first];
			mean.assign(cells,0.);
			for(std::size_t i=0;i<cells;++i)
				for(std::size_t z=0;z<zones;++z)mean[i]+=item.second[i*zones+z]/double(zones);
			if(zones>1)cell_fields["last_zone_"+item.first]=last_zone(item.first);
		}
		WriteLegacyVtk(path,organ,cell_fields);
	}
	json.Close();
	return json.Text();
}

std::string Describe(const Organ& organ,const std::filesystem::path& geometry_file)
{
	{
		std::ofstream file(geometry_file,std::ios::binary|std::ios::trunc);
		if(!file)throw std::runtime_error("cannot write geometry "+geometry_file.string());
		for(int axis=0;axis<3;++axis)
			for(const auto& c:organ.centroid)file.write(reinterpret_cast<const char*>(&c[axis]),sizeof(double));
		file.write(reinterpret_cast<const char*>(organ.volume.data()),
			static_cast<std::streamsize>(organ.volume.size()*sizeof(double)));
		std::vector<double> ids;
		for(const auto& cell:organ.mesh.cells)ids.push_back(static_cast<double>(cell.id));
		file.write(reinterpret_cast<const char*>(ids.data()),static_cast<std::streamsize>(ids.size()*sizeof(double)));
		if(!file)throw std::runtime_error("geometry write failed");
	}
	Json json;
	json.Open().Entry("ok",true).Entry("reaction_model",organ.reaction_model)
		.Entry("cells",organ.mesh.cells.size()).Entry("nodes",organ.mesh.points.size())
		.Entry("lobular_zones",organ.zones)
		.Entry("volume_m3",organ.total_volume)
		.Entry("case_sha256",organ.case_hash).Entry("mesh_sha256",organ.mesh_hash)
		.Entry("territories_sha256",organ.territories_hash);
	json.Key("geometry_file").Open().Entry("file",geometry_file.string()).Key("names").OpenArray()
		.Value("centroid_x_m").Value("centroid_y_m").Value("centroid_z_m").Value("volume_m3")
		.Value("gmsh_cell_id").CloseArray().Close();
	json.Key("species").OpenArray();
	for(const auto* name:iga::PerfusedOrganSpeciesNames())json.Value(name);
	json.CloseArray();
	json.Key("trees").Open();
	for(std::size_t t=0;t<organ.trees.size();++t){
		const auto& tree=organ.trees[t];
		const auto volumes=iga::PerfusedOrganTerritoryVolumes(organ.mesh,tree);
		json.Key(tree.id).Open().Entry("role",tree.inflow?"inflow":"outflow")
			.Entry("skeleton_sha256",organ.skeleton_hashes[t])
			.Entry("nodes",tree.network.nodes.size()).Entry("segments",tree.network.segments.size())
			.Entry("terminals_with_territory",tree.territories.size());
		json.Key("territory_volume_m3").Open();
		std::size_t index=0;
		for(const auto& item:tree.territories)
			json.Entry(std::to_string(organ.terminal_swc_id[t].at(item.first)),volumes[index++]);
		json.Close().Close();
	}
	json.Close().Close();
	return json.Text();
}

void Usage()
{
	std::cerr<<"usage: native_tet_perfused_organ CASE.json --check-input\n"
		"       native_tet_perfused_organ CASE.json --describe GEOMETRY.f64\n"
		"       native_tet_perfused_organ CASE.json --request REQUEST.json [--output RESPONSE.json]\n"
		"       native_tet_perfused_organ CASE.json --serve   (one request JSON per stdin line)\n";
}

int Run(int argc,char** argv)
{
	int ranks=1;
	MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	if(ranks!=1)throw std::invalid_argument("native_tet_perfused_organ runs on one MPI rank");
	if(argc<3){Usage();return 2;}
	const std::filesystem::path case_path=argv[1];
	const std::string mode=argv[2];
	const auto organ=LoadOrgan(case_path);
	if(mode=="--check-input"){
		std::cout<<Json().Open().Entry("ok",true).Entry("cells",organ.mesh.cells.size())
			.Entry("trees",organ.trees.size()).Close().Text()<<'\n';
		return 0;
	}
	if(mode=="--describe"&&argc>=4){
		std::cout<<Describe(organ,argv[3])<<'\n';
		return 0;
	}
	if(mode=="--request"&&argc>=4){
		const std::filesystem::path request_path=argv[3];
		const auto text=ReadText(request_path,"request");
		const auto response=Solve(organ,text,std::filesystem::absolute(request_path).parent_path());
		if(argc>=6&&std::string(argv[4])=="--output"){
			const std::string temporary=std::string(argv[5])+".partial";
			{
				std::ofstream file(temporary);
				file<<response<<'\n';
				if(!file)throw std::runtime_error("cannot write response");
			}
			std::filesystem::rename(temporary,argv[5]);
		}else std::cout<<response<<'\n';
		return 0;
	}
	if(mode=="--serve"){
		std::cout<<Json().Open().Entry("ok",true).Entry("ready",true)
			.Entry("cells",organ.mesh.cells.size()).Close().Text()<<std::endl;
		for(std::string line;std::getline(std::cin,line);){
			if(line.find_first_not_of(" \t\r")==std::string::npos)continue;
			try{
				if(line.find("\"command\"")!=std::string::npos&&line.find("\"exit\"")!=std::string::npos)break;
				std::cout<<Solve(organ,line,std::filesystem::current_path())<<std::endl;
			}catch(const std::exception& error){
				std::cout<<Json().Open().Entry("ok",false).Entry("error",std::string(error.what()))
					.Close().Text()<<std::endl;
			}
		}
		return 0;
	}
	Usage();
	return 2;
}

} // namespace

int main(int argc,char** argv)
{
	if(PetscInitialize(&argc,&argv,nullptr,nullptr))return 1;
	// Our own --flags are not PETSc options; keep stdout pure JSON.
	PetscOptionsSetValue(nullptr,"-options_left","no");
	int status=0;
	try{status=Run(argc,argv);}
	catch(const std::exception& error){
		std::cerr<<"native_tet_perfused_organ: "<<error.what()<<'\n';
		status=1;
	}
	PetscFinalize();
	return status;
}
