// Perfused-organ step solver: steady tree/Darcy hemodynamics plus quasi-steady
// multispecies tissue transport with hepatic kinetics. See
// docs/PERFUSED_ORGAN.md for the case, request and field-file contracts.

#include "CaseConfig.hpp"
#include "CheckedText.hpp"
#include "PerfusedOrganHemodynamics.hpp"
#include "PerfusedOrganKinetics.hpp"
#include "PerfusedOrganPlant.hpp"
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
	iga::PerfusedOrganModel model;
	std::vector<std::map<int,int>> terminal_swc_id; // per tree: node index -> SWC id
	std::vector<std::array<double,3>> centroid;
	double total_volume=0.;
	std::string case_hash,mesh_hash,territories_hash;
	std::vector<std::string> skeleton_hashes;
	std::string reaction_model;
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
		organ.model.zones=static_cast<std::size_t>(RequireInteger(*zones,"lobular_zones"));
	if(organ.model.zones<1||organ.model.zones>16)throw std::invalid_argument("lobular_zones must be 1..16");
	const auto base=organ.case_path.parent_path();
	const auto mesh_path=base/RequireString(Field(root,"tissue_mesh_file","perfused organ case"),"tissue_mesh_file");
	const double mesh_scale=NumberOr(root,"mesh_length_scale_to_m",1.,"perfused organ case");
	if(!(mesh_scale>0.))throw std::invalid_argument("mesh_length_scale_to_m must be positive");
	const auto mesh_text=ReadText(mesh_path,"tissue mesh");
	organ.mesh_hash=Hash(mesh_text);
	{
		std::istringstream input(mesh_text);
		organ.model.mesh=iga::ReadNativeTetMeshGmsh41(input);
	}
	for(auto& point:organ.model.mesh.points)for(auto& x:point)x*=mesh_scale;
	if(organ.model.mesh.cells.empty())throw std::invalid_argument("tissue mesh has no tetrahedra");
	std::unordered_map<std::uint64_t,std::size_t> cell_index;
	organ.model.cell_volume.resize(organ.model.mesh.cells.size());
	organ.centroid.resize(organ.model.mesh.cells.size());
	for(std::size_t i=0;i<organ.model.mesh.cells.size();++i){
		if(!cell_index.emplace(organ.model.mesh.cells[i].id,i).second)
			throw std::invalid_argument("tissue mesh repeats a cell id");
		organ.model.cell_volume[i]=iga::PerfusedOrganCellVolume(organ.model.mesh,i);
		organ.total_volume+=organ.model.cell_volume[i];
		for(const auto node:organ.model.mesh.cells[i].nodes)
			for(int axis=0;axis<3;++axis)organ.centroid[i][axis]+=0.25*organ.model.mesh.points[node][axis];
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
				// Range-check before the cast: converting a negative or huge double is undefined.
				if(!(id>=1.&&id<18446744073709551616.)||std::floor(id)!=id)
					throw std::invalid_argument("territory cell id is not in the tissue mesh");
				const auto cell=cell_index.find(static_cast<std::uint64_t>(id));
				if(cell==cell_index.end())
					throw std::invalid_argument("territory cell id is not in the tissue mesh");
				list.push_back({cell->second,RequireNumber(weights[k],"weights")});
			}
			swc_ids[node->second]=swc_id;
		}
		iga::ValidatePerfusedOrganTree(organ.model.mesh,tree);
		organ.model.trees.push_back(std::move(tree));
		organ.terminal_swc_id.push_back(swc_ids);
	}
	iga::FinalizePerfusedOrganModel(organ.model);
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
	iga::config_detail::RequireKnownKeys(object,{"vmax_o2_mol_m3_s","km_o2_mmhg","o2_hill_coefficient","atp_per_o2",
		"glycolysis_capacity_mol_m3_s","glycogen_half_saturation_mol_m3","km_glucose_mol_m3",
		"lactate_vmax_mol_m3_s","lactate_km_mol_m3","gluconeogenic_fraction",
		"glycogenolysis_vmax_mol_m3_s","glucose_uptake_vmax_mol_m3_s","insulin_half_effect_mU_L",
		"insulin_clearance_s_inv","respiratory_quotient","aerobic_glucose_fraction","glycolysis_ph50",
		"glycolysis_ph_hill","q10","temperature_c"},"kinetics");
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
	k.o2_hill_coefficient=NumberOr(object,"o2_hill_coefficient",k.o2_hill_coefficient,c);
	k.atp_per_o2=NumberOr(object,"atp_per_o2",k.atp_per_o2,c);
	k.glycogen_half_saturation_mol_m3=NumberOr(object,"glycogen_half_saturation_mol_m3",
		k.glycogen_half_saturation_mol_m3,c);
	k.km_glucose_mol_m3=NumberOr(object,"km_glucose_mol_m3",k.km_glucose_mol_m3,c);
	k.lactate_km_mol_m3=NumberOr(object,"lactate_km_mol_m3",k.lactate_km_mol_m3,c);
	k.gluconeogenic_fraction=NumberOr(object,"gluconeogenic_fraction",k.gluconeogenic_fraction,c);
	k.insulin_half_effect_mU_L=NumberOr(object,"insulin_half_effect_mU_L",k.insulin_half_effect_mU_L,c);
	k.respiratory_quotient=NumberOr(object,"respiratory_quotient",k.respiratory_quotient,c);
	k.aerobic_glucose_fraction=NumberOr(object,"aerobic_glucose_fraction",k.aerobic_glucose_fraction,c);
	k.glycolysis_ph50=NumberOr(object,"glycolysis_ph50",k.glycolysis_ph50,c);
	k.glycolysis_ph_hill=NumberOr(object,"glycolysis_ph_hill",k.glycolysis_ph_hill,c);
	k.q10=NumberOr(object,"q10",k.q10,c);
	return k;
}

iga::PerfusateChemistry ParseChemistry(const Object& object,bool require_hemoglobin=true)
{
	iga::config_detail::RequireKnownKeys(object,{"hemoglobin_g_dl","oxygen_solubility_ml_dl_mmhg",
		"p50_mmhg","hill_exponent","bohr_coefficient","co2_solubility_mol_m3_mmhg",
		"bicarbonate_pk","non_bicarbonate_buffer_mol_m3_ph","buffer_anion_at_7_4_mol_m3"},"chemistry");
	iga::PerfusateChemistry chemistry;
	const std::string c="chemistry";
	auto& oxygen=chemistry.oxygen;
	oxygen.hemoglobin_g_dl=require_hemoglobin?Number(object,"hemoglobin_g_dl",c)
		:NumberOr(object,"hemoglobin_g_dl",oxygen.hemoglobin_g_dl,c);
	oxygen.oxygen_solubility_ml_dl_mmhg=NumberOr(object,"oxygen_solubility_ml_dl_mmhg",
		oxygen.oxygen_solubility_ml_dl_mmhg,c);
	oxygen.p50_mmhg=NumberOr(object,"p50_mmhg",oxygen.p50_mmhg,c);
	oxygen.hill_exponent=NumberOr(object,"hill_exponent",oxygen.hill_exponent,c);
	chemistry.bohr_coefficient=NumberOr(object,"bohr_coefficient",chemistry.bohr_coefficient,c);
	chemistry.co2_solubility_mol_m3_mmhg=NumberOr(object,"co2_solubility_mol_m3_mmhg",
		chemistry.co2_solubility_mol_m3_mmhg,c);
	chemistry.bicarbonate_pk=NumberOr(object,"bicarbonate_pk",chemistry.bicarbonate_pk,c);
	chemistry.non_bicarbonate_buffer_mol_m3_ph=NumberOr(object,"non_bicarbonate_buffer_mol_m3_ph",
		chemistry.non_bicarbonate_buffer_mol_m3_ph,c);
	chemistry.buffer_anion_at_7_4_mol_m3=NumberOr(object,"buffer_anion_at_7_4_mol_m3",
		chemistry.buffer_anion_at_7_4_mol_m3,c);
	iga::ValidateOxygenCapacity(oxygen);
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
	out<<"POINTS "<<organ.model.mesh.points.size()<<" double\n";
	for(const auto& p:organ.model.mesh.points)out<<p[0]<<' '<<p[1]<<' '<<p[2]<<'\n';
	const auto cells=organ.model.mesh.cells.size();
	out<<"CELLS "<<cells<<' '<<cells*5<<'\n';
	for(const auto& cell:organ.model.mesh.cells)
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
		"field_input","field_output","vtk_file","hypoxic_po2_mmhg","tag","resolution","command"},"request");
	const std::size_t cells=organ.model.mesh.cells.size(),zones=organ.model.zones,nodes=cells*zones;
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
	for(const auto& tree:organ.model.trees){
		flow.root_pressure_pa.push_back(Number(pressures,tree.id,"root_pressure_pa"));
		flow.specific_conductance_inv_pa_s.push_back(Number(conductances,tree.id,
			"specific_conductance_inv_pa_s"));
	}
	if(pressures.size()!=organ.model.trees.size()||conductances.size()!=organ.model.trees.size())
		throw std::invalid_argument("root pressures and conductances must name every tree once");
	const double mobility=Number(request,"mobility_m2_pa_s","request");
	const auto scale=field_or("mobility_scale",1.);
	flow.mobility_by_cell_m2_pa_s.resize(cells);
	for(std::size_t i=0;i<cells;++i)flow.mobility_by_cell_m2_pa_s[i]=mobility*scale[i];
	if(fields.values.count("exchange_scale"))flow.exchange_scale_by_cell=fields.values.at("exchange_scale");

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
	const std::string resolution=Find(request,"resolution")
		?RequireString(*Find(request,"resolution"),"resolution"):"cell";
	const auto step_started=std::chrono::steady_clock::now();
	const auto step=iga::StepPerfusedOrgan(organ.model,flow,transport,resolution);
	const auto& hemodynamics=step.hemodynamics;
	const auto& species=step.species;
	const auto& node_volume=step.node_volume;
	const std::size_t groups=step.pools;
	const auto step_finished=std::chrono::steady_clock::now();

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
			perfusion[i]=hemodynamics.cell_exchange_in_m3_s[i]/organ.model.cell_volume[i];
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
	for(std::size_t t=0;t<organ.model.trees.size();++t){
		const auto& tree=organ.model.trees[t];
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
		WriteQuantiles(json,name,VolumeQuantiles(out.at(name),organ.model.cell_volume));
	WriteQuantiles(json,"last_zone_po2_mmHg",VolumeQuantiles(last_zone("po2_mmHg"),organ.model.cell_volume));
	WriteQuantiles(json,"last_zone_energy_deficit",VolumeQuantiles(last_zone("energy_deficit"),organ.model.cell_volume));
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
		.Entry("gauss_seidel_unconverged",species.gauss_seidel_unconverged)
		.Entry("newton_failures",species.newton_failures)
		.Entry("step_seconds",seconds(step_started,step_finished))
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


// ---- Plant advance ---------------------------------------------------------

std::map<std::string,double> NumberMap(const Object& object,const std::string& key,const std::string& context)
{
	std::map<std::string,double> result;
	if(const auto* value=iga::config_detail::Find(object,key))
		for(const auto& item:iga::config_detail::RequireObject(*value,context+"."+key))
			result[item.first]=iga::config_detail::RequireNumber(item.second,context+"."+key+"."+item.first);
	return result;
}

void RequireSpeciesKeys(const std::map<std::string,double>& values,const std::string& context)
{
	const auto& names=iga::PerfusedOrganSpeciesNames();
	for(const auto& item:values)
		if(std::find_if(names.begin(),names.end(),[&](const char* n){return item.first==n;})==names.end())
			throw std::invalid_argument(context+" names unknown species "+item.first);
}

iga::OrganPlantParameters ParsePlantParameters(const Object& object,const Organ& organ)
{
	using namespace iga::config_detail;
	RequireKnownKeys(object,{"viscosity_pa_s","outflow_pressure_pa","specific_conductance_inv_pa_s",
		"mobility_m2_pa_s","distribution_volume_fraction","kinetics","chemistry","hemoglobin_per_hematocrit",
		"atmospheric_pressure_mmhg","water_vapour_mmhg","carbogen_o2_fraction","carbogen_co2_fraction",
		"air_o2_fraction","oxygenator_efficiency","injury_deficit_threshold","injury_rate_s_inv",
		"oedema_deficit_rate_s_inv","oedema_pressure_threshold_pa","oedema_pressure_rate_inv_pa_s",
		"oedema_resolution_s_inv","oedema_maximum","oedema_exponent","vasodilator_elimination_s_inv",
		"vasodilator_emax","vasodilator_ec50_mg_l","hyperoxic_constriction_per_100mmhg",
		"hyperoxic_constriction_threshold_mmhg","oxidative_injury_rate_s_inv","oxidative_injury_threshold_mmhg",
		"shear_injury_rate_inv_pa_s","shear_injury_threshold_pa","sweep_gas_temperature_k","injury_ph50",
		"injury_ph_hill","injury_glucose_threshold_mol_m3","injury_glucose_gain_per_mol_m3",
		"vasodilator_respiration_ic50_mg_l","vasodilator_respiration_imax"},"parameters");
	iga::OrganPlantParameters p;
	const std::string c="parameters";
	p.viscosity_pa_s=Number(object,"viscosity_pa_s",c);
	p.outflow_pressure_pa=NumberMap(object,"outflow_pressure_pa",c);
	p.specific_conductance_inv_pa_s=NumberMap(object,"specific_conductance_inv_pa_s",c);
	for(const auto& tree:organ.model.trees){
		if(!p.specific_conductance_inv_pa_s.count(tree.id))
			throw std::invalid_argument("parameters.specific_conductance_inv_pa_s lacks tree "+tree.id);
		if(!tree.inflow&&!p.outflow_pressure_pa.count(tree.id))
			throw std::invalid_argument("parameters.outflow_pressure_pa lacks tree "+tree.id);
	}
	p.mobility_m2_pa_s=Number(object,"mobility_m2_pa_s",c);
	p.distribution_volume_fraction=SpeciesObject(RequireObject(Field(object,"distribution_volume_fraction",c),
		"distribution_volume_fraction"),"distribution_volume_fraction");
	p.kinetics=ParseKinetics(RequireObject(Field(object,"kinetics",c),"kinetics"));
	p.chemistry=ParseChemistry(RequireObject(Field(object,"chemistry",c),"chemistry"),false);
	const auto optional=[&](const char* key,double& value){value=NumberOr(object,key,value,c);};
	optional("hemoglobin_per_hematocrit",p.hemoglobin_per_hematocrit);
	optional("atmospheric_pressure_mmhg",p.atmospheric_pressure_mmhg);
	optional("water_vapour_mmhg",p.water_vapour_mmhg);
	optional("carbogen_o2_fraction",p.carbogen_o2_fraction);
	optional("carbogen_co2_fraction",p.carbogen_co2_fraction);
	optional("air_o2_fraction",p.air_o2_fraction);
	optional("oxygenator_efficiency",p.oxygenator_efficiency);
	optional("injury_deficit_threshold",p.injury_deficit_threshold);
	optional("injury_rate_s_inv",p.injury_rate_s_inv);
	optional("oedema_deficit_rate_s_inv",p.oedema_deficit_rate_s_inv);
	optional("oedema_pressure_threshold_pa",p.oedema_pressure_threshold_pa);
	optional("oedema_pressure_rate_inv_pa_s",p.oedema_pressure_rate_inv_pa_s);
	optional("oedema_resolution_s_inv",p.oedema_resolution_s_inv);
	optional("oedema_maximum",p.oedema_maximum);
	optional("oedema_exponent",p.oedema_exponent);
	optional("vasodilator_elimination_s_inv",p.vasodilator_elimination_s_inv);
	optional("vasodilator_emax",p.vasodilator_emax);
	optional("vasodilator_ec50_mg_l",p.vasodilator_ec50_mg_l);
	optional("hyperoxic_constriction_per_100mmhg",p.hyperoxic_constriction_per_100mmhg);
	optional("hyperoxic_constriction_threshold_mmhg",p.hyperoxic_constriction_threshold_mmhg);
	optional("oxidative_injury_rate_s_inv",p.oxidative_injury_rate_s_inv);
	optional("oxidative_injury_threshold_mmhg",p.oxidative_injury_threshold_mmhg);
	optional("shear_injury_rate_inv_pa_s",p.shear_injury_rate_inv_pa_s);
	optional("shear_injury_threshold_pa",p.shear_injury_threshold_pa);
	optional("sweep_gas_temperature_k",p.sweep_gas_temperature_k);
	optional("injury_ph50",p.injury_ph50);
	optional("injury_ph_hill",p.injury_ph_hill);
	optional("injury_glucose_threshold_mol_m3",p.injury_glucose_threshold_mol_m3);
	optional("injury_glucose_gain_per_mol_m3",p.injury_glucose_gain_per_mol_m3);
	optional("vasodilator_respiration_ic50_mg_l",p.vasodilator_respiration_ic50_mg_l);
	optional("vasodilator_respiration_imax",p.vasodilator_respiration_imax);
	return p;
}

iga::OrganPlantControls ParseControls(const Object& object)
{
	using namespace iga::config_detail;
	RequireKnownKeys(object,{"inlet_pressure_pa","fio2","dialysis_in_m3_s","dialysis_out_m3_s",
		"dialysate","infusion_rates","bolus","bolus_volume_m3","vasodilator_bolus_mg","temperature_c",
		"gas_co2_fraction","sweep_gas_flow_m3_s"},"controls");
	iga::OrganPlantControls controls;
	const std::string c="controls";
	controls.inlet_pressure_pa=Number(object,"inlet_pressure_pa",c);
	controls.fio2=Number(object,"fio2",c);
	controls.dialysis_in_m3_s=NumberOr(object,"dialysis_in_m3_s",0.,c);
	controls.dialysis_out_m3_s=NumberOr(object,"dialysis_out_m3_s",0.,c);
	controls.dialysate=NumberMap(object,"dialysate",c);
	controls.infusion_rates=NumberMap(object,"infusion_rates",c);
	controls.bolus=NumberMap(object,"bolus",c);
	RequireSpeciesKeys(controls.dialysate,"controls.dialysate");
	RequireSpeciesKeys(controls.infusion_rates,"controls.infusion_rates");
	RequireSpeciesKeys(controls.bolus,"controls.bolus");
	controls.bolus_volume_m3=NumberOr(object,"bolus_volume_m3",0.,c);
	controls.vasodilator_bolus_mg=NumberOr(object,"vasodilator_bolus_mg",0.,c);
	controls.temperature_c=NumberOr(object,"temperature_c",37.,c);
	controls.gas_co2_fraction=NumberOr(object,"gas_co2_fraction",controls.gas_co2_fraction,c);
	controls.sweep_gas_flow_m3_s=NumberOr(object,"sweep_gas_flow_m3_s",controls.sweep_gas_flow_m3_s,c);
	if(!std::isnan(controls.sweep_gas_flow_m3_s)&&!(controls.sweep_gas_flow_m3_s>0.&&std::isfinite(controls.sweep_gas_flow_m3_s)))
		throw std::invalid_argument("controls.sweep_gas_flow_m3_s must be positive");
	if(controls.dialysis_in_m3_s>0.){
		const auto& names=iga::PerfusedOrganSpeciesNames();
		for(const auto* name:names)
			if(!controls.dialysate.count(name))
				throw std::invalid_argument(std::string("controls.dialysate lacks ")+name);
	}
	return controls;
}

std::filesystem::path Resolve(const std::filesystem::path& path,const std::filesystem::path& base)
{
	return path.is_relative()?base/path:path;
}

// State files: a JSON header with scalars and the reservoir, and a float64
// field file with, in order, the six tissue species (per node), viability and
// glycogen (per node) and oedema (per cell).
iga::OrganPlantState LoadPlantState(const Object& spec,const Organ& organ,const std::filesystem::path& base)
{
	using namespace iga::config_detail;
	RequireKnownKeys(spec,{"json","fields"},"state_in");
	const auto header_text=ReadText(Resolve(RequireString(Field(spec,"json","state_in"),"state_in.json"),base),"state");
	const auto parsed=JsonParser(header_text).Parse();
	const auto& header=RequireObject(parsed,"state");
	const std::size_t cells=organ.model.mesh.cells.size(),nodes=cells*organ.model.zones;
	if(static_cast<std::size_t>(RequireInteger(Field(header,"cells","state"),"cells"))!=cells
		||static_cast<std::size_t>(RequireInteger(Field(header,"lobular_zones","state"),"lobular_zones"))!=organ.model.zones)
		throw std::invalid_argument("state does not match this organ model");
	iga::OrganPlantState state;
	state.time_s=Number(header,"time_s","state");
	state.vasodilator_mg_l=Number(header,"vasodilator_mg_l","state");
	const auto& reservoir=RequireObject(Field(header,"reservoir","state"),"state.reservoir");
	state.reservoir.volume_m3=Number(reservoir,"volume_m3","state.reservoir");
	state.reservoir.hematocrit_percent=Number(reservoir,"hematocrit_percent","state.reservoir");
	state.reservoir.temperature_c=NumberOr(reservoir,"temperature_c",37.,"state.reservoir");
	state.reservoir.species=NumberMap(reservoir,"species","state.reservoir");
	std::ifstream input(Resolve(RequireString(Field(spec,"fields","state_in"),"state_in.fields"),base),std::ios::binary);
	if(!input)throw std::runtime_error("cannot open state fields");
	const auto read=[&](std::vector<double>& values,std::size_t count){
		values.resize(count);
		input.read(reinterpret_cast<char*>(values.data()),static_cast<std::streamsize>(count*sizeof(double)));
		if(!input)throw std::runtime_error("state fields are truncated");
	};
	state.tissue.assign(nodes,iga::OrganConcentration{});
	std::vector<double> values;
	for(std::size_t s=0;s<iga::OrganSpeciesCount;++s){
		read(values,nodes);
		for(std::size_t n=0;n<nodes;++n)state.tissue[n][s]=values[n];
	}
	read(state.viability,nodes);
	read(state.glycogen_mol_m3,nodes);
	read(state.oedema,cells);
	char extra=0;
	if(input.read(&extra,1))throw std::invalid_argument("state fields are longer than expected");
	return state;
}

iga::OrganPlantState InitialPlantState(const Object& spec,const Organ& organ)
{
	using namespace iga::config_detail;
	RequireKnownKeys(spec,{"reservoir_volume_m3","hematocrit_percent","reservoir","viability",
		"glycogen_mol_m3","oedema","vasodilator_mg_l"},"initial_state");
	iga::OrganPlantState state;
	state.reservoir.volume_m3=Number(spec,"reservoir_volume_m3","initial_state");
	state.reservoir.hematocrit_percent=Number(spec,"hematocrit_percent","initial_state");
	const auto species=SpeciesObject(RequireObject(Field(spec,"reservoir","initial_state"),"reservoir"),"reservoir");
	const auto& names=iga::PerfusedOrganSpeciesNames();
	for(std::size_t s=0;s<iga::OrganSpeciesCount;++s)state.reservoir.species[names[s]]=species[s];
	const std::size_t cells=organ.model.mesh.cells.size(),nodes=cells*organ.model.zones;
	state.tissue.assign(nodes,species);
	state.viability.assign(nodes,NumberOr(spec,"viability",1.,"initial_state"));
	state.glycogen_mol_m3.assign(nodes,NumberOr(spec,"glycogen_mol_m3",0.,"initial_state"));
	state.oedema.assign(cells,NumberOr(spec,"oedema",0.,"initial_state"));
	state.vasodilator_mg_l=NumberOr(spec,"vasodilator_mg_l",0.,"initial_state");
	return state;
}

void SavePlantState(const Object& spec,const iga::OrganPlantState& state,const Organ& organ,
	const std::filesystem::path& base)
{
	using namespace iga::config_detail;
	RequireKnownKeys(spec,{"json","fields"},"state_out");
	const auto json_path=Resolve(RequireString(Field(spec,"json","state_out"),"state_out.json"),base);
	const auto fields_path=Resolve(RequireString(Field(spec,"fields","state_out"),"state_out.fields"),base);
	{
		std::ofstream out(fields_path.string()+".partial",std::ios::binary|std::ios::trunc);
		if(!out)throw std::runtime_error("cannot write state fields");
		const std::size_t nodes=state.tissue.size();
		std::vector<double> values(nodes);
		for(std::size_t s=0;s<iga::OrganSpeciesCount;++s){
			for(std::size_t n=0;n<nodes;++n)values[n]=state.tissue[n][s];
			out.write(reinterpret_cast<const char*>(values.data()),static_cast<std::streamsize>(nodes*sizeof(double)));
		}
		for(const auto* field:{&state.viability,&state.glycogen_mol_m3,&state.oedema})
			out.write(reinterpret_cast<const char*>(field->data()),static_cast<std::streamsize>(field->size()*sizeof(double)));
		if(!out)throw std::runtime_error("state field write failed");
	}
	Json json;
	json.Open().Entry("schema_version",1L).Entry("cells",organ.model.mesh.cells.size())
		.Entry("lobular_zones",organ.model.zones).Entry("time_s",state.time_s)
		.Entry("vasodilator_mg_l",state.vasodilator_mg_l);
	json.Key("reservoir").Open().Entry("volume_m3",state.reservoir.volume_m3)
		.Entry("hematocrit_percent",state.reservoir.hematocrit_percent)
		.Entry("temperature_c",state.reservoir.temperature_c).Key("species").Open();
	for(const auto& item:state.reservoir.species)json.Entry(item.first,item.second);
	json.Close().Close().Close();
	{
		std::ofstream out(json_path.string()+".partial",std::ios::trunc);
		out<<json.Text()<<'\n';
		if(!out)throw std::runtime_error("cannot write state header");
	}
	std::filesystem::rename(fields_path.string()+".partial",fields_path);
	std::filesystem::rename(json_path.string()+".partial",json_path);
}

void WriteReading(Json& json,const std::string& key,const iga::PerfusateReading& reading)
{
	const auto& names=iga::PerfusedOrganSpeciesNames();
	json.Key(key).Open().Entry("po2_mmHg",reading.po2_mmhg).Entry("saturation",reading.saturation)
		.Entry("ph",reading.ph).Entry("pco2_mmHg",reading.pco2_mmhg)
		.Entry("bicarbonate_mol_m3",reading.bicarbonate_mol_m3);
	for(std::size_t s=0;s<iga::OrganSpeciesCount;++s)json.Entry(names[s],reading.concentration[s]);
	json.Close();
}

std::string Advance(const Organ& organ,const Object& request,const std::filesystem::path& base)
{
	using namespace iga::config_detail;
	const auto started=std::chrono::steady_clock::now();
	RequireKnownKeys(request,{"command","tag","duration_s","substeps","resolution","parameters","controls",
		"state_in","initial_state","state_out","heterogeneity"},"advance request");
	const auto parameters=ParsePlantParameters(RequireObject(Field(request,"parameters","advance"),"parameters"),organ);
	const auto controls=ParseControls(RequireObject(Field(request,"controls","advance"),"controls"));
	iga::OrganPlantState state;
	if(const auto* in=Find(request,"state_in"))state=LoadPlantState(RequireObject(*in,"state_in"),organ,base);
	else state=InitialPlantState(RequireObject(Field(request,"initial_state","advance"),"initial_state"),organ);
	const std::size_t cells=organ.model.mesh.cells.size();
	Fields heterogeneity;
	if(const auto* spec=Find(request,"heterogeneity")){
		heterogeneity=ReadFields(RequireObject(*spec,"heterogeneity"),cells,organ.model.zones,base);
		for(const auto& item:heterogeneity.values)
			if(item.first!="mobility_scale"&&item.first!="exchange_scale")
				throw std::invalid_argument("heterogeneity accepts mobility_scale and exchange_scale only");
	}
	const auto field=[&](const char* name){
		const auto found=heterogeneity.values.find(name);
		return found==heterogeneity.values.end()?std::vector<double>{}:found->second;
	};
	const std::string resolution=Find(request,"resolution")?RequireString(*Find(request,"resolution"),"resolution"):"cell";
	const double duration=Number(request,"duration_s","advance");
	const int substeps=RequireInteger(Field(request,"substeps","advance"),"substeps");
	if(substeps<1)throw std::invalid_argument("advance substeps must be a positive integer");
	const auto records=iga::AdvancePerfusedOrganPlant(organ.model,parameters,controls,field("mobility_scale"),
		field("exchange_scale"),state,duration,static_cast<std::size_t>(substeps),resolution);
	if(const auto* out=Find(request,"state_out"))SavePlantState(RequireObject(*out,"state_out"),state,organ,base);
	Json json;
	json.Open().Entry("ok",true).Entry("command","advance");
	if(const auto* tag=Find(request,"tag"))json.Entry("tag",RequireString(*tag,"tag"));
	json.Entry("resolution",resolution).Entry("time_s",state.time_s);
	json.Key("records").OpenArray();
	const auto& names=iga::PerfusedOrganSpeciesNames();
	for(const auto& r:records){
		json.Open().Entry("time_s",r.time_s).Entry("flow_m3_s",r.flow_m3_s)
			.Entry("inlet_pressure_pa",r.inlet_pressure_pa).Entry("reservoir_volume_m3",r.reservoir_volume_m3)
			.Entry("hematocrit_percent",r.hematocrit_percent).Entry("vasodilator_mg_l",r.vasodilator_mg_l);
		WriteReading(json,"arterial",r.arterial);
		WriteReading(json,"venous",r.venous);
		json.Key("oxygenator").Open().Entry("exhaust_po2_mmHg",r.exhaust_po2_mmhg)
			.Entry("exhaust_pco2_mmHg",r.exhaust_pco2_mmhg).Entry("co2_removal_mol_s",r.oxygenator_co2_removal_mol_s).Close();
		json.Key("tissue").Open().Entry("hypoxic_volume_fraction",r.hypoxic_volume_fraction)
			.Entry("deficit_volume_fraction",r.deficit_volume_fraction)
			.Entry("mean_viability",r.mean_viability).Entry("mean_oedema",r.mean_oedema)
			.Entry("last_zone_po2_p05_mmHg",r.last_zone_po2_p05).Entry("energy_deficit_mean",r.energy_deficit_mean)
			.Entry("energy_deficit_p95",r.energy_deficit_p95).Entry("o2_uptake_mol_s",r.o2_uptake_mol_s).Close();
		json.Key("organ_uptake").Open();
		for(std::size_t s=0;s<iga::OrganSpeciesCount;++s)json.Entry(names[s],r.organ_uptake[s]);
		json.Close();
		json.Key("diagnostics").Open().Entry("newton_failures",r.newton_failures)
			.Entry("gauss_seidel_sweeps",r.gauss_seidel_sweeps)
			.Entry("gauss_seidel_unconverged",r.gauss_seidel_unconverged)
			.Entry("flow_residual_norm",r.flow_residual_norm)
			.Entry("exchange_imbalance_m3_s",r.exchange_imbalance_m3_s).Key("balance_defect").Open();
		for(std::size_t s=0;s<iga::OrganSpeciesCount;++s)json.Entry(names[s],r.balance_defect[s]);
		json.Close().Close().Close();
	}
	json.CloseArray();
	std::size_t newton_failures=0,unconverged=0;
	for(const auto& r:records){newton_failures+=r.newton_failures;unconverged+=r.gauss_seidel_unconverged;}
	json.Entry("newton_failures",newton_failures).Entry("gauss_seidel_unconverged",unconverged);
	json.Key("reservoir").Open().Entry("volume_m3",state.reservoir.volume_m3)
		.Entry("hematocrit_percent",state.reservoir.hematocrit_percent).Key("species").Open();
	for(const auto& item:state.reservoir.species)json.Entry(item.first,item.second);
	json.Close().Close();
	json.Entry("total_seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count());
	json.Close();
	return json.Text();
}

std::string Dispatch(const Organ& organ,const std::string& text,const std::filesystem::path& base)
{
	const auto parsed=iga::config_detail::JsonParser(text).Parse();
	const auto& request=iga::config_detail::RequireObject(parsed,"request");
	if(const auto* command=iga::config_detail::Find(request,"command")){
		const auto name=iga::config_detail::RequireString(*command,"command");
		if(name=="advance")return Advance(organ,request,base);
		if(name!="solve")throw std::invalid_argument("command must be solve or advance");
	}
	return Solve(organ,text,base);
}

std::string Describe(const Organ& organ,const std::filesystem::path& geometry_file)
{
	{
		std::ofstream file(geometry_file,std::ios::binary|std::ios::trunc);
		if(!file)throw std::runtime_error("cannot write geometry "+geometry_file.string());
		for(int axis=0;axis<3;++axis)
			for(const auto& c:organ.centroid)file.write(reinterpret_cast<const char*>(&c[axis]),sizeof(double));
		file.write(reinterpret_cast<const char*>(organ.model.cell_volume.data()),
			static_cast<std::streamsize>(organ.model.cell_volume.size()*sizeof(double)));
		std::vector<double> ids;
		for(const auto& cell:organ.model.mesh.cells)ids.push_back(static_cast<double>(cell.id));
		file.write(reinterpret_cast<const char*>(ids.data()),static_cast<std::streamsize>(ids.size()*sizeof(double)));
		if(!file)throw std::runtime_error("geometry write failed");
	}
	Json json;
	json.Open().Entry("ok",true).Entry("reaction_model",organ.reaction_model)
		.Entry("cells",organ.model.mesh.cells.size()).Entry("nodes",organ.model.mesh.points.size())
		.Entry("lobular_zones",organ.model.zones)
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
	for(std::size_t t=0;t<organ.model.trees.size();++t){
		const auto& tree=organ.model.trees[t];
		const auto volumes=iga::PerfusedOrganTerritoryVolumes(organ.model.mesh,tree);
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
		std::cout<<Json().Open().Entry("ok",true).Entry("cells",organ.model.mesh.cells.size())
			.Entry("trees",organ.model.trees.size()).Close().Text()<<'\n';
		return 0;
	}
	if(mode=="--describe"&&argc>=4){
		std::cout<<Describe(organ,argv[3])<<'\n';
		return 0;
	}
	if(mode=="--request"&&argc>=4){
		const std::filesystem::path request_path=argv[3];
		const auto text=ReadText(request_path,"request");
		const auto response=Dispatch(organ,text,std::filesystem::absolute(request_path).parent_path());
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
			.Entry("cells",organ.model.mesh.cells.size()).Close().Text()<<std::endl;
		for(std::string line;std::getline(std::cin,line);){
			if(line.find_first_not_of(" \t\r")==std::string::npos)continue;
			try{
				// Only a request whose command is "exit" stops the server.
				const auto parsed=iga::config_detail::JsonParser(line).Parse();
				if(parsed.type==iga::config_detail::JsonValue::Type::Object){
					const auto* command=iga::config_detail::Find(parsed.object,"command");
					if(command&&command->type==iga::config_detail::JsonValue::Type::String
						&&command->string=="exit")break;
				}
				std::cout<<Dispatch(organ,line,std::filesystem::current_path())<<std::endl;
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
