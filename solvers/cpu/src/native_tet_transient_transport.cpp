// Transient native Tet P2/P1 flow driven by a prescribed inflow waveform, with
// one-way passive species transported by each step's P2 velocity.
//
// Flow: backward-Euler Taylor-Hood Navier-Stokes on a rigid mesh. The inlet
// receives a parabolic profile scaled to Q(t) times an optional cosine start-up
// ramp; outlets have a prescribed gauge pressure (natural traction); walls are
// no-slip. Species: P1 transport with backward Euler or BDF2 (backward Euler on
// the first step), weak or strong Dirichlet inlet concentration c_in(t),
// constant source, first-order decay and linear wall exchange. PETSc supplies
// only algebra and MPI.
#include "CaseConfig.hpp"
#include "CheckedText.hpp"
#include "NativeTetAlePetscRuntime.hpp"
#include "NativeTetBoundaryFlow.hpp"
#include "NativeTetHydraulicVisualization.hpp"
#include "NativeTetMovingSpeciesPetscRuntime.hpp"
#include "ParallelVtkOutput.hpp"
#include "Sha256.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Json=iga::config_detail::JsonValue;
using Object=std::map<std::string,Json>;
namespace fs=std::filesystem;

const Json& Member(const Object& object,const std::string& key,const std::string& context)
{
	const auto* value=iga::config_detail::Find(object,key);
	if(!value)throw std::invalid_argument(context+"."+key+" is required");
	return *value;
}

const Object& Section(const Object& object,const std::string& key)
{
	return iga::config_detail::RequireObject(Member(object,key,"case"),key);
}

double Number(const Object& object,const std::string& key,const std::string& context)
{
	return iga::config_detail::RequireNumber(Member(object,key,context),context+"."+key);
}

double OptionalNumber(const Object& object,const std::string& key,double fallback)
{
	const auto* value=iga::config_detail::Find(object,key);
	return value?iga::config_detail::RequireNumber(*value,key):fallback;
}

int Integer(const Object& object,const std::string& key,const std::string& context)
{
	return iga::config_detail::RequireInteger(Member(object,key,context),context+"."+key);
}

// Piecewise-linear function of time, constant outside its samples. A JSON
// number is a constant; an array holds [time_s, value] pairs, times increasing.
struct Table
{
	std::vector<double> time,value;
	double operator()(double t) const
	{
		if(time.size()==1||t<=time.front())return value.front();
		if(t>=time.back())return value.back();
		const auto upper=std::upper_bound(time.begin(),time.end(),t);
		const auto i=static_cast<std::size_t>(upper-time.begin());
		const double w=(t-time[i-1])/(time[i]-time[i-1]);
		return (1.-w)*value[i-1]+w*value[i];
	}
};

Table ParseTable(const Json& json,const std::string& context)
{
	Table table;
	if(json.type==Json::Type::Number){
		table.time={0.};
		table.value={json.number};
	}else for(const auto& item:iga::config_detail::RequireArray(json,context)){
		const auto& pair=iga::config_detail::RequireArray(item,context+"[]");
		if(pair.size()!=2)throw std::invalid_argument(context+" entries are [time_s, value]");
		const double t=iga::config_detail::RequireNumber(pair[0],context);
		const double v=iga::config_detail::RequireNumber(pair[1],context);
		if(!table.time.empty()&&!(t>table.time.back()))
			throw std::invalid_argument(context+" times must increase");
		table.time.push_back(t);
		table.value.push_back(v);
	}
	if(table.time.empty())throw std::invalid_argument(context+" is empty");
	for(const double v:table.value)
		if(!std::isfinite(v))throw std::invalid_argument(context+" has a nonfinite value");
	return table;
}

struct SpeciesCase
{
	std::string id;
	double initial=0.,diffusivity=0.,source=0.,decay=0.,outlet_backflow=0.;
	bool monotone=false,bdf2=false,dirichlet=false;
	Table inflow;
	std::map<int,iga::NativeTetWallExchange> wall_exchange;
};

struct Case
{
	fs::path mesh_file;
	int inlet=-1;
	std::vector<int> outlets;
	std::set<int> walls;
	double density=0.,viscosity=0.,tolerance=0.;
	int iterations=0;
	Table flow;
	double ramp=0.,outlet_pressure=0.,backflow_beta=0.;
	double dt=0.;
	int steps=0,output_every=1;
	bool float32=false;
	std::set<int> visualization_steps;
	std::vector<SpeciesCase> species;
};

Case ParseCase(const fs::path& file)
{
	using namespace iga::config_detail;
	std::ifstream stream(file);
	const auto parsed=JsonParser(iga::ReadCheckedText(stream)).Parse();
	const auto& root=RequireObject(parsed,"case");
	RequireKnownKeys(root,{"schema_version","mesh_file","boundaries","fluid","inflow","outlets",
		"time","species","output"},"case");
	if(Integer(root,"schema_version","case")!=1)
		throw std::invalid_argument("unsupported transient transport case schema");
	Case value;
	value.mesh_file=file.parent_path()/RequireString(Member(root,"mesh_file","case"),"mesh_file");
	const auto& boundaries=Section(root,"boundaries");
	RequireKnownKeys(boundaries,{"inlet_label","outlet_labels","wall_labels"},"boundaries");
	value.inlet=Integer(boundaries,"inlet_label","boundaries");
	std::set<int> seen{value.inlet};
	for(const auto& item:RequireArray(Member(boundaries,"outlet_labels","boundaries"),"outlet_labels")){
		value.outlets.push_back(RequireInteger(item,"outlet_labels[]"));
		if(!seen.insert(value.outlets.back()).second)
			throw std::invalid_argument("a boundary label is assigned twice");
	}
	for(const auto& item:RequireArray(Member(boundaries,"wall_labels","boundaries"),"wall_labels")){
		const int label=RequireInteger(item,"wall_labels[]");
		if(!seen.insert(label).second)throw std::invalid_argument("a boundary label is assigned twice");
		value.walls.insert(label);
	}
	if(value.outlets.empty())throw std::invalid_argument("at least one outlet label is required");
	const auto& fluid=Section(root,"fluid");
	RequireKnownKeys(fluid,{"density_kg_m3","dynamic_viscosity_pa_s","nonlinear_tolerance",
		"maximum_iterations"},"fluid");
	value.density=Number(fluid,"density_kg_m3","fluid");
	value.viscosity=Number(fluid,"dynamic_viscosity_pa_s","fluid");
	value.tolerance=Number(fluid,"nonlinear_tolerance","fluid");
	value.iterations=Integer(fluid,"maximum_iterations","fluid");
	const auto& inflow=Section(root,"inflow");
	RequireKnownKeys(inflow,{"profile","flow_m3_s","startup_ramp_s"},"inflow");
	if(RequireString(Member(inflow,"profile","inflow"),"inflow.profile")!="parabolic")
		throw std::invalid_argument("inflow.profile must be parabolic");
	value.flow=ParseTable(Member(inflow,"flow_m3_s","inflow"),"inflow.flow_m3_s");
	value.ramp=OptionalNumber(inflow,"startup_ramp_s",0.);
	const auto& outlets=Section(root,"outlets");
	RequireKnownKeys(outlets,{"pressure_pa","backflow_stabilization_beta"},"outlets");
	value.outlet_pressure=Number(outlets,"pressure_pa","outlets");
	value.backflow_beta=OptionalNumber(outlets,"backflow_stabilization_beta",0.);
	const auto& time=Section(root,"time");
	RequireKnownKeys(time,{"dt_s","steps"},"time");
	value.dt=Number(time,"dt_s","time");
	value.steps=Integer(time,"steps","time");
	if(const auto* output=Find(root,"output")){
		const auto& object=RequireObject(*output,"output");
		RequireKnownKeys(object,{"every_steps","precision","visualization_steps"},"output");
		if(const auto* every=Find(object,"every_steps"))
			value.output_every=RequireInteger(*every,"output.every_steps");
		if(const auto* precision=Find(object,"precision")){
			const auto name=RequireString(*precision,"output.precision");
			if(name!="float32"&&name!="float64")
				throw std::invalid_argument("output.precision must be float32 or float64");
			value.float32=name=="float32";
		}
		if(const auto* steps=Find(object,"visualization_steps"))
			for(const auto& item:RequireArray(*steps,"output.visualization_steps"))
				value.visualization_steps.insert(RequireInteger(item,"output.visualization_steps[]"));
	}
	if(const auto* species=Find(root,"species")){
		std::set<std::string> ids;
		for(const auto& item:RequireArray(*species,"species")){
			const auto& object=RequireObject(item,"species[]");
			RequireKnownKeys(object,{"id","initial_concentration_mol_m3","diffusivity_m2_s",
				"time_integration","inlet_condition","inflow_concentration_mol_m3",
				"outlet_backflow_concentration_mol_m3","source_mol_m3_s",
				"first_order_decay_rate_s_inv","monotone","wall_exchange_by_label"},"species[]");
			SpeciesCase model;
			model.id=RequireString(Member(object,"id","species[]"),"species.id");
			if(model.id.empty()||!ids.insert(model.id).second
				||model.id.find_first_not_of("abcdefghijklmnopqrstuvwxyz"
					"ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_")!=std::string::npos)
				throw std::invalid_argument("species id must be a unique ASCII identifier");
			model.initial=Number(object,"initial_concentration_mol_m3","species[]");
			model.diffusivity=Number(object,"diffusivity_m2_s","species[]");
			model.source=OptionalNumber(object,"source_mol_m3_s",0.);
			model.decay=OptionalNumber(object,"first_order_decay_rate_s_inv",0.);
			model.outlet_backflow=OptionalNumber(object,"outlet_backflow_concentration_mol_m3",0.);
			if(const auto* monotone=Find(object,"monotone"))
				model.monotone=RequireBoolean(*monotone,"species.monotone");
			const std::string scheme=Find(object,"time_integration")
				?RequireString(Member(object,"time_integration","species[]"),"time_integration")
				:"backward_euler";
			if(scheme!="backward_euler"&&scheme!="bdf2")
				throw std::invalid_argument("species.time_integration must be backward_euler or bdf2");
			model.bdf2=scheme=="bdf2";
			const std::string inlet=Find(object,"inlet_condition")
				?RequireString(Member(object,"inlet_condition","species[]"),"inlet_condition")
				:"inflow_flux";
			if(inlet!="inflow_flux"&&inlet!="dirichlet")
				throw std::invalid_argument("species.inlet_condition must be inflow_flux or dirichlet");
			model.dirichlet=inlet=="dirichlet";
			model.inflow=ParseTable(Member(object,"inflow_concentration_mol_m3","species[]"),
				"species.inflow_concentration_mol_m3");
			if(const auto* exchange=Find(object,"wall_exchange_by_label"))
				for(const auto& entry:RequireObject(*exchange,"species.wall_exchange_by_label")){
					const auto& condition=RequireObject(entry.second,"species wall exchange");
					model.wall_exchange.emplace(std::stoi(entry.first),iga::NativeTetWallExchange{
						Number(condition,"transfer_coefficient_m_s","species wall exchange"),
						Number(condition,"external_concentration_mol_m3","species wall exchange")});
				}
			for(const auto& item:model.wall_exchange)
				if(!value.walls.count(item.first))
					throw std::invalid_argument("species wall exchange must use a wall label");
			value.species.push_back(std::move(model));
		}
	}
	if(!(value.density>0.)||!(value.viscosity>0.)||!(value.tolerance>0.)||value.iterations<1
		||!(value.dt>0.)||value.steps<1||value.output_every<1||value.ramp<0.
		||value.backflow_beta<0.)
		throw std::invalid_argument("transient transport case parameters are invalid");
	return value;
}

std::array<double,3> Coordinate(const iga::NativeTetMesh& mesh,
	const iga::NativeTaylorHoodTopology& topology,std::uint32_t node)
{
	if(node<mesh.points.size())return mesh.points[node];
	const auto& edge=topology.edges.at(node-mesh.points.size());
	std::array<double,3> point{};
	for(int c=0;c<3;++c)point[c]=0.5*(mesh.points[edge[0]][c]+mesh.points[edge[1]][c]);
	return point;
}

// Unit-inflow parabolic profile on the planar inlet cap, from the cap's own
// area-weighted centroid and normal; the rim radius is the largest in-plane
// vertex distance. Walls override shared rim nodes. Normalized to 1 m^3/s.
std::vector<double> InletProfile(const iga::NativeTetMesh& mesh,
	const iga::NativeTaylorHoodTopology& topology,int inlet,const std::set<int>& walls,
	std::array<double,3>& center,std::array<double,3>& normal,double& radius)
{
	double area=0.;
	center={};normal={};
	std::set<std::uint32_t> vertices;
	for(const auto& triangle:mesh.boundary_triangles){
		if(triangle.boundary_label!=inlet)continue;
		const auto& a=mesh.points[triangle.nodes[0]];
		const auto& b=mesh.points[triangle.nodes[1]];
		const auto& c=mesh.points[triangle.nodes[2]];
		std::array<double,3> u{},v{},n{};
		for(int k=0;k<3;++k){u[k]=b[k]-a[k];v[k]=c[k]-a[k];}
		n={{0.5*(u[1]*v[2]-u[2]*v[1]),0.5*(u[2]*v[0]-u[0]*v[2]),0.5*(u[0]*v[1]-u[1]*v[0])}};
		const double face=std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
		// Orient consistently with the first face before summing.
		if(area>0.&&n[0]*normal[0]+n[1]*normal[1]+n[2]*normal[2]<0.)for(auto& x:n)x=-x;
		for(int k=0;k<3;++k){
			center[k]+=face*(a[k]+b[k]+c[k])/3.;
			normal[k]+=n[k];
		}
		area+=face;
		vertices.insert(triangle.nodes.begin(),triangle.nodes.end());
	}
	if(!(area>0.))throw std::invalid_argument("inlet label has no boundary faces");
	const double length=std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
	for(int k=0;k<3;++k){center[k]/=area;normal[k]/=length;}
	radius=0.;
	double axial_max=0.;
	for(const auto vertex:vertices){
		double axial=0.,square=0.;
		for(int k=0;k<3;++k){
			const double d=mesh.points[vertex][k]-center[k];
			axial+=d*normal[k];square+=d*d;
		}
		axial_max=std::max(axial_max,std::abs(axial));
		radius=std::max(radius,std::sqrt(std::max(0.,square-axial*axial)));
	}
	if(!(radius>0.)||axial_max>1e-6*radius)
		throw std::invalid_argument("inlet cap is not planar");
	const auto nodes=mesh.points.size()+topology.edges.size();
	std::vector<double> profile(3*nodes+mesh.points.size(),0.);
	for(const auto node:topology.boundary_velocity_nodes.at(inlet)){
		const auto p=Coordinate(mesh,topology,node);
		double axial=0.,square=0.;
		for(int k=0;k<3;++k){
			const double d=p[k]-center[k];
			axial+=d*normal[k];square+=d*d;
		}
		const double shape=std::max(0.,2.*(1.-(square-axial*axial)/(radius*radius)));
		for(int k=0;k<3;++k)profile[3*node+k]=normal[k]*shape;
	}
	for(const int wall:walls)
		if(topology.boundary_velocity_nodes.count(wall))
			for(const auto node:topology.boundary_velocity_nodes.at(wall))
				for(int k=0;k<3;++k)profile[3*node+k]=0.;
	const auto operators=iga::BuildNativeTetAleBoundaryFluxOperators(mesh,topology);
	double inward=-operators.at(inlet).OutwardFlowM3S(profile);
	if(inward<0.){
		for(auto& x:profile)x=-x;
		for(auto& x:normal)x=-x;
		inward=-inward;
	}
	if(!(inward>0.))throw std::runtime_error("inlet profile carries no flow");
	for(auto& x:profile)x/=inward;
	for(auto& x:normal)x=-x;  // report the outward normal
	return profile;
}

template<class T>
void WriteArray(const fs::path& path,const std::vector<double>& values,std::size_t count)
{
	std::vector<T> data(values.begin(),values.begin()+static_cast<std::ptrdiff_t>(count));
	const auto temporary=path.string()+".tmp";
	{
		std::ofstream stream(temporary,std::ios::binary);
		stream.exceptions(std::ios::badbit|std::ios::failbit);
		stream.write(reinterpret_cast<const char*>(data.data()),
			static_cast<std::streamsize>(data.size()*sizeof(T)));
	}
	fs::rename(temporary,path);
}

std::vector<double> ReadDoubles(const fs::path& path,std::size_t count)
{
	std::ifstream stream(path,std::ios::binary);
	if(!stream)throw std::runtime_error("cannot read "+path.string());
	std::vector<double> values(count);
	stream.read(reinterpret_cast<char*>(values.data()),static_cast<std::streamsize>(count*sizeof(double)));
	if(stream.gcount()!=static_cast<std::streamsize>(count*sizeof(double))
		||stream.peek()!=std::char_traits<char>::eof())
		throw std::runtime_error("checkpoint size mismatch: "+path.string());
	return values;
}

// Restart identity: SHA-256 of the case text and the mesh file contents.
std::string Digest(const std::string& case_text,const std::string& mesh_text)
{
	iga::Sha256 hash;
	for(const auto* text:{&case_text,&mesh_text}){
		hash.AppendLittleEndian64(text->size());
		hash.Append(text->data(),text->size());
	}
	return hash.Hex();
}

std::string StepName(const std::string& prefix,int step,const std::string& extension)
{
	std::ostringstream name;
	name<<prefix<<std::setw(6)<<std::setfill('0')<<step<<extension;
	return name.str();
}

}  // namespace

int main(int argc,char** argv)
{
	// This driver's own arguments are removed before PETSc sees the rest.
	fs::path case_file,output;
	bool check_input=false,restart=false;
	int stop_after=-1;
	std::vector<char*> petsc_arguments{argv[0]};
	for(int i=1;i<argc;++i){
		const std::string arg=argv[i];
		if(arg=="--output-dir"&&i+1<argc)output=argv[++i];
		else if(arg=="--check-input")check_input=true;
		else if(arg=="--restart")restart=true;
		else if(arg=="--stop-after-step"&&i+1<argc)stop_after=std::stoi(argv[++i]);
		else if(case_file.empty()&&arg.rfind("-",0)!=0)case_file=arg;
		else petsc_arguments.push_back(argv[i]);
	}
	int petsc_count=static_cast<int>(petsc_arguments.size());
	char** petsc_values=petsc_arguments.data();
	if(PetscInitialize(&petsc_count,&petsc_values,nullptr,nullptr))return 1;
	int rank=0,ranks=1;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	try{
		if(case_file.empty()||(!check_input&&output.empty()))
			throw std::invalid_argument("usage: native_tet_transient_transport case.json "
				"(--check-input | --output-dir DIR [--restart] [--stop-after-step N]) [PETSc options]");
		std::ifstream case_stream(case_file);
		const std::string case_text=iga::ReadCheckedText(case_stream);
		const auto scenario=ParseCase(case_file);
		std::ifstream mesh_stream(scenario.mesh_file,std::ios::binary);
		if(!mesh_stream)throw std::invalid_argument("cannot open mesh "+scenario.mesh_file.string());
		const std::string mesh_text=iga::ReadCheckedText(mesh_stream);
		std::istringstream mesh_input(mesh_text);
		const auto mesh=iga::ReadNativeTetMeshGmsh41(mesh_input);
		const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
		for(const auto& cell:mesh.cells)iga::EvaluateNativeTetGeometry(mesh,cell);
		std::set<int> labels,assigned(scenario.walls);
		for(const auto& triangle:mesh.boundary_triangles)labels.insert(triangle.boundary_label);
		assigned.insert(scenario.inlet);
		assigned.insert(scenario.outlets.begin(),scenario.outlets.end());
		if(labels!=assigned)
			throw std::invalid_argument("mesh boundary labels differ from the case assignment");
		std::array<double,3> center{},normal{};
		double radius=0.;
		const auto profile=InletProfile(mesh,topology,scenario.inlet,scenario.walls,center,normal,radius);
		const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
		const std::size_t dofs=3*velocity_nodes+mesh.points.size();
		if(rank==0)
			std::cout<<"native transient transport input: vertices="<<mesh.points.size()
				<<" cells="<<mesh.cells.size()<<" p2_nodes="<<velocity_nodes<<" dofs="<<dofs
				<<" inlet_radius_m="<<radius<<" species="<<scenario.species.size()<<std::endl;
		if(check_input){PetscFinalize();return 0;}

		iga::NativeTetAleBoundaryConditions pressure;
		for(const int label:scenario.outlets){
			pressure.prescribed_pressure_pa[label]=scenario.outlet_pressure;
			if(scenario.backflow_beta>0.)pressure.backflow_stabilization_beta[label]=scenario.backflow_beta;
		}
		const std::vector<std::array<double,3>> grid(mesh.points.size(),{{0.,0.,0.}});
		std::vector<double> committed(dofs,0.);
		std::vector<std::vector<double>> concentration,older;
		for(const auto& model:scenario.species){
			concentration.emplace_back(mesh.points.size(),model.initial);
			older.emplace_back();
		}
		const std::string identity=Digest(case_text,mesh_text)+"-"+std::to_string(ranks);
		int first_step=1;
		if(restart){
			std::ifstream marker(output/"checkpoint.json");
			const auto parsed=iga::config_detail::JsonParser(iga::ReadCheckedText(marker)).Parse();
			const auto& object=iga::config_detail::RequireObject(parsed,"checkpoint");
			if(iga::config_detail::RequireString(Member(object,"identity","checkpoint"),"identity")!=identity)
				throw std::invalid_argument("restart case or rank count differs from the checkpoint");
			const int step=Integer(object,"step","checkpoint");
			std::size_t count=dofs;
			for(const auto& model:scenario.species)count+=mesh.points.size()*(model.bdf2?2:1);
			const auto values=ReadDoubles(output/"checkpoint.f64",count);
			std::size_t offset=0;
			std::copy_n(values.begin(),dofs,committed.begin());offset+=dofs;
			for(std::size_t s=0;s<scenario.species.size();++s){
				std::copy_n(values.begin()+static_cast<std::ptrdiff_t>(offset),mesh.points.size(),
					concentration[s].begin());
				offset+=mesh.points.size();
				if(scenario.species[s].bdf2){
					older[s].assign(values.begin()+static_cast<std::ptrdiff_t>(offset),
						values.begin()+static_cast<std::ptrdiff_t>(offset+mesh.points.size()));
					offset+=mesh.points.size();
				}
			}
			first_step=step+1;
			// A step can be logged before its checkpoint is written; drop such
			// lines so a resumed log has each accepted step exactly once.
			if(rank==0){
				std::ifstream log(output/"steps.jsonl");
				std::ostringstream kept;
				for(std::string line;std::getline(log,line);){
					// A line cut short by an interruption belongs to an unaccepted step.
					int logged=std::numeric_limits<int>::max();
					try{
						const auto row=iga::config_detail::JsonParser(line).Parse();
						logged=Integer(iga::config_detail::RequireObject(row,"steps.jsonl"),"step","steps.jsonl");
					}catch(const std::exception&){}
					if(logged<=step)kept<<line<<'\n';
				}
				log.close();
				{
					std::ofstream rewritten(output/"steps.jsonl.tmp",std::ios::trunc);
					rewritten<<kept.str();
					if(!rewritten.flush())throw std::runtime_error("cannot rewrite steps.jsonl");
				}
				fs::rename(output/"steps.jsonl.tmp",output/"steps.jsonl");
			}
		}else if(rank==0){
			if(fs::exists(output/"steps.jsonl"))
				throw std::runtime_error("output already contains a run; use --restart");
			fs::create_directories(output);
		}
		MPI_Barrier(PETSC_COMM_WORLD);
		const double pi=std::acos(-1.);
		const auto start=std::chrono::steady_clock::now();
		// Keep snapshots written before a restart in flow.pvd.
		std::vector<std::pair<double,fs::path>> visualization;
		for(const int step:scenario.visualization_steps)
			if(step<first_step&&fs::exists(output/StepName("vtk_step_",step,"")/"snapshot.pvtu"))
				visualization.push_back({step*scenario.dt,output/StepName("vtk_step_",step,"")/"snapshot.pvtu"});
		for(int step=first_step;step<=scenario.steps;++step){
			const double t=step*scenario.dt;
			const double ramp=scenario.ramp>0.&&t<scenario.ramp?0.5*(1.-std::cos(pi*t/scenario.ramp)):1.;
			const double flow=scenario.flow(t)*ramp;
			std::map<std::uint32_t,std::array<double,3>> boundary;
			for(const auto node:topology.boundary_velocity_nodes.at(scenario.inlet))
				for(int k=0;k<3;++k)boundary[node][k]=profile[3*node+k]*flow;
			for(const int wall:scenario.walls)
				if(topology.boundary_velocity_nodes.count(wall))
					for(const auto node:topology.boundary_velocity_nodes.at(wall))boundary[node]={{0.,0.,0.}};
			const auto begin=std::chrono::steady_clock::now();
			const auto solved=iga::SolveNativeTetAlePetscTransient(mesh,grid,committed,committed,boundary,
				std::numeric_limits<std::uint32_t>::max(),{scenario.density,scenario.viscosity},scenario.dt,
				scenario.tolerance,static_cast<std::size_t>(scenario.iterations),pressure);
			committed=solved.replicated_state;
			for(const double value:committed)
				if(!std::isfinite(value))throw std::runtime_error("nonfinite flow state");
			const auto flows=iga::EvaluateNativeTetBoundaryFlows(mesh,topology,committed);
			double net=0.;
			for(const auto& item:flows)net+=item.second.outward_flow_m3_s;
			const double scale=std::max(std::abs(flow),1e-300);
			const double imbalance=std::abs(net)/scale;
			const double inlet_error=std::abs(flows.at(scenario.inlet).outward_flow_m3_s+flow)/scale;
			if(std::abs(flow)>0.&&(imbalance>1e-5||inlet_error>1e-8))
				throw std::runtime_error("flow conservation gate failed at step "+std::to_string(step));
			std::vector<std::array<double,3>> velocity(velocity_nodes);
			for(std::size_t n=0;n<velocity_nodes;++n)
				for(int k=0;k<3;++k)velocity[n][k]=committed[3*n+k];
			std::ostringstream species_log;
			for(std::size_t s=0;s<scenario.species.size();++s){
				const auto& model=scenario.species[s];
				std::map<int,double> inflow{{scenario.inlet,model.inflow(t)}};
				for(const int label:scenario.outlets)inflow[label]=model.outlet_backflow;
				iga::NativeTetMovingSpeciesStepOptions options;
				options.allow_negative_concentration=!model.monotone;
				if(model.dirichlet)options.dirichlet_labels={scenario.inlet};
				if(model.bdf2&&!older[s].empty()){
					options.time_integration=iga::NativeTetSpeciesTimeIntegration::Bdf2;
					options.older_concentration_mol_m3=older[s];
				}
				const auto result=iga::SolveNativeTetMovingSpeciesPetscStep(mesh,mesh,velocity,grid,
					concentration[s],inflow,model.diffusivity,model.source,scenario.dt,model.monotone,
					model.decay,model.wall_exchange,{},{},{},options);
				const auto& values=result.step.concentration_mol_m3;
				const auto [low,high]=std::minmax_element(values.begin(),values.end());
				const double budget=std::max({1e-30,std::abs(result.step.inventory_rate_mol_s),
					std::abs(result.step.outward_advective_flux_mol_s),std::abs(result.step.source_mol_s)});
				if(std::abs(result.step.balance_defect_mol_s)>1e-6*budget
					&&std::abs(result.step.balance_defect_mol_s)>1e-15)
					throw std::runtime_error("species budget gate failed for "+model.id);
				species_log<<(s?",":"")<<'"'<<model.id<<"\":{\"inventory_mol\":"<<result.step.current_inventory_mol
					<<",\"outward_transfer_mol_s\":"<<result.step.outward_advective_flux_mol_s
					<<",\"wall_exchange_mol_s\":"<<result.step.outward_wall_exchange_mol_s
					<<",\"source_mol_s\":"<<result.step.source_mol_s
					<<",\"balance_defect_mol_s\":"<<result.step.balance_defect_mol_s
					<<",\"minimum\":"<<*low<<",\"maximum\":"<<*high
					<<",\"inlet_value\":"<<model.inflow(t)
					<<",\"scheme\":\""<<(options.time_integration==iga::NativeTetSpeciesTimeIntegration::Bdf2
						?"bdf2":"backward_euler")<<"\",\"linear_iterations\":"<<result.linear_iterations<<'}';
				older[s]=concentration[s];
				concentration[s]=values;
			}
			if(rank==0){
				std::ostringstream line;
				line<<std::setprecision(17)<<"{\"step\":"<<step<<",\"time_s\":"<<t<<",\"inflow_m3_s\":"<<flow
					<<",\"relative_mass_imbalance\":"<<imbalance<<",\"relative_inlet_error\":"<<inlet_error
					<<",\"newton_iterations\":"<<solved.newton_iterations<<",\"residual_l2\":"<<solved.final_residual_l2
					<<",\"linear_iterations\":"<<solved.total_linear_iterations
					<<",\"step_wall_s\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()
					<<",\"elapsed_s\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
					<<",\"boundaries\":{";
				bool first=true;
				for(const auto& item:flows){
					line<<(first?"":",")<<'"'<<item.first<<"\":{\"outward_flow_m3_s\":"<<item.second.outward_flow_m3_s
						<<",\"average_pressure_pa\":"<<item.second.average_pressure_pa<<'}';
					first=false;
				}
				line<<"},\"species\":{"<<species_log.str()<<"}}";
				{std::ofstream log(output/"steps.jsonl",std::ios::app);log<<line.str()<<'\n';}
				if(step%scenario.output_every==0||step==scenario.steps){
					const auto write=[&](const fs::path& path,const std::vector<double>& values){
						if(scenario.float32)WriteArray<float>(path,values,values.size());
						else WriteArray<double>(path,values,values.size());};
					const std::string extension=scenario.float32?".f32":".f64";
					write(output/StepName("flow_",step,extension),committed);
					for(std::size_t s=0;s<scenario.species.size();++s)
						write(output/StepName("species_"+scenario.species[s].id+"_",step,extension),concentration[s]);
				}
				std::vector<double> checkpoint(committed);
				for(std::size_t s=0;s<scenario.species.size();++s){
					checkpoint.insert(checkpoint.end(),concentration[s].begin(),concentration[s].end());
					if(scenario.species[s].bdf2)checkpoint.insert(checkpoint.end(),older[s].begin(),older[s].end());
				}
				WriteArray<double>(output/"checkpoint.f64",checkpoint,checkpoint.size());
				std::ofstream marker(output/"checkpoint.json.tmp");
				marker<<"{\"identity\":\""<<identity<<"\",\"step\":"<<step<<",\"time_s\":"<<std::setprecision(17)<<t<<"}\n";
				marker.close();
				fs::rename(output/"checkpoint.json.tmp",output/"checkpoint.json");
			}
			if(scenario.visualization_steps.count(step)){
				// The solved P1 pressure is already the gauge pressure (outlets carry
				// outlet_pressure_pa as traction), the same field as flow_NNNNNN.
				auto piece=iga::BuildNativeTetHydraulicVtkPartition(mesh,mesh,committed,scenario.viscosity,rank,ranks);
				const auto snapshot=output/StepName("vtk_step_",step,"");
				iga::WriteParallelVtkSnapshot(PETSC_COMM_WORLD,snapshot,piece,t);
				visualization.push_back({t,snapshot/"snapshot.pvtu"});
				iga::WriteParallelVtkSeries(PETSC_COMM_WORLD,output/"flow.pvd",visualization);
			}
			MPI_Barrier(PETSC_COMM_WORLD);
			if(step==stop_after&&step<scenario.steps){
				if(rank==0)std::cout<<"stopped after accepted step "<<step<<std::endl;
				PetscFinalize();
				return 0;
			}
		}
		if(rank==0){
			std::ofstream series(output/"series.json");
			series<<std::setprecision(17)<<"{\n  \"schema\": \"native-tet-transient-transport-output-v1\",\n"
				<<"  \"mesh_file\": \""<<scenario.mesh_file.string()<<"\",\n"
				<<"  \"vertices\": "<<mesh.points.size()<<",\n  \"p2_nodes\": "<<velocity_nodes<<",\n"
				<<"  \"dt_s\": "<<scenario.dt<<",\n  \"steps\": "<<scenario.steps<<",\n"
				<<"  \"output_every_steps\": "<<scenario.output_every<<",\n"
				<<"  \"precision\": \""<<(scenario.float32?"float32":"float64")<<"\",\n"
				<<"  \"flow_layout\": \"interleaved P2 velocity xyz (m/s) in native vertex then edge order, "
				<<"followed by P1 gauge pressure (Pa); outlets carry outlet_pressure_pa\",\n"
				<<"  \"species_layout\": \"P1 concentration (mol/m^3) per vertex\",\n"
				<<"  \"outlet_pressure_pa\": "<<scenario.outlet_pressure<<",\n"
				<<"  \"inlet\": {\"center_m\": ["<<center[0]<<", "<<center[1]<<", "<<center[2]<<"], "
				<<"\"outward_normal\": ["<<normal[0]<<", "<<normal[1]<<", "<<normal[2]<<"], \"radius_m\": "<<radius<<"},\n"
				<<"  \"species\": [";
			for(std::size_t s=0;s<scenario.species.size();++s)
				series<<(s?", ":"")<<'"'<<scenario.species[s].id<<'"';
			series<<"],\n  \"ranks\": "<<ranks<<"\n}\n";
			std::ofstream done(output/"completed.json");
			done<<"{\"completed\": true, \"steps\": "<<scenario.steps<<", \"identity\": \""<<identity<<"\"}\n";
		}
	}catch(const std::exception& error){
		std::cerr<<"rank="<<rank<<" native transient transport failed: "<<error.what()<<std::endl;
		MPI_Abort(PETSC_COMM_WORLD,2);
		return 2;
	}
	PetscFinalize();
	return 0;
}
