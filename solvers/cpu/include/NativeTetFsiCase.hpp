#ifndef IGA_NATIVE_TET_FSI_CASE_HPP
#define IGA_NATIVE_TET_FSI_CASE_HPP

// Native tetrahedral FSI: case files, interface matching, and the strongly
// coupled step loop shared by native_tet_fsi_steps (one combined case file)
// and native_tet_fsi_graph (fluid and solid domains of a multidomain graph).

#include "CaseConfig.hpp"
#include "NativeTetAleConservation.hpp"
#include "NativeTetAleFsiRuntime.hpp"
#include "NativeTetAlePetscRuntime.hpp"
#include "NativeTetHydraulicVisualization.hpp"
#include "NativeTetLinearElasticVisualization.hpp"
#include "NativeTetMovingSpeciesPetscRuntime.hpp"
#include "NativeTetPetscLinearSolve.hpp"
#include "NativeTetSpeciesVisualization.hpp"
#include "NativeTetSolidFsiRuntime.hpp"
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
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

struct NativeTetFsiSpecies
{
	std::string id;
	double initial=0.,diffusivity=0.,source=0.,decay=0.;
	std::map<int,double> inflow;
	std::map<int,NativeTetWallExchange> wall_exchange;
};

struct NativeTetFsiFluidCase
{
	std::filesystem::path mesh;
	double density=0.,viscosity=0.,initial_velocity=0.;
	std::map<int,std::array<double,3>> prescribed_velocity;
	std::set<int> natural_labels;
	std::vector<NativeTetFsiSpecies> species;
};

struct NativeTetFsiSolidCase
{
	std::filesystem::path mesh;
	double young=0.,poisson=0.,density=0.;
	std::vector<std::size_t> fixed_nodes;
};

struct NativeTetFsiCase
{
	NativeTetFsiFluidCase fluid;
	NativeTetFsiSolidCase solid;
	int interface_label=0,steps=0;
	double dt_s=0.;
};

namespace native_tet_fsi_detail {

template<class Object>
const config_detail::JsonValue& Field(const Object& object,const char* name)
{
	const auto* value=config_detail::Find(object,name);
	if(!value)throw std::invalid_argument(std::string("FSI case requires ")+name);
	return *value;
}

template<class Object>
const auto& Section(const Object& object,const char* name)
{
	return config_detail::RequireObject(Field(object,name),name);
}

template<class Object>
double Number(const Object& object,const char* name)
{
	return config_detail::RequireNumber(Field(object,name),name);
}

template<class Object>
int Integer(const Object& object,const char* name)
{
	return config_detail::RequireInteger(Field(object,name),name);
}

inline config_detail::JsonValue ParseFile(const std::filesystem::path& file,const char* what)
{
	std::ifstream input(file);
	if(!input)throw std::runtime_error(std::string("cannot open ")+what);
	std::ostringstream stream;stream<<input.rdbuf();
	return config_detail::JsonParser(stream.str()).Parse();
}

template<class Object>
void ParseFluidSection(const Object& fluid,NativeTetFsiFluidCase& result)
{
	using namespace config_detail;
	result.density=Number(fluid,"density_kg_m3");
	result.viscosity=Number(fluid,"dynamic_viscosity_pa_s");
	result.initial_velocity=Number(fluid,"initial_edge_velocity_x_m_s");
	if(const auto* prescribed=Find(fluid,"boundary_velocity_by_label_m_s"))
		for(const auto& entry:RequireObject(*prescribed,"fluid boundary velocity")){
			const auto& values=RequireArray(entry.second,"boundary velocity vector");
			if(values.size()!=3)throw std::invalid_argument("boundary velocity needs three components");
			std::array<double,3> vector{};
			for(int axis=0;axis<3;++axis)
				vector[axis]=RequireNumber(values[axis],"boundary velocity component");
			result.prescribed_velocity.emplace(std::stoi(entry.first),vector);
		}
	if(const auto* natural=Find(fluid,"natural_boundary_labels"))
		for(const auto& item:RequireArray(*natural,"natural boundary labels"))
			result.natural_labels.insert(RequireInteger(item,"natural boundary label"));
}

template<class Object>
void ParseSpecies(const Object& root,NativeTetFsiFluidCase& result)
{
	using namespace config_detail;
	if(const auto* species=Find(root,"species"))
		for(const auto& entry:RequireArray(*species,"FSI species")){
			const auto& data=RequireObject(entry,"FSI species entry");
			NativeTetFsiSpecies model;
			model.id=RequireString(Field(data,"id"),"species id");
			model.initial=Number(data,"initial_concentration_mol_m3");
			model.diffusivity=Number(data,"diffusivity_m2_s");
			model.source=Number(data,"source_mol_m3_s");
			if(const auto* decay=Find(data,"first_order_decay_rate_s_inv"))
				model.decay=RequireNumber(*decay,"species decay");
			for(const auto& donor:RequireObject(Field(data,
				"inflow_concentration_by_label_mol_m3"),"species inflow"))
				model.inflow.emplace(std::stoi(donor.first),
					RequireNumber(donor.second,"species inflow concentration"));
			if(const auto* exchange=Find(data,"wall_exchange_by_label"))
				for(const auto& wall:RequireObject(*exchange,"species wall exchange")){
					const auto& condition=RequireObject(wall.second,"wall exchange");
					model.wall_exchange.emplace(std::stoi(wall.first),
						NativeTetWallExchange{
							Number(condition,"transfer_coefficient_m_s"),
							Number(condition,"external_concentration_mol_m3")});
				}
			result.species.push_back(std::move(model));
		}
}

template<class Object>
void ParseSolidSection(const Object& solid,NativeTetFsiSolidCase& result)
{
	using namespace config_detail;
	result.young=Number(solid,"young_modulus_pa");
	result.poisson=Number(solid,"poisson_ratio");
	result.density=Number(solid,"density_kg_m3");
	for(const auto& item:RequireArray(Field(solid,"fixed_nodes"),"fixed_nodes"))
		result.fixed_nodes.push_back(RequireInteger(item,"fixed node"));
}

} // namespace native_tet_fsi_detail

inline void ValidateNativeTetFsiFluidCase(const NativeTetFsiFluidCase& fluid)
{
	if(!(fluid.density>0.)||!(fluid.viscosity>0.))
		throw std::invalid_argument("native FSI fluid material is invalid");
}

inline void ValidateNativeTetFsiSolidCase(const NativeTetFsiSolidCase& solid)
{
	if(!(solid.young>0.)||!(solid.density>0.)||solid.poisson<=-1.||solid.poisson>=.5)
		throw std::invalid_argument("native FSI solid material is invalid");
}

inline void ValidateNativeTetFsiCase(const NativeTetFsiCase& scenario)
{
	if(!(scenario.fluid.density>0.)||!(scenario.fluid.viscosity>0.)||!(scenario.solid.young>0.)
		||!(scenario.solid.density>0.)||!(scenario.dt_s>0.)||scenario.steps<1
		||scenario.solid.poisson<=-1.||scenario.solid.poisson>=.5)
		throw std::invalid_argument("native FSI material or time is invalid");
}

// The combined case of native_tet_fsi_steps: fluid, solid, interface, and time.
inline NativeTetFsiCase ReadNativeTetFsiCase(const std::filesystem::path& file)
{
	using namespace config_detail;
	using namespace native_tet_fsi_detail;
	const auto parsed=ParseFile(file,"native FSI case");
	const auto& root=RequireObject(parsed,"native FSI case");
	if(Integer(root,"schema_version")!=1)
		throw std::invalid_argument("native FSI schema differs");
	NativeTetFsiCase result;
	result.fluid.mesh=file.parent_path()/RequireString(Field(root,"fluid_mesh_file"),
		"fluid_mesh_file");
	result.solid.mesh=file.parent_path()/RequireString(Field(root,"solid_mesh_file"),
		"solid_mesh_file");
	result.interface_label=Integer(root,"interface_label");
	ParseFluidSection(Section(root,"fluid"),result.fluid);
	ParseSolidSection(Section(root,"solid"),result.solid);
	const auto& time=Section(root,"time");
	result.dt_s=Number(time,"dt_s");result.steps=Integer(time,"steps");
	ParseSpecies(root,result.fluid);
	ValidateNativeTetFsiCase(result);
	return result;
}

// Graph domain cases: fluid.json holds mesh_file, fluid, and optional species;
// solid.json holds mesh_file and solid. Interface and time come from the graph.
inline NativeTetFsiFluidCase ReadNativeTetFsiFluidCase(const std::filesystem::path& file)
{
	using namespace config_detail;
	using namespace native_tet_fsi_detail;
	const auto parsed=ParseFile(file,"native FSI fluid case");
	const auto& root=RequireObject(parsed,"native FSI fluid case");
	if(Integer(root,"schema_version")!=1)
		throw std::invalid_argument("native FSI fluid schema differs");
	NativeTetFsiFluidCase result;
	result.mesh=file.parent_path()/RequireString(Field(root,"mesh_file"),"mesh_file");
	ParseFluidSection(Section(root,"fluid"),result);
	ParseSpecies(root,result);
	ValidateNativeTetFsiFluidCase(result);
	return result;
}

inline NativeTetFsiSolidCase ReadNativeTetFsiSolidCase(const std::filesystem::path& file)
{
	using namespace config_detail;
	using namespace native_tet_fsi_detail;
	const auto parsed=ParseFile(file,"native FSI solid case");
	const auto& root=RequireObject(parsed,"native FSI solid case");
	if(Integer(root,"schema_version")!=1)
		throw std::invalid_argument("native FSI solid schema differs");
	NativeTetFsiSolidCase result;
	result.mesh=file.parent_path()/RequireString(Field(root,"mesh_file"),"mesh_file");
	ParseSolidSection(Section(root,"solid"),result);
	ValidateNativeTetFsiSolidCase(result);
	return result;
}

inline NativeTetMesh ReadNativeTetFsiMesh(const std::filesystem::path& file)
{
	std::ifstream input(file);
	if(!input)throw std::runtime_error("cannot open native FSI tetra mesh");
	return ReadNativeTetMeshGmsh41(input);
}

// Matching interface: both meshes get the shared interface vertices first, in
// the coordinate order of `common`, so fluid and solid interface node IDs agree.
inline std::map<std::array<double,3>,std::uint32_t> NativeTetFsiInterfaceNodes(
	const NativeTetMesh& mesh,int label)
{
	std::map<std::array<double,3>,std::uint32_t> result;
	for(const auto& face:mesh.boundary_triangles)
		if(face.boundary_label==label)
			for(auto node:face.nodes)
				result.emplace(mesh.points[node],0);
	std::uint32_t index=0;
	for(auto& item:result)item.second=index++;
	return result;
}

inline std::vector<std::uint32_t> RenumberNativeTetFsiInterface(NativeTetMesh& mesh,int label,
	const std::map<std::array<double,3>,std::uint32_t>& common)
{
	std::vector<std::uint32_t> old_to_new(mesh.points.size(),
		std::numeric_limits<std::uint32_t>::max());
	std::set<std::uint32_t> interface_nodes;
	for(const auto& face:mesh.boundary_triangles)
		if(face.boundary_label==label)
			for(auto node:face.nodes)interface_nodes.insert(node);
	if(interface_nodes.size()!=common.size())
		throw std::invalid_argument("fluid and solid interface node counts differ");
	for(auto node:interface_nodes)
		old_to_new[node]=common.at(mesh.points[node]);
	std::uint32_t next=static_cast<std::uint32_t>(common.size());
	for(std::size_t node=0;node<mesh.points.size();++node)
		if(old_to_new[node]==std::numeric_limits<std::uint32_t>::max())
			old_to_new[node]=next++;
	std::vector<std::array<double,3>> reordered(mesh.points.size());
	for(std::size_t node=0;node<mesh.points.size();++node)
		reordered[old_to_new[node]]=mesh.points[node];
	mesh.points=std::move(reordered);
	for(auto& cell:mesh.cells)
		for(auto& node:cell.nodes)node=old_to_new[node];
	for(auto& face:mesh.boundary_triangles)
		for(auto& node:face.nodes)node=old_to_new[node];
	return old_to_new;
}

inline std::string NativeTetFsiInterfaceIdentity(const NativeTetMesh& mesh,int label)
{
	Sha256 hash;
	for(const auto& face:mesh.boundary_triangles)
		if(face.boundary_label==label){
			std::array<std::array<double,3>,3> points{{mesh.points.at(face.nodes[0]),
				mesh.points.at(face.nodes[1]),mesh.points.at(face.nodes[2])}};
			std::sort(points.begin(),points.end());
			for(const auto& point:points)
				for(double component:point)hash.AppendNormalizedDouble(component);
		}
	return hash.Hex();
}

inline DistributedSurfaceLayout BuildNativeTetFsiSurfaceLayout(const NativeTetMesh& mesh,
	int label,const std::string& identity)
{
	DistributedSurfaceLayout result;
	result.reference_mesh_identity_sha256=identity;
	std::set<std::uint64_t> ids;
	for(const auto& triangle:mesh.boundary_triangles)
		if(triangle.boundary_label==label){
			result.reference_triangles.push_back({{triangle.nodes[0],triangle.nodes[1],
				triangle.nodes[2]}});
			for(auto node:triangle.nodes)ids.insert(node);
		}
	result.global_node_count=ids.size();
	result.partition_count=1;result.partition_rank=0;
	result.owned_global_node_ids.assign(ids.begin(),ids.end());
	for(auto node:ids)
		result.reference_positions.push_back({node,mesh.points[node]});
	result.owned_reference_lumped_areas_m2.assign(ids.size(),0.);
	for(const auto& triangle:result.reference_triangles){
		const auto& a=mesh.points[triangle[0]];
		const auto& b=mesh.points[triangle[1]];
		const auto& c=mesh.points[triangle[2]];
		std::array<double,3> ab{},ac{};
		for(int axis=0;axis<3;++axis){ab[axis]=b[axis]-a[axis];ac[axis]=c[axis]-a[axis];}
		const double area=.5*NativeTetInterfaceNorm(NativeTetInterfaceCross(ab,ac));
		for(auto node:triangle){
			const auto row=std::lower_bound(result.owned_global_node_ids.begin(),
				result.owned_global_node_ids.end(),node)-
				result.owned_global_node_ids.begin();
			result.owned_reference_lumped_areas_m2[row]+=area/3.;
		}
	}
	result.layout_identity_sha256=BuildDistributedSurfaceLayoutIdentitySha256(result);
	return result;
}

inline DistributedSurfaceInterface MakeNativeTetFsiSurface(const std::string& domain,
	const std::string& subsystem,const std::string& interface_id,int label,
	const std::string& identity,bool fluid)
{
	DistributedSurfaceInterface value;
	value.id={domain,subsystem,interface_id};value.subsystem_id=subsystem;
	value.boundary_labels={label};value.reference_mesh_identity_sha256=identity;
	if(fluid){
		value.provides={SurfaceFieldQuantity::TractionOnStructure};
		value.requires={SurfaceFieldQuantity::Displacement,SurfaceFieldQuantity::Velocity};
	}else{
		value.provides={SurfaceFieldQuantity::Displacement,SurfaceFieldQuantity::Velocity};
		value.requires={SurfaceFieldQuantity::TractionOnStructure};
	}
	return value;
}

inline NativeTetMesh NativeTetFsiCurrentMesh(const NativeTetMesh& reference,
	const std::vector<std::array<double,3>>& displacement)
{
	auto current=reference;
	for(std::size_t node=0;node<current.points.size();++node)
		for(int axis=0;axis<3;++axis)
			current.points[node][axis]+=displacement[node][axis];
	return current;
}

inline std::vector<std::array<double,3>> NativeTetFsiVertexVelocity(
	const std::vector<double>& state,std::size_t points)
{
	std::vector<std::array<double,3>> result(points);
	for(std::size_t node=0;node<points;++node)
		for(int axis=0;axis<3;++axis)result[node][axis]=state[3*node+axis];
	return result;
}

// The fluid, mesh motion, and solid solve with PETSc on all ranks (each rank
// holds the whole mesh); dense_solvers keeps the single-rank dense solvers
// for comparison. Options prefixed fsi_mesh_ and fsi_solid_ tune the
// mesh-motion (CG/Jacobi) and solid (LU) solves; unprefixed KSP options tune
// the fluid, as for the other native ALE solvers.
struct NativeTetFsiSolverOptions
{
	bool dense_solvers=false;
	double fluid_tolerance=1e-6;       // fluid Newton absolute residual
	double solid_relative_tolerance=1e-5;
};

inline NativeTetFsiSolverOptions ReadNativeTetFsiSolverOptions(int ranks)
{
	NativeTetFsiSolverOptions options;
	PetscBool dense=PETSC_FALSE;
	PetscOptionsGetBool(nullptr,nullptr,"-fsi_dense_solvers",&dense,nullptr);
	options.dense_solvers=dense;
	if(options.dense_solvers&&ranks!=1)
		throw std::invalid_argument("-fsi_dense_solvers runs on one MPI rank");
	PetscReal fluid_tolerance=options.fluid_tolerance,solid_tolerance=options.solid_relative_tolerance;
	PetscOptionsGetReal(nullptr,nullptr,"-fsi_fluid_tolerance",&fluid_tolerance,nullptr);
	PetscOptionsGetReal(nullptr,nullptr,"-fsi_solid_rtol",&solid_tolerance,nullptr);
	options.fluid_tolerance=fluid_tolerance;
	options.solid_relative_tolerance=solid_tolerance;
	return options;
}

// Domain, subsystem, interface, and edge identifiers of the coupled pair.
struct NativeTetFsiNames
{
	std::string fluid_domain="fluid",fluid_subsystem="native-ale",fluid_interface="wall";
	std::string solid_domain="solid",solid_subsystem="native-solid",solid_interface="wall";
	std::string edge="wall";
};

// One strongly coupled fluid/solid pair: meshes with a matching interface,
// the surface contract, the two FSI runtimes, and the coordinator. RunStep
// advances one macro step, writes the snapshots, and prints the step ledger.
class NativeTetFsiSession
{
public:
	NativeTetFsiSession(NativeTetFsiCase scenario,const NativeTetFsiSolverOptions& options,
		NativeTetFsiNames names={})
		:scenario_(std::move(scenario)),names_(std::move(names))
	{
		ValidateNativeTetFsiCase(scenario_);
		fluid_mesh_=ReadNativeTetFsiMesh(scenario_.fluid.mesh);
		solid_mesh_=ReadNativeTetFsiMesh(scenario_.solid.mesh);
		const int label=scenario_.interface_label;
		const auto interface_coordinates=NativeTetFsiInterfaceNodes(fluid_mesh_,label);
		RenumberNativeTetFsiInterface(fluid_mesh_,label,interface_coordinates);
		const auto solid_node_map=RenumberNativeTetFsiInterface(solid_mesh_,label,
			interface_coordinates);
		identity_=NativeTetFsiInterfaceIdentity(fluid_mesh_,label);
		layout_=BuildNativeTetFsiSurfaceLayout(fluid_mesh_,label,identity_);
		fluid_surface_=MakeNativeTetFsiSurface(names_.fluid_domain,names_.fluid_subsystem,
			names_.fluid_interface,label,identity_,true);
		solid_surface_=MakeNativeTetFsiSurface(names_.solid_domain,names_.solid_subsystem,
			names_.solid_interface,label,identity_,false);
		edge_=std::make_unique<FsiCouplingEdge>(names_.edge,fluid_surface_.id,solid_surface_.id,
			FsiCouplingLaw::FluidStructureTractionKinematics);
		topology_=BuildNativeTaylorHoodTopology(fluid_mesh_);
		velocity_nodes_=fluid_mesh_.points.size()+topology_.edges.size();
		std::vector<double> initial(3*velocity_nodes_+fluid_mesh_.points.size(),0.);
		std::set<std::uint32_t> boundary_nodes;
		std::set<std::uint32_t> interface_nodes;
		for(const auto& face:fluid_mesh_.boundary_triangles)
			for(auto node:face.nodes){
				boundary_nodes.insert(node);
				if(face.boundary_label==label)
					interface_nodes.insert(node);
			}
		for(auto node:topology_.boundary_velocity_nodes.at(label))
			interface_nodes.insert(node);
		for(std::size_t edge_index=0;edge_index<topology_.edges.size();++edge_index){
			const auto edge_nodes=topology_.edges[edge_index];
			if(!boundary_nodes.count(edge_nodes[0])
				||!boundary_nodes.count(edge_nodes[1]))
				initial[3*(fluid_mesh_.points.size()+edge_index)]=scenario_.fluid.initial_velocity;
		}
		std::map<std::uint32_t,std::array<double,3>> essential_velocity;
		for(const auto& by_label:topology_.boundary_velocity_nodes)
			if(by_label.first!=label
				&&!scenario_.fluid.natural_labels.count(by_label.first)
				&&!scenario_.fluid.prescribed_velocity.count(by_label.first))
				for(auto node:by_label.second)
					essential_velocity.emplace(node,std::array<double,3>{{0.,0.,0.}});
		for(const auto& boundary:scenario_.fluid.prescribed_velocity)
			for(auto node:topology_.boundary_velocity_nodes.at(boundary.first))
				essential_velocity.emplace(node,interface_nodes.count(node)
					?std::array<double,3>{{0.,0.,0.}}:boundary.second);
		NativeTetAleFsiSolvers solvers;
		if(!options.dense_solvers){
			solvers.fluid=[](const NativeTetMesh& mesh,
				const std::vector<std::array<double,3>>& mesh_velocity,
				const std::vector<double>& committed,const std::vector<double>& trial,
				const std::map<std::uint32_t,std::array<double,3>>& prescribed,
				std::uint32_t gauge,const NativeNavierStokesParameters& parameters,
				double dt,double tolerance,std::size_t iterations){
				auto solve=SolveNativeTetAlePetscTransient(mesh,mesh_velocity,committed,trial,
					prescribed,gauge,parameters,dt,tolerance,iterations);
				NativeTetAleDenseSolveResult result;
				result.state=std::move(solve.replicated_state);
				result.newton_iterations=solve.newton_iterations;
				result.final_free_residual_l2=solve.final_residual_l2;
				return result;};
			NativeTetPetscLinearSolveDefaults motion;
			motion.ksp_type=KSPCG;motion.pc_type=PCJACOBI;motion.factor_package.clear();
			solvers.mesh_motion=MakeNativeTetPetscLinearSolve("fsi_mesh_",motion);
		}
		fluid_=std::make_unique<NativeTetAleFsiRuntime>(names_.fluid_domain,names_.fluid_subsystem,
			*edge_,fluid_surface_,solid_surface_,layout_,layout_,fluid_mesh_,label,
			NativeNavierStokesParameters{scenario_.fluid.density,scenario_.fluid.viscosity},0u,
			initial,0.,options.fluid_tolerance,30,essential_velocity,scenario_.fluid.natural_labels,
			std::vector<std::array<double,3>>{},std::vector<std::array<double,3>>{},std::move(solvers));
		std::map<std::size_t,double> constraints;
		for(auto node:scenario_.solid.fixed_nodes)
			for(int axis=0;axis<3;++axis)
				constraints.emplace(3*solid_node_map.at(node)+axis,0.);
		NativeTetSolidStaticOptions solid_options;
		solid_options.relative_tolerance=options.solid_relative_tolerance;
		solid_options.absolute_tolerance_n=1e-8;
		solid_options.maximum_iterations=40;
		if(!options.dense_solvers)solid_options.linear_solve=MakeNativeTetPetscLinearSolve("fsi_solid_");
		solid_=std::make_unique<NativeTetSolidFsiRuntime>(names_.solid_domain,names_.solid_subsystem,
			*edge_,solid_surface_,fluid_surface_,layout_,layout_,solid_mesh_,label,
			NativeTetSolidMaterial{scenario_.solid.young,scenario_.solid.poisson,scenario_.solid.density},
			constraints,solid_options);
		StrongFluidStructureCouplingOptions controls;
		controls.maximum_iterations=12;
		controls.absolute_displacement_tolerance_m=2e-5;
		controls.relative_displacement_tolerance=2e-2;
		controls.absolute_nodal_force_residual_tolerance_n=2e-2;
		controls.absolute_traction_residual_tolerance_pa=1e-1;
		controls.fluid_residual_tolerance=1e-6;
		controls.structure_residual_tolerance_n=1e-7;
		controls.interface_power_defect_tolerance_w=1e-10;
		controls.aitken_controls.initial_relaxation=.8;
		coordinator_=std::make_unique<Coordinator>(*fluid_,*solid_,*edge_,layout_,controls);
		for(const auto& species:scenario_.fluid.species)
			species_concentration_.emplace(species.id,
				std::vector<double>(fluid_mesh_.points.size(),species.initial));
	}

	const NativeTetFsiCase& Scenario()const{return scenario_;}
	const NativeTetMesh& FluidMesh()const{return fluid_mesh_;}
	const NativeTetMesh& SolidMesh()const{return solid_mesh_;}
	const std::string& InterfaceIdentity()const{return identity_;}
	const DistributedSurfaceLayout& Layout()const{return layout_;}
	const DistributedSurfaceInterface& FluidSurface()const{return fluid_surface_;}
	const DistributedSurfaceInterface& SolidSurface()const{return solid_surface_;}
	const FsiCouplingEdge& Edge()const{return *edge_;}

	void RunStep(int step,const std::filesystem::path& output,int rank,int ranks)
	{
		const auto& scenario=scenario_;
		auto& fluid=*fluid_;auto& solid=*solid_;
		const auto previous=NativeTetFsiCurrentMesh(fluid_mesh_,fluid.CommittedAleDisplacementM());
		const auto before=NativeTetFsiVertexVelocity(fluid.CommittedFluidState(),fluid_mesh_.points.size());
		const auto result=coordinator_->Execute({step,(step-1)*scenario.dt_s,scenario.dt_s});
		if(!result.converged)throw std::runtime_error("native FSI step did not converge");
		const auto current=NativeTetFsiCurrentMesh(fluid_mesh_,fluid.CommittedAleDisplacementM());
		const auto after=NativeTetFsiVertexVelocity(fluid.CommittedFluidState(),fluid_mesh_.points.size());
		std::vector<std::array<double,3>> mesh_velocity(fluid_mesh_.points.size());
		for(std::size_t node=0;node<mesh_velocity.size();++node)
			for(int axis=0;axis<3;++axis)
				mesh_velocity[node][axis]=(current.points[node][axis]
					-previous.points[node][axis])/scenario.dt_s;
		const auto conservation=EvaluateNativeTetAleConservation(
			previous,current,mesh_velocity,before,after,scenario.dt_s);
		std::map<std::uint64_t,std::array<std::uint32_t,3>> interface_faces;
		for(const auto& face:current.boundary_triangles)
			if(face.boundary_label==scenario.interface_label)
				interface_faces.emplace(face.id,face.nodes);
		std::array<double,3> interface_force{};
		double interface_power=0.;
		for(const auto& traction:EvaluateNativeTetFluidTractionOnMatchingStructure(
			current,topology_,fluid.CommittedFluidState(),scenario.fluid.viscosity,
			scenario.interface_label)){
			const auto& nodes=interface_faces.at(traction.triangle_id);
			const auto& a=current.points[nodes[0]];
			const auto& b=current.points[nodes[1]];
			const auto& c=current.points[nodes[2]];
			std::array<double,3> ab{},ac{},mean_velocity{};
			for(int axis=0;axis<3;++axis){
				ab[axis]=b[axis]-a[axis];ac[axis]=c[axis]-a[axis];
				mean_velocity[axis]=(mesh_velocity[nodes[0]][axis]
					+mesh_velocity[nodes[1]][axis]
					+mesh_velocity[nodes[2]][axis])/3.;
			}
			const double area=.5*NativeTetInterfaceNorm(NativeTetInterfaceCross(ab,ac));
			for(int axis=0;axis<3;++axis){
				const double force=area*traction.traction_on_structure_pa[axis];
				interface_force[axis]+=force;
				interface_power+=force*mean_velocity[axis];
			}
		}
		const double force_norm=NativeTetInterfaceNorm(interface_force);
		double maximum_wall_displacement=0.;
		for(double value:solid.CommittedDisplacementM())
			maximum_wall_displacement=std::max(maximum_wall_displacement,std::abs(value));
		NativeTetLinearElasticResult solid_field;
		solid_field.displacement_m.resize(solid_mesh_.points.size());
		for(std::size_t node=0;node<solid_mesh_.points.size();++node)
			for(int axis=0;axis<3;++axis)
				solid_field.displacement_m[node][axis]=solid.CommittedDisplacementM()[3*node+axis];
		const auto flow_path=output/("flow_step_"+std::to_string(step));
		const auto solid_path=output/("wall_step_"+std::to_string(step));
		WriteParallelVtkSnapshot(PETSC_COMM_WORLD,flow_path,
			BuildNativeTetHydraulicVtkPartition(fluid_mesh_,current,
				fluid.CommittedFluidState(),scenario.fluid.viscosity,rank,ranks),
			step*scenario.dt_s);
		WriteParallelVtkSnapshot(PETSC_COMM_WORLD,solid_path,
			BuildNativeTetLinearElasticVtkPartition(solid_mesh_,solid_field,rank,ranks),
			step*scenario.dt_s);
		std::vector<std::array<double,3>> fluid_velocity(velocity_nodes_);
		for(std::size_t node=0;node<velocity_nodes_;++node)
			for(int axis=0;axis<3;++axis)
				fluid_velocity[node][axis]=fluid.CommittedFluidState()[3*node+axis];
		for(const auto& species:scenario.fluid.species){
			auto& concentration=species_concentration_.at(species.id);
			const auto transport=SolveNativeTetMovingSpeciesPetscStep(
				previous,current,fluid_velocity,mesh_velocity,concentration,
				species.inflow,species.diffusivity,species.source,
				scenario.dt_s,false,species.decay,species.wall_exchange);
			concentration=transport.step.concentration_mol_m3;
			const auto species_path=output/("species_"+species.id+"_step_"+std::to_string(step));
			WriteParallelVtkSnapshot(PETSC_COMM_WORLD,species_path,
				BuildNativeTetSpeciesVtkPartition(fluid_mesh_,current,concentration,rank,ranks),
				step*scenario.dt_s);
			species_series_[species.id].push_back({step*scenario.dt_s,species_path/"snapshot.pvtu"});
			if(rank==0)std::cout<<std::setprecision(17)
				<<"native_fsi_species_step id="<<species.id
				<<" step="<<step
				<<" inventory_mol="<<transport.step.current_inventory_mol
				<<" balance_defect_mol_s="<<transport.step.balance_defect_mol_s
				<<" linear_iterations="<<transport.linear_iterations<<'\n';
		}
		flow_series_.push_back({step*scenario.dt_s,flow_path/"snapshot.pvtu"});
		solid_series_.push_back({step*scenario.dt_s,solid_path/"snapshot.pvtu"});
		if(rank==0)std::cout<<std::setprecision(17)
			<<"native_fsi_step step="<<step
			<<" iterations="<<result.iterations
			<<" gcl_residual_m3_s="<<conservation.gcl_residual_m3_s
				<<" volume_rate_m3_s="<<conservation.volume_rate_m3_s
				<<" mesh_boundary_flux_m3_s="<<conservation.mesh_boundary_flux_m3_s
				<<" relative_boundary_flux_m3_s="<<conservation.relative_boundary_flux_m3_s
				<<" moving_balance_residual_m3_s="<<conservation.moving_domain_balance_residual_m3_s
				<<" interface_force_n="<<force_norm
				<<" interface_power_w="<<interface_power
				<<" wall_max_displacement_m="<<maximum_wall_displacement
			<<" interface_power_defect_w="<<*result.history.back().interface_power_defect_w<<'\n';
	}

	void WriteSeries(const std::filesystem::path& output)const
	{
		WriteParallelVtkSeries(PETSC_COMM_WORLD,output/"flow.pvd",flow_series_);
		WriteParallelVtkSeries(PETSC_COMM_WORLD,output/"wall.pvd",solid_series_);
		for(const auto& series:species_series_)
			WriteParallelVtkSeries(PETSC_COMM_WORLD,output/("species_"+series.first+".pvd"),series.second);
	}

private:
	using Coordinator=StrongFluidStructureCoupling<NativeTetAleFsiRuntime,NativeTetSolidFsiRuntime>;

	NativeTetFsiCase scenario_;
	NativeTetFsiNames names_;
	NativeTetMesh fluid_mesh_,solid_mesh_;
	std::string identity_;
	DistributedSurfaceLayout layout_;
	DistributedSurfaceInterface fluid_surface_,solid_surface_;
	std::unique_ptr<FsiCouplingEdge> edge_;
	NativeTaylorHoodTopology topology_;
	std::size_t velocity_nodes_=0;
	std::unique_ptr<NativeTetAleFsiRuntime> fluid_;
	std::unique_ptr<NativeTetSolidFsiRuntime> solid_;
	std::unique_ptr<Coordinator> coordinator_;
	std::vector<std::pair<double,std::filesystem::path>> flow_series_,solid_series_;
	std::map<std::string,std::vector<double>> species_concentration_;
	std::map<std::string,std::vector<std::pair<double,std::filesystem::path>>> species_series_;
};

} // namespace iga

#endif
