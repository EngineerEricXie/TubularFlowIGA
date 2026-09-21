#include "CaseConfig.hpp"
#include "NativeTetAleConservation.hpp"
#include "NativeTetAleFsiRuntime.hpp"
#include "NativeTetHydraulicVisualization.hpp"
#include "NativeTetLinearElasticVisualization.hpp"
#include "NativeTetMovingSpeciesPetscRuntime.hpp"
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
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using Object=std::map<std::string,iga::config_detail::JsonValue>;

const iga::config_detail::JsonValue& Field(const Object& object,const char* name)
{
	const auto* value=iga::config_detail::Find(object,name);
	if(!value)throw std::invalid_argument(std::string("FSI case requires ")+name);
	return *value;
}

const Object& Section(const Object& object,const char* name)
{
	return iga::config_detail::RequireObject(Field(object,name),name);
}

double Number(const Object& object,const char* name)
{
	return iga::config_detail::RequireNumber(Field(object,name),name);
}

int Integer(const Object& object,const char* name)
{
	return iga::config_detail::RequireInteger(Field(object,name),name);
}

struct Case
{
	struct Species
	{
		std::string id;
		double initial=0.,diffusivity=0.,source=0.,decay=0.;
		std::map<int,double> inflow;
		std::map<int,iga::NativeTetWallExchange> wall_exchange;
	};
	std::filesystem::path fluid_mesh,solid_mesh;
	int interface_label=0,steps=0;
	double dt_s=0.,density=0.,viscosity=0.,initial_velocity=0.;
	double young=0.,poisson=0.,solid_density=0.;
	std::vector<std::size_t> fixed_nodes;
	std::map<int,std::array<double,3>> prescribed_velocity;
	std::set<int> natural_labels;
	std::vector<Species> species;
};

Case Parse(const std::filesystem::path& file)
{
	using namespace iga::config_detail;
	std::ifstream input(file);
	if(!input)throw std::runtime_error("cannot open native FSI case");
	std::ostringstream stream;stream<<input.rdbuf();
	const auto parsed=JsonParser(stream.str()).Parse();
	const auto& root=RequireObject(parsed,"native FSI case");
	if(Integer(root,"schema_version")!=1)
		throw std::invalid_argument("native FSI schema differs");
	Case result;
	result.fluid_mesh=file.parent_path()/RequireString(Field(root,"fluid_mesh_file"),
		"fluid_mesh_file");
	result.solid_mesh=file.parent_path()/RequireString(Field(root,"solid_mesh_file"),
		"solid_mesh_file");
	result.interface_label=Integer(root,"interface_label");
	const auto& fluid=Section(root,"fluid");
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
	const auto& solid=Section(root,"solid");
	result.young=Number(solid,"young_modulus_pa");
	result.poisson=Number(solid,"poisson_ratio");
	result.solid_density=Number(solid,"density_kg_m3");
	for(const auto& item:RequireArray(Field(solid,"fixed_nodes"),"fixed_nodes"))
		result.fixed_nodes.push_back(RequireInteger(item,"fixed node"));
	const auto& time=Section(root,"time");
	result.dt_s=Number(time,"dt_s");result.steps=Integer(time,"steps");
	if(const auto* species=Find(root,"species"))
		for(const auto& entry:RequireArray(*species,"FSI species")){
			const auto& data=RequireObject(entry,"FSI species entry");
			Case::Species model;
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
						iga::NativeTetWallExchange{
							Number(condition,"transfer_coefficient_m_s"),
							Number(condition,"external_concentration_mol_m3")});
				}
			result.species.push_back(std::move(model));
		}
	if(!(result.density>0.)||!(result.viscosity>0.)||!(result.young>0.)
		||!(result.solid_density>0.)||!(result.dt_s>0.)||result.steps<1
		||result.poisson<=-1.||result.poisson>=.5)
		throw std::invalid_argument("native FSI material or time is invalid");
	return result;
}

iga::NativeTetMesh ReadMesh(const std::filesystem::path& file)
{
	std::ifstream input(file);
	if(!input)throw std::runtime_error("cannot open native FSI tetra mesh");
	return iga::ReadNativeTetMeshGmsh41(input);
}

std::vector<std::uint32_t> RenumberInterface(iga::NativeTetMesh& mesh,int label,
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

std::map<std::array<double,3>,std::uint32_t> InterfaceNodes(
	const iga::NativeTetMesh& mesh,int label)
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

std::string InterfaceIdentity(const iga::NativeTetMesh& mesh,int label)
{
	iga::Sha256 hash;
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

iga::DistributedSurfaceLayout Layout(const iga::NativeTetMesh& mesh,int label,
	const std::string& identity)
{
	iga::DistributedSurfaceLayout result;
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
		const double area=.5*iga::NativeTetInterfaceNorm(
			iga::NativeTetInterfaceCross(ab,ac));
		for(auto node:triangle){
			const auto row=std::lower_bound(result.owned_global_node_ids.begin(),
				result.owned_global_node_ids.end(),node)-
				result.owned_global_node_ids.begin();
			result.owned_reference_lumped_areas_m2[row]+=area/3.;
		}
	}
	result.layout_identity_sha256=iga::BuildDistributedSurfaceLayoutIdentitySha256(result);
	return result;
}

iga::DistributedSurfaceInterface Surface(const std::string& domain,
	const std::string& subsystem,int label,const std::string& identity,bool fluid)
{
	iga::DistributedSurfaceInterface value;
	value.id={domain,subsystem,"wall"};value.subsystem_id=subsystem;
	value.boundary_labels={label};value.reference_mesh_identity_sha256=identity;
	if(fluid){
		value.provides={iga::SurfaceFieldQuantity::TractionOnStructure};
		value.requires={iga::SurfaceFieldQuantity::Displacement,
			iga::SurfaceFieldQuantity::Velocity};
	}else{
		value.provides={iga::SurfaceFieldQuantity::Displacement,
			iga::SurfaceFieldQuantity::Velocity};
		value.requires={iga::SurfaceFieldQuantity::TractionOnStructure};
	}
	return value;
}

iga::NativeTetMesh CurrentMesh(const iga::NativeTetMesh& reference,
	const std::vector<std::array<double,3>>& displacement)
{
	auto current=reference;
	for(std::size_t node=0;node<current.points.size();++node)
		for(int axis=0;axis<3;++axis)
			current.points[node][axis]+=displacement[node][axis];
	return current;
}

std::vector<std::array<double,3>> VertexVelocity(const std::vector<double>& state,
	std::size_t points)
{
	std::vector<std::array<double,3>> result(points);
	for(std::size_t node=0;node<points;++node)
		for(int axis=0;axis<3;++axis)result[node][axis]=state[3*node+axis];
	return result;
}

} // namespace

int main(int argc,char** argv)
{
	if(PetscInitialize(&argc,&argv,nullptr,nullptr))return 1;
	int rank=0,ranks=1,status=0;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	try{
		if(argc!=3||ranks!=1)throw std::invalid_argument(
			"usage: native_tet_fsi_steps case.json output_dir (one MPI rank)");
		const auto scenario=Parse(argv[1]);
		auto fluid_mesh=ReadMesh(scenario.fluid_mesh);
		auto solid_mesh=ReadMesh(scenario.solid_mesh);
		const auto interface_coordinates=InterfaceNodes(fluid_mesh,scenario.interface_label);
		RenumberInterface(fluid_mesh,scenario.interface_label,interface_coordinates);
		const auto solid_node_map=RenumberInterface(solid_mesh,
			scenario.interface_label,interface_coordinates);
		const auto identity=InterfaceIdentity(fluid_mesh,scenario.interface_label);
		const auto layout=Layout(fluid_mesh,scenario.interface_label,identity);
		const auto fluid_surface=Surface("fluid","native-ale",scenario.interface_label,
			identity,true);
		const auto solid_surface=Surface("solid","native-solid",scenario.interface_label,
			identity,false);
		const iga::FsiCouplingEdge edge("wall",fluid_surface.id,solid_surface.id,
			iga::FsiCouplingLaw::FluidStructureTractionKinematics);
		const auto topology=iga::BuildNativeTaylorHoodTopology(fluid_mesh);
		const std::size_t velocity_nodes=fluid_mesh.points.size()+topology.edges.size();
		std::vector<double> initial(3*velocity_nodes+fluid_mesh.points.size(),0.);
		std::set<std::uint32_t> boundary_nodes;
		std::set<std::uint32_t> interface_nodes;
		for(const auto& face:fluid_mesh.boundary_triangles)
			for(auto node:face.nodes){
				boundary_nodes.insert(node);
				if(face.boundary_label==scenario.interface_label)
					interface_nodes.insert(node);
			}
		for(auto node:topology.boundary_velocity_nodes.at(scenario.interface_label))
			interface_nodes.insert(node);
		for(std::size_t edge_index=0;edge_index<topology.edges.size();++edge_index){
			const auto edge_nodes=topology.edges[edge_index];
			if(!boundary_nodes.count(edge_nodes[0])
				||!boundary_nodes.count(edge_nodes[1]))
				initial[3*(fluid_mesh.points.size()+edge_index)]=scenario.initial_velocity;
		}
		std::map<std::uint32_t,std::array<double,3>> essential_velocity;
		for(const auto& by_label:topology.boundary_velocity_nodes)
			if(by_label.first!=scenario.interface_label
				&&!scenario.natural_labels.count(by_label.first)
				&&!scenario.prescribed_velocity.count(by_label.first))
				for(auto node:by_label.second)
					essential_velocity.emplace(node,std::array<double,3>{{0.,0.,0.}});
		for(const auto& boundary:scenario.prescribed_velocity)
			for(auto node:topology.boundary_velocity_nodes.at(boundary.first))
				essential_velocity.emplace(node,interface_nodes.count(node)
					?std::array<double,3>{{0.,0.,0.}}:boundary.second);
		iga::NativeTetAleFsiRuntime fluid("fluid","native-ale",edge,fluid_surface,
			solid_surface,layout,layout,fluid_mesh,scenario.interface_label,
			{scenario.density,scenario.viscosity},0,initial,0.,1e-6,30,
			essential_velocity,scenario.natural_labels);
		std::map<std::size_t,double> constraints;
		for(auto node:scenario.fixed_nodes)
			for(int axis=0;axis<3;++axis)
				constraints.emplace(3*solid_node_map.at(node)+axis,0.);
		iga::NativeTetSolidStaticOptions solid_options;
		solid_options.relative_tolerance=1e-5;
		solid_options.absolute_tolerance_n=1e-8;
		solid_options.maximum_iterations=40;
		iga::NativeTetSolidFsiRuntime solid("solid","native-solid",edge,solid_surface,
			fluid_surface,layout,layout,solid_mesh,scenario.interface_label,
			{scenario.young,scenario.poisson,scenario.solid_density},constraints,
			solid_options);
		iga::StrongFluidStructureCouplingOptions controls;
		controls.maximum_iterations=12;
		controls.absolute_displacement_tolerance_m=2e-5;
		controls.relative_displacement_tolerance=2e-2;
		controls.absolute_nodal_force_residual_tolerance_n=2e-2;
		controls.absolute_traction_residual_tolerance_pa=1e-1;
		controls.fluid_residual_tolerance=1e-6;
		controls.structure_residual_tolerance_n=1e-7;
		controls.interface_power_defect_tolerance_w=1e-10;
		controls.aitken_controls.initial_relaxation=.8;
		iga::StrongFluidStructureCoupling<iga::NativeTetAleFsiRuntime,
			iga::NativeTetSolidFsiRuntime> coordinator(fluid,solid,edge,layout,controls);
		const std::filesystem::path output(argv[2]);
		std::vector<std::pair<double,std::filesystem::path>> flow_series,solid_series;
		std::map<std::string,std::vector<double>> species_concentration;
		std::map<std::string,std::vector<std::pair<double,std::filesystem::path>>>
			species_series;
		for(const auto& species:scenario.species)
			species_concentration.emplace(species.id,
				std::vector<double>(fluid_mesh.points.size(),species.initial));
		for(int step=1;step<=scenario.steps;++step){
			const auto previous=CurrentMesh(fluid_mesh,fluid.CommittedAleDisplacementM());
			const auto before=VertexVelocity(fluid.CommittedFluidState(),fluid_mesh.points.size());
			const auto result=coordinator.Execute({step,
				(step-1)*scenario.dt_s,scenario.dt_s});
			if(!result.converged)throw std::runtime_error("native FSI step did not converge");
			const auto current=CurrentMesh(fluid_mesh,fluid.CommittedAleDisplacementM());
			const auto after=VertexVelocity(fluid.CommittedFluidState(),fluid_mesh.points.size());
			std::vector<std::array<double,3>> mesh_velocity(fluid_mesh.points.size());
			for(std::size_t node=0;node<mesh_velocity.size();++node)
				for(int axis=0;axis<3;++axis)
					mesh_velocity[node][axis]=(current.points[node][axis]
						-previous.points[node][axis])/scenario.dt_s;
			const auto conservation=iga::EvaluateNativeTetAleConservation(
				previous,current,mesh_velocity,before,after,scenario.dt_s);
			std::map<std::uint64_t,std::array<std::uint32_t,3>> interface_faces;
			for(const auto& face:current.boundary_triangles)
				if(face.boundary_label==scenario.interface_label)
					interface_faces.emplace(face.id,face.nodes);
			std::array<double,3> interface_force{};
			double interface_power=0.;
			for(const auto& traction:iga::EvaluateNativeTetFluidTractionOnMatchingStructure(
				current,topology,fluid.CommittedFluidState(),scenario.viscosity,
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
				const double area=.5*iga::NativeTetInterfaceNorm(
					iga::NativeTetInterfaceCross(ab,ac));
				for(int axis=0;axis<3;++axis){
					const double force=area*traction.traction_on_structure_pa[axis];
					interface_force[axis]+=force;
					interface_power+=force*mean_velocity[axis];
				}
			}
			const double force_norm=iga::NativeTetInterfaceNorm(interface_force);
			double maximum_wall_displacement=0.;
			for(double value:solid.CommittedDisplacementM())
				maximum_wall_displacement=std::max(maximum_wall_displacement,
					std::abs(value));
			iga::NativeTetLinearElasticResult solid_field;
			solid_field.displacement_m.resize(solid_mesh.points.size());
			for(std::size_t node=0;node<solid_mesh.points.size();++node)
				for(int axis=0;axis<3;++axis)
					solid_field.displacement_m[node][axis]=
						solid.CommittedDisplacementM()[3*node+axis];
			const auto flow_path=output/("flow_step_"+std::to_string(step));
			const auto solid_path=output/("wall_step_"+std::to_string(step));
			iga::WriteParallelVtkSnapshot(PETSC_COMM_WORLD,flow_path,
				iga::BuildNativeTetHydraulicVtkPartition(fluid_mesh,current,
					fluid.CommittedFluidState(),scenario.viscosity,rank,ranks),
				step*scenario.dt_s);
			iga::WriteParallelVtkSnapshot(PETSC_COMM_WORLD,solid_path,
				iga::BuildNativeTetLinearElasticVtkPartition(solid_mesh,solid_field,
					rank,ranks),step*scenario.dt_s);
			std::vector<std::array<double,3>> fluid_velocity(velocity_nodes);
			for(std::size_t node=0;node<velocity_nodes;++node)
				for(int axis=0;axis<3;++axis)
					fluid_velocity[node][axis]=fluid.CommittedFluidState()[3*node+axis];
			for(const auto& species:scenario.species){
				auto& concentration=species_concentration.at(species.id);
				const auto transport=iga::SolveNativeTetMovingSpeciesPetscStep(
					previous,current,fluid_velocity,mesh_velocity,concentration,
					species.inflow,species.diffusivity,species.source,
					scenario.dt_s,false,species.decay,species.wall_exchange);
				concentration=transport.step.concentration_mol_m3;
				const auto species_path=output/("species_"+species.id+"_step_"
					+std::to_string(step));
				iga::WriteParallelVtkSnapshot(PETSC_COMM_WORLD,species_path,
					iga::BuildNativeTetSpeciesVtkPartition(fluid_mesh,current,
						concentration,rank,ranks),step*scenario.dt_s);
				species_series[species.id].push_back({step*scenario.dt_s,
					species_path/"snapshot.pvtu"});
				if(rank==0)std::cout<<std::setprecision(17)
					<<"native_fsi_species_step id="<<species.id
					<<" step="<<step
					<<" inventory_mol="<<transport.step.current_inventory_mol
					<<" balance_defect_mol_s="<<transport.step.balance_defect_mol_s
					<<" linear_iterations="<<transport.linear_iterations<<'\n';
			}
			flow_series.push_back({step*scenario.dt_s,flow_path/"snapshot.pvtu"});
			solid_series.push_back({step*scenario.dt_s,solid_path/"snapshot.pvtu"});
		if(rank==0)std::cout<<std::setprecision(17)
			<<"native_fsi_step step="<<step
			<<" iterations="<<result.iterations
			<<" gcl_residual_m3_s="<<conservation.gcl_residual_m3_s
				<<" volume_rate_m3_s="<<conservation.volume_rate_m3_s
				<<" mesh_boundary_flux_m3_s="
					<<conservation.mesh_boundary_flux_m3_s
				<<" relative_boundary_flux_m3_s="
					<<conservation.relative_boundary_flux_m3_s
				<<" moving_balance_residual_m3_s="
					<<conservation.moving_domain_balance_residual_m3_s
				<<" interface_force_n="<<force_norm
				<<" interface_power_w="<<interface_power
				<<" wall_max_displacement_m="<<maximum_wall_displacement
			<<" interface_power_defect_w="
				<<*result.history.back().interface_power_defect_w<<'\n';
		}
		iga::WriteParallelVtkSeries(PETSC_COMM_WORLD,output/"flow.pvd",flow_series);
		iga::WriteParallelVtkSeries(PETSC_COMM_WORLD,output/"wall.pvd",solid_series);
		for(const auto& series:species_series)
			iga::WriteParallelVtkSeries(PETSC_COMM_WORLD,
				output/("species_"+series.first+".pvd"),series.second);
		if(rank==0)std::cout<<"native_fsi_complete accepted_steps="<<scenario.steps<<'\n';
	}catch(const std::exception& error){
		if(rank==0)std::cerr<<"native_tet_fsi_steps: "<<error.what()<<'\n';
		status=1;
	}
	PetscFinalize();return status;
}
