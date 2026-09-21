#include "CaseConfig.hpp"
#include "NativeTetAlePetscRuntime.hpp"
#include "NativeTetBoundaryFlow.hpp"
#include "NativeTetDarcyPetsc.hpp"
#include "NativeTetDarcyVisualization.hpp"
#include "NativeTetHydraulicVisualization.hpp"
#include "NativeTetMatchingInterfaceSource.hpp"
#include "NativeTetMovingSpeciesPetscRuntime.hpp"
#include "NativeTetSpeciesVisualization.hpp"
#include "ParallelVtkOutput.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <cmath>
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
#include <vector>

namespace {

using Object=std::map<std::string,iga::config_detail::JsonValue>;

const iga::config_detail::JsonValue& Field(const Object& object,const char* name)
{
	const auto* found=iga::config_detail::Find(object,name);
	if(!found)throw std::invalid_argument(std::string("chain requires ")+name);
	return *found;
}

const Object& Section(const Object& object,const char* name)
{
	return iga::config_detail::RequireObject(Field(object,name),name);
}

int Label(const Object& object,const char* name)
{
	return iga::config_detail::RequireInteger(Field(object,name),name);
}

double Number(const Object& object,const char* name)
{
	return iga::config_detail::RequireNumber(Field(object,name),name);
}

struct Case
{
	struct Species
	{
		std::string id;
		double initial=0.,diffusivity=0.,inlet=0.;
	};
	std::filesystem::path artery,tissue,vein;
	int artery_inlet=0,artery_wall=0,vein_wall=0,vein_outlet=0;
	std::array<double,3> inlet_velocity{};
	double density=0.,viscosity=0.,mobility=0.,venous_pressure=0.;
	std::vector<iga::NativeTetMatchingInterfacePort> arterial_ports;
	std::vector<std::pair<int,int>> venous_ports;
	std::vector<Species> species;
	double dt_s=0.;
	int steps=0;
};

Case Parse(const std::filesystem::path& file)
{
	using namespace iga::config_detail;
	std::ifstream input(file);
	if(!input)throw std::runtime_error("cannot open flow-Darcy case");
	std::ostringstream stream;stream<<input.rdbuf();
	const auto parsed=JsonParser(stream.str()).Parse();
	const auto& root=RequireObject(parsed,"flow-Darcy case");
	if(Label(root,"schema_version")!=1)
		throw std::invalid_argument("flow-Darcy case schema differs");
	Case result;
	const auto& meshes=Section(root,"meshes");
	result.artery=file.parent_path()/RequireString(Field(meshes,"artery"),"artery mesh");
	result.tissue=file.parent_path()/RequireString(Field(meshes,"tissue"),"tissue mesh");
	result.vein=file.parent_path()/RequireString(Field(meshes,"vein"),"vein mesh");
	const auto& artery=Section(root,"artery");
	result.artery_inlet=Label(artery,"inlet_label");
	result.artery_wall=Label(artery,"wall_label");
	const auto& velocity=RequireArray(Field(artery,"inlet_velocity_m_s"),"inlet velocity");
	if(velocity.size()!=3)throw std::invalid_argument("inlet velocity needs three components");
	for(int axis=0;axis<3;++axis)
		result.inlet_velocity[axis]=RequireNumber(velocity[axis],"inlet velocity");
	const auto& vein=Section(root,"vein");
	result.vein_wall=Label(vein,"wall_label");
	result.vein_outlet=Label(vein,"outlet_label");
	result.venous_pressure=Number(vein,"outlet_pressure_pa");
	const auto& material=Section(root,"materials");
	result.density=Number(material,"density_kg_m3");
	result.viscosity=Number(material,"dynamic_viscosity_pa_s");
	result.mobility=Number(material,"darcy_mobility_m2_pa_s");
	if(!(result.density>0.)||!(result.viscosity>0.)||!(result.mobility>0.))
		throw std::invalid_argument("flow-Darcy material must be positive");
	for(const auto& item:RequireArray(Field(root,"artery_to_tissue"),"artery_to_tissue")){
		const auto& port=RequireObject(item,"arterial exchange port");
		result.arterial_ports.push_back({RequireString(Field(port,"name"),"port name"),
			Label(port,"artery_label"),Label(port,"tissue_label")});
	}
	for(const auto& item:RequireArray(Field(root,"tissue_to_vein"),"tissue_to_vein")){
		const auto& port=RequireObject(item,"venous exchange port");
		result.venous_ports.push_back({Label(port,"tissue_label"),
			Label(port,"vein_label")});
	}
	if(result.arterial_ports.empty()||result.venous_ports.empty())
		throw std::invalid_argument("flow-Darcy chain requires both exchange sets");
	if(const auto* species=Find(root,"species")){
		const auto& time=Section(root,"time");
		result.dt_s=Number(time,"dt_s");result.steps=Label(time,"steps");
		if(!(result.dt_s>0.)||result.steps<1)
			throw std::invalid_argument("chain species time is invalid");
		std::set<std::string> ids;
		for(const auto& item:RequireArray(*species,"species")){
			const auto& object=RequireObject(item,"species[]");
			Case::Species entry;
			entry.id=RequireString(Field(object,"id"),"species id");
			entry.initial=Number(object,"initial_concentration_mol_m3");
			entry.diffusivity=Number(object,"diffusivity_m2_s");
			entry.inlet=Number(object,"inlet_concentration_mol_m3");
			if(entry.id.empty()||!ids.insert(entry.id).second
				||entry.id.find_first_not_of("abcdefghijklmnopqrstuvwxyz"
					"ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_")!=std::string::npos
				||entry.initial<0.||entry.diffusivity<0.||entry.inlet<0.)
				throw std::invalid_argument("chain species definition is invalid");
			result.species.push_back(std::move(entry));
		}
	}
	return result;
}

iga::NativeTetMesh ReadMesh(const std::filesystem::path& file)
{
	std::ifstream input(file);
	if(!input)throw std::runtime_error("cannot open tetra mesh: "+file.string());
	return iga::ReadNativeTetMeshGmsh41(input);
}

std::vector<double> ZeroState(const iga::NativeTetMesh& mesh,
	const iga::NativeTaylorHoodTopology& topology)
{
	return std::vector<double>(3*(mesh.points.size()+topology.edges.size())
		+mesh.points.size(),0.);
}

std::map<std::uint32_t,std::array<double,3>> VelocityBoundary(
	const iga::NativeTaylorHoodTopology& topology,int wall,int inlet,
	const std::array<double,3>& velocity)
{
	std::map<std::uint32_t,std::array<double,3>> result;
	for(auto node:topology.boundary_velocity_nodes.at(wall))
		result.emplace(node,std::array<double,3>{{0.,0.,0.}});
	if(inlet>=0)
		for(auto node:topology.boundary_velocity_nodes.at(inlet))
			result.emplace(node,velocity);
	return result;
}

double BoundaryFlow(const std::map<int,iga::NativeTetBoundaryFlowValue>& values,
	int label)
{
	return values.at(label).outward_flow_m3_s;
}

void Snapshot(MPI_Comm comm,const std::filesystem::path& output,const char* name,
	const iga::VtkPartition& piece)
{
	const auto path=output/name;
	iga::WriteParallelVtkSnapshot(comm,path,piece,0.);
	iga::WriteParallelVtkSeries(comm,output/(std::string(name)+".pvd"),
		{{0.,path/"snapshot.pvtu"}});
}

} // namespace

int main(int argc,char** argv)
{
	if(PetscInitialize(&argc,&argv,nullptr,nullptr))return 1;
	int rank=0,ranks=1,status=0;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	try{
		if(argc!=3)throw std::invalid_argument(
			"usage: native_tet_flow_darcy_chain case.json output_dir");
		const auto scenario=Parse(argv[1]);
		const auto artery=ReadMesh(scenario.artery);
		const auto tissue=ReadMesh(scenario.tissue);
		const auto vein=ReadMesh(scenario.vein);
		const auto artery_topology=iga::BuildNativeTaylorHoodTopology(artery);
		const auto vein_topology=iga::BuildNativeTaylorHoodTopology(vein);
		const std::vector<std::array<double,3>> artery_mesh_velocity(
			artery.points.size(),{{0.,0.,0.}});
		const std::vector<std::array<double,3>> vein_mesh_velocity(
			vein.points.size(),{{0.,0.,0.}});
		iga::NativeTetAleBoundaryConditions arterial_pressure;
		for(const auto& port:scenario.arterial_ports)
			arterial_pressure.prescribed_pressure_pa.emplace(
				port.vessel_boundary_label,0.);
		const auto arterial=iga::SolveNativeTetAlePetscSteady(artery,
			artery_mesh_velocity,ZeroState(artery,artery_topology),
			VelocityBoundary(artery_topology,scenario.artery_wall,
				scenario.artery_inlet,scenario.inlet_velocity),
			std::numeric_limits<std::uint32_t>::max(),
			{scenario.density,scenario.viscosity},1e-12,15,arterial_pressure);
		if(arterial.final_linear_converged_reason<=0)
			throw std::runtime_error("arterial flow did not converge");
		const auto artery_flows=iga::EvaluateNativeTetBoundaryFlows(artery,
			artery_topology,arterial.replicated_state);
		const auto mapped=iga::MapNativeTetMatchingInterfaceToTissueSource(
			artery,artery_topology,arterial.replicated_state,tissue,
			scenario.arterial_ports);
		std::map<int,double> tissue_pressure;
		for(const auto& port:scenario.venous_ports)
			tissue_pressure.emplace(port.first,scenario.venous_pressure);
		const std::vector<double> mobility(tissue.cells.size(),scenario.mobility);
		const auto darcy=iga::SolveNativeTetDarcyPetsc(tissue,mobility,
			mapped.source.tissue_source_s_inv,tissue_pressure);
		if(darcy.converged_reason<=0)
			throw std::runtime_error("Darcy solve did not converge");
		iga::NativeTetAleBoundaryConditions venous_ports;
		venous_ports.prescribed_pressure_pa.emplace(scenario.vein_outlet,
			scenario.venous_pressure);
		double transferred=0.;
		for(const auto& port:scenario.venous_ports){
			const double flow=darcy.conservative_outward_boundary_flow_m3_s.at(
				port.first);
			venous_ports.flow_rate_controls.push_back({port.second,-flow});
			transferred+=flow;
		}
		const auto venous=iga::SolveNativeTetAlePetscSteady(vein,
			vein_mesh_velocity,ZeroState(vein,vein_topology),
			VelocityBoundary(vein_topology,scenario.vein_wall,-1,{{0.,0.,0.}}),
			std::numeric_limits<std::uint32_t>::max(),
			{scenario.density,scenario.viscosity},1e-12,15,venous_ports);
		if(venous.final_linear_converged_reason<=0)
			throw std::runtime_error("venous flow did not converge");
		const auto vein_flows=iga::EvaluateNativeTetBoundaryFlows(vein,
			vein_topology,venous.replicated_state);
		const double inlet=BoundaryFlow(artery_flows,scenario.artery_inlet);
		const double outlet=BoundaryFlow(vein_flows,scenario.vein_outlet);
		const double scale=std::max(std::abs(inlet),1e-18);
		for(const auto& port:scenario.arterial_ports)
			if(std::abs(BoundaryFlow(artery_flows,port.vessel_boundary_label)
				-mapped.source.vessel_outward_flow_m3_s.at(port.name))>1e-6*scale)
				throw std::runtime_error("arterial interface flow differs from Darcy source");
		for(const auto& port:scenario.venous_ports)
			if(std::abs(BoundaryFlow(vein_flows,port.second)
				+darcy.conservative_outward_boundary_flow_m3_s.at(port.first))>1e-6*scale)
				throw std::runtime_error("venous interface flow differs from Darcy outlet");
		if(std::abs(inlet+mapped.source.total_vessel_outward_flow_m3_s)>1e-5*scale
			||std::abs(darcy.volume_source_m3_s-transferred)>1e-6*scale
			||std::abs(outlet-transferred)>1e-5*scale
			||darcy.maximum_cell_balance_defect_m3_s>1e-6*scale)
			throw std::runtime_error("flow-Darcy chain mass balance failed");
		const std::filesystem::path output(argv[2]);
		Snapshot(PETSC_COMM_WORLD,output,"artery",
			iga::BuildNativeTetHydraulicVtkPartition(artery,artery,
				arterial.replicated_state,scenario.viscosity,rank,ranks));
		Snapshot(PETSC_COMM_WORLD,output,"tissue",
			iga::BuildNativeTetDarcyVtkPartition(tissue,darcy,mobility,
				mapped.source.tissue_source_s_inv,rank,ranks));
		Snapshot(PETSC_COMM_WORLD,output,"vein",
			iga::BuildNativeTetHydraulicVtkPartition(vein,vein,
				venous.replicated_state,scenario.viscosity,rank,ranks));
		if(!scenario.species.empty()){
			for(const auto& port:scenario.arterial_ports)
				if(!(mapped.source.vessel_outward_flow_m3_s.at(port.name)>0.))
					throw std::runtime_error("species requires forward arterial exchange");
			for(const auto& transfer:mapped.facet_transfers)
				if(!(transfer.vessel_outward_flow_m3_s>0.))
					throw std::runtime_error("species requires forward arterial facets");
			for(const auto& port:scenario.venous_ports)
				if(!(darcy.conservative_outward_boundary_flow_m3_s.at(port.first)>0.))
					throw std::runtime_error("species requires forward tissue-to-vein exchange");
			const auto velocity=[&](const std::vector<double>& state,std::size_t count){
				std::vector<std::array<double,3>> values(count);
				for(std::size_t node=0;node<count;++node)
					for(int axis=0;axis<3;++axis)
						values[node][axis]=state[3*node+axis];
				return values;
			};
			const auto arterial_velocity=velocity(arterial.replicated_state,
				artery.points.size()+artery_topology.edges.size());
			const auto venous_velocity=velocity(venous.replicated_state,
				vein.points.size()+vein_topology.edges.size());
			const std::vector<std::array<double,3>> artery_grid(
				artery.points.size(),{{0.,0.,0.}}),
				tissue_grid(tissue.points.size(),{{0.,0.,0.}}),
				vein_grid(vein.points.size(),{{0.,0.,0.}});
			std::map<std::uint64_t,std::size_t> tissue_cells;
			for(std::size_t cell=0;cell<tissue.cells.size();++cell)
				tissue_cells.emplace(tissue.cells[cell].id,cell);
			for(const auto& species:scenario.species){
				std::vector<double> arterial_concentration(artery.points.size(),species.initial);
				std::vector<double> tissue_concentration(tissue.points.size(),species.initial);
				std::vector<double> venous_concentration(vein.points.size(),species.initial);
				std::array<std::vector<std::pair<double,std::filesystem::path>>,3> series;
				for(int step=1;step<=scenario.steps;++step){
					const auto arterial_step=iga::SolveNativeTetMovingSpeciesPetscStep(
						artery,artery,arterial_velocity,artery_grid,
						arterial_concentration,{{scenario.artery_inlet,species.inlet}},
						species.diffusivity,0.,scenario.dt_s,true);
					std::map<std::string,double> arterial_moles;
					for(const auto& port:scenario.arterial_ports)
						arterial_moles.emplace(port.name,
							arterial_step.step.outward_advective_flux_by_label_mol_s.at(
								port.vessel_boundary_label));
					for(const auto& item:arterial_moles)
						if(item.second<0.)
							throw std::runtime_error("arterial species exchange reversed");
					std::vector<double> cell_source(tissue.cells.size(),0.);
					for(const auto& transfer:mapped.facet_transfers){
						const auto cell=tissue_cells.at(transfer.tissue_cell_id);
						const double volume=iga::EvaluateNativeTetGeometry(tissue,
							tissue.cells[cell]).determinant/6.;
						cell_source[cell]+=transfer.vessel_outward_flow_m3_s*
							(arterial_moles.at(transfer.port_name)/
							mapped.source.vessel_outward_flow_m3_s.at(transfer.port_name))/volume;
					}
					const auto tissue_step=iga::SolveNativeTetMovingSpeciesPetscStep(
						tissue,tissue,{},tissue_grid,tissue_concentration,{},
						species.diffusivity,0.,scenario.dt_s,true,0.,{},
						darcy.conservative_face_flow_m3_s,cell_source);
					std::map<int,double> vein_donor;
					double arterial_mol=0.,tissue_mol=0.;
					for(const auto& item:arterial_moles)arterial_mol+=item.second;
					for(const auto& port:scenario.venous_ports){
						const double amount=tissue_step.step.outward_advective_flux_by_label_mol_s.at(
							port.first);
						if(amount<0.)
							throw std::runtime_error("tissue species exchange reversed");
						tissue_mol+=amount;
						vein_donor.emplace(port.second,amount/
							darcy.conservative_outward_boundary_flow_m3_s.at(port.first));
					}
					const auto venous_step=iga::SolveNativeTetMovingSpeciesPetscStep(
						vein,vein,venous_velocity,vein_grid,venous_concentration,
						vein_donor,species.diffusivity,0.,scenario.dt_s,true);
					double vein_inward=0.;
					for(const auto& port:scenario.venous_ports)
						vein_inward-=venous_step.step.outward_advective_flux_by_label_mol_s.at(
							port.second);
					const double defect=arterial_step.step.balance_defect_mol_s+
						tissue_step.step.balance_defect_mol_s+
						venous_step.step.balance_defect_mol_s+
						(tissue_step.step.source_mol_s-arterial_mol)+
						(vein_inward-tissue_mol);
					if(std::abs(defect)>1e-14+1e-7*std::abs(inlet*species.inlet))
						throw std::runtime_error("flow-Darcy species balance failed");
					arterial_concentration=arterial_step.step.concentration_mol_m3;
					tissue_concentration=tissue_step.step.concentration_mol_m3;
					venous_concentration=venous_step.step.concentration_mol_m3;
					const std::array<const char*,3> regions{{"artery","tissue","vein"}};
					const std::array<const iga::NativeTetMesh*,3> meshes{{&artery,&tissue,&vein}};
					const std::array<const std::vector<double>*,3> values{{
						&arterial_concentration,&tissue_concentration,&venous_concentration}};
					for(int region=0;region<3;++region){
						const auto path=output/("species_"+species.id+"_"+regions[region]
							+"_step_"+std::to_string(step));
						iga::WriteParallelVtkSnapshot(PETSC_COMM_WORLD,path,
							iga::BuildNativeTetSpeciesVtkPartition(*meshes[region],
								*meshes[region],*values[region],rank,ranks),step*scenario.dt_s);
						series[region].push_back({step*scenario.dt_s,path/"snapshot.pvtu"});
					}
					if(rank==0)std::cout<<std::setprecision(17)
						<<"native_flow_darcy_species_step id="<<species.id
						<<" step="<<step<<" global_defect_mol_s="<<defect<<'\n';
				}
				for(int region=0;region<3;++region){
					const std::array<const char*,3> names{{"artery","tissue","vein"}};
					iga::WriteParallelVtkSeries(PETSC_COMM_WORLD,
						output/("species_"+species.id+"_"+names[region]+".pvd"),
						series[region]);
				}
			}
		}
		if(rank==0)std::cout<<std::setprecision(17)
			<<"native_flow_darcy_chain: PASS inlet_m3_s="<<inlet
			<<" artery_to_tissue_m3_s="
			<<mapped.source.total_vessel_outward_flow_m3_s
			<<" tissue_to_vein_m3_s="<<transferred
			<<" outlet_m3_s="<<outlet
			<<" darcy_cell_defect_m3_s="
			<<darcy.maximum_cell_balance_defect_m3_s<<'\n';
	}catch(const std::exception& error){
		if(rank==0)std::cerr<<"native_tet_flow_darcy_chain: "<<error.what()<<'\n';
		status=1;
	}
	PetscFinalize();return status;
}
