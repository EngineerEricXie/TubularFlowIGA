#include "NativeTetAlePetscRuntime.hpp"
#include "NativeTetBoundaryFlow.hpp"
#include "NativeTetMatchingInterfaceSource.hpp"
#include "NativeTetMeshComponents.hpp"
#include "NativeTetMovingSpeciesPetscRuntime.hpp"
#include "NativeTetSpeciesVisualization.hpp"
#include "ParallelVtkOutput.hpp"

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
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Point=std::array<double,3>;

iga::NativeTetMesh ReadMesh(const char* path)
{
	std::ifstream input(path);
	if(!input)throw std::runtime_error(std::string("cannot open tetra mesh: ")+path);
	return iga::ReadNativeTetMeshGmsh41(input);
}

int Integer(const char* text,const char* name)
{
	std::size_t used=0;
	const int value=std::stoi(text,&used);
	if(used!=std::string(text).size()||value<0)
		throw std::invalid_argument(std::string(name)+" is invalid");
	return value;
}

double Number(const char* text,const char* name)
{
	std::size_t used=0;
	const double value=std::stod(text,&used);
	if(used!=std::string(text).size()||!std::isfinite(value))
		throw std::invalid_argument(std::string(name)+" is invalid");
	return value;
}

struct FluidState
{
	iga::NativeTetMesh current;
	iga::NativeTaylorHoodTopology topology;
	std::vector<Point> velocity_m_s;
	std::vector<double> state;
};

FluidState ReadFluidState(const char* path,const iga::NativeTetMesh& reference)
{
	FluidState result;
	result.current=reference;
	result.topology=iga::BuildNativeTaylorHoodTopology(reference);
	const std::size_t vertices=reference.points.size();
	const std::size_t count=vertices+result.topology.edges.size();
	std::ifstream input(path);
	std::size_t rows=0;
	if(!input||!(input>>rows)||rows!=count)
		throw std::invalid_argument("frozen fluid field count differs from native P2 mesh");
	result.velocity_m_s.resize(count);
	result.state.assign(3*count+vertices,0.);
	std::vector<Point> current(count),reported_reference(count);
	std::vector<bool> seen(count,false);
	for(std::size_t row=0;row<count;++row){
		std::size_t id=0;
		Point velocity{},base{},position{};
		if(!(input>>id>>velocity[0]>>velocity[1]>>velocity[2]
			>>base[0]>>base[1]>>base[2]
			>>position[0]>>position[1]>>position[2])||id>=count||seen[id])
			throw std::invalid_argument("frozen fluid field row or global ID is invalid");
		seen[id]=true;
		for(int axis=0;axis<3;++axis){
			if(!std::isfinite(velocity[axis])||!std::isfinite(base[axis])
				||!std::isfinite(position[axis]))
				throw std::invalid_argument("frozen fluid field is nonfinite");
			result.state[3*id+axis]=velocity[axis];
		}
		result.velocity_m_s[id]=velocity;
		reported_reference[id]=base;
		current[id]=position;
	}
	std::string trailing;
	if(input>>trailing)
		throw std::invalid_argument("frozen fluid field has trailing data");
	for(std::size_t id=0;id<count;++id){
		const Point expected=id<vertices?reference.points[id]:Point{{
			0.5*(reference.points[result.topology.edges[id-vertices][0]][0]
				+reference.points[result.topology.edges[id-vertices][1]][0]),
			0.5*(reference.points[result.topology.edges[id-vertices][0]][1]
				+reference.points[result.topology.edges[id-vertices][1]][1]),
			0.5*(reference.points[result.topology.edges[id-vertices][0]][2]
				+reference.points[result.topology.edges[id-vertices][1]][2])}};
		for(int axis=0;axis<3;++axis)
			if(std::abs(reported_reference[id][axis]-expected[axis])>1e-12)
				throw std::invalid_argument("frozen fluid field reference mesh differs");
		if(id<vertices)result.current.points[id]=current[id];
	}
	for(std::size_t id=vertices;id<count;++id){
		const auto& edge=result.topology.edges[id-vertices];
		for(int axis=0;axis<3;++axis)
			if(std::abs(current[id][axis]-0.5*(result.current.points[edge[0]][axis]
				+result.current.points[edge[1]][axis]))>1e-12)
				throw std::invalid_argument("frozen P2 edge coordinate differs from moved tetra");
	}
	for(const auto& cell:result.current.cells)
		iga::EvaluateNativeTetGeometry(result.current,cell);
	return result;
}

std::vector<std::array<double,4>> ReadTissueFlux(const char* path,
	const iga::NativeTetMesh& tissue)
{
	std::ifstream input(path);
	std::size_t rows=0;
	if(!input||!(input>>rows)||rows!=tissue.cells.size())
		throw std::invalid_argument("frozen Darcy RT0 face-flow count differs");
	std::map<std::uint64_t,std::size_t> cells;
	for(std::size_t index=0;index<tissue.cells.size();++index)
		if(!cells.emplace(tissue.cells[index].id,index).second)
			throw std::invalid_argument("native tissue cell IDs are duplicate");
	std::vector<std::array<double,4>> result(rows);
	std::set<std::uint64_t> seen;
	for(std::size_t row=0;row<rows;++row){
		std::uint64_t id=0;
		std::array<double,4> flow{};
		if(!(input>>id>>flow[0]>>flow[1]>>flow[2]>>flow[3])
			||!cells.count(id)||!seen.insert(id).second)
			throw std::invalid_argument("frozen Darcy RT0 cell row or ID is invalid");
		for(const double value:flow)
			if(!std::isfinite(value))
				throw std::invalid_argument("frozen Darcy RT0 face flow is nonfinite");
		result[cells.at(id)]=flow;
	}
	std::string trailing;
	if(input>>trailing)
		throw std::invalid_argument("frozen Darcy RT0 field has trailing data");
	return result;
}

double BoundaryFlow(const iga::NativeTetMesh& mesh,
	const std::vector<std::array<double,4>>& faces,int label)
{
	using namespace iga::native_tet_matching_detail;
	const auto owners=BoundaryOwners(mesh);
	double total=0.;std::size_t count=0;
	for(const auto& triangle:mesh.boundary_triangles){
		if(triangle.boundary_label!=label)continue;
		auto nodes=triangle.nodes;
		std::sort(nodes.begin(),nodes.end());
		const auto found=owners.find(nodes);
		if(found==owners.end())
			throw std::invalid_argument("Darcy labelled face has no unique owner");
		total+=faces.at(found->second.cell).at(found->second.opposite);
		++count;
	}
	if(count==0||!std::isfinite(total))
		throw std::invalid_argument("Darcy boundary flow label is absent or invalid");
	return total;
}

double LabelFlux(const iga::NativeTetMovingSpeciesStep& step,int label)
{
	const auto found=step.outward_advective_flux_by_label_mol_s.find(label);
	if(found==step.outward_advective_flux_by_label_mol_s.end()
		||!std::isfinite(found->second))
		throw std::runtime_error("species boundary flux label is absent or nonfinite");
	return found->second;
}

void GateStep(const iga::NativeTetMovingSpeciesPetscResult& result,
	const char* region,double scale)
{
	if(result.converged_reason<=0||!std::isfinite(result.residual_norm)
		||std::abs(result.step.balance_defect_mol_s)>1e-14+1e-7*scale)
		throw std::runtime_error(std::string(region)+" species solve/balance failed");
	for(const double concentration:result.step.concentration_mol_m3)
		if(!std::isfinite(concentration)||concentration<-1e-12)
			throw std::runtime_error(std::string(region)+" species concentration is invalid");
}

void Near(double actual,double expected,double scale,const char* description)
{
	if(!std::isfinite(actual)||!std::isfinite(expected)
		||std::abs(actual-expected)>1e-14+1e-6*scale)
		throw std::runtime_error(std::string(description)+" conservation gate failed");
}

} // namespace

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int rank=0,ranks=1,status=0;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	try{
		if(argc!=24)throw std::invalid_argument(
			"usage: native_tet_cube_oxygen artery.msh tissue.msh vein.msh "
			"artery_state.txt tissue_rt0.txt vein_state.txt artery_inlet "
			"artery_tip0..3 vein_tip0..3 vein_outlet diffusivity_m2_s "
			"inlet_concentration_mol_m3 dt_s min_steps max_steps "
			"breakthrough_fraction output_dir");
		const auto artery=ReadMesh(argv[1]),tissue=ReadMesh(argv[2]),vein=ReadMesh(argv[3]);
		iga::RequireNativeTetSingleFaceComponent(artery,"oxygen arterial lumen");
		iga::RequireNativeTetSingleFaceComponent(tissue,"oxygen fixed tissue");
		iga::RequireNativeTetSingleFaceComponent(vein,"oxygen venous lumen");
		const auto a=ReadFluidState(argv[4],artery);
		const auto face_flow=ReadTissueFlux(argv[5],tissue);
		const auto v=ReadFluidState(argv[6],vein);
		const int artery_inlet=Integer(argv[7],"artery inlet label");
		std::array<int,4> artery_tips{},vein_tips{};
		for(int i=0;i<4;++i)artery_tips[i]=Integer(argv[8+i],"artery tip label");
		for(int i=0;i<4;++i)vein_tips[i]=Integer(argv[12+i],"vein tip label");
		const int vein_outlet=Integer(argv[16],"vein outlet label");
		std::set<int> distinct{artery_inlet,vein_outlet};
		for(int label:artery_tips)distinct.insert(label);
		for(int label:vein_tips)distinct.insert(label);
		if(distinct.size()!=10)
			throw std::invalid_argument("oxygen port labels are not distinct");
		const double diffusion=Number(argv[17],"diffusivity");
		const double inlet_concentration=Number(argv[18],"inlet concentration");
		const double dt=Number(argv[19],"time step");
		const int min_steps=Integer(argv[20],"minimum steps");
		const int max_steps=Integer(argv[21],"maximum steps");
		const double breakthrough=Number(argv[22],"breakthrough fraction");
		if(diffusion<0.||!(inlet_concentration>0.)||!(dt>0.)||min_steps<2
			||max_steps<min_steps||max_steps>1000
			||!(breakthrough>0.&&breakthrough<1.))
			throw std::invalid_argument("oxygen functional material/time contract is invalid");
		const std::filesystem::path output(argv[23]);
		const std::vector<Point> artery_grid(artery.points.size(),Point{{0,0,0}});
		const std::vector<Point> tissue_grid(tissue.points.size(),Point{{0,0,0}});
		const std::vector<Point> vein_grid(vein.points.size(),Point{{0,0,0}});
		std::vector<double> artery_concentration(artery.points.size(),0.);
		std::vector<double> tissue_concentration(tissue.points.size(),0.);
		std::vector<double> vein_concentration(vein.points.size(),0.);
		std::vector<iga::NativeTetMatchingInterfacePort> ports;
		for(int i=0;i<4;++i)ports.push_back({"tip_"+std::to_string(i),
			artery_tips[i],artery_tips[i]});
		const auto mapped=iga::MapNativeTetMatchingInterfaceToTissueSource(
			a.current,a.topology,a.state,tissue,ports);
		std::map<std::uint64_t,std::size_t> tissue_cell_ids;
		for(std::size_t index=0;index<tissue.cells.size();++index)
			tissue_cell_ids.emplace(tissue.cells[index].id,index);
		std::array<double,4> artery_water{},tissue_water{},vein_water{};
		const auto artery_boundary=iga::EvaluateNativeTetBoundaryFlows(
			a.current,a.topology,a.state);
		const auto vein_boundary=iga::EvaluateNativeTetBoundaryFlows(
			v.current,v.topology,v.state);
		for(int i=0;i<4;++i){
			const auto key="tip_"+std::to_string(i);
			artery_water[i]=mapped.source.vessel_outward_flow_m3_s.at(key);
			tissue_water[i]=BoundaryFlow(tissue,face_flow,vein_tips[i]);
			vein_water[i]=-vein_boundary.at(vein_tips[i]).outward_flow_m3_s;
			if(!(artery_water[i]>0.)||!(tissue_water[i]>0.)||!(vein_water[i]>0.))
				throw std::runtime_error("oxygen terminal water flow is not forward");
			Near(vein_water[i],tissue_water[i],tissue_water[i],
				"oxygen tissue/vein hydraulic port");
		}
		const double inlet_water=-artery_boundary.at(artery_inlet).outward_flow_m3_s;
		const double outlet_water=vein_boundary.at(vein_outlet).outward_flow_m3_s;
		if(!(inlet_water>0.)||!(outlet_water>0.))
			throw std::runtime_error("oxygen root water flow is not forward");
		Near(outlet_water,inlet_water,inlet_water,"oxygen root hydraulic");
		for(int label:artery_tips)
			Near(BoundaryFlow(tissue,face_flow,label),0.,inlet_water,
				"oxygen tissue arterial cap no-flux");
		iga::CollectiveLocalStage(PETSC_COMM_WORLD,"oxygen output directory",[&]{
			if(rank==0&&!std::filesystem::create_directory(output))
				throw std::runtime_error("oxygen output already exists");
		});
		std::ofstream ledger;
		if(rank==0){
			ledger.open(output/"ledger.csv");
			if(!ledger)throw std::runtime_error("cannot create oxygen ledger");
			ledger<<"step,time_s,artery_inventory_mol,tissue_inventory_mol,vein_inventory_mol,"
				<<"artery_to_tissue_mol_s,tissue_to_vein_mol_s,vein_inward_mol_s,"
				<<"venous_outlet_mol_s,venous_outlet_concentration_mol_m3,"
				<<"global_balance_defect_mol_s,artery_reason,tissue_reason,vein_reason\n";
		}
		std::array<std::vector<std::pair<double,std::filesystem::path>>,3> series;
		auto write_step=[&](int step){
			const double time=step*dt;
			const std::array<std::string,3> names{{"artery","tissue","vein"}};
			const std::array<const iga::NativeTetMesh*,3> references{{&artery,&tissue,&vein}};
			const std::array<const iga::NativeTetMesh*,3> currents{{
				&a.current,&tissue,&v.current}};
			const std::array<const std::vector<double>*,3> values{{
				&artery_concentration,&tissue_concentration,&vein_concentration}};
			for(int domain=0;domain<3;++domain){
				const auto snapshot=output/(names[domain]+"_step_"+std::to_string(step));
				iga::WriteParallelVtkSnapshot(PETSC_COMM_WORLD,snapshot,
					iga::BuildNativeTetSpeciesVtkPartition(*references[domain],
						*currents[domain],*values[domain],rank,ranks),time);
				series[domain].push_back({time,snapshot/"snapshot.pvtu"});
			}
		};
		write_step(0);
		int accepted=0;double final_outlet_concentration=0.;
		for(int step=1;step<=max_steps;++step){
			const auto arterial=iga::SolveNativeTetMovingSpeciesPetscStep(
				a.current,a.current,a.velocity_m_s,artery_grid,artery_concentration,
				{{artery_inlet,inlet_concentration}},diffusion,0.,dt,true);
			GateStep(arterial,"arterial",inlet_water*inlet_concentration);
			std::array<double,4> artery_mol{};
			for(int i=0;i<4;++i){
				artery_mol[i]=LabelFlux(arterial.step,artery_tips[i]);
				if(artery_mol[i]<0.)
					throw std::runtime_error("oxygen arterial terminal species reversed");
			}
			std::vector<double> cell_source(tissue.cells.size(),0.);
			for(const auto& transfer:mapped.facet_transfers){
				const int index=std::stoi(transfer.port_name.substr(4));
				const double face_water=transfer.vessel_outward_flow_m3_s;
				if(face_water<0.)
					throw std::runtime_error("oxygen arterial terminal facet reverses");
				const std::size_t cell=tissue_cell_ids.at(transfer.tissue_cell_id);
				const double volume=iga::EvaluateNativeTetGeometry(tissue,
					tissue.cells[cell]).determinant/6.;
				cell_source[cell]+=face_water*(artery_mol[index]/artery_water[index])/volume;
			}
			const auto tissue_step=iga::SolveNativeTetMovingSpeciesPetscStep(
				tissue,tissue,{},tissue_grid,tissue_concentration,{},
				diffusion,0.,dt,true,0.,{},face_flow,cell_source);
			GateStep(tissue_step,"tissue",inlet_water*inlet_concentration);
			double artery_transfer=0.,tissue_transfer=0.;
			std::map<int,double> vein_donor;
			for(int i=0;i<4;++i){
				artery_transfer+=artery_mol[i];
				const double mol=LabelFlux(tissue_step.step,vein_tips[i]);
				if(mol<0.)
					throw std::runtime_error("oxygen tissue venous terminal species reversed");
				tissue_transfer+=mol;
				vein_donor.emplace(vein_tips[i],mol/tissue_water[i]);
			}
			Near(tissue_step.step.source_mol_s,artery_transfer,
				inlet_water*inlet_concentration,"oxygen artery/tissue species");
			const auto venous=iga::SolveNativeTetMovingSpeciesPetscStep(
				v.current,v.current,v.velocity_m_s,vein_grid,vein_concentration,
				vein_donor,diffusion,0.,dt,true);
			GateStep(venous,"venous",inlet_water*inlet_concentration);
			double vein_inward=0.;
			for(int label:vein_tips)vein_inward-=LabelFlux(venous.step,label);
			Near(vein_inward,tissue_transfer,inlet_water*inlet_concentration,
				"oxygen tissue/vein species");
			const double outlet_mol=LabelFlux(venous.step,vein_outlet);
			if(outlet_mol<0.)
				throw std::runtime_error("oxygen venous outlet species reversed");
			const double global_defect=arterial.step.balance_defect_mol_s+
				tissue_step.step.balance_defect_mol_s+venous.step.balance_defect_mol_s+
				(tissue_step.step.source_mol_s-artery_transfer)+
				(vein_inward-tissue_transfer);
			if(std::abs(global_defect)>1e-14+1e-7*inlet_water*inlet_concentration)
				throw std::runtime_error("oxygen three-domain global balance failed");
			artery_concentration=arterial.step.concentration_mol_m3;
			tissue_concentration=tissue_step.step.concentration_mol_m3;
			vein_concentration=venous.step.concentration_mol_m3;
			write_step(step);
			accepted=step;
			final_outlet_concentration=outlet_mol/outlet_water;
			if(rank==0){
				ledger<<std::setprecision(17)<<step<<','<<step*dt<<','
					<<arterial.step.current_inventory_mol<<','
					<<tissue_step.step.current_inventory_mol<<','
					<<venous.step.current_inventory_mol<<','
					<<artery_transfer<<','<<tissue_transfer<<','<<vein_inward<<','
					<<outlet_mol<<','<<final_outlet_concentration<<','
					<<global_defect<<','<<arterial.converged_reason<<','
					<<tissue_step.converged_reason<<','<<venous.converged_reason<<'\n';
				ledger.flush();
				std::cout<<std::setprecision(12)<<"oxygen_step="<<step
					<<" tissue_inventory_mol="<<tissue_step.step.current_inventory_mol
					<<" venous_outlet_concentration_mol_m3="
					<<final_outlet_concentration<<std::endl;
			}
			if(step>=min_steps&&final_outlet_concentration>=
				breakthrough*inlet_concentration)break;
		}
		if(accepted<min_steps||final_outlet_concentration<
			breakthrough*inlet_concentration)
			throw std::runtime_error("oxygen venous outlet breakthrough not reached");
		for(int domain=0;domain<3;++domain)
			iga::WriteParallelVtkSeries(PETSC_COMM_WORLD,
				output/(std::array<std::string,3>{{"artery","tissue","vein"}}[domain]
					+"_oxygen.pvd"),series[domain]);
		if(rank==0){
			std::ofstream report(output/"summary.json");
			if(!report)throw std::runtime_error("cannot write oxygen summary");
			report<<std::setprecision(17)
				<<"{\n  \"kind\": \"idealized_cube_passive_oxygen_functional_only\",\n"
				<<"  \"physiological_validation\": false,\n"
				<<"  \"frozen_post_fsi_velocity\": true,\n"
				<<"  \"tissue_porosity_assumed\": 1.0,\n"
				<<"  \"tissue_consumption_mol_m3_s\": 0.0,\n"
				<<"  \"ranks\": "<<ranks<<",\n"
				<<"  \"accepted_steps\": "<<accepted<<",\n"
				<<"  \"final_time_s\": "<<accepted*dt<<",\n"
				<<"  \"inlet_concentration_mol_m3\": "<<inlet_concentration<<",\n"
				<<"  \"diffusivity_m2_s\": "<<diffusion<<",\n"
				<<"  \"venous_outlet_concentration_mol_m3\": "
				<<final_outlet_concentration<<",\n"
				<<"  \"breakthrough_fraction\": "<<breakthrough<<"\n}\n";
			std::cout<<"native_cube_oxygen: PASS steps="<<accepted
				<<" outlet_concentration_mol_m3="<<final_outlet_concentration
				<<std::endl;
		}
	}catch(const std::exception& error){
		if(rank==0)std::cerr<<"native_cube_oxygen: ERROR: "<<error.what()<<'\n';
		status=2;
	}
	PetscFinalize();return status;
}
