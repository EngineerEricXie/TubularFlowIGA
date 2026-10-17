#include "PerfusedOrganHemodynamics.hpp"
#include "PerfusedOrganKinetics.hpp"
#include "PerfusedOrganTransport.hpp"

#include <petscsys.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

iga::NativeTetMesh Cube(unsigned subdivisions)
{
	iga::NativeTetMesh mesh;
	const auto node=[&](unsigned x,unsigned y,unsigned z){
		return static_cast<std::uint32_t>((z*(subdivisions+1)+y)*(subdivisions+1)+x);
	};
	for(unsigned z=0;z<=subdivisions;++z)
		for(unsigned y=0;y<=subdivisions;++y)
			for(unsigned x=0;x<=subdivisions;++x)
				mesh.points.push_back({{double(x)/subdivisions,double(y)/subdivisions,
					double(z)/subdivisions}});
	const std::array<std::array<unsigned,4>,6> pattern{{
		{{0,1,3,7}},{{0,3,2,7}},{{0,2,6,7}},{{0,6,4,7}},{{0,4,5,7}},{{0,5,1,7}}}};
	std::map<std::array<std::uint32_t,3>,unsigned> uses;
	std::uint64_t id=1;
	for(unsigned z=0;z<subdivisions;++z)
		for(unsigned y=0;y<subdivisions;++y)
			for(unsigned x=0;x<subdivisions;++x){
				const std::array<std::uint32_t,8> corners{{node(x,y,z),node(x+1,y,z),
					node(x,y+1,z),node(x+1,y+1,z),node(x,y,z+1),node(x+1,y,z+1),
					node(x,y+1,z+1),node(x+1,y+1,z+1)}};
				for(const auto& local:pattern){
					iga::NativeTetCell cell;
					cell.id=id++;
					for(std::size_t i=0;i<4;++i)cell.nodes[i]=corners[local[i]];
					mesh.cells.push_back(cell);
					for(std::size_t opposite=0;opposite<4;++opposite){
						std::array<std::uint32_t,3> face{};
						std::size_t entry=0;
						for(std::size_t k=0;k<4;++k)if(k!=opposite)face[entry++]=cell.nodes[k];
						std::sort(face.begin(),face.end());
						++uses[face];
					}
				}
			}
	for(const auto& item:uses)
		if(item.second==1)mesh.boundary_triangles.push_back({id++,item.first,1});
	return mesh;
}

void Close(double actual,double expected,double tolerance,const char* what)
{
	const bool ok=std::isfinite(actual)
		&&std::abs(actual-expected)<=tolerance*std::max(1e-300,std::abs(expected));
	if(!ok)std::cerr<<std::setprecision(15)<<what<<": actual="<<actual<<" expected="<<expected<<"\n";
	assert(ok);
}

std::filesystem::path WriteText(const std::string& name,const std::string& text)
{
	const auto path=std::filesystem::temp_directory_path()/name;
	std::ofstream(path)<<text;
	return path;
}

std::array<double,3> Centroid(const iga::NativeTetMesh& mesh,std::size_t cell)
{
	std::array<double,3> c{};
	for(const auto node:mesh.cells[cell].nodes)
		for(int axis=0;axis<3;++axis)c[axis]+=0.25*mesh.points[node][axis];
	return c;
}

// Inflow tree: root -> branch point -> terminals in the left and right halves.
// Outflow tree: one terminal draining every cell.
std::vector<iga::PerfusedOrganTree> CubeTrees(const iga::NativeTetMesh& mesh,double viscosity)
{
	const auto inflow_path=WriteText("perfused_organ_inflow.swc",
		"1 1 -0.2 0.5 0.5 0.01 -1\n2 3 0.0 0.5 0.5 0.01 1\n"
		"3 3 0.25 0.5 0.5 0.005 2\n4 3 0.75 0.5 0.5 0.005 2\n");
	const auto outflow_path=WriteText("perfused_organ_outflow.swc",
		"1 1 1.2 0.5 0.5 0.01 -1\n2 3 0.5 0.5 0.9 0.01 1\n");
	iga::PerfusedOrganTree inflow{"inflow",true,iga::ReadOneDNetwork(inflow_path,1.,1,viscosity),{}};
	iga::PerfusedOrganTree outflow{"outflow",false,iga::ReadOneDNetwork(outflow_path,1.,1,viscosity),{}};
	const int left=inflow.network.node_index.at(3),right=inflow.network.node_index.at(4);
	const int drain=outflow.network.node_index.at(2);
	double left_volume=0.,right_volume=0.,total=0.;
	for(std::size_t i=0;i<mesh.cells.size();++i){
		const double v=iga::PerfusedOrganCellVolume(mesh,i);
		(Centroid(mesh,i)[0]<0.5?left_volume:right_volume)+=v;
		total+=v;
	}
	for(std::size_t i=0;i<mesh.cells.size();++i){
		const double v=iga::PerfusedOrganCellVolume(mesh,i);
		if(Centroid(mesh,i)[0]<0.5)inflow.territories[left].push_back({i,v/left_volume});
		else inflow.territories[right].push_back({i,v/right_volume});
		outflow.territories[drain].push_back({i,v/total});
	}
	return {inflow,outflow};
}

iga::HepaticKinetics Liver()
{
	iga::HepaticKinetics k;
	k.vmax_o2_mol_m3_s=0.02;
	k.glycolysis_capacity_mol_m3_s=0.02;
	k.lactate_vmax_mol_m3_s=0.004;
	k.glycogenolysis_vmax_mol_m3_s=0.003;
	k.glucose_uptake_vmax_mol_m3_s=0.003;
	k.insulin_clearance_s_inv=1e-3;
	return k;
}

iga::PerfusateChemistry Perfusate()
{
	iga::PerfusateChemistry chemistry;
	chemistry.haemoglobin_g_dl=4.6;
	return chemistry;
}

iga::OrganConcentration Inlet()
{
	return {{3.0,8.0,1.0,24.0,1.2,100.0}};
}

iga::PerfusedOrganFlowResult CubeFlow(const iga::NativeTetMesh& mesh,double mobility)
{
	const auto trees=CubeTrees(mesh,1e-3);
	iga::PerfusedOrganFlowInput input;
	input.viscosity_pa_s=1e-3;
	input.root_pressure_pa={1000.,0.};
	input.specific_conductance_inv_pa_s={2e-5,4e-5};
	input.mobility_by_cell_m2_pa_s.assign(mesh.cells.size(),mobility);
	return iga::SolvePerfusedOrganFlow(mesh,trees,input);
}

void TestFlowConservation()
{
	// Tissue conductance K L comparable to the terminal conductances, as in an
	// organ whose tissue resistance matters.
	const auto mesh=Cube(4);
	const auto result=CubeFlow(mesh,1e-5);
	const double q_in=result.root_inflow_m3_s[0],q_out=-result.root_inflow_m3_s[1];
	Close(q_out,q_in,1e-9,"organ outflow equals inflow");
	assert(std::abs(result.exchange_imbalance_m3_s)<=1e-10*q_in);
	double terminals=0.;
	for(const auto& item:result.terminal_flow_m3_s[0])terminals+=item.second;
	Close(terminals,q_in,1e-9,"terminal flows sum to the root flow");
	const auto& flows=result.terminal_flow_m3_s[0];
	// The near terminal (shorter branch) receives more flow.
	assert(std::next(flows.begin())->second>0.);
	assert(flows.begin()->second>std::next(flows.begin())->second);
	assert(result.darcy.maximum_cell_balance_defect_m3_s<=1e-12*q_in);
	double in=0.,out=0.;
	for(std::size_t i=0;i<mesh.cells.size();++i){
		in+=result.cell_exchange_in_m3_s[i];
		out+=result.cell_exchange_out_m3_s[i];
	}
	Close(in,q_in,1e-9,"cell exchange inflow");
	Close(out,q_in,1e-9,"cell exchange outflow");
}

void TestSeriesLimit()
{
	// With near-infinite tissue mobility the organ is a resistor network. The
	// 1e6 conductance ratio limits conservation to about 1e-6 by round-off in
	// the Darcy rows, so this check is loose by construction.
	const auto mesh=Cube(4);
	const double viscosity=1e-3;
	const auto result=CubeFlow(mesh,10.);
	const auto conductance=[&](double radius,double length){
		return M_PI*std::pow(radius,4)/(8.*viscosity*length);};
	const auto series=[](double a,double b){return a*b/(a+b);};
	const double trunk=conductance(0.01,0.2);
	const double halves=series(conductance(0.0075,0.25),2e-5*0.5)
		+series(conductance(0.0075,0.75),2e-5*0.5);
	const double drain_vessel=conductance(0.01,std::hypot(0.7,0.4));
	const double total=1./(1./trunk+1./halves+1./(4e-5*1.)+1./drain_vessel);
	Close(result.root_inflow_m3_s[0],total*1000.,1e-4,"series resistance limit");
	Close(-result.root_inflow_m3_s[1],total*1000.,1e-4,"series resistance outflow");
}

void TestOxygenInversion()
{
	auto chemistry=Perfusate();
	for(const double ph:{7.0,7.4,7.7})
		for(const double po2:{0.5,5.,26.8,90.,600.}){
			const double content=iga::OxygenContent(po2,ph,chemistry);
			Close(iga::PartialPressureO2(content,ph,chemistry),po2,1e-9,"pO2 inversion");
		}
	assert(iga::PartialPressureO2(0.,7.4,chemistry)==0.);
}

void TestKineticsLimits()
{
	const auto chemistry=Perfusate();
	auto k=Liver();
	iga::HepaticCellState state{1.,100.};
	auto c=Inlet();
	c[iga::OrganOxygen]=iga::OxygenContent(500.,7.4,chemistry);
	const auto oxygenated=iga::EvaluateHepaticRates(c,state,k,chemistry);
	assert(oxygenated.energy_deficit<1e-3);
	assert(oxygenated.lactate_production<1e-2*k.vmax_o2_mol_m3_s*k.atp_per_o2);
	c[iga::OrganOxygen]=0.;
	k.glycolysis_capacity_mol_m3_s=0.01;
	const auto anoxic=iga::EvaluateHepaticRates(c,state,k,chemistry);
	assert(anoxic.o2_uptake==0.&&anoxic.lactate_clearance==0.);
	const double fuel=100./120.+(20./120.)*8./9.;
	Close(anoxic.lactate_production,0.01*fuel,1e-12,"glycolysis capacity");
	Close(anoxic.energy_deficit,1.-0.01*fuel/(0.02*5.),1e-12,"energy deficit");
	Close(anoxic.reaction[iga::OrganBicarbonate],-anoxic.lactate_production,1e-12,"acid load");
	state.viability=0.;
	const auto dead=iga::EvaluateHepaticRates(Inlet(),state,k,chemistry);
	for(const double value:dead.reaction)assert(value==0.);
}

iga::OrganTransportGraph SingleCell(double volume,double flow)
{
	iga::OrganTransportGraph graph;
	graph.volume_m3={volume};
	graph.face_outflow_m3_s={0.};
	graph.boundary_outflow_m3_s={0.};
	graph.exchange_in_m3_s={flow};
	graph.exchange_out_m3_s={flow};
	graph.upstream={{}};
	return graph;
}

iga::OrganTransportInput TransportInput(std::size_t cells,double dt)
{
	iga::OrganTransportInput input;
	input.dt_s=dt;
	input.inlet=Inlet();
	input.distribution_volume_fraction={{0.2,0.6,0.6,0.3,0.3,0.2}};
	input.previous.assign(cells,Inlet());
	input.state.assign(cells,{1.,100.});
	input.kinetics=Liver();
	input.chemistry=Perfusate();
	return input;
}

void TestSingleCellOxygenAgainstScalarRoot()
{
	const double volume=1e-3,flow=2e-6;
	auto input=TransportInput(1,1e9);
	input.chemistry.bohr_coefficient=0.;
	const auto result=iga::SolveOrganTransport(SingleCell(volume,flow),input);
	assert(result.newton_failures==0);
	// Q (C_in - C) = V vmax p/(Km + p) with C = content(p), solved by bisection.
	const auto f=[&](double p){
		return flow*(input.inlet[iga::OrganOxygen]-iga::OxygenContent(p,7.4,input.chemistry))
			-volume*input.kinetics.vmax_o2_mol_m3_s*p/(input.kinetics.km_o2_mmhg+p);};
	double low=0.,high=iga::PartialPressureO2(input.inlet[iga::OrganOxygen],7.4,input.chemistry);
	for(int i=0;i<200;++i){const double mid=0.5*(low+high);(f(mid)>0.?low:high)=mid;}
	// Content converges to ~1e-10; pO2 amplifies that near anoxia.
	Close(result.rates[0].po2_mmhg,0.5*(low+high),1e-5,"single-cell pO2");
	for(std::size_t s=0;s<iga::OrganSpeciesCount;++s)
		assert(std::abs(result.balance_defect[s])<=1e-9*std::max(1e-12,std::abs(result.inflow[s])));
}

void TestStagnantCellRelaxes()
{
	auto input=TransportInput(1,100.);
	input.kinetics=iga::HepaticKinetics{};
	input.kinetics.insulin_clearance_s_inv=1e-3;
	const auto result=iga::SolveOrganTransport(SingleCell(1e-3,0.),input);
	Close(result.concentration[0][iga::OrganInsulin],100./(1.+100.*1e-3/0.2),1e-9,
		"backward-Euler insulin clearance");
	Close(result.concentration[0][iga::OrganGlucose],8.,1e-12,"inert glucose");
}

void TestCubeTransport()
{
	const auto mesh=Cube(4);
	const auto hemodynamics=CubeFlow(mesh,1e-5);
	const auto graph=iga::BuildTetTransportGraph(mesh,hemodynamics.darcy.conservative_face_flow_m3_s,
		hemodynamics.cell_exchange_in_m3_s,hemodynamics.cell_exchange_out_m3_s);
	// Inert species stay at the inlet value everywhere in the steady limit.
	auto inert=TransportInput(mesh.cells.size(),1e9);
	inert.previous.assign(mesh.cells.size(),iga::OrganConcentration{});
	inert.kinetics=iga::HepaticKinetics{};
	const auto passive=iga::SolveOrganTransport(graph,inert);
	for(const auto& c:passive.concentration)
		for(std::size_t s=0;s<iga::OrganSpeciesCount;++s)
			Close(c[s],inert.inlet[s],1e-6,"inert steady concentration");
	auto living=TransportInput(mesh.cells.size(),600.);
	const auto result=iga::SolveOrganTransport(graph,living);
	assert(result.newton_failures==0);
	for(std::size_t s=0;s<iga::OrganSpeciesCount;++s){
		const double scale=std::max({std::abs(result.inflow[s]),std::abs(result.reaction[s]),1e-30});
		assert(std::abs(result.balance_defect[s])<=1e-6*scale);
	}
	for(const auto& c:result.concentration)for(const double value:c)assert(value>=0.);
	assert(result.reaction[iga::OrganOxygen]<0.);
	std::cout<<"cube: Q="<<hemodynamics.root_inflow_m3_s[0]<<" components="<<result.components
		<<" largest="<<result.largest_component<<" sweeps="<<result.gauss_seidel_sweeps<<'\n';
}

} // namespace

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	TestFlowConservation();
	TestSeriesLimit();
	TestOxygenInversion();
	TestKineticsLimits();
	TestSingleCellOxygenAgainstScalarRoot();
	TestStagnantCellRelaxes();
	TestCubeTransport();
	std::cout<<"perfused organ tests passed\n";
	PetscFinalize();
	return 0;
}
