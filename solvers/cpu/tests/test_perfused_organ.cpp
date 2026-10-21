#include "PerfusedOrganHemodynamics.hpp"
#include "PerfusedOrganKinetics.hpp"
#include "PerfusedOrganPlant.hpp"
#include "PerfusedOrganTransport.hpp"

#include <petscsys.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
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
	chemistry.oxygen.hemoglobin_g_dl=4.6;
	return chemistry;
}

// HCO3- 24 mmol/L and pCO2 40 mmHg at Hb 4.6 g/dL: pH from Henderson-
// Hasselbalch, buffer base = HCO3- + A0 + beta (pH - 7.4).
iga::OrganConcentration Inlet()
{
	const double ph=6.1+std::log10(24.0/(0.0307*40.0));
	const double buffer_base=24.0+17.7+0.42*4.6+(1.43*4.6+7.7)*(ph-7.4);
	return {{3.0,8.0,1.0,buffer_base,24.0+0.0307*40.0,100.0}};
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

void TestAcidBase()
{
	const auto chemistry=Perfusate();
	const auto inlet=Inlet();
	const auto normal=iga::SolveAcidBase(inlet[iga::OrganBufferBase],inlet[iga::OrganCarbonDioxide],chemistry);
	Close(normal.bicarbonate_mol_m3,24.,1e-9,"normal bicarbonate");
	Close(normal.pco2_mmhg,40.,1e-9,"normal pCO2");
	Close(normal.ph,6.1+std::log10(24./(0.0307*40.)),1e-12,"Henderson-Hasselbalch");
	// Adding CO2 at fixed buffer base: most becomes HCO3- through the
	// non-bicarbonate buffer, so pCO2 rises far less than without it.
	const auto loaded=iga::SolveAcidBase(inlet[iga::OrganBufferBase],inlet[iga::OrganCarbonDioxide]+2.,chemistry);
	assert(loaded.ph<normal.ph&&loaded.bicarbonate_mol_m3>24.);
	assert(loaded.pco2_mmhg<40.+2./0.0307);
	// A strong acid lowers buffer base one-for-one and lowers pH.
	const auto acid=iga::SolveAcidBase(inlet[iga::OrganBufferBase]-5.,inlet[iga::OrganCarbonDioxide],chemistry);
	assert(acid.ph<normal.ph);
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
	Close(anoxic.reaction[iga::OrganBufferBase],-anoxic.lactate_production,1e-12,"acid load");
	// A steep O2 dependence removes the aerobic shortfall at normoxic pO2.
	{
		auto steep=Liver();
		steep.km_o2_mmhg=0.5;
		steep.o2_hill_coefficient=3.;
		auto normoxic=Inlet();
		normoxic[iga::OrganOxygen]=iga::OxygenContent(20.,7.4,chemistry);
		const auto rates=iga::EvaluateHepaticRates(normoxic,state,steep,chemistry);
		Close(rates.o2_uptake,steep.vmax_o2_mol_m3_s/(1.+std::pow(0.5/rates.po2_mmhg,3.)),1e-12,"Hill O2 uptake");
		assert(rates.lactate_production<1e-4*steep.vmax_o2_mol_m3_s);
	}
	// Aerobic glucose oxidation: six O2 per glucose, RQ 1 for that part.
	{
		auto oxidiser=Liver();
		oxidiser.aerobic_glucose_fraction=1.;
		oxidiser.glycolysis_capacity_mol_m3_s=0.;
		oxidiser.glycogenolysis_vmax_mol_m3_s=0.;
		oxidiser.glucose_uptake_vmax_mol_m3_s=0.;
		oxidiser.lactate_vmax_mol_m3_s=0.;
		auto rich=Inlet();
		rich[iga::OrganOxygen]=iga::OxygenContent(500.,7.4,chemistry);
		const iga::HepaticCellState no_glycogen{1.,0.};
		const auto burning=iga::EvaluateHepaticRates(rich,no_glycogen,oxidiser,chemistry);
		const double fuel=8./(1.+8.);
		Close(-burning.reaction[iga::OrganGlucose],fuel*burning.o2_uptake/6.,1e-12,"glucose oxidation");
		Close(burning.reaction[iga::OrganCarbonDioxide],(fuel+(1.-fuel)*oxidiser.respiratory_quotient)
			*burning.o2_uptake,1e-12,"oxidation CO2");
	}
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

// An anoxic, glycogen-fuelled cell with no buffer-base supply consumes buffer
// base faster than inflow and storage can provide: the species balance has no
// nonnegative root, the clip blocks every Newton step, and the cell must be
// reported as a Newton failure with a visible balance defect rather than
// counted as converged.
void TestBlockedCellIsReportedAsFailure()
{
	const double volume=1e-3,flow=2e-6;
	// Hypoxic lactate production consumes more buffer base than the inlet
	// supplies, so the clip pins buffer base at zero. With the Bohr shift off
	// the other species converge, the clipped step stops moving, and the cell
	// must still be counted as a Newton failure instead of converged.
	auto input=TransportInput(1,1e9);
	input.chemistry.bohr_coefficient=0.;
	input.inlet[iga::OrganOxygen]*=0.05;
	input.inlet[iga::OrganBufferBase]=0.1;
	input.previous[0]=input.inlet;
	input.state[0].glycogen_mol_m3=500.;
	const auto result=iga::SolveOrganTransport(SingleCell(volume,flow),input);
	const auto& rates=result.rates[0];
	assert(rates.lactate_production>rates.lactate_clearance);
	assert(result.concentration[0][iga::OrganBufferBase]==0.);
	assert(result.newton_failures==1);
	assert(std::abs(result.balance_defect[iga::OrganBufferBase])
		>0.5*volume*(rates.lactate_production-rates.lactate_clearance));
}

// A recirculating two-cell loop: Gauss-Seidel contracts by about F/(F+X) per
// sweep, so a large recirculation cannot settle within the sweep cap and must
// be reported rather than returned as converged.
iga::OrganTransportResult SolveTwoCellLoop(double recirculation)
{
	const double exchange=1e-6;
	iga::OrganTransportGraph graph;
	graph.volume_m3={1e-3,1e-3};
	graph.face_outflow_m3_s={recirculation+exchange,recirculation};
	graph.boundary_outflow_m3_s={0.,0.};
	graph.exchange_in_m3_s={exchange,0.};
	graph.exchange_out_m3_s={0.,exchange};
	graph.upstream={{{1,recirculation}},{{0,recirculation+exchange}}};
	// A long step makes storage negligible against the 1e-6 m^3/s exchange.
	auto input=TransportInput(2,1e15);
	input.kinetics=iga::HepaticKinetics{};
	input.previous.assign(2,iga::OrganConcentration{});
	return iga::SolveOrganTransport(graph,input);
}

void TestGaussSeidelNonConvergenceIsReported()
{
	const auto settled=SolveTwoCellLoop(1e-6);
	assert(settled.largest_component==2&&settled.gauss_seidel_unconverged==0);
	Close(settled.concentration[1][iga::OrganGlucose],8.,1e-9,"settled loop outflow");
	const auto stalled=SolveTwoCellLoop(1e-1);
	assert(stalled.gauss_seidel_sweeps==500&&stalled.gauss_seidel_unconverged==1);
}

template<typename Function>
bool Rejects(Function function)
{
	try{function();}catch(const std::invalid_argument&){return true;}
	return false;
}

void TestGlycolysisAcidInhibition()
{
	// Anoxic tissue at a lowered buffer base: glycolysis runs at its cap.
	auto c=Inlet();
	c[iga::OrganOxygen]=0.;
	c[iga::OrganBufferBase]-=12.;
	const iga::HepaticCellState state{1.,100.};
	const auto chemistry=Perfusate();
	auto k=Liver();
	const auto free=iga::EvaluateHepaticRates(c,state,k,chemistry);
	k.glycolysis_ph50=std::nan("");
	const auto disabled=iga::EvaluateHepaticRates(c,state,k,chemistry);
	assert(disabled.lactate_production==free.lactate_production&&disabled.reaction==free.reaction);
	k.glycolysis_ph50=6.8;
	k.glycolysis_ph_hill=4.;
	const auto inhibited=iga::EvaluateHepaticRates(c,state,k,chemistry);
	Close(inhibited.ph,free.ph,0.,"inhibition leaves the acid-base state");
	Close(inhibited.lactate_production,free.lactate_production/(1.+std::pow(10.,4.*(6.8-free.ph))),1e-14,
		"acid-inhibited lactate production");
	Close(inhibited.reaction[iga::OrganBufferBase],-inhibited.lactate_production,1e-14,
		"buffer base follows the inhibited production");
	assert(inhibited.energy_deficit>free.energy_deficit);
	for(const auto& bad:std::vector<std::pair<double,double>>{{5.,4.},{9.5,4.},{7.,0.},{7.,-1.}}){
		k.glycolysis_ph50=bad.first;
		k.glycolysis_ph_hill=bad.second;
		assert(Rejects([&]{iga::ValidateHepaticKinetics(k,chemistry);}));
	}
}

void TestAcidInhibitionKeepsBufferBase()
{
	// The clip-blocked cell of TestBlockedCellIsReportedAsFailure: with acid
	// inhibition, production falls with pH before buffer base runs out.
	const double volume=1e-3,flow=2e-6;
	auto input=TransportInput(1,1e9);
	input.chemistry.bohr_coefficient=0.;
	input.inlet[iga::OrganOxygen]*=0.05;
	input.inlet[iga::OrganBufferBase]=0.1;
	input.previous[0]=input.inlet;
	input.state[0].glycogen_mol_m3=500.;
	input.kinetics.glycolysis_ph50=6.8;
	input.kinetics.glycolysis_ph_hill=4.;
	const auto result=iga::SolveOrganTransport(SingleCell(volume,flow),input);
	assert(result.newton_failures==0);
	assert(result.concentration[0][iga::OrganBufferBase]>0.);
	for(std::size_t s=0;s<iga::OrganSpeciesCount;++s)
		assert(std::abs(result.balance_defect[s])<=1e-8*flow*std::max(input.inlet[s],1e-3));
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

void TestExchangeScaleMatchesConductance()
{
	// Doubling every cell's exchange multiplier equals doubling both trees'
	// specific conductances.
	const auto mesh=Cube(3);
	const auto trees=CubeTrees(mesh,1e-3);
	iga::PerfusedOrganFlowInput input;
	input.viscosity_pa_s=1e-3;
	input.root_pressure_pa={1000.,0.};
	input.specific_conductance_inv_pa_s={4e-5,8e-5};
	input.mobility_by_cell_m2_pa_s.assign(mesh.cells.size(),1e-5);
	const auto doubled=iga::SolvePerfusedOrganFlow(mesh,trees,input);
	input.specific_conductance_inv_pa_s={2e-5,4e-5};
	input.exchange_scale_by_cell.assign(mesh.cells.size(),2.);
	const auto scaled=iga::SolvePerfusedOrganFlow(mesh,trees,input);
	Close(scaled.root_inflow_m3_s[0],doubled.root_inflow_m3_s[0],1e-10,"exchange scale");
}

void TestLobularZonesAreTanksInSeries()
{
	// Steady first-order insulin clearance through N equal zones:
	// c_out = c_in / (1 + k V / (N Q))^N.
	const double volume=1e-3,flow=1e-6,k=2e-3;
	for(const std::size_t zones:{std::size_t(1),std::size_t(3)}){
		auto input=TransportInput(zones,1e12);
		input.kinetics=iga::HepaticKinetics{};
		input.kinetics.insulin_clearance_s_inv=k;
		const auto graph=iga::SplitIntoLobularZones(SingleCell(volume,flow),zones);
		const auto result=iga::SolveOrganTransport(graph,input);
		const double expected=100./std::pow(1.+k*volume/(zones*flow),double(zones));
		Close(result.concentration.back()[iga::OrganInsulin],expected,1e-8,"tanks in series");
		for(std::size_t z=1;z<zones;++z)
			assert(result.concentration[z][iga::OrganInsulin]<result.concentration[z-1][iga::OrganInsulin]);
		Close(result.outflow[iga::OrganInsulin],flow*expected,1e-8,"zone outflow");
	}
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
	// Three lobular zones: same balance, oxygen falls from zone 1 to zone 3.
	const auto zoned_graph=iga::SplitIntoLobularZones(graph,3);
	auto zoned_input=TransportInput(3*mesh.cells.size(),600.);
	const auto zoned=iga::SolveOrganTransport(zoned_graph,zoned_input);
	assert(zoned.newton_failures==0);
	for(std::size_t s=0;s<iga::OrganSpeciesCount;++s){
		const double scale=std::max({std::abs(zoned.inflow[s]),std::abs(zoned.reaction[s]),1e-30});
		assert(std::abs(zoned.balance_defect[s])<=1e-6*scale);
	}
	for(std::size_t i=0;i<mesh.cells.size();++i)
		if(graph.exchange_in_m3_s[i]>0.)
			assert(zoned.rates[3*i+2].po2_mmhg<=zoned.rates[3*i].po2_mmhg+1e-9);
	// Lumping into two pools (x halves) and into one pool conserves species.
	std::vector<std::size_t> halves(mesh.cells.size());
	for(std::size_t i=0;i<mesh.cells.size();++i)halves[i]=Centroid(mesh,i)[0]<0.5?0:1;
	for(const std::size_t groups:{std::size_t(1),std::size_t(2)}){
		std::vector<std::size_t> group=halves;
		if(groups==1)std::fill(group.begin(),group.end(),0);
		const auto pooled_graph=iga::SplitIntoLobularZones(iga::AggregateTransportGraph(graph,group,groups),3);
		double volume=0.;
		for(const double v:pooled_graph.volume_m3)volume+=v;
		Close(volume,1.,1e-12,"pooled volume");
		auto pooled_input=TransportInput(3*groups,600.);
		const auto pooled=iga::SolveOrganTransport(pooled_graph,pooled_input);
		assert(pooled.newton_failures==0);
		for(std::size_t s=0;s<iga::OrganSpeciesCount;++s){
			const double scale=std::max({std::abs(pooled.inflow[s]),std::abs(pooled.reaction[s]),1e-30});
			assert(std::abs(pooled.balance_defect[s])<=1e-6*scale);
		}
	}
	std::cout<<"cube: Q="<<hemodynamics.root_inflow_m3_s[0]<<" components="<<result.components
		<<" largest="<<result.largest_component<<" sweeps="<<result.gauss_seidel_sweeps<<'\n';
}

iga::CouplingDefinition ReservoirCircuit(double volume)
{
	iga::CouplingDefinition definition;
	definition.mode=iga::SimulationScopeMode::VcaClosedLoop;
	definition.external_circuit.reservoir.enabled=true;
	definition.external_circuit.reservoir.volume_m3=volume;
	definition.external_circuit.reservoir.species={{"x",2.0},{"y",1.0}};
	definition.external_circuit.pump.mode="pressure_control";
	definition.external_circuit.pump.pressure_pa=1000.;
	definition.perfusate.oxygen.hematocrit_percent=20.;
	return definition;
}

void TestCircuitExtensions()
{
	const double volume=2e-3,k=1e-6,dt=3600.;
	{
		// First-order clearance U = k c is integrated exactly for any step.
		iga::VcaExternalCircuit circuit(ReservoirCircuit(volume));
		const auto inlet=circuit.InletState(0.);
		assert(inlet.has_pressure&&!inlet.has_flow&&inlet.pressure_pa==1000.);
		circuit.AdvanceWithUptake({{"x",k*2.0}},{},0.,dt,0.);
		Close(circuit.State().species.at("x"),2.0*std::exp(-k*dt/volume),1e-12,"exact clearance");
		Close(circuit.State().species.at("y"),1.0,1e-15,"untouched species");
	}
	{
		// Dialysate inflow (zero solute) dilutes: c V is conserved.
		auto definition=ReservoirCircuit(volume);
		definition.external_circuit.dialyzer.enabled=true;
		definition.external_circuit.dialyzer.inflow_m3_s=1e-7;
		definition.external_circuit.dialyzer.dialysate={{"x",0.},{"y",0.}};
		iga::VcaExternalCircuit circuit(definition);
		circuit.AdvanceWithUptake({},{},0.,dt,0.);
		const double v1=volume+1e-7*dt;
		Close(circuit.State().volume_m3,v1,1e-12,"dialysis volume");
		Close(circuit.State().species.at("x"),2.0*volume/v1,1e-12,"exact dilution");
		Close(circuit.State().hematocrit_percent,20.*volume/v1,1e-12,"haematocrit dilution");
	}
	{
		// Outflow removes perfusate at reservoir composition.
		auto definition=ReservoirCircuit(volume);
		definition.external_circuit.dialyzer.enabled=true;
		definition.external_circuit.dialyzer.outflow_m3_s=1e-7;
		iga::VcaExternalCircuit circuit(definition);
		circuit.AdvanceWithUptake({},{},0.,dt,0.);
		Close(circuit.State().volume_m3,volume-1e-7*dt,1e-12,"outflow volume");
		Close(circuit.State().species.at("x"),2.0,1e-12,"outflow keeps composition");
	}
	{
		// A gas-reset species relaxes to arterial - U/q at the circuit flow.
		iga::VcaExternalCircuit circuit(ReservoirCircuit(volume));
		const double q=2e-5,arterial=4.0,uptake=1e-5;
		circuit.AdvanceWithUptake({{"y",-uptake}},{{"y",arterial}},q,1e5,0.);
		Close(circuit.State().species.at("y"),arterial+uptake/q,1e-9,"reset species steady state");
		iga::VcaExternalCircuit bolus(ReservoirCircuit(volume));
		bolus.ApplyBolus({{"x",1e-3}},volume);
		Close(bolus.State().species.at("x"),(2.0*volume+1e-3)/(2.*volume),1e-12,"bolus mass");
		Close(bolus.State().volume_m3,2.*volume,1e-15,"bolus volume");
		// A rejected bolus must leave the reservoir untouched.
		iga::VcaExternalCircuit rejected(ReservoirCircuit(volume));
		bool threw=false;
		try{rejected.ApplyBolus({{"x",1e-3},{"z",-1.}},volume);}catch(const std::runtime_error&){threw=true;}
		assert(threw);
		Close(rejected.State().species.at("x"),2.0,1e-15,"rejected bolus keeps composition");
		Close(rejected.State().volume_m3,volume,1e-15,"rejected bolus keeps volume");
	}
}

struct PlantCase
{
	iga::PerfusedOrganModel model;
	iga::OrganPlantParameters parameters;
	iga::OrganPlantState state;
	iga::OrganPlantControls controls;
};

// The unit-cube plant. `perfusion` scales both exchange conductances and the
// mobility; `metabolism` scales every Vmax of the liver kinetics.
PlantCase MakePlantCase(double perfusion,double metabolism)
{
	PlantCase plant;
	auto& model=plant.model;
	model.mesh=Cube(3);
	model.trees=CubeTrees(model.mesh,3e-3);
	model.zones=3;
	iga::FinalizePerfusedOrganModel(model);
	const std::size_t cells=model.mesh.cells.size(),nodes=3*cells;
	auto& p=plant.parameters;
	p.viscosity_pa_s=3e-3;
	p.outflow_pressure_pa={{"outflow",0.}};
	p.specific_conductance_inv_pa_s={{"inflow",2e-8*perfusion},{"outflow",6e-8*perfusion}};
	p.mobility_m2_pa_s=1e-8*perfusion;
	p.distribution_volume_fraction={{0.2,0.7,0.7,0.3,0.5,0.25}};
	p.kinetics=Liver();
	for(double* rate:{&p.kinetics.vmax_o2_mol_m3_s,&p.kinetics.glycolysis_capacity_mol_m3_s,
		&p.kinetics.lactate_vmax_mol_m3_s,&p.kinetics.glycogenolysis_vmax_mol_m3_s,
		&p.kinetics.glucose_uptake_vmax_mol_m3_s})*rate*=metabolism;
	p.chemistry=Perfusate();
	auto& state=plant.state;
	state.reservoir.volume_m3=2.0;
	state.reservoir.hematocrit_percent=18.;
	const auto& names=iga::PerfusedOrganSpeciesNames();
	const auto inlet=Inlet();
	for(std::size_t s=0;s<iga::OrganSpeciesCount;++s)state.reservoir.species[names[s]]=inlet[s];
	state.tissue.assign(nodes,inlet);
	state.viability.assign(nodes,1.);
	state.glycogen_mol_m3.assign(nodes,100.);
	state.oedema.assign(cells,0.);
	auto& controls=plant.controls;
	controls.inlet_pressure_pa=1000.;
	controls.fio2=0.6;
	controls.bolus={{names[iga::OrganGlucose],5.0}};
	controls.dialysis_in_m3_s=1e-8;
	return plant;
}

void TestPlantAdvance()
{
	// The baseline cube is anoxic: lactate production exhausts buffer base,
	// which the clip pins at zero, so every substep must report Newton
	// failures instead of committing silently.
	const auto plant=MakePlantCase(1.,1.);
	const auto& names=iga::PerfusedOrganSpeciesNames();
	const auto inlet=Inlet();
	for(const std::string resolution:{"cell","portal_territory","organ"}){
		auto copy=plant.state;
		const auto records=iga::AdvancePerfusedOrganPlant(plant.model,plant.parameters,plant.controls,{},{},
			copy,3600.,3,resolution);
		assert(records.size()==3);
		Close(copy.time_s,3600.,1e-12,"plant time");
		Close(copy.reservoir.volume_m3,2.0+1e-8*3600.,1e-9,"plant volume");
		for(const auto& record:records){
			assert(record.flow_m3_s>0.&&record.arterial.po2_mmhg>record.venous.po2_mmhg);
			assert(record.arterial.pco2_mmhg>0.&&record.mean_viability<=1.);
			assert(record.newton_failures>0&&record.gauss_seidel_unconverged==0);
			assert(std::isfinite(record.flow_residual_norm)&&record.flow_residual_norm<1e-12);
		}
		assert(copy.reservoir.species.at(names[iga::OrganGlucose])>inlet[iga::OrganGlucose]);
		for(const double v:copy.viability)assert(v>0.&&v<=1.);
	}
	const auto advance=[&](const iga::OrganPlantParameters& parameters,iga::OrganPlantState trial){
		iga::AdvancePerfusedOrganPlant(plant.model,parameters,plant.controls,{},{},trial,60.,1,"organ");};
	for(const auto& change:std::vector<std::function<void(iga::OrganPlantParameters&)>>{
		[](iga::OrganPlantParameters& q){q.oedema_maximum=1.;},
		[](iga::OrganPlantParameters& q){q.oxygenator_efficiency=1.5;},
		[](iga::OrganPlantParameters& q){q.carbogen_o2_fraction=q.air_o2_fraction;},
		[](iga::OrganPlantParameters& q){q.injury_rate_s_inv=-1.;},
		[](iga::OrganPlantParameters& q){q.vasodilator_ec50_mg_l=0.;}}){
		auto bad=plant.parameters;
		change(bad);
		assert(Rejects([&]{advance(bad,plant.state);}));
	}
	for(const auto& change:std::vector<std::function<void(iga::OrganPlantState&)>>{
		[](iga::OrganPlantState& q){q.viability[0]=1.5;},
		[](iga::OrganPlantState& q){q.glycogen_mol_m3[0]=-1.;},
		[](iga::OrganPlantState& q){q.oedema[0]=1.;},
		[](iga::OrganPlantState& q){q.tissue[0][iga::OrganLactate]=std::nan("");},
		[](iga::OrganPlantState& q){q.vasodilator_mg_l=-1.;}}){
		auto bad=plant.state;
		change(bad);
		assert(Rejects([&]{advance(plant.parameters,bad);}));
	}
	assert(!Rejects([&]{advance(plant.parameters,plant.state);}));
}

void TestAcidInhibitedPlantConserves()
{
	// The anoxic cube of TestPlantAdvance with acid inhibition of glycolysis:
	// buffer base stays positive, and every substep converges and conserves.
	const auto inlet=Inlet();
	for(const auto& setting:std::vector<std::pair<double,double>>{{6.8,4.},{6.6,6.}}){
		auto plant=MakePlantCase(1.,1.);
		plant.parameters.kinetics.glycolysis_ph50=setting.first;
		plant.parameters.kinetics.glycolysis_ph_hill=setting.second;
		for(const std::string resolution:{"cell","portal_territory","organ"}){
			auto copy=plant.state;
			const auto records=iga::AdvancePerfusedOrganPlant(plant.model,plant.parameters,plant.controls,{},{},
				copy,3600.,3,resolution);
			for(const auto& record:records){
				assert(record.newton_failures==0&&record.gauss_seidel_unconverged==0);
				for(std::size_t s=0;s<iga::OrganSpeciesCount;++s)
					assert(std::abs(record.balance_defect[s])<=1e-8*record.flow_m3_s*inlet[s]);
			}
			for(const auto& c:copy.tissue)assert(c[iga::OrganBufferBase]>0.);
		}
	}
}

void TestOxygenatedPlantConserves()
{
	// Thirtyfold perfusion and a hundredth of the metabolism keep every node
	// oxygenated; then every substep converges and conserves every species.
	const auto plant=MakePlantCase(30.,0.01);
	const auto inlet=Inlet();
	for(const std::string resolution:{"cell","portal_territory","organ"}){
		auto copy=plant.state;
		const auto records=iga::AdvancePerfusedOrganPlant(plant.model,plant.parameters,plant.controls,{},{},
			copy,3600.,3,resolution);
		for(const auto& record:records){
			assert(record.hypoxic_volume_fraction==0.&&record.venous.ph>7.3);
			assert(record.newton_failures==0&&record.gauss_seidel_unconverged==0);
			for(std::size_t s=0;s<iga::OrganSpeciesCount;++s)
				assert(std::abs(record.balance_defect[s])<=1e-8*record.flow_m3_s*inlet[s]);
		}
	}
}

void TestPlantControlCosts()
{
	// Costs of aggressive control on the oxygenated cube: hyperoxic
	// vasoconstriction lowers flow and oxidative and shear injury lower
	// viability; all are off by default.
	auto plant=MakePlantCase(30.,0.01);
	plant.controls.fio2=0.95;
	plant.controls.inlet_pressure_pa=2000.;
	auto baseline=plant.state,costly=plant.state;
	const auto plain=iga::AdvancePerfusedOrganPlant(plant.model,plant.parameters,plant.controls,{},{},
		baseline,3600.,2,"cell");
	auto harsh=plant.parameters;
	harsh.hyperoxic_constriction_per_100mmhg=0.2;
	harsh.oxidative_injury_rate_s_inv=1e-5;
	harsh.shear_injury_rate_inv_pa_s=1e-8;
	harsh.shear_injury_threshold_pa=100.;
	const auto penalised=iga::AdvancePerfusedOrganPlant(plant.model,harsh,plant.controls,{},{},
		costly,3600.,2,"cell");
	for(const auto& record:penalised)assert(record.newton_failures==0&&record.gauss_seidel_unconverged==0);
	assert(penalised.back().flow_m3_s<plain.back().flow_m3_s);
	assert(penalised.back().mean_viability<plain.back().mean_viability);
	for(const auto& change:std::vector<std::function<void(iga::OrganPlantParameters&)>>{
		[](iga::OrganPlantParameters& q){q.hyperoxic_constriction_per_100mmhg=-1.;},
		[](iga::OrganPlantParameters& q){q.oxidative_injury_rate_s_inv=-1.;},
		[](iga::OrganPlantParameters& q){q.shear_injury_rate_inv_pa_s=std::nan("");},
		[](iga::OrganPlantParameters& q){q.oxidative_injury_threshold_mmhg=std::nan("");}}){
		auto bad=plant.parameters;
		change(bad);
		auto trial=plant.state;
		assert(Rejects([&]{iga::AdvancePerfusedOrganPlant(plant.model,bad,plant.controls,{},{},trial,60.,1,
			"organ");}));
	}
}

void TestDecoupledGasCo2()
{
	// An independent sweep-gas CO2 decouples arterial pCO2 from FiO2.
	const auto plant=MakePlantCase(1.,1.);
	const auto chemistry=Perfusate();
	const auto venous=Inlet();
	const auto low=iga::Oxygenate(venous,0.3,plant.parameters,chemistry,0.05);
	const auto high=iga::Oxygenate(venous,0.9,plant.parameters,chemistry,0.05);
	Close(iga::ReadPerfusate(low,chemistry).pco2_mmhg,iga::ReadPerfusate(high,chemistry).pco2_mmhg,1e-9,
		"decoupled gas CO2");
	const auto blended=iga::Oxygenate(venous,0.3,plant.parameters,chemistry);
	assert(iga::ReadPerfusate(blended,chemistry).pco2_mmhg<iga::ReadPerfusate(low,chemistry).pco2_mmhg);
	assert(Rejects([&]{iga::Oxygenate(venous,0.3,plant.parameters,chemistry,0.5);}));
	// The dry sweep gas cannot hold more than all of itself.
	assert(Rejects([&]{iga::Oxygenate(venous,0.9,plant.parameters,chemistry,0.15);}));
	assert(!Rejects([&]{iga::Oxygenate(venous,0.9,plant.parameters,chemistry,0.1);}));
}

void TestSweepOxygenator()
{
	// A finite sweep gas flow: the exhaust gas closes both gas-side balances,
	// CO2 removal rises with sweep, and a very large sweep recovers the
	// infinite-sweep exchanger.
	const auto plant=MakePlantCase(1.,1.);
	const auto& p=plant.parameters;
	const auto chemistry=Perfusate();
	auto venous=Inlet();
	venous[iga::OrganOxygen]=1.5;
	venous[iga::OrganCarbonDioxide]+=3.;
	const double q=1e-5,fio2=0.21,dry=p.atmospheric_pressure_mmhg-p.water_vapour_mmhg;
	const auto molar=[&](double sweep){
		return sweep*p.atmospheric_pressure_mmhg*133.322387415/(8.314462618*p.sweep_gas_temperature_k);};
	double previous_pco2=1e9;
	for(const double sweep:{2e-6,8e-6,3e-5,1e-4}){
		const auto out=iga::OxygenateWithSweep(venous,fio2,sweep,q,p,chemistry);
		const auto inlet=iga::InletGasPartialPressures(fio2,p);
		const double co2_to_gas=molar(sweep)*(out.exhaust_pco2_mmhg-inlet.second)/dry;
		Close(co2_to_gas,q*(venous[iga::OrganCarbonDioxide]-out.arterial[iga::OrganCarbonDioxide]),
			1e-9,"sweep CO2 balance");
		Close(out.co2_removal_mol_s,co2_to_gas,1e-9,"reported CO2 removal");
		const double o2_from_gas=molar(sweep)*(inlet.first-out.exhaust_po2_mmhg)/dry;
		Close(o2_from_gas,q*(out.arterial[iga::OrganOxygen]-venous[iga::OrganOxygen]),
			1e-9,"sweep O2 balance");
		const double pco2=iga::ReadPerfusate(out.arterial,chemistry).pco2_mmhg;
		assert(pco2<previous_pco2);
		previous_pco2=pco2;
	}
	// Compared on a carbogen blend: with CO2-free gas the equilibrium content
	// is singular at pCO2 = 0.
	const auto large=iga::OxygenateWithSweep(venous,0.6,1e3,q,p,chemistry);
	const auto infinite=iga::Oxygenate(venous,0.6,p,chemistry);
	for(const auto s:{iga::OrganOxygen,iga::OrganCarbonDioxide})
		Close(large.arterial[s],infinite[s],1e-6,"large sweep recovers the equilibrium exchanger");
	assert(Rejects([&]{iga::OxygenateWithSweep(venous,fio2,0.,q,p,chemistry);}));
	// Through the plant: a lower sweep keeps more CO2 in the arterial perfusate.
	auto low=plant,high=plant;
	low.controls.fio2=high.controls.fio2=0.21;
	low.controls.sweep_gas_flow_m3_s=1e-7;
	high.controls.sweep_gas_flow_m3_s=1e-4;
	auto low_state=plant.state,high_state=plant.state;
	const auto low_records=iga::AdvancePerfusedOrganPlant(low.model,low.parameters,low.controls,{},{},
		low_state,600.,1,"organ");
	const auto high_records=iga::AdvancePerfusedOrganPlant(high.model,high.parameters,high.controls,{},{},
		high_state,600.,1,"organ");
	assert(low_records[0].arterial.pco2_mmhg>high_records[0].arterial.pco2_mmhg);
	assert(low_records[0].exhaust_pco2_mmhg>high_records[0].exhaust_pco2_mmhg);
}

void TestInjuryModifiers()
{
	// The anoxic cube with acid inhibition: the pH paradox term protects the
	// acidotic energy-depleted tissue, and hyperglycaemia amplifies its injury.
	auto plant=MakePlantCase(1.,1.);
	plant.parameters.kinetics.glycolysis_ph50=6.8;
	plant.parameters.kinetics.glycolysis_ph_hill=4.;
	const auto viability=[&](const iga::OrganPlantParameters& parameters){
		auto state=plant.state;
		return iga::AdvancePerfusedOrganPlant(plant.model,parameters,plant.controls,{},{},state,3600.,3,
			"organ").back().mean_viability;};
	const double base=viability(plant.parameters);
	auto acid=plant.parameters;
	acid.injury_ph50=7.2;
	assert(viability(acid)>base);
	auto sugar=plant.parameters;
	sugar.injury_glucose_threshold_mol_m3=0.;
	sugar.injury_glucose_gain_per_mol_m3=0.14;
	assert(viability(sugar)<base);
	for(const auto& change:std::vector<std::function<void(iga::OrganPlantParameters&)>>{
		[](iga::OrganPlantParameters& q){q.injury_ph50=10.;},
		[](iga::OrganPlantParameters& q){q.injury_ph50=7.2;q.injury_ph_hill=0.;},
		[](iga::OrganPlantParameters& q){q.injury_glucose_gain_per_mol_m3=-1.;},
		[](iga::OrganPlantParameters& q){q.sweep_gas_temperature_k=0.;}}){
		auto bad=plant.parameters;
		change(bad);
		assert(Rejects([&]{viability(bad);}));
	}
}

void TestRespiratoryInhibition()
{
	// Inhibited respiration lowers O2 uptake and lactate clearance by its factor
	// and leaves the ATP demand, so glycolysis and the deficit rise.
	const auto chemistry=Perfusate();
	auto k=Liver();
	const iga::HepaticCellState state{1.,100.};
	auto c=Inlet();
	c[iga::OrganOxygen]=iga::OxygenContent(500.,7.4,chemistry);
	const auto intact=iga::EvaluateHepaticRates(c,state,k,chemistry);
	k.respiration_capacity_scale=0.5;
	const auto inhibited=iga::EvaluateHepaticRates(c,state,k,chemistry);
	Close(inhibited.o2_uptake,0.5*intact.o2_uptake,1e-12,"inhibited O2 uptake");
	Close(inhibited.lactate_clearance,0.5*intact.lactate_clearance,1e-12,"inhibited lactate clearance");
	assert(inhibited.lactate_production>intact.lactate_production);
	k.respiration_capacity_scale=1.5;
	assert(Rejects([&]{iga::ValidateHepaticKinetics(k,chemistry);}));
	// Through the plant: a vasodilator bolus with a finite respiratory IC50 lowers O2 uptake.
	auto plant=MakePlantCase(30.,0.01);
	plant.controls.vasodilator_bolus_mg=20.;
	auto plain_state=plant.state,inhibited_state=plant.state;
	const auto plain=iga::AdvancePerfusedOrganPlant(plant.model,plant.parameters,plant.controls,{},{},
		plain_state,600.,1,"organ");
	auto parameters=plant.parameters;
	parameters.vasodilator_respiration_ic50_mg_l=3.4;
	const auto with_effect=iga::AdvancePerfusedOrganPlant(plant.model,parameters,plant.controls,{},{},
		inhibited_state,600.,1,"organ");
	assert(with_effect.back().o2_uptake_mol_s<plain.back().o2_uptake_mol_s);
	assert(with_effect.back().deficit_volume_fraction>=plain.back().deficit_volume_fraction);
	assert(plain.back().deficit_volume_fraction>=0.&&plain.back().deficit_volume_fraction<=1.);
	// A smaller Imax leaves more respiration at the same concentration.
	auto partial=parameters;
	partial.vasodilator_respiration_imax=0.5;
	auto partial_state=plant.state;
	const auto partial_records=iga::AdvancePerfusedOrganPlant(plant.model,partial,plant.controls,{},{},
		partial_state,600.,1,"organ");
	assert(partial_records.back().o2_uptake_mol_s>with_effect.back().o2_uptake_mol_s);
	partial.vasodilator_respiration_imax=1.5;
	auto bad_imax=plant.state;
	assert(Rejects([&]{iga::AdvancePerfusedOrganPlant(plant.model,partial,plant.controls,{},{},
		bad_imax,600.,1,"organ");}));
	parameters.vasodilator_respiration_ic50_mg_l=0.;
	auto rejected=plant.state;
	assert(Rejects([&]{iga::AdvancePerfusedOrganPlant(plant.model,parameters,plant.controls,{},{},
		rejected,600.,1,"organ");}));
}

} // namespace

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	TestFlowConservation();
	TestSeriesLimit();
	TestOxygenInversion();
	TestAcidBase();
	TestKineticsLimits();
	TestSingleCellOxygenAgainstScalarRoot();
	TestBlockedCellIsReportedAsFailure();
	TestGaussSeidelNonConvergenceIsReported();
	TestGlycolysisAcidInhibition();
	TestAcidInhibitionKeepsBufferBase();
	TestStagnantCellRelaxes();
	TestExchangeScaleMatchesConductance();
	TestLobularZonesAreTanksInSeries();
	TestCubeTransport();
	TestCircuitExtensions();
	TestPlantAdvance();
	TestAcidInhibitedPlantConserves();
	TestOxygenatedPlantConserves();
	TestPlantControlCosts();
	TestDecoupledGasCo2();
	TestSweepOxygenator();
	TestInjuryModifiers();
	TestRespiratoryInhibition();
	std::cout<<"perfused organ tests passed\n";
	PetscFinalize();
	return 0;
}
