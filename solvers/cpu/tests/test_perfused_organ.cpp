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
	}
}

void TestPlantAdvance()
{
	iga::PerfusedOrganModel model;
	model.mesh=Cube(3);
	model.trees=CubeTrees(model.mesh,3e-3);
	model.zones=3;
	iga::FinalizePerfusedOrganModel(model);
	const std::size_t cells=model.mesh.cells.size(),nodes=3*cells;
	iga::OrganPlantParameters p;
	p.viscosity_pa_s=3e-3;
	p.outflow_pressure_pa={{"outflow",0.}};
	p.specific_conductance_inv_pa_s={{"inflow",2e-8},{"outflow",6e-8}};
	p.mobility_m2_pa_s=1e-8;
	p.distribution_volume_fraction={{0.2,0.7,0.7,0.3,0.5,0.25}};
	p.kinetics=Liver();
	p.chemistry=Perfusate();
	iga::OrganPlantState state;
	state.reservoir.volume_m3=2.0;
	state.reservoir.hematocrit_percent=18.;
	const auto& names=iga::PerfusedOrganSpeciesNames();
	const auto inlet=Inlet();
	for(std::size_t s=0;s<iga::OrganSpeciesCount;++s)state.reservoir.species[names[s]]=inlet[s];
	state.tissue.assign(nodes,inlet);
	state.viability.assign(nodes,1.);
	state.glycogen_mol_m3.assign(nodes,100.);
	state.oedema.assign(cells,0.);
	iga::OrganPlantControls controls;
	controls.inlet_pressure_pa=1000.;
	controls.fio2=0.6;
	controls.bolus={{names[iga::OrganGlucose],5.0}};
	controls.dialysis_in_m3_s=1e-8;
	for(const std::string resolution:{"cell","portal_territory","organ"}){
		auto copy=state;
		const auto records=iga::AdvancePerfusedOrganPlant(model,p,controls,{},{},copy,3600.,3,resolution);
		assert(records.size()==3);
		Close(copy.time_s,3600.,1e-12,"plant time");
		Close(copy.reservoir.volume_m3,2.0+1e-8*3600.,1e-9,"plant volume");
		for(const auto& record:records){
			assert(record.flow_m3_s>0.&&record.arterial.po2_mmhg>record.venous.po2_mmhg);
			assert(record.arterial.pco2_mmhg>0.&&record.mean_viability<=1.);
		}
		assert(copy.reservoir.species.at(names[iga::OrganGlucose])>inlet[iga::OrganGlucose]);
		for(const double v:copy.viability)assert(v>0.&&v<=1.);
	}
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
	TestStagnantCellRelaxes();
	TestExchangeScaleMatchesConductance();
	TestLobularZonesAreTanksInSeries();
	TestCubeTransport();
	TestCircuitExtensions();
	TestPlantAdvance();
	std::cout<<"perfused organ tests passed\n";
	PetscFinalize();
	return 0;
}
