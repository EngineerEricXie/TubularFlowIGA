#include "NativeTetMovingSpeciesPorts.hpp"
#include "SimulationGraph.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

namespace {

using Vector=std::array<double,3>;

iga::NativeTetMesh Mesh()
{
	iga::NativeTetMesh mesh;
	mesh.points={{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}},{{0,0,-1}}};
	mesh.cells={{1,{{0,1,2,3}}},{2,{{0,2,1,4}}}};
	using Face=std::array<std::uint32_t,3>;
	std::map<Face,unsigned> uses;
	for(const auto& cell:mesh.cells)
		for(std::size_t opposite=0;opposite<4;++opposite){
			Face face{};std::size_t entry=0;
			for(std::size_t local=0;local<4;++local)
				if(local!=opposite)face[entry++]=cell.nodes[local];
			std::sort(face.begin(),face.end());++uses[face];
		}
	std::uint64_t id=1;
	for(const auto& item:uses)
		if(item.second==1){
			const bool inlet=std::all_of(item.first.begin(),item.first.end(),
				[&](std::uint32_t node){return mesh.points[node][0]==0.;});
			mesh.boundary_triangles.push_back({id++,item.first,inlet?2:3});
		}
	return mesh;
}

iga::CouplingPort Port(const char* id,int label)
{
	iga::CouplingPort value;
	value.id=id;value.subsystem_id="tet";
	value.locator_kind="boundary_label";value.locator=std::to_string(label);
	value.provides={iga::PortQuantity::Area,iga::PortQuantity::FlowRate,
		iga::PortQuantity::MeanPressure,
		iga::PortQuantity::SpeciesConcentration,iga::PortQuantity::SpeciesFlux};
	value.requires={iga::PortQuantity::MeanPressure,
		iga::PortQuantity::SpeciesConcentration,iga::PortQuantity::SpeciesFlux};
	value.species={"tracer"};return value;
}

void Close(double actual,double expected)
{
	assert(std::abs(actual-expected)<=1e-11*std::max(1.,std::abs(expected)));
}

template<class Function> void Reject(Function action)
{
	bool rejected=false;
	try{action();}catch(const std::exception&){rejected=true;}
	assert(rejected);
}

} // namespace

int main()
{
	const auto mesh=Mesh();
	const std::vector<iga::CouplingPort> ports{Port("inlet",2),Port("outlet",3)};
	iga::CouplingPort peer;
	peer.id="root";peer.subsystem_id="one_d";
	peer.locator_kind="runtime_port";peer.locator="root";
	peer.provides={iga::PortQuantity::FlowRate,iga::PortQuantity::MeanPressure,
		iga::PortQuantity::SpeciesConcentration,iga::PortQuantity::SpeciesFlux};
	peer.requires={iga::PortQuantity::FlowRate,
		iga::PortQuantity::SpeciesConcentration,iga::PortQuantity::SpeciesFlux};
	peer.species={"tracer"};
	const iga::SimulationGraph graph(
		{{"tet",iga::DomainKind::ThreeDBodyFittedFlow,ports,{{"tracer","tracer"}}},
		{"one_d",iga::DomainKind::OneDFlow,{peer},{{"tracer","tracer"}}}},
		{{"exchange",{"tet","inlet"},{"one_d","root"},
			iga::CouplingLaw::PressureFlow,{"tracer"}}},
		{{"tracer","mol/m^3"}});
	assert(graph.Edges().size()==1);
	const auto plan=iga::MakeAcyclicPressureFlowPlan(graph,"tet");
	assert(plan.interfaces.size()==1);
	const iga::NativeTetMovingSpeciesPorts bridge("tet","tracer",mesh,ports);
	iga::PortBoundaryData inlet;inlet.time_s=0.1;inlet.concentration={{"tracer",2.}};
	const auto inflow=bridge.InflowByLabel({{"inlet",-1.},{"outlet",1.}},
		{{"inlet",inlet}},0.1);
	assert(inflow.size()==1&&inflow.at(2)==2.);
	Reject([&]{bridge.InflowByLabel({{"inlet",-1.},{"outlet",1.}},
		{},0.1);});
	Reject([&]{bridge.InflowByLabel({{"inlet",1.},{"outlet",-1.}},
		{{"inlet",inlet}},0.1);});
	Reject([&]{bridge.InflowByLabel({{"inlet",-1.},{"outlet",1.}},
		{{"unknown",inlet}},0.1);});
	const std::vector<Vector> zero(mesh.points.size(),Vector{{0,0,0}});
	const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
	const std::vector<Vector> flow(mesh.points.size()+topology.edges.size(),
		Vector{{1,0,0}});
	const std::vector<double> uniform(mesh.points.size(),2.);
	const auto step=iga::SolveNativeTetMovingSpeciesDenseStep(mesh,mesh,
		flow,zero,uniform,inflow,0.,0.,0.1);
	const auto observed=bridge.Observe(mesh,step,0.1);
	Close(*observed.at("inlet").outward_flow_m3_s,-1.);
	Close(*observed.at("outlet").outward_flow_m3_s,1.);
	Close(observed.at("inlet").concentration.at("tracer"),2.);
	Close(observed.at("outlet").concentration.at("tracer"),2.);
	Close(observed.at("inlet").outward_species_flux.at("tracer"),-2.);
	Close(observed.at("outlet").outward_species_flux.at("tracer"),2.);
	iga::PortState peer_out;
	peer_out.time_s=0.1;peer_out.outward_flow_m3_s=1.;
	peer_out.concentration={{"tracer",2.}};
	peer_out.outward_species_flux={{"tracer",2.}};
	const auto edge_residual=iga::ConservativeEdgeResidual(
		observed.at("inlet"),peer_out);
	Close(edge_residual.outward_flow_m3_s,0.);
	Close(edge_residual.outward_species_flux.at("tracer"),0.);
	const auto accounting=bridge.Accounting(step,0.1);
	Close(accounting.initial_mass,2./3.);
	Close(accounting.final_mass,2./3.);
	Close(accounting.outward_port_amount.at("inlet"),-0.2);
	Close(accounting.outward_port_amount.at("outlet"),0.2);
	Close(accounting.residual,0.);
	iga::PortBoundaryData reverse_input;
	reverse_input.time_s=0.1;reverse_input.concentration={{"tracer",2.}};
	const auto reverse_inflow=bridge.InflowByLabel(
		{{"inlet",1.},{"outlet",-1.}},{{"outlet",reverse_input}},0.1);
	assert(reverse_inflow.size()==1&&reverse_inflow.at(3)==2.);
	const std::vector<Vector> reverse_flow(flow.size(),Vector{{-1,0,0}});
	const auto reverse_step=iga::SolveNativeTetMovingSpeciesDenseStep(mesh,mesh,
		reverse_flow,zero,uniform,reverse_inflow,0.,0.,0.1);
	const auto reverse_observed=bridge.Observe(mesh,reverse_step,0.1);
	Close(*reverse_observed.at("inlet").outward_flow_m3_s,1.);
	Close(*reverse_observed.at("outlet").outward_flow_m3_s,-1.);
	Close(reverse_observed.at("inlet").outward_species_flux.at("tracer"),2.);
	Close(reverse_observed.at("outlet").outward_species_flux.at("tracer"),-2.);
	const auto reverse_edge_residual=iga::ConservativeEdgeResidual(
		reverse_observed.at("outlet"),peer_out);
	Close(reverse_edge_residual.outward_flow_m3_s,0.);
	Close(reverse_edge_residual.outward_species_flux.at("tracer"),0.);
	const auto reverse_accounting=bridge.Accounting(reverse_step,0.1);
	Close(reverse_accounting.residual,0.);
	const std::vector<Vector> stopped(flow.size(),Vector{{0,0,0}});
	const auto stopped_step=iga::SolveNativeTetMovingSpeciesDenseStep(mesh,mesh,
		stopped,zero,uniform,{},0.,0.,0.1);
	const auto stopped_ports=bridge.Observe(mesh,stopped_step,0.1);
	assert(bridge.InflowByLabel({{"inlet",0.},{"outlet",0.}},
		{{"inlet",inlet}},0.1).empty());
	for(const auto& item:stopped_ports){
		Close(*item.second.outward_flow_m3_s,0.);
		Close(item.second.outward_species_flux.at("tracer"),0.);
		Close(item.second.concentration.at("tracer"),2.);
	}
	const std::vector<double> varying{1.,2.,3.,4.,5.};
	const auto stopped_varying=iga::SolveNativeTetMovingSpeciesDenseStep(mesh,mesh,
		stopped,zero,varying,{},0.,0.,0.1);
	const auto varying_ports=bridge.Observe(mesh,stopped_varying,0.1);
	Close(varying_ports.at("inlet").concentration.at("tracer"),17./6.);
	auto mixed=step;
	mixed.negative_relative_flow_by_label_m3_s[3]=-0.1;
	Reject([&]{bridge.Observe(mesh,mixed,0.1);});
	auto negative=step;negative.concentration_mol_m3[0]=-0.1;
	Reject([&]{bridge.Observe(mesh,negative,0.1);});
	auto unported=step;
	unported.outward_advective_flux_by_label_mol_s[4]=0.1;
	Reject([&]{bridge.Accounting(unported,0.1);});
	auto duplicate=ports;duplicate[1].locator="2";
	Reject([&]{iga::NativeTetMovingSpeciesPorts("tet","tracer",mesh,duplicate);});
	auto extra_species=ports;extra_species[0].species.insert("other");
	Reject([&]{iga::NativeTetMovingSpeciesPorts("tet","tracer",mesh,extra_species);});
	std::cout<<"native moving species graph-port contract passed"
		<<" inlet_mol="<<accounting.outward_port_amount.at("inlet")
		<<" outlet_mol="<<accounting.outward_port_amount.at("outlet")<<'\n';
}
