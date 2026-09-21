#include "NativeTetAleGraphPorts.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

using Vector=std::array<double,3>;

iga::NativeTetMesh Mesh()
{
	iga::NativeTetMesh mesh;
	mesh.points={{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}}};
	mesh.cells={{1,{{0,1,2,3}}}};
	mesh.boundary_triangles={{1,{{0,2,3}},2},
		{2,{{1,2,3}},3},{3,{{0,1,3}},3},{4,{{0,1,2}},3}};
	return mesh;
}

iga::CouplingPort Port(const char* id,int label,iga::PortQuantity input)
{
	iga::CouplingPort port;
	port.id=id;port.subsystem_id="tet";
	port.locator_kind="boundary_label";port.locator=std::to_string(label);
	port.provides={iga::PortQuantity::Area,iga::PortQuantity::FlowRate,
		iga::PortQuantity::MeanPressure};
	port.requires={input};return port;
}

std::vector<double> State(const iga::NativeTetMesh& mesh,double ux,double pressure)
{
	const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
	const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
	std::vector<double> state(3*velocity_nodes+mesh.points.size(),0.);
	for(std::size_t node=0;node<velocity_nodes;++node)state[3*node]=ux;
	for(std::size_t node=0;node<mesh.points.size();++node)
		state[3*velocity_nodes+node]=pressure;
	return state;
}

void Close(double actual,double expected)
{
	assert(std::isfinite(actual));
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
	const iga::NativeTetAleGraphPorts bridge("tet",mesh,
		{Port("inlet",2,iga::PortQuantity::FlowRate),
			Port("outlet",3,iga::PortQuantity::MeanPressure)});
	const std::vector<Vector> zero(mesh.points.size(),Vector{{0,0,0}});
	iga::PortBoundaryData inlet;inlet.time_s=0.1;inlet.outward_flow_m3_s=-0.5;
	iga::PortBoundaryData outlet;outlet.time_s=0.1;outlet.mean_pressure_pa=7.;
	const auto fixed_bc=bridge.BoundaryConditions(mesh,zero,
		{{"inlet",inlet},{"outlet",outlet}},0.1);
	assert(fixed_bc.flow_rate_controls.size()==1);
	assert(fixed_bc.flow_rate_controls[0].boundary_label==2);
	Close(fixed_bc.flow_rate_controls[0].target_outward_flow_m3_s,-0.5);
	Close(fixed_bc.prescribed_pressure_pa.at(3),7.);
	const auto fixed=bridge.Observe(mesh,zero,State(mesh,1.,7.),0.1);
	Close(fixed.at("inlet").area_m2.value(),0.5);
	Close(fixed.at("inlet").outward_flow_m3_s.value(),-0.5);
	Close(fixed.at("outlet").outward_flow_m3_s.value(),0.5);
	Close(fixed.at("outlet").mean_pressure_pa.value(),7.);
	auto moved=mesh;
	for(auto& point:moved.points)point[0]+=0.1;
	const std::vector<Vector> grid(mesh.points.size(),Vector{{1,0,0}});
	const auto moved_bc=bridge.BoundaryConditions(moved,grid,
		{{"inlet",inlet},{"outlet",outlet}},0.1);
	Close(moved_bc.flow_rate_controls[0].target_outward_flow_m3_s,-1.);
	const auto moving=bridge.Observe(moved,grid,State(moved,2.,7.),0.1);
	Close(moving.at("inlet").outward_flow_m3_s.value(),-0.5);
	Close(moving.at("outlet").outward_flow_m3_s.value(),0.5);
	Close(moving.at("inlet").mean_pressure_pa.value(),7.);
	Reject([&]{bridge.BoundaryConditions(mesh,zero,{{"inlet",inlet}},0.1);});
	Reject([&]{bridge.BoundaryConditions(mesh,zero,
		{{"inlet",outlet},{"outlet",outlet}},0.1);});
	Reject([&]{bridge.BoundaryConditions(mesh,zero,
		{{"inlet",inlet},{"outlet",outlet}},0.2);});
	Reject([&]{bridge.Observe(mesh,zero,std::vector<double>{1.,2.},0.1);});
	std::cout<<"native tetra ALE graph hydraulic ports passed"
		<<" relative_inlet_m3_s="<<moving.at("inlet").outward_flow_m3_s.value()
		<<" absolute_controller_m3_s="
		<<moved_bc.flow_rate_controls[0].target_outward_flow_m3_s<<'\n';
}
