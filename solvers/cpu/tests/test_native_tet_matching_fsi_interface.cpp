#include "NativeTetMatchingFsiInterface.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
bool Close(double a,double b){return std::abs(a-b)<=1e-12*std::max({1.0,std::abs(a),std::abs(b)});}
template<class F>void Reject(F&&f){bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}assert(rejected);}
}

int main()
{
	iga::NativeTetMesh mesh;mesh.points={{{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}}}};
	mesh.cells={iga::NativeTetCell{1,{{0,1,2,3}}}};
	mesh.boundary_triangles={iga::NativeTetTriangle{10,{{0,1,2}},7},iga::NativeTetTriangle{11,{{0,3,1}},0},iga::NativeTetTriangle{12,{{0,2,3}},0},iga::NativeTetTriangle{13,{{1,3,2}},0}};
	const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
	std::vector<double> displacement(12),velocity(12);
	for(std::size_t node=0;node<4;++node){
		displacement[3*node]=0.1*mesh.points[node][0];displacement[3*node+1]=0.2*mesh.points[node][1];
		velocity[3*node]=1.0+mesh.points[node][0];velocity[3*node+1]=2.0*mesh.points[node][1];velocity[3*node+2]=-0.5;
	}
	const std::array<double,3> traction{{2,-1,3}};
	const auto transfer=iga::BuildNativeTetMatchingFsiTransfer(mesh,topology,7,displacement,velocity,{{10,traction}});
	const double area=0.5*1.1*1.2;
	assert(transfer.triangles==1&&Close(transfer.current_area_m2,area));
	for(int component=0;component<3;++component)assert(Close(transfer.total_force_n[component],area*traction[component]));
	const std::array<double,3> centroid{{1.1/3.0,1.2/3.0,0}};
	const auto expected_moment=iga::NativeTetInterfaceCross(centroid,transfer.total_force_n);
	for(int component=0;component<3;++component)assert(Close(transfer.total_moment_n_m[component],expected_moment[component]));
	assert(Close(transfer.nodal_power_w,transfer.quadrature_power_w));
	for(const auto node:{0u,1u,2u})for(int component=0;component<3;++component)
		assert(Close(transfer.solid_nodal_force_n[3*node+component],area*traction[component]/3.0));
	std::array<std::uint32_t,2> edge{{0,1}};const auto edge_position=std::lower_bound(topology.edges.begin(),topology.edges.end(),edge);
	const auto midpoint=static_cast<std::uint32_t>(mesh.points.size()+(edge_position-topology.edges.begin()));
	for(int component=0;component<3;++component)assert(Close(transfer.fluid_no_slip_velocity_m_s.at(midpoint)[component],0.5*(velocity[component]+velocity[3+component])));
	Reject([&]{(void)iga::BuildNativeTetMatchingFsiTransfer(mesh,topology,99,displacement,velocity,{});});
	Reject([&]{(void)iga::BuildNativeTetMatchingFsiTransfer(mesh,topology,7,displacement,velocity,{{11,traction}});});
	Reject([&]{(void)iga::BuildNativeTetMatchingFsiTransfer(mesh,topology,7,displacement,velocity,{{10,traction},{10,traction}});});
	const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();std::vector<double> pressure_state(3*velocity_nodes+mesh.points.size(),0.0);
	for(std::size_t node=0;node<mesh.points.size();++node)pressure_state[3*velocity_nodes+node]=2.0;
	const auto fluid_traction=iga::EvaluateNativeTetFluidTractionOnMatchingStructure(mesh,topology,pressure_state,0.5,7);
	assert(fluid_traction.size()==1&&fluid_traction[0].triangle_id==10);
	assert(Close(fluid_traction[0].traction_on_structure_pa[0],0.0));assert(Close(fluid_traction[0].traction_on_structure_pa[1],0.0));
	// z=0 has fluid outward normal -z: sigma*n=+2z, hence fluid-on-structure is -2z.
	assert(Close(fluid_traction[0].traction_on_structure_pa[2],-2.0));
	std::cout<<"native matching tetrahedral FSI interface test passed power_error="
		<<std::abs(transfer.nodal_power_w-transfer.quadrature_power_w)<<'\n';
	return 0;
}
