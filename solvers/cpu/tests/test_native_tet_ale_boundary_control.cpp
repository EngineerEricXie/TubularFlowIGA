#include "NativeTetAleBoundaryControl.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

bool Rejects(const iga::NativeTetAleBoundaryConditions& conditions,
	const std::map<int,iga::NativeTetAleBoundaryFluxOperator>& operators)
{
	try { iga::ValidateNativeTetAleBoundaryConditions(conditions,operators); }
	catch(const std::invalid_argument&) { return true; }
	return false;
}

}

int main()
{
	iga::NativeTetMesh mesh;
	mesh.points={{{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}}}};
	mesh.cells={{1,{{0,1,2,3}}}};
	mesh.boundary_triangles={iga::NativeTetTriangle{1,{{1,2,3}},7},
		iga::NativeTetTriangle{2,{{0,3,2}},8},
		iga::NativeTetTriangle{3,{{0,1,3}},8},
		iga::NativeTetTriangle{4,{{0,2,1}},8}};
	const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
	const auto operators=iga::BuildNativeTetAleBoundaryFluxOperators(mesh,topology);
	assert(operators.size()==2);
	const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
	const std::array<double,3> velocity{{0.2,-0.1,0.04}};
	std::vector<double> state(3*velocity_nodes+mesh.points.size(),0.0);
	for(std::size_t node=0;node<velocity_nodes;++node)
		for(int component=0;component<3;++component)
			state[3*node+component]=velocity[component];
	for(const auto& labelled:operators){
		double independent=0.0;
		for(const auto& triangle:mesh.boundary_triangles)
			if(triangle.boundary_label==labelled.first){
				const auto& a=mesh.points[triangle.nodes[0]];
				const auto& b=mesh.points[triangle.nodes[1]];
				const auto& c=mesh.points[triangle.nodes[2]];
				std::array<double,3> first{},second{},area_vector{};
				for(int component=0;component<3;++component){
					first[component]=b[component]-a[component];
					second[component]=c[component]-a[component];
				}
				area_vector={{0.5*(first[1]*second[2]-first[2]*second[1]),
					0.5*(first[2]*second[0]-first[0]*second[2]),
					0.5*(first[0]*second[1]-first[1]*second[0])}};
				const auto& opposite=mesh.points[labelled.first==7?0:1];
				double inward=0.0;
				for(int component=0;component<3;++component)
					inward+=area_vector[component]*(opposite[component]-a[component]);
				if(inward>0.0)for(double& value:area_vector)value=-value;
				for(int component=0;component<3;++component)
					independent+=velocity[component]*area_vector[component];
			}
		assert(std::abs(labelled.second.OutwardFlowM3S(state)-independent)<1e-15);
	}
	iga::NativeTetAleBoundaryConditions valid;
	valid.prescribed_pressure_pa[8]=0.0;valid.flow_rate_controls.push_back({7,1e-6});
	iga::ValidateNativeTetAleBoundaryConditions(valid,operators);
	auto duplicate=valid;duplicate.flow_rate_controls.push_back({7,2e-6});
	assert(Rejects(duplicate,operators));
	auto overlap=valid;overlap.prescribed_pressure_pa[7]=2.0;
	assert(Rejects(overlap,operators));
	auto absent=valid;absent.prescribed_pressure_pa[9]=0.0;
	assert(Rejects(absent,operators));
	auto nonfinite=valid;nonfinite.flow_rate_controls[0].target_outward_flow_m3_s=
		std::numeric_limits<double>::quiet_NaN();
	assert(Rejects(nonfinite,operators));
	auto stabilized=valid;stabilized.backflow_stabilization_beta[8]=0.5;
	iga::ValidateNativeTetAleBoundaryConditions(stabilized,operators);
	auto wrong_port=valid;wrong_port.backflow_stabilization_beta[7]=0.5;
	assert(Rejects(wrong_port,operators));
	auto negative_beta=valid;negative_beta.backflow_stabilization_beta[8]=-0.1;
	assert(Rejects(negative_beta,operators));
	std::cout<<"native tetrahedral ALE pressure/flow boundary operator tests passed\n";
}
