#include "NativeTetAleConservation.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

namespace {

iga::NativeTetMesh ReferenceMesh()
{
	iga::NativeTetMesh mesh;
	mesh.points={{{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}}}};
	mesh.cells={iga::NativeTetCell{1,{{0,1,2,3}}}};
	mesh.boundary_triangles={iga::NativeTetTriangle{1,{{1,2,3}},0},
		iga::NativeTetTriangle{2,{{0,3,2}},0},iga::NativeTetTriangle{3,{{0,1,3}},0},
		iga::NativeTetTriangle{4,{{0,2,1}},0}};
	return mesh;
}

} // namespace

int main()
{
	const auto reference=ReferenceMesh();
	const double dt=0.2;
	auto translated=reference;
	std::vector<std::array<double,3>> translation(reference.points.size(),{{0.3,-0.2,0.1}});
	for(std::size_t node=0;node<translated.points.size();++node)
		for(int component=0;component<3;++component)
			translated.points[node][component]+=dt*translation[node][component];
	const auto rigid=iga::EvaluateNativeTetAleConservation(reference,translated,translation,
		translation,translation,dt);
	assert(std::abs(rigid.gcl_residual_m3_s)<1e-14);
	assert(std::abs(rigid.relative_boundary_flux_m3_s)<1e-14);
	assert(std::abs(rigid.moving_domain_balance_residual_m3_s)<1e-14);

	auto expanded=reference;
	std::vector<std::array<double,3>> expansion(reference.points.size());
	for(std::size_t node=0;node<expanded.points.size();++node)
		for(int component=0;component<3;++component) {
			expansion[node][component]=0.5*reference.points[node][component];
			expanded.points[node][component]+=dt*expansion[node][component];
		}
	const std::vector<std::array<double,3>> zero(reference.points.size(),{{0,0,0}});
	const auto dilation=iga::EvaluateNativeTetAleConservation(reference,expanded,expansion,
		zero,zero,dt);
	assert(std::abs(dilation.gcl_residual_m3_s)<1e-14);
	assert(std::abs(dilation.moving_domain_balance_residual_m3_s)<1e-14);
	assert(std::abs(dilation.current_volume_m3-std::pow(1.1,3)/6.0)<1e-14);
	auto contracted=reference;
	std::vector<std::array<double,3>> contraction(reference.points.size());
	for(std::size_t node=0;node<contracted.points.size();++node)
		for(int component=0;component<3;++component) {
			contraction[node][component]=-0.5*reference.points[node][component];
			contracted.points[node][component]+=dt*contraction[node][component];
		}
	const auto shrinkage=iga::EvaluateNativeTetAleConservation(reference,contracted,contraction,
		zero,zero,dt);
	assert(std::abs(shrinkage.gcl_residual_m3_s)<1e-14);
	assert(std::abs(shrinkage.moving_domain_balance_residual_m3_s)<1e-14);
	assert(std::abs(shrinkage.current_volume_m3-std::pow(0.9,3)/6.0)<1e-14);
	std::cout<<"native tetrahedral ALE GCL and moving-domain balance tests passed\n";
	return 0;
}
