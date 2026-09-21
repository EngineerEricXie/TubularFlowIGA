#include "NativeTetVelocityField.hpp"

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

int main()
{
	iga::NativeTetMesh mesh;
	mesh.points={{{0.,0.,0.}},{{1.,0.,0.}},{{0.,1.,0.}},{{0.,0.,1.}}};
	mesh.cells.push_back({1,{{0,1,2,3}}});
	const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
	const auto velocity_nodes=mesh.points.size()+topology.edges.size();
	std::vector<double> state(3*velocity_nodes+mesh.points.size(),0.0);
	for(std::size_t node=0;node<velocity_nodes;++node){
		const auto point=node<mesh.points.size()?mesh.points[node]:
			std::array<double,3>{{
				.5*(mesh.points[topology.edges[node-mesh.points.size()][0]][0]
					+mesh.points[topology.edges[node-mesh.points.size()][1]][0]),
				.5*(mesh.points[topology.edges[node-mesh.points.size()][0]][1]
					+mesh.points[topology.edges[node-mesh.points.size()][1]][1]),
				.5*(mesh.points[topology.edges[node-mesh.points.size()][0]][2]
					+mesh.points[topology.edges[node-mesh.points.size()][1]][2])}};
		state[3*node]=2.*point[0]-point[1]+.5;
		state[3*node+1]=point[1]+3.*point[2];
		state[3*node+2]=-point[0]+.25*point[2];
	}
	const iga::NativeTetVelocityField field(mesh,topology,state);
	for(const auto& point:std::vector<std::array<double,3>>{
		{{.1,.2,.3}},{{.25,.25,.25}},{{0.,0.,0.}},{{.5,.5,0.}}}){
		const auto actual=field.At(point);
		if(!actual||std::abs((*actual)[0]-(2.*point[0]-point[1]+.5))>1e-12
			||std::abs((*actual)[1]-(point[1]+3.*point[2]))>1e-12
			||std::abs((*actual)[2]-(-point[0]+.25*point[2]))>1e-12)
			throw std::runtime_error("native tetrahedral P2 field failed affine reproduction");
	}
	if(field.At({{.5,.5,.5}}))
		throw std::runtime_error("native tetrahedral P2 field extrapolated outside its mesh");
	return 0;
}
