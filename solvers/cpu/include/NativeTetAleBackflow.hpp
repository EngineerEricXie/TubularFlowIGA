#ifndef IGA_NATIVE_TET_ALE_BACKFLOW_HPP
#define IGA_NATIVE_TET_ALE_BACKFLOW_HPP

#include "NativeTetFem.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <vector>

namespace iga {

struct NativeTetAleBackflowFaceSystem
{
	std::array<double,30> residual{};
	std::array<double,900> jacobian{};
	std::array<std::uint32_t,10> velocity_nodes{};
	double inward_relative_flux_m3_s = 0.0;
};

inline std::vector<std::size_t> NativeTetAleBoundaryFaceOwners(const NativeTetMesh& mesh)
{
	using Face=std::array<std::uint32_t,3>;
	std::map<Face,std::pair<std::size_t,unsigned>> uses;
	for(std::size_t cell=0;cell<mesh.cells.size();++cell)
		for(unsigned opposite=0;opposite<4;++opposite){
			Face face{};std::size_t entry=0;
			for(unsigned local=0;local<4;++local)
				if(local!=opposite)face[entry++]=mesh.cells[cell].nodes[local];
			std::sort(face.begin(),face.end());
			const auto found=uses.find(face);
			if(found==uses.end())uses.emplace(face,std::make_pair(cell,1));
			else ++found->second.second;
		}
	std::vector<std::size_t> result;
	result.reserve(mesh.boundary_triangles.size());
	for(const auto& triangle:mesh.boundary_triangles){
		auto face=triangle.nodes;std::sort(face.begin(),face.end());
		const auto found=uses.find(face);
		if(found==uses.end()||found->second.second!=1)
			throw std::runtime_error("native ALE backflow face has no unique tetra owner");
		result.push_back(found->second.first);
	}
	return result;
}

// Adds beta*rho*max(-(u-w).n,0)*u to the weak momentum residual. This is
// optional and changes the open-boundary model; it must be selected explicitly.
inline NativeTetAleBackflowFaceSystem BuildNativeTetAleBackflowFaceSystem(
	const NativeTetMesh& mesh,const NativeTaylorHoodTopology& topology,
	std::size_t cell_index,const NativeTetTriangle& triangle,
	const std::array<double,30>& local_velocity,
	const std::array<std::array<double,3>,4>& local_mesh_velocity,
	double density_kg_m3,double beta)
{
	if(cell_index>=mesh.cells.size()||!(density_kg_m3>0.0)
		||!std::isfinite(density_kg_m3)||!(beta>=0.0)||!std::isfinite(beta))
		throw std::invalid_argument("native ALE backflow input is invalid");
	const auto& cell=mesh.cells[cell_index];
	NativeTetAleBackflowFaceSystem result;
	result.velocity_nodes=topology.cell_velocity_nodes.at(cell_index);
	const auto& a=mesh.points.at(triangle.nodes[0]);
	const auto& b=mesh.points.at(triangle.nodes[1]);
	const auto& c=mesh.points.at(triangle.nodes[2]);
	std::array<double,3> first{},second{},area_vector{};
	for(int component=0;component<3;++component){
		first[component]=b[component]-a[component];
		second[component]=c[component]-a[component];
	}
	area_vector={{0.5*(first[1]*second[2]-first[2]*second[1]),
		0.5*(first[2]*second[0]-first[0]*second[2]),
		0.5*(first[0]*second[1]-first[1]*second[0])}};
	std::array<int,3> local_vertex{{-1,-1,-1}};
	for(std::size_t vertex=0;vertex<3;++vertex){
		for(std::size_t local=0;local<4;++local)
			if(cell.nodes[local]==triangle.nodes[vertex])local_vertex[vertex]=static_cast<int>(local);
		if(local_vertex[vertex]<0)
			throw std::invalid_argument("native ALE backflow triangle is absent from owner");
	}
	std::size_t opposite=0;
	while(opposite<4&&(cell.nodes[opposite]==triangle.nodes[0]
		||cell.nodes[opposite]==triangle.nodes[1]
		||cell.nodes[opposite]==triangle.nodes[2]))++opposite;
	if(opposite==4)throw std::invalid_argument("native ALE backflow owner has no opposite vertex");
	double inward=0.0;
	for(int component=0;component<3;++component)
		inward+=area_vector[component]*(mesh.points[cell.nodes[opposite]][component]-a[component]);
	if(inward>0.0)for(double& component:area_vector)component=-component;
	const double area=std::sqrt(area_vector[0]*area_vector[0]
		+area_vector[1]*area_vector[1]+area_vector[2]*area_vector[2]);
	if(!(area>0.0)||!std::isfinite(area))
		throw std::runtime_error("native ALE backflow face is degenerate");
	std::array<double,3> normal{};
	for(int component=0;component<3;++component)
		normal[component]=area_vector[component]/area;
	const std::array<std::array<double,3>,3> quadrature{{
		{{2.0/3.0,1.0/6.0,1.0/6.0}},{{1.0/6.0,2.0/3.0,1.0/6.0}},
		{{1.0/6.0,1.0/6.0,2.0/3.0}}}};
	for(const auto& face_lambda:quadrature){
		std::array<double,4> lambda{{0,0,0,0}};
		for(std::size_t vertex=0;vertex<3;++vertex)
			lambda[static_cast<std::size_t>(local_vertex[vertex])]=face_lambda[vertex];
		const auto basis=EvaluateNativeTaylorHoodBasis(lambda[1],lambda[2],lambda[3]);
		std::array<double,3> velocity{},grid_velocity{};
		for(std::size_t node=0;node<10;++node)
			for(int component=0;component<3;++component)
				velocity[component]+=basis.velocity[node]*local_velocity[3*node+component];
		for(std::size_t node=0;node<4;++node)
			for(int component=0;component<3;++component)
				grid_velocity[component]+=basis.pressure[node]*local_mesh_velocity[node][component];
		double relative_normal=0.0;
		for(int component=0;component<3;++component)
			relative_normal+=(velocity[component]-grid_velocity[component])*normal[component];
		const double weight=area/3.0;
		if(relative_normal>=0.0||beta==0.0)continue;
		result.inward_relative_flux_m3_s+=-relative_normal*weight;
		const double coefficient=beta*density_kg_m3*weight;
		for(std::size_t test=0;test<10;++test)
			for(int component=0;component<3;++component){
				const std::size_t row=3*test+component;
				result.residual[row]+=coefficient*basis.velocity[test]
					*(-relative_normal)*velocity[component];
				for(std::size_t unknown=0;unknown<10;++unknown)
					for(int derivative=0;derivative<3;++derivative){
						const std::size_t column=3*unknown+derivative;
						result.jacobian[30*row+column]+=coefficient*basis.velocity[test]
							*basis.velocity[unknown]*((component==derivative?-relative_normal:0.0)
								-normal[derivative]*velocity[component]);
					}
			}
	}
	return result;
}

} // namespace iga

#endif
