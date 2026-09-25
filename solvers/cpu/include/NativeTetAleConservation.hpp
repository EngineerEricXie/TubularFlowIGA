#ifndef IGA_NATIVE_TET_ALE_CONSERVATION_HPP
#define IGA_NATIVE_TET_ALE_CONSERVATION_HPP

#include "NativeTetAleKinematics.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <vector>

namespace iga {

struct NativeTetAleConservationDiagnostics
{
	double reference_volume_m3 = 0.0;
	double current_volume_m3 = 0.0;
	double volume_rate_m3_s = 0.0;
	double mesh_boundary_flux_m3_s = 0.0;
	double relative_boundary_flux_m3_s = 0.0;
	double gcl_residual_m3_s = 0.0;
	double moving_domain_balance_residual_m3_s = 0.0;
};

namespace native_tet_ale_conservation_detail {

inline std::array<std::uint32_t,3> SortedFace(std::array<std::uint32_t,3> face)
{
	std::sort(face.begin(),face.end());return face;
}

inline double Volume(const NativeTetMesh& mesh)
{
	double volume=0.0;
	for (const auto& cell : mesh.cells)
		volume += EvaluateNativeTetGeometry(mesh,cell).determinant/6.0;
	return volume;
}

inline double BoundaryFlux(const NativeTetMesh& mesh,
	const std::vector<std::array<double,3>>& velocity)
{
	if (velocity.size()!=mesh.points.size())
		throw std::invalid_argument("native ALE boundary-flux velocity size is invalid");
	struct Owner { std::uint32_t opposite=0; };
	std::map<std::array<std::uint32_t,3>,Owner> owners;
	for (const auto& cell : mesh.cells)
		for (std::size_t opposite=0;opposite<4;++opposite) {
			std::array<std::uint32_t,3> face{};std::size_t entry=0;
			for (std::size_t local=0;local<4;++local)
				if(local!=opposite) face[entry++]=cell.nodes[local];
			const auto key=SortedFace(face);
			const auto found=owners.find(key);
			if(found==owners.end()) owners.emplace(key,Owner{cell.nodes[opposite]});
			else owners.erase(found);
		}
	double flux=0.0;
	for (const auto& triangle : mesh.boundary_triangles) {
		const auto owner=owners.find(SortedFace(triangle.nodes));
		if(owner==owners.end())
			throw std::runtime_error("native ALE boundary triangle has no unique volume owner");
		const auto& p0=mesh.points[triangle.nodes[0]];
		const auto& p1=mesh.points[triangle.nodes[1]];
		const auto& p2=mesh.points[triangle.nodes[2]];
		const auto& interior=mesh.points[owner->second.opposite];
		std::array<double,3> first{},second{},area_vector{},average_velocity{};
		for(int component=0;component<3;++component) {
			first[component]=p1[component]-p0[component];
			second[component]=p2[component]-p0[component];
			average_velocity[component]=(velocity[triangle.nodes[0]][component]
				+velocity[triangle.nodes[1]][component]
				+velocity[triangle.nodes[2]][component])/3.0;
		}
		area_vector={{0.5*(first[1]*second[2]-first[2]*second[1]),
			0.5*(first[2]*second[0]-first[0]*second[2]),
			0.5*(first[0]*second[1]-first[1]*second[0])}};
		double inward=0.0;
		for(int component=0;component<3;++component)
			inward+=area_vector[component]*(interior[component]-p0[component]);
		if(inward>0.0) for(auto& component:area_vector) component=-component;
		for(int component=0;component<3;++component)
			flux+=average_velocity[component]*area_vector[component];
	}
	return flux;
}

inline NativeTetMesh InterpolateMesh(const NativeTetMesh& first,
	const NativeTetMesh& second,double alpha)
{
	if(first.points.size()!=second.points.size()||first.cells.size()!=second.cells.size()
		||first.boundary_triangles.size()!=second.boundary_triangles.size())
		throw std::invalid_argument("native ALE conservation meshes have different topology sizes");
	NativeTetMesh result=first;
	for(std::size_t node=0;node<result.points.size();++node)
		for(int component=0;component<3;++component)
			result.points[node][component]=(1.0-alpha)*first.points[node][component]
				+alpha*second.points[node][component];
	return result;
}

inline std::vector<std::array<double,3>> InterpolateVelocity(
	const std::vector<std::array<double,3>>& first,
	const std::vector<std::array<double,3>>& second,double alpha)
{
	if(first.size()!=second.size())
		throw std::invalid_argument("native ALE conservation velocity sizes differ");
	auto result=first;
	for(std::size_t node=0;node<result.size();++node)
		for(int component=0;component<3;++component)
			result[node][component]=(1.0-alpha)*first[node][component]
				+alpha*second[node][component];
	return result;
}

} // namespace native_tet_ale_conservation_detail

// Simpson time integration is exact for the quadratic boundary-flux history
// induced by linearly moving P1 tetrahedral nodes.
inline NativeTetAleConservationDiagnostics EvaluateNativeTetAleConservation(
	const NativeTetMesh& reference_mesh,const NativeTetMesh& current_mesh,
	const std::vector<std::array<double,3>>& mesh_velocity_m_s,
	const std::vector<std::array<double,3>>& fluid_velocity_reference_m_s,
	const std::vector<std::array<double,3>>& fluid_velocity_current_m_s,double dt_s)
{
	using namespace native_tet_ale_conservation_detail;
	if(!(dt_s>0.0)||!std::isfinite(dt_s)
		||mesh_velocity_m_s.size()!=reference_mesh.points.size())
		throw std::invalid_argument("native ALE conservation step is invalid");
	const auto midpoint_mesh=InterpolateMesh(reference_mesh,current_mesh,0.5);
	const auto midpoint_fluid=InterpolateVelocity(fluid_velocity_reference_m_s,
		fluid_velocity_current_m_s,0.5);
	auto relative_reference=fluid_velocity_reference_m_s;
	auto relative_midpoint=midpoint_fluid;
	auto relative_current=fluid_velocity_current_m_s;
	if(relative_reference.size()!=mesh_velocity_m_s.size()
		||relative_current.size()!=mesh_velocity_m_s.size())
		throw std::invalid_argument("native ALE conservation fluid velocity size is invalid");
	for(std::size_t node=0;node<mesh_velocity_m_s.size();++node)
		for(int component=0;component<3;++component) {
			relative_reference[node][component]-=mesh_velocity_m_s[node][component];
			relative_midpoint[node][component]-=mesh_velocity_m_s[node][component];
			relative_current[node][component]-=mesh_velocity_m_s[node][component];
		}
	NativeTetAleConservationDiagnostics result;
	result.reference_volume_m3=Volume(reference_mesh);
	result.current_volume_m3=Volume(current_mesh);
	result.volume_rate_m3_s=(result.current_volume_m3-result.reference_volume_m3)/dt_s;
	const double mesh_first=BoundaryFlux(reference_mesh,mesh_velocity_m_s);
	const double mesh_mid=BoundaryFlux(midpoint_mesh,mesh_velocity_m_s);
	const double mesh_last=BoundaryFlux(current_mesh,mesh_velocity_m_s);
	result.mesh_boundary_flux_m3_s=(mesh_first+4.0*mesh_mid+mesh_last)/6.0;
	const double relative_first=BoundaryFlux(reference_mesh,relative_reference);
	const double relative_mid=BoundaryFlux(midpoint_mesh,relative_midpoint);
	const double relative_last=BoundaryFlux(current_mesh,relative_current);
	result.relative_boundary_flux_m3_s=(relative_first+4.0*relative_mid+relative_last)/6.0;
	result.gcl_residual_m3_s=result.volume_rate_m3_s-result.mesh_boundary_flux_m3_s;
	result.moving_domain_balance_residual_m3_s=result.volume_rate_m3_s
		+result.relative_boundary_flux_m3_s;
	return result;
}

} // namespace iga

#endif
