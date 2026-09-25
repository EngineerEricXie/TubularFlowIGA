#ifndef IGA_NATIVE_TET_MATCHING_FSI_INTERFACE_HPP
#define IGA_NATIVE_TET_MATCHING_FSI_INTERFACE_HPP

#include "NativeTetFem.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

namespace iga {

struct NativeTetInterfaceTriangleTraction
{
	std::uint64_t triangle_id=0;
	std::array<double,3> traction_on_structure_pa{};
};

struct NativeTetMatchingFsiTransfer
{
	NativeTetMesh current_mesh;
	std::vector<std::array<double,3>> ale_mesh_velocity_m_s;
	std::map<std::uint32_t,std::array<double,3>> fluid_no_slip_velocity_m_s;
	std::vector<double> solid_nodal_force_n;
	std::array<double,3> total_force_n{};
	std::array<double,3> total_moment_n_m{};
	double nodal_power_w=0.0;
	double quadrature_power_w=0.0;
	double current_area_m2=0.0;
	std::size_t triangles=0;
};

inline std::array<double,3> NativeTetInterfaceCross(const std::array<double,3>& a,
	const std::array<double,3>& b);
inline double NativeTetInterfaceNorm(const std::array<double,3>& a);

// Area-average the native fluid Cauchy traction over each selected current
// triangle and flip sigma*n_fluid to the force exerted on the structure.
inline std::vector<NativeTetInterfaceTriangleTraction>
EvaluateNativeTetFluidTractionOnMatchingStructure(const NativeTetMesh& current_mesh,
	const NativeTaylorHoodTopology& topology,const std::vector<double>& fluid_state,
	double dynamic_viscosity_pa_s,int interface_label)
{
	if(!(dynamic_viscosity_pa_s>0.0)||!std::isfinite(dynamic_viscosity_pa_s))
		throw std::invalid_argument("native matching FSI fluid viscosity is invalid");
	const std::size_t velocity_nodes=current_mesh.points.size()+topology.edges.size();
	if(fluid_state.size()!=3*velocity_nodes+current_mesh.points.size())
		throw std::invalid_argument("native matching FSI fluid state size mismatch");
	using Face=std::array<std::uint32_t,3>;
	struct Owner{std::size_t cell=0,opposite=0;};std::map<Face,Owner> owners;
	for(std::size_t cell=0;cell<current_mesh.cells.size();++cell)for(std::size_t opposite=0;opposite<4;++opposite){
		Face face{};std::size_t entry=0;for(std::size_t local=0;local<4;++local)if(local!=opposite)face[entry++]=current_mesh.cells[cell].nodes[local];
		std::sort(face.begin(),face.end());const auto found=owners.find(face);if(found==owners.end())owners.emplace(face,Owner{cell,opposite});else owners.erase(found);
	}
	const std::array<std::array<double,3>,3> triangle_points{{{{2.0/3,1.0/6,1.0/6}},{{1.0/6,2.0/3,1.0/6}},{{1.0/6,1.0/6,2.0/3}}}};
	std::vector<NativeTetInterfaceTriangleTraction> result;
	for(const auto& triangle:current_mesh.boundary_triangles){
		if(triangle.boundary_label!=interface_label)continue;
		auto face=triangle.nodes;std::sort(face.begin(),face.end());const auto owner=owners.find(face);
		if(owner==owners.end())throw std::runtime_error("native matching FSI fluid triangle has no unique owner");
		const auto cell=owner->second.cell;const auto geometry=EvaluateNativeTetGeometry(current_mesh,current_mesh.cells[cell]);
		const auto&p0=current_mesh.points[triangle.nodes[0]];const auto&p1=current_mesh.points[triangle.nodes[1]];const auto&p2=current_mesh.points[triangle.nodes[2]];
		std::array<double,3>a{},b{},area_vector{},normal{};for(int c=0;c<3;++c){a[c]=p1[c]-p0[c];b[c]=p2[c]-p0[c];}
		area_vector={{0.5*(a[1]*b[2]-a[2]*b[1]),0.5*(a[2]*b[0]-a[0]*b[2]),0.5*(a[0]*b[1]-a[1]*b[0])}};
		const auto& interior=current_mesh.points[current_mesh.cells[cell].nodes[owner->second.opposite]];double inward=0.0;for(int c=0;c<3;++c)inward+=area_vector[c]*(interior[c]-p0[c]);if(inward>0.0)for(auto&v:area_vector)v=-v;
		const double area=NativeTetInterfaceNorm(area_vector);if(!(area>0.0))throw std::runtime_error("native matching FSI fluid triangle area is nonpositive");for(int c=0;c<3;++c)normal[c]=area_vector[c]/area;
		std::array<int,3> local_vertex{};for(std::size_t vertex=0;vertex<3;++vertex){local_vertex[vertex]=-1;for(std::size_t local=0;local<4;++local)if(current_mesh.cells[cell].nodes[local]==triangle.nodes[vertex])local_vertex[vertex]=static_cast<int>(local);if(local_vertex[vertex]<0)throw std::runtime_error("native matching FSI fluid face is absent from owner");}
		NativeTetInterfaceTriangleTraction value;value.triangle_id=triangle.id;
		for(const auto& face_lambda:triangle_points){
			std::array<double,4> lambda{{0,0,0,0}};for(std::size_t vertex=0;vertex<3;++vertex)lambda[static_cast<std::size_t>(local_vertex[vertex])]=face_lambda[vertex];
			const auto basis=EvaluateNativeTaylorHoodPhysicalBasis(geometry,{{lambda[1],lambda[2],lambda[3]}});
			std::array<std::array<double,3>,3> gradient{};double pressure=0.0;
			for(std::size_t local=0;local<10;++local){const auto global=topology.cell_velocity_nodes[cell][local];for(int component=0;component<3;++component)for(int derivative=0;derivative<3;++derivative)gradient[component][derivative]+=basis.velocity_gradients[local][derivative]*fluid_state[3*global+component];}
			for(std::size_t local=0;local<4;++local)pressure+=basis.pressure[local]*fluid_state[3*velocity_nodes+current_mesh.cells[cell].nodes[local]];
			for(int row=0;row<3;++row)for(int column=0;column<3;++column)value.traction_on_structure_pa[row]-=((row==column?-pressure:0.0)+dynamic_viscosity_pa_s*(gradient[row][column]+gradient[column][row]))*normal[column]/3.0;
		}
		result.push_back(value);
	}
	if(result.empty())throw std::invalid_argument("native matching FSI fluid label has no interface triangles");
	return result;
}

inline std::array<double,3> NativeTetInterfaceCross(const std::array<double,3>& a,
	const std::array<double,3>& b)
{
	return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}};
}

inline double NativeTetInterfaceNorm(const std::array<double,3>& a)
{
	return std::sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

inline NativeTetMatchingFsiTransfer BuildNativeTetMatchingFsiTransfer(
	const NativeTetMesh& reference_mesh,const NativeTaylorHoodTopology& fluid_topology,
	int interface_label,const std::vector<double>& solid_displacement_m,
	const std::vector<double>& solid_velocity_m_s,
	const std::vector<NativeTetInterfaceTriangleTraction>& triangle_tractions)
{
	const std::size_t nodes=reference_mesh.points.size();
	if(solid_displacement_m.size()!=3*nodes||solid_velocity_m_s.size()!=3*nodes)
		throw std::invalid_argument("native matching FSI solid state size mismatch");
	for(double value:solid_displacement_m)if(!std::isfinite(value))
		throw std::invalid_argument("native matching FSI displacement is nonfinite");
	for(double value:solid_velocity_m_s)if(!std::isfinite(value))
		throw std::invalid_argument("native matching FSI velocity is nonfinite");
	std::map<std::uint64_t,std::array<double,3>> traction;
	for(const auto& item:triangle_tractions){
		for(double value:item.traction_on_structure_pa)if(!std::isfinite(value))
			throw std::invalid_argument("native matching FSI traction is nonfinite");
		if(!traction.emplace(item.triangle_id,item.traction_on_structure_pa).second)
			throw std::invalid_argument("native matching FSI has duplicate triangle traction");
	}
	NativeTetMatchingFsiTransfer result;result.current_mesh=reference_mesh;
	result.ale_mesh_velocity_m_s.assign(nodes,{{0,0,0}});
	result.solid_nodal_force_n.assign(3*nodes,0.0);
	std::set<std::uint64_t> consumed;
	for(const auto& triangle:reference_mesh.boundary_triangles){
		if(triangle.boundary_label!=interface_label)continue;
		++result.triangles;
		const auto found=traction.find(triangle.id);
		if(found==traction.end())throw std::invalid_argument("native matching FSI interface triangle lacks traction");
		consumed.insert(triangle.id);
		std::array<std::array<double,3>,3> current{},velocity{};
		for(std::size_t a=0;a<3;++a){
			const auto node=triangle.nodes[a];if(node>=nodes)throw std::out_of_range("native matching FSI triangle node is invalid");
			for(int component=0;component<3;++component){
				current[a][component]=reference_mesh.points[node][component]+solid_displacement_m[3*node+component];
				velocity[a][component]=solid_velocity_m_s[3*node+component];
				result.current_mesh.points[node][component]=current[a][component];
				result.ale_mesh_velocity_m_s[node][component]=velocity[a][component];
			}
			result.fluid_no_slip_velocity_m_s[node]=velocity[a];
		}
		std::array<double,3> edge0{},edge1{},centroid{},centroid_velocity{};
		for(int component=0;component<3;++component){
			edge0[component]=current[1][component]-current[0][component];
			edge1[component]=current[2][component]-current[0][component];
			centroid[component]=(current[0][component]+current[1][component]+current[2][component])/3.0;
			centroid_velocity[component]=(velocity[0][component]+velocity[1][component]+velocity[2][component])/3.0;
		}
		const double area=0.5*NativeTetInterfaceNorm(NativeTetInterfaceCross(edge0,edge1));
		if(!(area>0.0)||!std::isfinite(area))throw std::runtime_error("native matching FSI current triangle area is nonpositive");
		result.current_area_m2+=area;
		std::array<double,3> triangle_force{};
		for(int component=0;component<3;++component){
			triangle_force[component]=area*found->second[component];result.total_force_n[component]+=triangle_force[component];
			for(std::size_t a=0;a<3;++a)result.solid_nodal_force_n[3*triangle.nodes[a]+component]+=triangle_force[component]/3.0;
		}
		const auto moment=NativeTetInterfaceCross(centroid,triangle_force);
		for(int component=0;component<3;++component)result.total_moment_n_m[component]+=moment[component];
		for(int component=0;component<3;++component)result.quadrature_power_w+=triangle_force[component]*centroid_velocity[component];
		for(const auto pair:{std::array<int,2>{{0,1}},std::array<int,2>{{0,2}},std::array<int,2>{{1,2}}}){
			std::array<std::uint32_t,2> edge{{triangle.nodes[pair[0]],triangle.nodes[pair[1]]}};std::sort(edge.begin(),edge.end());
			// Topology edge storage follows first cell encounter so it is stable for
			// DOF numbering, but is not globally lexicographically sorted on a
			// multi-tetrahedron mesh.
			const auto edge_position=std::find(fluid_topology.edges.begin(),fluid_topology.edges.end(),edge);
			if(edge_position==fluid_topology.edges.end()||*edge_position!=edge)
				throw std::runtime_error("native matching FSI P2 boundary edge is absent from fluid topology");
			const auto velocity_node=static_cast<std::uint32_t>(nodes+(edge_position-fluid_topology.edges.begin()));
			std::array<double,3> midpoint{};for(int component=0;component<3;++component)midpoint[component]=0.5*(solid_velocity_m_s[3*edge[0]+component]+solid_velocity_m_s[3*edge[1]+component]);
			result.fluid_no_slip_velocity_m_s[velocity_node]=midpoint;
		}
	}
	if(result.triangles==0)throw std::invalid_argument("native matching FSI label has no interface triangles");
	if(consumed.size()!=traction.size())throw std::invalid_argument("native matching FSI received extraneous triangle traction");
	for(std::size_t node=0;node<nodes;++node)for(int component=0;component<3;++component)
		result.nodal_power_w+=result.solid_nodal_force_n[3*node+component]*solid_velocity_m_s[3*node+component];
	return result;
}

} // namespace iga

#endif
