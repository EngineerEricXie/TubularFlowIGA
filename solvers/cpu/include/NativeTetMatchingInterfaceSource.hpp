#ifndef IGA_NATIVE_TET_MATCHING_INTERFACE_SOURCE_HPP
#define IGA_NATIVE_TET_MATCHING_INTERFACE_SOURCE_HPP

#include "NativeTetVesselTissueSourceMap.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

// This is deliberately an exact, matching-facet contract. A nonconforming
// interface needs its own conservative intersection/quadrature algorithm.
struct NativeTetMatchingInterfacePort
{
	std::string name;
	int vessel_boundary_label=0;
	int tissue_boundary_label=0;
};

struct NativeTetMatchingInterfaceSourceResult
{
	NativeTetVesselTissueSourceResult source;
	std::map<std::string,std::size_t> matched_facets;
	struct FacetTransfer
	{
		std::string port_name;
		std::uint64_t vessel_cell_id=0;
		std::uint64_t tissue_cell_id=0;
		std::array<std::array<double,3>,3> triangle_m{};
		double vessel_outward_flow_m3_s=0.;
	};
	std::vector<FacetTransfer> facet_transfers;
};

namespace native_tet_matching_detail {

using Face=std::array<std::uint32_t,3>;
using Point=std::array<double,3>;
using GeometricFace=std::array<Point,3>;

struct Owner
{
	std::size_t cell=0;
	std::size_t opposite=0;
};

struct Facet
{
	std::size_t cell=0;
	Face nodes{};
	Point outward_area_vector{};
};

inline std::map<Face,Owner> BoundaryOwners(const NativeTetMesh& mesh)
{
	std::map<Face,std::vector<Owner>> uses;
	for(std::size_t cell=0;cell<mesh.cells.size();++cell) {
		EvaluateNativeTetGeometry(mesh,mesh.cells[cell]);
		for(std::size_t opposite=0;opposite<4;++opposite) {
			Face face{};std::size_t entry=0;
			for(std::size_t local=0;local<4;++local)
				if(local!=opposite) face[entry++]=mesh.cells[cell].nodes[local];
			std::sort(face.begin(),face.end());
			uses[face].push_back({cell,opposite});
		}
	}
	std::map<Face,Owner> result;
	for(const auto& item:uses) {
		if(item.second.size()>2)
			throw std::invalid_argument("matching interface has a nonmanifold tetrahedral face");
		if(item.second.size()==1) result.emplace(item.first,item.second[0]);
	}
	return result;
}

inline GeometricFace Coordinates(const NativeTetMesh& mesh,const Face& nodes)
{
	GeometricFace result{};
	for(std::size_t vertex=0;vertex<3;++vertex) {
		result[vertex]=mesh.points.at(nodes[vertex]);
		for(const double value:result[vertex])
			if(!std::isfinite(value))
				throw std::invalid_argument("matching interface coordinate is nonfinite");
	}
	std::sort(result.begin(),result.end());
	return result;
}

inline Point OutwardAreaVector(const NativeTetMesh& mesh,const Face& nodes,
	const Owner& owner)
{
	const auto& a=mesh.points.at(nodes[0]);
	const auto& b=mesh.points.at(nodes[1]);
	const auto& c=mesh.points.at(nodes[2]);
	const auto& interior=mesh.points.at(mesh.cells[owner.cell].nodes[owner.opposite]);
	Point value{{
		0.5*((b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1])),
		0.5*((b[2]-a[2])*(c[0]-a[0])-(b[0]-a[0])*(c[2]-a[2])),
		0.5*((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]))}};
	double toward_interior=0.,squared_area=0.;
	for(int axis=0;axis<3;++axis) {
		toward_interior+=value[axis]*(interior[axis]-a[axis]);
		squared_area+=value[axis]*value[axis];
	}
	if(!(squared_area>0.)||!std::isfinite(squared_area)||!std::isfinite(toward_interior)
			||toward_interior==0.)
		throw std::invalid_argument("matching interface facet has invalid geometry");
	if(toward_interior>0.) for(double& component:value) component=-component;
	return value;
}

inline std::map<GeometricFace,Facet> LabelledFacets(const NativeTetMesh& mesh,int label,
	const std::map<Face,Owner>& owners)
{
	std::map<GeometricFace,Facet> result;
	for(const auto& triangle:mesh.boundary_triangles) {
		if(triangle.boundary_label!=label) continue;
		Face nodes=triangle.nodes;
		std::sort(nodes.begin(),nodes.end());
		const auto owner=owners.find(nodes);
		if(owner==owners.end())
			throw std::invalid_argument("matching interface triangle lacks one boundary owner");
		const auto coordinates=Coordinates(mesh,nodes);
		const Facet facet{owner->second.cell,nodes,
			OutwardAreaVector(mesh,nodes,owner->second)};
		if(!result.emplace(coordinates,facet).second)
			throw std::invalid_argument("matching interface has a duplicate geometric facet");
	}
	if(result.empty())
		throw std::invalid_argument("matching interface boundary label is absent");
	return result;
}

inline double VesselFacetFlow(const NativeTetMesh& mesh,const std::vector<double>& state,
	const Facet& facet,const std::map<std::array<std::uint32_t,2>,std::size_t>& edges)
{
	Point mid_edge_velocity_sum{{0.,0.,0.}};
	for(const auto pair:std::array<std::array<int,2>,3>{{{{0,1}},{{1,2}},{{2,0}}}}) {
		std::array<std::uint32_t,2> edge{{facet.nodes[pair[0]],facet.nodes[pair[1]]}};
		std::sort(edge.begin(),edge.end());
		const auto found=edges.find(edge);
		if(found==edges.end())
			throw std::invalid_argument("matching interface P2 edge is absent");
		const std::size_t node=mesh.points.size()+found->second;
		for(int axis=0;axis<3;++axis)
			mid_edge_velocity_sum[axis]+=state[3*node+axis];
	}
	// P2 vertex basis integrals vanish; each P2 mid-edge basis integrates to A/3.
	double flow=0.;
	for(int axis=0;axis<3;++axis)
		flow+=mid_edge_velocity_sum[axis]*facet.outward_area_vector[axis]/3.;
	if(!std::isfinite(flow))
		throw std::invalid_argument("matching interface facet flow is nonfinite");
	return flow;
}

} // namespace native_tet_matching_detail

inline NativeTetMatchingInterfaceSourceResult MapNativeTetMatchingInterfaceToTissueSource(
	const NativeTetMesh& vessel_mesh,const NativeTaylorHoodTopology& topology,
	const std::vector<double>& vessel_state,const NativeTetMesh& tissue_mesh,
	const std::vector<NativeTetMatchingInterfacePort>& ports)
{
	using namespace native_tet_matching_detail;
	if(ports.empty()||tissue_mesh.cells.empty())
		throw std::invalid_argument("matching interface requires ports and tissue cells");
	const std::size_t velocity_nodes=vessel_mesh.points.size()+topology.edges.size();
	if(vessel_state.size()!=3*velocity_nodes+vessel_mesh.points.size())
		throw std::invalid_argument("matching interface vessel state size is invalid");
	for(const double value:vessel_state)
		if(!std::isfinite(value))
			throw std::invalid_argument("matching interface vessel state is nonfinite");
	std::map<std::array<std::uint32_t,2>,std::size_t> edges;
	for(std::size_t index=0;index<topology.edges.size();++index) {
		auto edge=topology.edges[index];std::sort(edge.begin(),edge.end());
		if(!edges.emplace(edge,index).second)
			throw std::invalid_argument("matching interface P2 topology has duplicate edges");
	}
	std::vector<double> tissue_volumes(tissue_mesh.cells.size());
	for(std::size_t cell=0;cell<tissue_mesh.cells.size();++cell)
		tissue_volumes[cell]=EvaluateNativeTetGeometry(tissue_mesh,
			tissue_mesh.cells[cell]).determinant/6.;
	const auto vessel_owners=BoundaryOwners(vessel_mesh);
	const auto tissue_owners=BoundaryOwners(tissue_mesh);
	NativeTetMatchingInterfaceSourceResult result;
	result.source.tissue_source_s_inv.assign(tissue_mesh.cells.size(),0.);
	std::vector<double> tissue_cell_flow_m3_s(tissue_mesh.cells.size(),0.);
	std::set<int> vessel_labels,tissue_labels;
	std::set<GeometricFace> used_geometric_facets;
	for(const auto& port:ports) {
		if(port.name.empty()||!vessel_labels.insert(port.vessel_boundary_label).second
				||!tissue_labels.insert(port.tissue_boundary_label).second
				||result.source.vessel_outward_flow_m3_s.count(port.name))
			throw std::invalid_argument("matching interface port name or label is duplicate");
		const auto vessel=LabelledFacets(vessel_mesh,port.vessel_boundary_label,
			vessel_owners);
		const auto tissue=LabelledFacets(tissue_mesh,port.tissue_boundary_label,
			tissue_owners);
		if(vessel.size()!=tissue.size())
			throw std::invalid_argument("matching interface facet counts differ");
		double port_flow=0.,absolute_face_flow_sum=0.;
		std::map<std::size_t,double> port_cell_flow;
		for(const auto& item:vessel) {
			if(!used_geometric_facets.insert(item.first).second)
				throw std::invalid_argument("matching interface geometric facet is reused by two ports");
			const auto found=tissue.find(item.first);
			if(found==tissue.end())
				throw std::invalid_argument("matching interface has a missing tissue facet");
			const auto& v=item.second.outward_area_vector;
			const auto& t=found->second.outward_area_vector;
			double dot=0.,v2=0.,t2=0.;
			for(int axis=0;axis<3;++axis) {
				dot+=v[axis]*t[axis];v2+=v[axis]*v[axis];t2+=t[axis]*t[axis];
			}
			if(!(dot<0.)||std::abs(dot+std::sqrt(v2*t2))>1e-12*std::sqrt(v2*t2))
				throw std::invalid_argument("matching interface outward normals are not opposed");
			const double flow=VesselFacetFlow(vessel_mesh,vessel_state,
				item.second,edges);
			result.facet_transfers.push_back({port.name,
				vessel_mesh.cells[item.second.cell].id,
				tissue_mesh.cells[found->second.cell].id,item.first,flow});
			port_cell_flow[found->second.cell]+=flow;
			port_flow+=flow;
			absolute_face_flow_sum+=std::abs(flow);
		}
		double distributed=0.;
		for(const auto& cell:port_cell_flow) {
			tissue_cell_flow_m3_s[cell.first]+=cell.second;
			distributed+=cell.second;
		}
		const double defect=std::abs(port_flow-distributed);
		if(!std::isfinite(port_flow)||!std::isfinite(distributed)
				||!std::isfinite(absolute_face_flow_sum)
				||defect>1e-12*std::max(1e-12,absolute_face_flow_sum))
			throw std::invalid_argument("matching interface port flow or balance is invalid");
		result.source.vessel_outward_flow_m3_s.emplace(port.name,port_flow);
		result.source.tissue_inward_flow_m3_s.emplace(port.name,distributed);
		result.source.total_vessel_outward_flow_m3_s+=port_flow;
		result.source.total_tissue_inward_flow_m3_s+=distributed;
		result.source.maximum_port_balance_defect_m3_s=std::max(
			result.source.maximum_port_balance_defect_m3_s,defect);
		result.matched_facets.emplace(port.name,vessel.size());
	}
	for(std::size_t cell=0;cell<tissue_cell_flow_m3_s.size();++cell) {
		result.source.tissue_source_s_inv[cell]=tissue_cell_flow_m3_s[cell]
			/tissue_volumes[cell];
		if(!std::isfinite(result.source.tissue_source_s_inv[cell]))
			throw std::invalid_argument("matching interface tissue source is nonfinite");
	}
	if(!std::isfinite(result.source.total_vessel_outward_flow_m3_s)
			||!std::isfinite(result.source.total_tissue_inward_flow_m3_s))
		throw std::invalid_argument("matching interface total flow is nonfinite");
	return result;
}

} // namespace iga

#endif
