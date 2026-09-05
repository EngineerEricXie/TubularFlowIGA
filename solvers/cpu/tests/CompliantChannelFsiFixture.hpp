#ifndef IGA_COMPLIANT_CHANNEL_FSI_FIXTURE_HPP
#define IGA_COMPLIANT_CHANNEL_FSI_FIXTURE_HPP

// Shared, deliberately small production FSI channel: a closed .6 m cube with
// a 3x3 top-wall patch (label 7), eight clamped edge nodes, and one free
// centre node.  Keep the material patch, surface layout, and moving-flow
// setup together so runtime and coupled benchmark coverage cannot drift.
#include "MaterialSurfacePatchMap.hpp"
#include "MovingImmersedTransientFlowRuntime.hpp"
#include "PrescribedSurfaceMotion.hpp"

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace iga {
namespace compliant_channel_fixture {

inline RawSurfaceTriangle Face(int a, int b, int c, std::uint32_t label)
{ RawSurfaceTriangle value; value.indices={{a,b,c}}; value.boundary_id=label; return value; }

inline MaterialSurfaceKinematics InitialMaterial()
{
	RawSurfaceSoup soup;
	std::vector<std::array<int,3>> lattice;
	auto node=[&](int x,int y,int z) {
		const std::array<int,3> key{{x,y,z}};
		for(std::size_t i=0;i<lattice.size();++i) if(lattice[i]==key) return static_cast<int>(i);
		lattice.push_back(key); soup.vertices.push_back({{.2+.3*x,.2+.3*y,.2+.3*z}});
		return static_cast<int>(lattice.size()-1);
	};
	auto add_face=[&](std::uint32_t label, auto coordinate) {
		int grid[3][3];
		for(int i=0;i<3;++i) for(int j=0;j<3;++j) { const auto p=coordinate(i,j); grid[i][j]=node(p[0],p[1],p[2]); }
		for(int i=0;i<2;++i) for(int j=0;j<2;++j) {
			const int a=grid[i][j],b=grid[i+1][j],c=grid[i][j+1],d=grid[i+1][j+1];
			auto add_outward=[&](int q0,int q1,int q2) {
				const auto& p0=soup.vertices[q0]; const auto& p1=soup.vertices[q1]; const auto& p2=soup.vertices[q2];
				const std::array<double,3> normal{{(p1[1]-p0[1])*(p2[2]-p0[2])-(p1[2]-p0[2])*(p2[1]-p0[1]),(p1[2]-p0[2])*(p2[0]-p0[0])-(p1[0]-p0[0])*(p2[2]-p0[2]),(p1[0]-p0[0])*(p2[1]-p0[1])-(p1[1]-p0[1])*(p2[0]-p0[0])}};
				const std::array<double,3> radial{{(p0[0]+p1[0]+p2[0])/3.-.5,(p0[1]+p1[1]+p2[1])/3.-.5,(p0[2]+p1[2]+p2[2])/3.-.5}};
				if(normal[0]*radial[0]+normal[1]*radial[1]+normal[2]*radial[2]<0.) std::swap(q1,q2);
				soup.triangles.push_back(Face(q0,q1,q2,label));
			};
			add_outward(a,b,d); add_outward(a,d,c);
		}
	};
	add_face(7,[](int x,int y){return std::array<int,3>{{x,y,2}};});
	add_face(9,[](int x,int y){return std::array<int,3>{{x,y,0}};});
	add_face(1,[](int y,int z){return std::array<int,3>{{0,y,z}};});
	add_face(2,[](int y,int z){return std::array<int,3>{{2,y,z}};});
	add_face(9,[](int x,int z){return std::array<int,3>{{x,0,z}};});
	add_face(9,[](int x,int z){return std::array<int,3>{{x,2,z}};});
	return PrescribedSurfaceMotion({{0.0,soup},{1.0,soup}}).Evaluate(1.0,0.0,1.0);
}

inline DistributedSurfaceInterface Structure(const std::string& ref)
{
	DistributedSurfaceInterface v; v.id={"structure","membrane","patch"}; v.subsystem_id="membrane";
	v.boundary_labels={7}; v.reference_mesh_identity_sha256=ref;
	v.provides={SurfaceFieldQuantity::Displacement,SurfaceFieldQuantity::Velocity};
	v.requires={SurfaceFieldQuantity::TractionOnStructure}; return v;
}
inline DistributedSurfaceInterface Fluid(const std::string& ref)
{
	DistributedSurfaceInterface v; v.id={"fluid","immersed","patch"}; v.subsystem_id="immersed";
	v.boundary_labels={7}; v.reference_mesh_identity_sha256=ref;
	v.provides={SurfaceFieldQuantity::TractionOnStructure};
	v.requires={SurfaceFieldQuantity::Displacement,SurfaceFieldQuantity::Velocity}; return v;
}
inline MaterialSurfacePatchMap PatchMap(const MaterialSurfaceKinematics& initial)
{
	DistributedSurfaceLayout layout; layout.global_node_count=9; layout.partition_count=1; layout.partition_rank=0;
	layout.owned_global_node_ids={11,12,13,14,15,16,17,18,19};
	for(std::size_t i=0;i<9;++i) layout.reference_positions.push_back({layout.owned_global_node_ids[i],initial.ReferenceMaterialVerticesM()[i]});
	for(std::uint64_t x=0;x<2;++x) for(std::uint64_t y=0;y<2;++y) { const auto a=11+3*x+y,b=a+3,c=a+1,d=b+1; layout.reference_triangles.push_back({{a,b,d}}); layout.reference_triangles.push_back({{a,d,c}}); }
	layout.owned_reference_lumped_areas_m2={.03,.045,.015,.045,.09,.045,.015,.045,.03};
	std::vector<MaterialSurfacePatchMap::GlobalToSourceVertex> vertices; for(std::size_t i=0;i<9;++i) vertices.push_back({layout.owned_global_node_ids[i],static_cast<std::uint32_t>(i)});
	const std::vector<std::uint32_t> triangles{0,1,2,3,4,5,6,7}; const std::vector<std::uint64_t> clamps{11,12,13,14,16,17,18,19};
	layout.reference_mesh_identity_sha256=MaterialSurfacePatchMap::BuildReferenceIdentitySha256(initial,7,vertices,layout.reference_triangles,triangles,clamps);
	layout.layout_identity_sha256=BuildDistributedSurfaceLayoutIdentitySha256(layout);
	return MaterialSurfacePatchMap::Create(Structure(layout.reference_mesh_identity_sha256),layout,initial,7,vertices,triangles,clamps);
}
inline MovingImmersedTransientFlowOptions FlowOptions()
{
	MovingImmersedTransientFlowOptions v; v.grid={{{0,0,0}},{{1,1,1}},{{3,3,3}}}; v.geometry.volume.max_depth=3;
	v.geometry.volume.max_nodes=v.geometry.volume.max_leaves=v.geometry.volume.max_points=1000000;
	v.flow.parameters={1.,1.,1.}; v.flow.wall_labels={7,9};
	v.flow.ports={{"inlet",1,ImmersedFlowPortControlMode::Pressure,.10},{"outlet",2,ImmersedFlowPortControlMode::Pressure,0.}};
	v.flow.ksp_relative_tolerance=1e-12; v.flow.lu_pivot_shift=1e-20; return v;
}

} // namespace compliant_channel_fixture
} // namespace iga

#endif
